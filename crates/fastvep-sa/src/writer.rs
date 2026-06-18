//! Writer for .osa position/allele-level annotation files.
//!
//! Records must be added in chromosome-sorted, position-sorted order.
//! The writer accumulates entries into blocks, compresses them, and writes
//! to the data file while building the index.

use crate::block::{BlockEntry, SaBlock};
use crate::common::{AnnotationRecord, DEFAULT_BLOCK_SIZE, OSA_MAGIC, SCHEMA_VERSION};
use crate::index::{BlockRef, IndexHeader, SaIndex};
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Blocks compressed per parallel window. Enough to keep every core busy without
/// holding the whole source's compressed output in memory at once.
const COMPRESS_WINDOW_BLOCKS: usize = 64;

/// A block partitioned off and ready to be compressed, retaining the index
/// metadata that must be recorded once its compressed length is known.
struct PendingBlock {
    chrom: String,
    start_pos: u32,
    end_pos: u32,
    block: SaBlock,
}

/// Push the current block onto `blocks` (capturing its index metadata) and
/// reset it to empty. No-op when the block holds no entries.
fn finalize_block(block: &mut SaBlock, chrom: &str, blocks: &mut Vec<PendingBlock>) {
    if block.is_empty() {
        return;
    }
    let start_pos = block.start_position().unwrap();
    let end_pos = block.end_position().unwrap();
    let finished = std::mem::replace(block, SaBlock::new(DEFAULT_BLOCK_SIZE));
    blocks.push(PendingBlock {
        chrom: chrom.to_string(),
        start_pos,
        end_pos,
        block: finished,
    });
}

/// Builds an .osa data file and its .osa.idx index file.
pub struct SaWriter {
    index: SaIndex,
    block: SaBlock,
    current_chrom: Option<String>,
    last_key: Option<(u16, u32)>,
    /// Chromosome name -> numeric index mapping.
    chrom_names: Vec<String>,
    data_offset: u64,
}

impl SaWriter {
    pub fn new(header: IndexHeader) -> Self {
        Self {
            index: SaIndex::new(header),
            block: SaBlock::new(DEFAULT_BLOCK_SIZE),
            current_chrom: None,
            last_key: None,
            chrom_names: Vec::new(),
            data_offset: 0,
        }
    }

    /// Build .osa and .osa.idx from an iterator of sorted annotation records.
    ///
    /// Records MUST be sorted by (chrom_idx, position).
    /// `chrom_map` maps chrom_idx -> chromosome name string.
    ///
    /// The records are partitioned into blocks sequentially (preserving order
    /// and the block-size/chromosome boundaries), the blocks are zstd-compressed
    /// in parallel via rayon (compression is CPU-bound and per-block
    /// independent), then written and indexed sequentially so the on-disk layout
    /// is byte-identical to a single-threaded build. All records are
    /// materialised into blocks first; callers that must stream without holding
    /// every record in memory should use [`Self::write_all_results`].
    pub fn write_all<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = AnnotationRecord>,
        chrom_map: &[String],
    ) -> Result<()> {
        self.chrom_names = chrom_map.to_vec();
        let blocks = Self::collect_blocks(records.map(Ok), chrom_map)?;
        self.write_compressed_blocks(data_writer, blocks, COMPRESS_WINDOW_BLOCKS)
    }

    /// Partition sorted records into blocks without compressing them, applying
    /// the same sort validation and block/chromosome boundary rules as the
    /// streaming path.
    fn collect_blocks(
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<Vec<PendingBlock>> {
        let mut blocks: Vec<PendingBlock> = Vec::new();
        let mut current_chrom: Option<String> = None;
        let mut last_key: Option<(u16, u32)> = None;
        let mut block = SaBlock::new(DEFAULT_BLOCK_SIZE);

        for record in records {
            let record = record?;
            if let Some((last_chrom, last_pos)) = last_key {
                if (record.chrom_idx, record.position) < (last_chrom, last_pos) {
                    anyhow::bail!(
                        "SA records are not sorted: previous chrom_idx={}, position={}; current chrom_idx={}, position={}. \
                         The streaming .osa builder requires input sorted by chromosome (chr1..chr22,X,Y,M) then position \
                         — sort the source file (e.g. `bcftools sort` / `sort -k1,1 -k2,2n`) and rebuild.",
                        last_chrom,
                        last_pos,
                        record.chrom_idx,
                        record.position
                    );
                }
            }
            last_key = Some((record.chrom_idx, record.position));

            let chrom_name = &chrom_map[record.chrom_idx as usize];
            if current_chrom.as_deref() != Some(chrom_name.as_str()) {
                finalize_block(
                    &mut block,
                    current_chrom.as_deref().unwrap_or(""),
                    &mut blocks,
                );
                current_chrom = Some(chrom_name.clone());
            }

            let entry = BlockEntry {
                position: record.position,
                ref_allele: record.ref_allele,
                alt_allele: record.alt_allele,
                json: record.json,
            };

            if !block.add(entry.clone()) {
                finalize_block(&mut block, current_chrom.as_deref().unwrap(), &mut blocks);
                assert!(block.add(entry), "Single entry exceeds block size");
            }
        }

        finalize_block(
            &mut block,
            current_chrom.as_deref().unwrap_or(""),
            &mut blocks,
        );
        Ok(blocks)
    }

    /// Compress all blocks in parallel, then write and index them sequentially.
    fn write_compressed_blocks<W: Write>(
        &mut self,
        data_writer: &mut W,
        blocks: Vec<PendingBlock>,
        window_blocks: usize,
    ) -> Result<()> {
        data_writer.write_all(OSA_MAGIC)?;
        data_writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        self.data_offset = (OSA_MAGIC.len() + 2) as u64;

        // Compress a bounded window at a time rather than the whole source at
        // once: holding every compressed block alongside every uncompressed one
        // doubles peak memory on a genome-scale source, and these callers have
        // already materialised the records.
        for window in blocks.chunks(window_blocks) {
            // Booting rayon's global pool spawns a thread per core, which is
            // not worth paying to run a single compression.
            let compressed: Vec<Vec<u8>> = match window {
                [only] => vec![only.block.compress()?],
                _ => window
                    .par_iter()
                    .map(|pending| pending.block.compress())
                    .collect::<Result<Vec<_>>>()?,
            };
            self.write_window(data_writer, window, compressed)?;
        }
        Ok(())
    }

    /// Write and index one already-compressed window, in order.
    fn write_window<W: Write>(
        &mut self,
        data_writer: &mut W,
        window: &[PendingBlock],
        compressed: Vec<Vec<u8>>,
    ) -> Result<()> {
        for (pending, data) in window.iter().zip(compressed) {
            let compressed_len = data.len() as u32;
            data_writer.write_all(&compressed_len.to_le_bytes())?;
            data_writer.write_all(&data)?;
            self.index.add_block(
                &pending.chrom,
                BlockRef {
                    start_pos: pending.start_pos,
                    end_pos: pending.end_pos,
                    file_offset: self.data_offset,
                    compressed_len,
                },
            );
            self.data_offset += 4 + compressed_len as u64;
        }
        Ok(())
    }

    /// Build .osa and .osa.idx from an iterator that can surface parse errors.
    ///
    /// Records MUST be sorted by (chrom_idx, position).
    pub fn write_all_results<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<()> {
        self.chrom_names = chrom_map.to_vec();

        data_writer.write_all(OSA_MAGIC)?;
        data_writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        self.data_offset = (OSA_MAGIC.len() + 2) as u64;

        for record in records {
            self.write_record(data_writer, record?, chrom_map)?;
        }

        self.flush_block(data_writer)?;
        Ok(())
    }

    fn write_record<W: Write>(
        &mut self,
        data_writer: &mut W,
        record: AnnotationRecord,
        chrom_map: &[String],
    ) -> Result<()> {
        if let Some((last_chrom, last_pos)) = self.last_key {
            if (record.chrom_idx, record.position) < (last_chrom, last_pos) {
                anyhow::bail!(
                    "SA records are not sorted: previous chrom_idx={}, position={}; current chrom_idx={}, position={}. \
                     The streaming .osa builder requires input sorted by chromosome (chr1..chr22,X,Y,M) then position \
                     — sort the source file (e.g. `bcftools sort` / `sort -k1,1 -k2,2n`) and rebuild.",
                    last_chrom,
                    last_pos,
                    record.chrom_idx,
                    record.position
                );
            }
        }
        self.last_key = Some((record.chrom_idx, record.position));

        let chrom_name = &chrom_map[record.chrom_idx as usize];

        // If we've moved to a new chromosome, flush the current block
        if self.current_chrom.as_ref() != Some(chrom_name) {
            self.flush_block(data_writer)?;
            self.current_chrom = Some(chrom_name.clone());
        }

        let entry = BlockEntry {
            position: record.position,
            ref_allele: record.ref_allele,
            alt_allele: record.alt_allele,
            json: record.json,
        };

        if !self.block.add(entry.clone()) {
            // Block is full, flush and retry
            self.flush_block(data_writer)?;
            assert!(self.block.add(entry), "Single entry exceeds block size");
        }

        Ok(())
    }

    fn flush_block<W: Write>(&mut self, writer: &mut W) -> Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }

        let chrom = self.current_chrom.as_ref().unwrap().clone();
        let start_pos = self.block.start_position().unwrap();
        let end_pos = self.block.end_position().unwrap();

        let compressed = self.block.compress()?;
        let compressed_len = compressed.len() as u32;

        // Write compressed block length prefix + data
        writer.write_all(&compressed_len.to_le_bytes())?;
        writer.write_all(&compressed)?;

        self.index.add_block(
            &chrom,
            BlockRef {
                start_pos,
                end_pos,
                file_offset: self.data_offset,
                compressed_len,
            },
        );

        self.data_offset += 4 + compressed_len as u64;
        self.block.clear();
        Ok(())
    }

    /// Write the index file.
    pub fn write_index<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.index.write_to(writer)
    }

    /// Convenience: write .osa and .osa.idx to files at the given base path.
    pub fn write_to_files(
        &mut self,
        base_path: &Path,
        records: impl Iterator<Item = AnnotationRecord>,
        chrom_map: &[String],
    ) -> Result<()> {
        let data_path = base_path.with_extension("osa");
        let idx_path = base_path.with_extension("osa.idx");

        let data_file = std::fs::File::create(&data_path).with_context(|| {
            format!(
                "Creating output file {} (does the output directory exist?)",
                data_path.display()
            )
        })?;
        let mut data_writer = BufWriter::new(data_file);
        self.write_all(&mut data_writer, records, chrom_map)?;
        data_writer.flush()?;

        let idx_file = std::fs::File::create(&idx_path)
            .with_context(|| format!("Creating index file {}", idx_path.display()))?;
        let mut idx_writer = BufWriter::new(idx_file);
        self.write_index(&mut idx_writer)?;
        idx_writer.flush()?;

        Ok(())
    }

    /// Convenience: write .osa and .osa.idx to files from fallible records.
    pub fn write_results_to_files(
        &mut self,
        base_path: &Path,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<()> {
        let data_path = base_path.with_extension("osa");
        let idx_path = base_path.with_extension("osa.idx");

        let data_file = std::fs::File::create(&data_path).with_context(|| {
            format!(
                "Creating output file {} (does the output directory exist?)",
                data_path.display()
            )
        })?;
        let mut data_writer = BufWriter::new(data_file);
        self.write_all_results(&mut data_writer, records, chrom_map)?;
        data_writer.flush()?;

        let idx_file = std::fs::File::create(&idx_path)
            .with_context(|| format!("Creating index file {}", idx_path.display()))?;
        let mut idx_writer = BufWriter::new(idx_file);
        self.write_index(&mut idx_writer)?;
        idx_writer.flush()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> IndexHeader {
        IndexHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test".into(),
            version: "test".into(),
            description: "test".into(),
            assembly: "GRCh38".into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        }
    }

    fn record(chrom_idx: u16, position: u32) -> AnnotationRecord {
        AnnotationRecord {
            chrom_idx,
            position,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        }
    }

    #[test]
    fn write_all_rejects_unsorted_records() {
        let mut writer = SaWriter::new(header());
        let mut out = Vec::new();
        let err = writer
            .write_all(
                &mut out,
                vec![record(0, 20), record(0, 10)].into_iter(),
                &["1".into()],
            )
            .unwrap_err();

        assert!(err.to_string().contains("SA records are not sorted"));
    }

    #[test]
    fn write_all_is_byte_identical_whatever_the_compression_window() {
        // Blocks are compressed a window at a time, so the on-disk layout must
        // not depend on that window or a build's output would change with a
        // tuning constant. One block at a time against the real window, over a
        // source spanning several windows.
        //
        // The data stream is compared byte for byte. The index is compared field
        // by field rather than through its serialized form: `SaIndex.chromosomes`
        // is a `HashMap`, so the encoded byte order varies between two maps in
        // the same process regardless of what was inserted.
        let chroms: Vec<String> = (0..COMPRESS_WINDOW_BLOCKS * 2 + 5)
            .map(|i| i.to_string())
            .collect();
        let records: Vec<AnnotationRecord> = (0..chroms.len())
            .map(|i| record(i as u16, 100 + i as u32))
            .collect();

        type Layout = Vec<(String, Vec<(u32, u32, u64, u32)>)>;
        let run = |window: usize| -> (Vec<u8>, Layout) {
            let mut writer = SaWriter::new(header());
            let mut out = Vec::new();
            let blocks =
                SaWriter::collect_blocks(records.clone().into_iter().map(Ok), &chroms).unwrap();
            writer
                .write_compressed_blocks(&mut out, blocks, window)
                .unwrap();
            let mut layout: Layout = writer
                .index
                .chromosomes
                .iter()
                .map(|(chrom, refs)| {
                    let refs = refs
                        .iter()
                        .map(|r| (r.start_pos, r.end_pos, r.file_offset, r.compressed_len))
                        .collect();
                    (chrom.clone(), refs)
                })
                .collect();
            layout.sort();
            (out, layout)
        };

        assert!(
            chroms.len() > COMPRESS_WINDOW_BLOCKS,
            "test must span more than one window"
        );
        let (sequential, seq_layout) = run(1);
        let (windowed, win_layout) = run(COMPRESS_WINDOW_BLOCKS);

        assert_eq!(sequential, windowed, "data stream differs by window size");
        assert_eq!(seq_layout, win_layout, "index differs by window size");
    }

    #[test]
    fn write_all_preserves_order_across_blocks() {
        // Each chromosome change starts a new block, so this produces three
        // blocks that are compressed in parallel and must be written back in
        // the original order. Decompressing the data stream must reproduce the
        // input records in order.
        use crate::block::SaBlock;

        let mut writer = SaWriter::new(header());
        let mut out = Vec::new();
        let records = vec![record(0, 10), record(0, 20), record(1, 5), record(2, 100)];
        writer
            .write_all(
                &mut out,
                records.into_iter(),
                &["1".into(), "2".into(), "3".into()],
            )
            .unwrap();

        // Walk the data stream: magic + schema version, then length-prefixed
        // compressed blocks.
        let mut cursor = OSA_MAGIC.len() + 2;
        let mut positions = Vec::new();
        let mut block_count = 0;
        while cursor < out.len() {
            let len = u32::from_le_bytes(out[cursor..cursor + 4].try_into().unwrap()) as usize;
            cursor += 4;
            let entries = SaBlock::decompress(&out[cursor..cursor + len]).unwrap();
            positions.extend(entries.iter().map(|e| e.position));
            cursor += len;
            block_count += 1;
        }

        assert_eq!(block_count, 3, "one block per chromosome");
        assert_eq!(positions, vec![10, 20, 5, 100]);
    }
}
