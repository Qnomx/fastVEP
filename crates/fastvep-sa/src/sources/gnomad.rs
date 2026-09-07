//! gnomAD VCF parser for building .osa annotation files.
//!
//! Parses gnomAD's sites-only VCF to extract allele frequencies
//! per population, allele counts, and homozygote counts.
//!
//! Auto-detects the INFO field-naming convention used by the input:
//!
//! | Release | Top-level | Per-population |
//! |---------|-----------|----------------|
//! | v2.x (exomes/genomes) | `AF` / `AC` / `AN` / `nhomalt` | `AF_afr` (underscore) |
//! | v3.x (genomes)        | `AF` / `AC` / `AN` / `nhomalt` | none at top level (only subset-qualified, e.g. `AF-non_neuro-afr`) |
//! | v4.0/v4.1 (exomes/genomes) | `AF` / `AC` / `AN` / `nhomalt` | `AF_afr` (underscore) |
//! | v4.1 joint            | `AF_joint` / `AC_joint` / `AN_joint` / `nhomalt_joint` | `AF_joint_afr` |
//!
//! Some locally-processed files expose simple per-population fields with a
//! hyphen separator (`AF-afr`). The parser detects and handles this case.
//!
//! We pick the scheme by scanning `##INFO=<ID=...>` header lines.

use crate::common::AnnotationRecord;
use crate::fields::{Field, FieldType};
use crate::writer_v2::{Osa2Metadata, Osa2Record};
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet, VecDeque};
use std::io::BufRead;

/// Population keys to extract from gnomAD VCF INFO field. Covers both
/// gnomAD v2.1 codes (`oth`) and the v4.1 codes (`mid`, `remaining`); a
/// missing key is silently skipped per VCF, so listing all is harmless.
const POPULATIONS: &[&str] = &[
    "afr",
    "amr",
    "asj",
    "eas",
    "fin",
    "mid",
    "nfe",
    "oth",
    "remaining",
    "sas",
];

// =============================================================================
// Extended quality-control and stratified columns
// =============================================================================
//
// Everything below is additive: a database built before these columns existed
// simply does not carry the keys, and every consumer treats an absent key the
// same way it behaved before the columns were introduced. Both the v1 JSON
// builder and the v2 u32 encoder walk the same descriptor tables, in the same
// order, so the two output formats cannot drift apart.

/// gnomAD FILTER-column entries, as (VCF filter name, output alias, description).
///
/// A record carrying any of these is not a trustworthy basis for benign
/// frequency evidence: `AC0` means every genotype backing the call was filtered
/// out, `AS_VQSR` that the site failed the variant-quality model, and
/// `InbreedingCoeff` that the genotypes are distributed in a way real
/// population data is not. Stored as one flag per filter rather than a single
/// "not PASS" bit so a reviewer reading the annotation can see *which* filter
/// fired.
const FILTER_FLAGS: &[(&str, &str, &str)] = &[
    (
        "AC0",
        "filterAc0",
        "FILTER=AC0: allele count is zero after filtering out low-confidence genotypes",
    ),
    (
        "AS_VQSR",
        "filterVqsr",
        "FILTER=AS_VQSR: site failed allele-specific VQSR thresholds",
    ),
    (
        "InbreedingCoeff",
        "filterInbreeding",
        "FILTER=InbreedingCoeff: inbreeding coefficient < -0.3",
    ),
];

/// Region-quality INFO flags, as (INFO key, output alias, description).
///
/// `lcr` and `segdup` mark the low-complexity and segmentally-duplicated
/// regions where short-read allele frequencies are systematically unreliable.
/// This is the same concern as the hand-curated homologous-gene list
/// (Mandelker 2016, PMID 27228465), but resolved per site rather than per gene.
/// `non_par` is required to interpret the XY-stratified counts: only outside a
/// pseudoautosomal region is an XY sample hemizygous.
const INFO_FLAGS: &[(&str, &str, &str)] = &[
    ("lcr", "lcr", "Variant falls within a low-complexity region"),
    (
        "segdup",
        "segdup",
        "Variant falls within a segmental duplication",
    ),
    (
        "non_par",
        "nonPar",
        "Sex-chromosome variant falls outside a pseudoautosomal region",
    ),
];

/// Per-allele (`Number=A`) XY-stratified counts.
///
/// On the non-PAR regions of chrX and chrY, gnomAD calls XY samples haploid, so
/// `AC_XY` is the count of hemizygous individuals. That is the observation ACMG
/// BS2 asks for on an X-linked gene, and it is invisible in `nhomalt`.
///
/// `nhomalt_XY` is deliberately *not* extracted. Because those calls are
/// haploid, gnomAD never records an XY sample as homozygous: across all 6,955
/// non-PAR chrX records in the IDS region it is zero, including at a site with
/// 109,916 XY carriers. Inside the PAR, XY samples are diploid and the global
/// `nhomalt` already counts them. The field therefore has no reading under
/// which it adds information.
const XY_ALLELE_INTS: &[(&str, &str, &str)] = &[(
    "AC_XY",
    "allAcXY",
    "Alternate allele count in XY samples (hemizygous individuals outside the PAR)",
)];

/// Per-site (`Number=1`) XY-stratified counts.
const XY_SITE_INTS: &[(&str, &str, &str)] =
    &[("AN_XY", "allAnXY", "Total allele number in XY samples")];

/// Filtering allele frequencies (Whiffin 2017, PMID 28518168).
///
/// The FAF is the lower bound of the 95 % confidence interval on the allele
/// frequency, computed per genetic-ancestry group. `fafmax_faf95_max` is the
/// maximum across groups, which is the statistic ACMG BA1/BS1 actually want:
/// it is robust to a founder variant that is common in one population and
/// absent elsewhere, and to a frequency estimated from very few alleles. Both
/// failure modes were raised in the round-2 medical-genetics review.
const FAF_FLOATS: &[(&str, &str, &str)] = &[
    (
        "faf95",
        "faf95",
        "Filtering allele frequency (95% CI lower bound), all samples",
    ),
    (
        "fafmax_faf95_max",
        "faf95Max",
        "Maximum filtering allele frequency (95% CI lower bound) across genetic-ancestry groups",
    ),
    (
        "fafmax_faf99_max",
        "faf99Max",
        "Maximum filtering allele frequency (99% CI lower bound) across genetic-ancestry groups",
    ),
];

/// The genetic-ancestry group a grpmax FAF was drawn from, as (INFO key,
/// output alias, description).
///
/// `faf95Max` alone says how common the variant is in the group where it is
/// commonest; this says which group that is. A founder variant common in one
/// ancestry and absent elsewhere reads very differently from one at the same
/// frequency everywhere, and only this field separates them.
///
/// Stored as a categorical column over [`POPULATIONS`], which is the
/// vocabulary gnomAD draws these codes from. A code outside it encodes as
/// absent: the frequency remains, only the attribution is lost.
const FAF_CATEGORICALS: &[(&str, &str, &str)] = &[(
    "fafmax_faf95_max_gen_anc",
    "faf95MaxGenAnc",
    "Genetic-ancestry group the maximum filtering allele frequency was drawn from",
)];

/// A decoded value for one extended column, before it is rendered into either
/// output format.
#[derive(Debug, Clone, Copy, PartialEq)]
enum ExtValue {
    Int(Option<i64>),
    Float(Option<f64>),
    Flag(bool),
    /// A [`POPULATIONS`] index, for a categorical column. `None` when the
    /// value is absent or is a code outside that vocabulary.
    Category(Option<u32>),
}

/// Resolve a genetic-ancestry code to its [`POPULATIONS`] index.
fn population_index(value: &str) -> Option<u32> {
    POPULATIONS
        .iter()
        .position(|p| *p == value)
        .map(|i| i as u32)
}

/// Extract every extended column for one alternate allele, in schema order.
///
/// Shared by both encoders so the v1 `.osa` JSON and the v2 `.osa2` columns are
/// populated from exactly the same logic. Must stay in the same order as
/// [`extended_fields`]; `test_extended_fields_match_extended_values` enforces it.
fn extended_values(
    info_map: &HashMap<String, String>,
    filter_column: &str,
    allele_idx: usize,
) -> Vec<(&'static str, ExtValue)> {
    let mut out = Vec::with_capacity(
        XY_ALLELE_INTS.len()
            + XY_SITE_INTS.len()
            + FAF_FLOATS.len()
            + FAF_CATEGORICALS.len()
            + INFO_FLAGS.len()
            + FILTER_FLAGS.len(),
    );

    for (key, alias, _) in XY_ALLELE_INTS {
        let vals = split_info_values(info_map.get(*key).map(|s| s.as_str()));
        out.push((
            *alias,
            ExtValue::Int(vals.get(allele_idx).and_then(|s| s.parse::<i64>().ok())),
        ));
    }
    for (key, alias, _) in XY_SITE_INTS {
        let vals = split_info_values(info_map.get(*key).map(|s| s.as_str()));
        out.push((
            *alias,
            ExtValue::Int(vals.first().and_then(|s| s.parse::<i64>().ok())),
        ));
    }
    for (key, alias, _) in FAF_FLOATS {
        let vals = split_info_values(info_map.get(*key).map(|s| s.as_str()));
        out.push((
            *alias,
            ExtValue::Float(vals.get(allele_idx).and_then(|s| s.parse::<f64>().ok())),
        ));
    }
    for (key, alias, _) in FAF_CATEGORICALS {
        let vals = split_info_values(info_map.get(*key).map(|s| s.as_str()));
        out.push((
            *alias,
            ExtValue::Category(
                vals.get(allele_idx)
                    .map(|s| s.as_str())
                    .and_then(population_index),
            ),
        ));
    }
    for (key, alias, _) in INFO_FLAGS {
        out.push((*alias, ExtValue::Flag(info_map.contains_key(*key))));
    }
    for (name, alias, _) in FILTER_FLAGS {
        out.push((*alias, ExtValue::Flag(filter_has(filter_column, name))));
    }

    out
}

/// Whether the VCF FILTER column lists `name`.
///
/// FILTER is a semicolon-separated list, or `PASS` when nothing fired.
fn filter_has(filter_column: &str, name: &str) -> bool {
    filter_column.split(';').any(|f| f == name)
}

/// The FILTER value that means no filter fired.
const FILTER_PASS: &str = "PASS";

/// Reject a FILTER column outside the vocabulary [`FILTER_FLAGS`] models.
///
/// The gnomAD v4.0/v4.1 VCF header declares exactly four FILTER IDs — `AC0`,
/// `AS_VQSR`, `InbreedingCoeff`, `PASS` — and [`FILTER_FLAGS`] carries the
/// three non-`PASS` ones. Anything else encodes as three unset flags, which is
/// indistinguishable from `PASS`, and `PASS` is what permits benign frequency
/// evidence. So an unmodelled verdict is not a value to pass through: a `.` or
/// empty column records no verdict at all, and an unrecognised name means the
/// file is a release whose vocabulary differs — v2.1 used `RF`, `LCR` and
/// `SEGDUP`. Either way the flags this encoder would write are wrong, and the
/// build stops rather than producing them.
fn validate_filter_column(filter_column: &str) -> Result<()> {
    if filter_column.is_empty() || filter_column == "." {
        anyhow::bail!(
            "gnomAD FILTER column is {filter_column:?}: a record with no QC verdict cannot be \
             distinguished from one that passed"
        );
    }
    for name in filter_column.split(';') {
        if name != FILTER_PASS && !FILTER_FLAGS.iter().any(|(f, _, _)| *f == name) {
            let known: Vec<&str> = std::iter::once(FILTER_PASS)
                .chain(FILTER_FLAGS.iter().map(|(f, _, _)| *f))
                .collect();
            anyhow::bail!(
                "Unrecognized gnomAD FILTER value {name:?} in {filter_column:?}: this encoder \
                 models {known:?}"
            );
        }
    }
    Ok(())
}

/// INFO field names for a particular gnomAD release flavor.
///
/// Built from the VCF header so we use whatever names the upstream file
/// actually exposes. See module docs for the v4.1-joint vs. v4.0 split.
#[derive(Debug, Clone)]
struct FieldNames {
    af: String,
    an: String,
    ac: String,
    nhomalt: String,
    /// Separator between a statistic's name and a population code: `_` in
    /// every gnomAD release, `-` in some locally-processed files. The release
    /// flavor is already carried by the four names above, so `AC_joint` and
    /// `_` give `AC_joint_afr`.
    pop_separator: String,
}

impl FieldNames {
    fn standard() -> Self {
        Self {
            af: "AF".into(),
            an: "AN".into(),
            ac: "AC".into(),
            nhomalt: "nhomalt".into(),
            pop_separator: "_".into(),
        }
    }

    fn joint() -> Self {
        Self {
            af: "AF_joint".into(),
            an: "AN_joint".into(),
            ac: "AC_joint".into(),
            nhomalt: "nhomalt_joint".into(),
            pop_separator: "_".into(),
        }
    }

    fn af_pop_key(&self, pop: &str) -> String {
        format!("{}{}{}", self.af, self.pop_separator, pop)
    }

    fn an_pop_key(&self, pop: &str) -> String {
        format!("{}{}{}", self.an, self.pop_separator, pop)
    }

    fn ac_pop_key(&self, pop: &str) -> String {
        format!("{}{}{}", self.ac, self.pop_separator, pop)
    }

    fn nhomalt_pop_key(&self, pop: &str) -> String {
        format!("{}{}{}", self.nhomalt, self.pop_separator, pop)
    }
}

/// The per-population count statistics, in the order both encoders emit them
/// after that population's AF.
///
/// A per-population AF says how common the variant is in that group; it does
/// not say how many alleles the estimate rests on. An AF of 1e-4 over 200
/// alleles and the same AF over 200,000 are different observations, and only
/// AC/AN separate them. `Hc` is the homozygote count, which is the observation
/// ACMG BS2 asks for.
#[derive(Debug, Clone, Copy, PartialEq)]
enum PopCount {
    Ac,
    An,
    Hc,
}

const POP_COUNTS: &[PopCount] = &[PopCount::Ac, PopCount::An, PopCount::Hc];

impl PopCount {
    /// The suffix appended to the population code to form the output alias,
    /// as in `nfeAc`.
    fn alias_suffix(self) -> &'static str {
        match self {
            PopCount::Ac => "Ac",
            PopCount::An => "An",
            PopCount::Hc => "Hc",
        }
    }

    /// The INFO key holding this statistic for `pop` in this release flavor.
    fn info_key(self, names: &FieldNames, pop: &str) -> String {
        match self {
            PopCount::Ac => names.ac_pop_key(pop),
            PopCount::An => names.an_pop_key(pop),
            PopCount::Hc => names.nhomalt_pop_key(pop),
        }
    }

    /// Whether gnomAD reports the statistic per alternate allele
    /// (`Number=A`). AN is a property of the population's sample, not of an
    /// allele, so it is reported once per site (`Number=1`).
    fn per_allele(self) -> bool {
        self != PopCount::An
    }

    /// The INFO key as it appears in the standard release, for the v2 field
    /// descriptor's `field` metadata.
    fn canonical_field(self, pop: &str) -> String {
        self.info_key(&FieldNames::standard(), pop)
    }

    fn description(self, pop: &str) -> String {
        let upper = pop.to_uppercase();
        match self {
            PopCount::Ac => format!("{upper} alternate allele count"),
            PopCount::An => format!("{upper} total allele number"),
            PopCount::Hc => format!("{upper} homozygote count"),
        }
    }
}

/// One population's value for a count statistic, sliced out of the INFO map.
fn pop_count_value(
    info_map: &HashMap<String, String>,
    field_names: &FieldNames,
    stat: PopCount,
    pop: &str,
    allele_idx: usize,
) -> Option<i64> {
    let vals = split_info_values(
        info_map
            .get(&stat.info_key(field_names, pop))
            .map(|s| s.as_str()),
    );
    let raw = if stat.per_allele() {
        vals.get(allele_idx)
    } else {
        vals.first()
    };
    raw.and_then(|s| s.parse::<i64>().ok())
}

/// Pick a field-naming scheme based on which INFO IDs the VCF header
/// declares. Prefer standard names when present; fall back to the joint
/// names if the standard `AF` is absent but `AF_joint` is declared.
///
/// Some locally-processed files may use a hyphen separator for per-population
/// fields (`AF-afr`) rather than the standard underscore (`AF_afr`). We
/// detect this by checking whether any `AF-<pop>` field appears in the header.
fn detect_field_names(info_ids: &HashSet<String>) -> FieldNames {
    if info_ids.contains("AF") {
        let mut names = FieldNames::standard();
        if POPULATIONS
            .iter()
            .any(|p| info_ids.contains(&format!("AF-{p}")))
        {
            names.pop_separator = "-".into();
        }
        names
    } else if info_ids.contains("AF_joint") {
        FieldNames::joint()
    } else {
        FieldNames::standard()
    }
}

/// Extract the `ID=` value from a `##INFO=<ID=...,Number=...,...>` header
/// line. Returns `None` for non-INFO lines or malformed entries.
///
/// Handles the two real-world quirks of VCF headers:
/// - The trailing `>` when `ID` is the last attribute
///   (`##INFO=<...,ID=AF>` must yield `"AF"`, not `"AF>"`).
/// - Commas inside quoted `Description="foo, bar"` values, which a naive
///   `split(',')` would split on. We walk the body tracking quote state
///   so attributes can appear in any order.
fn parse_info_id(line: &str) -> Option<&str> {
    let body = line.strip_prefix("##INFO=<")?.strip_suffix('>')?;
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        // Skip leading whitespace between attributes.
        while i < bytes.len() && bytes[i] == b' ' {
            i += 1;
        }
        let key_start = i;
        while i < bytes.len() && bytes[i] != b'=' && bytes[i] != b',' {
            i += 1;
        }
        let key = &body[key_start..i];
        // Bare attribute with no value: skip the separator and continue.
        if i >= bytes.len() || bytes[i] == b',' {
            if i < bytes.len() {
                i += 1;
            }
            continue;
        }
        // Consume '=' and read the value, respecting double-quoted strings.
        i += 1;
        let value_start = i;
        let mut in_quotes = false;
        while i < bytes.len() {
            match bytes[i] {
                b'"' => in_quotes = !in_quotes,
                b',' if !in_quotes => break,
                _ => {}
            }
            i += 1;
        }
        if key == "ID" {
            let mut value = &body[value_start..i];
            if value.len() >= 2 && value.starts_with('"') && value.ends_with('"') {
                value = &value[1..value.len() - 1];
            }
            return Some(value);
        }
        if i < bytes.len() {
            i += 1; // skip the ',' separator
        }
    }
    None
}

/// Parse a gnomAD sites-only VCF and produce sorted `AnnotationRecord`s.
///
/// Collects all records into memory and sorts. For genome-scale VCFs use
/// `iter_gnomad_vcf` instead, which streams one block at a time.
pub fn parse_gnomad_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_gnomad_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });
    Ok(records)
}

/// Stream a coordinate-sorted gnomAD VCF as `AnnotationRecord`s without
/// buffering the whole file in memory.
///
/// The input must already be sorted by chromosome and position (all standard
/// gnomAD releases are). The writer will detect and error on out-of-order
/// records.
pub fn iter_gnomad_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> GnomadRecordIter<'a, R> {
    GnomadRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        pending: VecDeque::new(),
        info_ids: HashSet::new(),
        field_names: None,
    }
}

pub struct GnomadRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<AnnotationRecord>,
    info_ids: HashSet<String>,
    field_names: Option<FieldNames>,
}

impl<R: BufRead> Iterator for GnomadRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.pending.pop_front() {
                return Some(Ok(record));
            }

            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading gnomAD VCF line")),
            };

            if line.starts_with('#') {
                if let Some(id) = parse_info_id(&line) {
                    self.info_ids.insert(id.to_string());
                }
                continue;
            }

            let field_names = self
                .field_names
                .get_or_insert_with(|| detect_field_names(&self.info_ids));

            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 8 {
                continue;
            }

            let chrom = normalize_chrom(fields[0]);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match fields[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };

            let ref_allele = fields[3].to_string();
            let alt_field = fields[4];
            let filter_column = fields[6];
            if let Err(e) = validate_filter_column(filter_column) {
                return Some(Err(e));
            }
            let info = fields[7];

            let info_map = parse_info(info);
            let alts: Vec<&str> = alt_field.split(',').collect();
            let all_afs = split_info_values(info_map.get(&field_names.af).map(|s| s.as_str()));
            let all_ans = split_info_values(info_map.get(&field_names.an).map(|s| s.as_str()));
            let all_acs = split_info_values(info_map.get(&field_names.ac).map(|s| s.as_str()));
            let all_nhomalt =
                split_info_values(info_map.get(&field_names.nhomalt).map(|s| s.as_str()));

            for (i, alt) in alts.iter().enumerate() {
                if *alt == "." || *alt == "*" {
                    continue;
                }
                let json = build_gnomad_json(
                    GlobalCounts {
                        af: all_afs.get(i).map(|s| s.as_str()),
                        an: all_ans.first().map(|s| s.as_str()),
                        ac: all_acs.get(i).map(|s| s.as_str()),
                        nhomalt: all_nhomalt.get(i).map(|s| s.as_str()),
                    },
                    &info_map,
                    filter_column,
                    i,
                    field_names,
                );
                self.pending.push_back(AnnotationRecord {
                    chrom_idx,
                    position: pos,
                    ref_allele: ref_allele.clone(),
                    alt_allele: alt.to_string(),
                    json,
                });
            }
        }
    }
}

/// The global AF/AN/AC/nhomalt for one alternate allele, already sliced out of
/// the INFO column by the caller. Bundled so they are split once per VCF line
/// rather than once per alternate allele.
struct GlobalCounts<'a> {
    af: Option<&'a str>,
    an: Option<&'a str>,
    ac: Option<&'a str>,
    nhomalt: Option<&'a str>,
}

fn build_gnomad_json(
    counts: GlobalCounts<'_>,
    info_map: &HashMap<String, String>,
    filter_column: &str,
    allele_idx: usize,
    field_names: &FieldNames,
) -> String {
    let GlobalCounts {
        af,
        an,
        ac,
        nhomalt,
    } = counts;
    let mut parts = Vec::new();

    if let Some(af_str) = af {
        if let Ok(f) = af_str.parse::<f64>() {
            parts.push(format!("\"allAf\":{:.6e}", f));
        }
    }

    // AN/AC/nhomalt are written unquoted, so each must be a validated
    // non-negative integer — raw INFO-field text (e.g. the "." missing-value
    // sentinel or other garbage) would otherwise land in the JSON as a bare,
    // unquoted token and break every downstream serde_json::from_str on this
    // record. Mirrors the CNT validation in sources/cosmic.rs.
    if let Some(an_str) = an {
        if let Ok(n) = an_str.parse::<u64>() {
            parts.push(format!("\"allAn\":{}", n));
        }
    }

    if let Some(ac_str) = ac {
        if let Ok(n) = ac_str.parse::<u64>() {
            parts.push(format!("\"allAc\":{}", n));
        }
    }

    if let Some(nh) = nhomalt {
        if let Ok(n) = nh.parse::<u64>() {
            parts.push(format!("\"allHc\":{}", n));
        }
    }

    // Per-population AF, then the counts the AF was computed from.
    for pop in POPULATIONS {
        let key = field_names.af_pop_key(pop);
        if let Some(val) = info_map.get(&key) {
            let vals = split_info_values(Some(val.as_str()));
            if let Some(af_str) = vals.get(allele_idx) {
                if let Ok(f) = af_str.parse::<f64>() {
                    parts.push(format!("\"{}Af\":{:.6e}", pop, f));
                }
            }
        }
        for stat in POP_COUNTS {
            if let Some(n) = pop_count_value(info_map, field_names, *stat, pop, allele_idx) {
                parts.push(format!("\"{}{}\":{}", pop, stat.alias_suffix(), n));
            }
        }
    }

    // Extended QC / stratified columns. Absent values and unset flags are
    // omitted entirely, so a site with nothing to report costs no bytes and an
    // older consumer that does not know these keys is unaffected.
    for (alias, value) in extended_values(info_map, filter_column, allele_idx) {
        match value {
            ExtValue::Int(Some(n)) => parts.push(format!("\"{}\":{}", alias, n)),
            ExtValue::Float(Some(f)) if f.is_finite() => {
                parts.push(format!("\"{}\":{:.6e}", alias, f))
            }
            ExtValue::Flag(true) => parts.push(format!("\"{}\":true", alias)),
            ExtValue::Category(Some(i)) => {
                parts.push(format!("\"{}\":\"{}\"", alias, POPULATIONS[i as usize]))
            }
            _ => {}
        }
    }

    format!("{{{}}}", parts.join(","))
}

/// Parse a VCF INFO column into a key -> value map.
///
/// Bare flags (`Number=0` entries such as `segdup`, which appear with no `=`)
/// are stored with an empty value, so `contains_key` answers "is this flag
/// set?" while `get` on a valued key is unaffected.
fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in info.split(';') {
        match pair.split_once('=') {
            Some((key, value)) => {
                map.insert(key.to_string(), value.to_string());
            }
            None if !pair.is_empty() && pair != "." => {
                map.insert(pair.to_string(), String::new());
            }
            None => {}
        }
    }
    map
}

fn split_info_values(value: Option<&str>) -> Vec<String> {
    match value {
        Some(v) => v.split(',').map(|s| s.to_string()).collect(),
        None => Vec::new(),
    }
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

// =============================================================================
// v2 (.osa2) encoding
// =============================================================================

/// Multiplier for allele-frequency floats stored as u32 (`stored = (af * M) as u32`).
/// This gives a fixed *absolute* resolution of `1 / M` = 5e-7 (not a fixed
/// number of significant digits — precision in relative terms degrades as AF
/// shrinks), while staying well under `u32::MAX` for any AF in 0..1. The v1
/// builder formats the raw VCF float, so a v1 and a v2 database agree for every
/// AF coarser than 5e-7 — which is all of gnomAD except the very rarest v4
/// singletons (AF ~6e-7 at ~1.6M alleles), where v2 truncates down to the
/// nearest 5e-7 step and v1 does not. AC/AN/nhomalt are integers and always exact.
const AF_MULTIPLIER: u32 = 2_000_000;

fn af_field(field: &str, alias: &str, description: &str) -> Field {
    Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Float,
        multiplier: AF_MULTIPLIER,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

fn count_field(field: &str, alias: &str, description: &str) -> Field {
    Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Integer,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

/// A column holding one of a fixed set of strings, stored as an index into
/// that set. The vocabulary is written into the archive beside the values, so
/// the reader renders the string rather than the index.
fn categorical_field(field: &str, alias: &str, description: &str) -> Field {
    Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Categorical,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

/// A boolean column. `missing_value` is deliberately *not* `u32::MAX`: a flag
/// has no missing state, only set (1) or unset (0), and leaving the sentinel at
/// `u32::MAX` would make every unset flag decode as `null` rather than `false`.
fn flag_field(field: &str, alias: &str, description: &str) -> Field {
    Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Flag,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

/// The extended QC / stratified columns, in the order [`extended_values`]
/// produces them.
fn extended_fields() -> Vec<Field> {
    let mut fields = Vec::new();
    for (key, alias, desc) in XY_ALLELE_INTS {
        fields.push(count_field(key, alias, desc));
    }
    for (key, alias, desc) in XY_SITE_INTS {
        fields.push(count_field(key, alias, desc));
    }
    for (key, alias, desc) in FAF_FLOATS {
        fields.push(af_field(key, alias, desc));
    }
    for (key, alias, desc) in FAF_CATEGORICALS {
        fields.push(categorical_field(key, alias, desc));
    }
    for (key, alias, desc) in INFO_FLAGS {
        fields.push(flag_field(key, alias, desc));
    }
    for (key, alias, desc) in FILTER_FLAGS {
        fields.push(flag_field(key, alias, desc));
    }
    fields
}

/// Encode one extended column's value into its u32 slot.
fn encode_ext(field: &Field, value: ExtValue) -> u32 {
    match value {
        ExtValue::Int(Some(n)) => field.encode_int(n),
        ExtValue::Float(Some(f)) => field.encode_float(f),
        ExtValue::Category(Some(i)) => i,
        ExtValue::Int(None) | ExtValue::Float(None) | ExtValue::Category(None) => {
            field.missing_value
        }
        ExtValue::Flag(b) => u32::from(b),
    }
}

/// Number of columns each population contributes: its AF, then one per entry
/// of [`POP_COUNTS`].
const COLUMNS_PER_POPULATION: usize = 1 + POP_COUNTS.len();

/// Canonical gnomAD v2 (`.osa2`) field schema. The value vector produced by
/// [`iter_gnomad_osa2`] is parallel to this list, in this exact order: global
/// AF / AN / AC / nhomalt, then AF / AC / AN / nhomalt for each entry in
/// [`POPULATIONS`], then the extended columns. The aliases match the JSON keys
/// the v1 builder emits (`allAf`, `allAn`, `allAc`, `allHc`, `<pop>Af`,
/// `<pop>Ac`, `<pop>An`, `<pop>Hc`). Counts are byte-identical to v1; AFs are
/// quantized to a fixed 5e-7 resolution (see [`AF_MULTIPLIER`]), so they match
/// v1 for every AF except the rarest v4 singletons below that floor.
pub fn gnomad_osa2_fields() -> Vec<Field> {
    let mut fields = vec![
        af_field("AF", "allAf", "Global allele frequency"),
        count_field("AN", "allAn", "Total allele number"),
        count_field("AC", "allAc", "Allele count"),
        count_field("nhomalt", "allHc", "Homozygote count"),
    ];
    for pop in POPULATIONS {
        fields.push(af_field(
            &format!("AF_{pop}"),
            &format!("{pop}Af"),
            &format!("{} allele frequency", pop.to_uppercase()),
        ));
        for stat in POP_COUNTS {
            fields.push(count_field(
                &stat.canonical_field(pop),
                &format!("{pop}{}", stat.alias_suffix()),
                &stat.description(pop),
            ));
        }
    }
    fields.extend(extended_fields());
    fields
}

/// The categorical vocabularies [`gnomad_osa2_fields`] declares, paired with
/// their field indices, ready for the v2 writer's `set_string_table`.
pub fn gnomad_string_tables() -> Vec<(usize, Vec<String>)> {
    let table: Vec<String> = POPULATIONS.iter().map(|p| p.to_string()).collect();
    gnomad_osa2_fields()
        .iter()
        .enumerate()
        .filter(|(_, f)| f.ftype == FieldType::Categorical)
        .map(|(i, _)| (i, table.clone()))
        .collect()
}

/// Index of the first extended column in [`gnomad_osa2_fields`].
fn extended_offset() -> usize {
    4 + POPULATIONS.len() * COLUMNS_PER_POPULATION
}

/// Standard gnomAD `.osa2` metadata. Mirrors the v1 header (`json_key =
/// "gnomad"`, allele-matched, non-positional) so the annotation pipeline
/// treats a v1 and v2 gnomAD database identically.
pub fn gnomad_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "gnomAD".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("gnomAD population frequencies for {assembly}"),
    }
}

/// Stream a coordinate-sorted gnomAD VCF as `.osa2` records with parallel u32
/// value arrays (parallel to [`gnomad_osa2_fields`]), without buffering the
/// whole file. Reuses the same INFO field-name detection as the v1 path, so
/// v2/v3/v4/joint releases and hyphen-separated local files are all handled.
pub fn iter_gnomad_osa2<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> GnomadOsa2Iter<'a, R> {
    GnomadOsa2Iter {
        lines: reader.lines(),
        chrom_to_idx,
        fields: gnomad_osa2_fields(),
        pending: VecDeque::new(),
        info_ids: HashSet::new(),
        field_names: None,
    }
}

pub struct GnomadOsa2Iter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    fields: Vec<Field>,
    pending: VecDeque<Osa2Record>,
    info_ids: HashSet<String>,
    field_names: Option<FieldNames>,
}

impl<R: BufRead> Iterator for GnomadOsa2Iter<'_, R> {
    type Item = Result<Osa2Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.pending.pop_front() {
                return Some(Ok(record));
            }

            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading gnomAD VCF line")),
            };

            if line.starts_with('#') {
                if let Some(id) = parse_info_id(&line) {
                    self.info_ids.insert(id.to_string());
                }
                continue;
            }

            let field_names = self
                .field_names
                .get_or_insert_with(|| detect_field_names(&self.info_ids));

            let cols: Vec<&str> = line.splitn(9, '\t').collect();
            if cols.len() < 8 {
                continue;
            }

            let chrom = normalize_chrom(cols[0]);
            if !self.chrom_to_idx.contains_key(&chrom) {
                continue;
            }
            let pos: u32 = match cols[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };
            let ref_allele = cols[3].as_bytes().to_vec();
            let alt_field = cols[4];
            let filter_column = cols[6];
            if let Err(e) = validate_filter_column(filter_column) {
                return Some(Err(e));
            }
            let info_map = parse_info(cols[7]);

            let all_afs = split_info_values(info_map.get(&field_names.af).map(|s| s.as_str()));
            // AN is a per-site (Number=1) value; AF/AC/nhomalt are per-allele.
            let an = all_from(&info_map, &field_names.an);
            let all_acs = split_info_values(info_map.get(&field_names.ac).map(|s| s.as_str()));
            let all_nh = split_info_values(info_map.get(&field_names.nhomalt).map(|s| s.as_str()));

            for (ai, alt) in alt_field.split(',').enumerate() {
                if alt == "." || alt == "*" {
                    continue;
                }
                // Value order MUST match `gnomad_osa2_fields()`.
                let mut values = Vec::with_capacity(self.fields.len());
                values.push(enc_float(
                    &self.fields[0],
                    all_afs.get(ai).map(|s| s.as_str()),
                ));
                values.push(enc_int(&self.fields[1], an.first().map(|s| s.as_str())));
                values.push(enc_int(
                    &self.fields[2],
                    all_acs.get(ai).map(|s| s.as_str()),
                ));
                values.push(enc_int(&self.fields[3], all_nh.get(ai).map(|s| s.as_str())));
                for (pi, pop) in POPULATIONS.iter().enumerate() {
                    let base = 4 + pi * COLUMNS_PER_POPULATION;
                    let key = field_names.af_pop_key(pop);
                    let vals = split_info_values(info_map.get(&key).map(|s| s.as_str()));
                    values.push(enc_float(
                        &self.fields[base],
                        vals.get(ai).map(|s| s.as_str()),
                    ));
                    for (ci, stat) in POP_COUNTS.iter().enumerate() {
                        let field = &self.fields[base + 1 + ci];
                        values.push(
                            match pop_count_value(&info_map, field_names, *stat, pop, ai) {
                                Some(n) => field.encode_int(n),
                                None => field.missing_value,
                            },
                        );
                    }
                }
                let ext_off = extended_offset();
                for (ei, (_, value)) in extended_values(&info_map, filter_column, ai)
                    .into_iter()
                    .enumerate()
                {
                    values.push(encode_ext(&self.fields[ext_off + ei], value));
                }

                self.pending.push_back(Osa2Record {
                    chrom: chrom.clone(),
                    position: pos,
                    ref_allele: ref_allele.clone(),
                    alt_allele: alt.as_bytes().to_vec(),
                    values,
                    json_blob: None,
                });
            }
        }
    }
}

fn all_from(info: &HashMap<String, String>, key: &str) -> Vec<String> {
    split_info_values(info.get(key).map(|s| s.as_str()))
}

fn enc_float(field: &Field, raw: Option<&str>) -> u32 {
    match raw.and_then(|s| s.parse::<f64>().ok()) {
        Some(f) => field.encode_float(f),
        None => field.missing_value,
    }
}

fn enc_int(field: &Field, raw: Option<&str>) -> u32 {
    match raw.and_then(|s| s.parse::<i64>().ok()) {
        Some(i) => field.encode_int(i),
        None => field.missing_value,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gnomad_vcf() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele number\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_afr=0.002;AF_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF=0.01,0.005;AN=140000;AC=1400,700;nhomalt=10,3;AF_eas=0.02,0.01
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3); // 1 SNV + 2 from multi-allelic

        // First record
        assert_eq!(records[0].position, 10001);
        assert!(records[0].json.contains("\"allAf\":"));
        assert!(records[0].json.contains("\"afrAf\":"));
        assert!(records[0].json.contains("\"nfeAf\":"));

        // Multi-allelic: second alt
        assert_eq!(records[2].position, 20000);
        assert_eq!(records[2].alt_allele, "A");
    }

    #[test]
    fn test_parse_gnomad_v41_joint_vcf() {
        // Regression for issue #39: v4.1 joint release uses *_joint suffixes
        // and FV_GNOMAD came out empty because the parser only looked for
        // the bare AF/AC/AN names.
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF_joint,Number=A,Type=Float,Description=\"Joint allele frequency\">
##INFO=<ID=AN_joint,Number=1,Type=Integer,Description=\"Joint total allele number\">
##INFO=<ID=AC_joint,Number=A,Type=Integer,Description=\"Joint allele count\">
##INFO=<ID=nhomalt_joint,Number=A,Type=Integer,Description=\"Joint homozygote count\">
##INFO=<ID=AF_joint_afr,Number=A,Type=Float,Description=\"Joint AF AFR\">
##INFO=<ID=AF_joint_nfe,Number=A,Type=Float,Description=\"Joint AF NFE\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF_joint=0.001;AN_joint=150000;AC_joint=150;nhomalt_joint=2;AF_joint_afr=0.002;AF_joint_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF_joint=0.01,0.005;AN_joint=140000;AC_joint=1400,700;nhomalt_joint=10,3;AF_joint_eas=0.02,0.01
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3);

        // First record: all frequency fields populated, not empty.
        assert_eq!(records[0].position, 10001);
        assert!(
            records[0].json.contains("\"allAf\":"),
            "missing allAf in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"allAn\":150000"),
            "missing allAn in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"allAc\":150"),
            "missing allAc in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"allHc\":2"),
            "missing allHc in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"afrAf\":"),
            "missing afrAf in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"nfeAf\":"),
            "missing nfeAf in: {}",
            records[0].json
        );

        // Multi-allelic: second alt also gets per-allele values.
        assert_eq!(records[2].position, 20000);
        assert_eq!(records[2].alt_allele, "A");
        assert!(
            records[2].json.contains("\"allAc\":700"),
            "second alt should pick second AC value: {}",
            records[2].json
        );
    }

    #[test]
    fn test_garbage_an_ac_nhomalt_omitted_not_emitted_unescaped() {
        // AN/AC/nhomalt are spliced into the JSON unquoted, so malformed
        // upstream text (e.g. the "." missing-value sentinel or other
        // garbage) must be dropped rather than emitted as a bare token that
        // would break serde_json::from_str on the whole record.
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele number\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=not_a_number;AC=garbage;nhomalt=.
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);
        // Malformed fields must be entirely absent, not emitted as garbage.
        assert!(!records[0].json.contains("allAn"));
        assert!(!records[0].json.contains("allAc"));
        assert!(!records[0].json.contains("allHc"));
        // The emitted JSON must still be valid.
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert!(v.get("allAf").is_some());
    }

    #[test]
    fn test_detect_field_names_prefers_standard() {
        let mut ids = HashSet::new();
        ids.insert("AF".into());
        ids.insert("AF_joint".into());
        let names = detect_field_names(&ids);
        assert_eq!(names.af, "AF");
    }

    #[test]
    fn test_detect_field_names_hyphen_separator() {
        // Some locally-processed files use AF-afr (hyphen) instead of the standard AF_afr.
        let mut ids = HashSet::new();
        ids.insert("AF".into());
        ids.insert("AN".into());
        ids.insert("AF-afr".into());
        ids.insert("AF-nfe".into());
        let names = detect_field_names(&ids);
        assert_eq!(names.af, "AF");
        assert_eq!(names.af_pop_key("afr"), "AF-afr");
        assert_eq!(names.af_pop_key("nfe"), "AF-nfe");
        // The separator applies to every per-population statistic, not just AF.
        assert_eq!(names.ac_pop_key("nfe"), "AC-nfe");
        assert_eq!(names.an_pop_key("nfe"), "AN-nfe");
        assert_eq!(names.nhomalt_pop_key("nfe"), "nhomalt-nfe");
    }

    #[test]
    fn test_detect_field_names_falls_back_to_joint() {
        let mut ids = HashSet::new();
        ids.insert("AF_joint".into());
        ids.insert("AC_joint".into());
        let names = detect_field_names(&ids);
        assert_eq!(names.af, "AF_joint");
        assert_eq!(names.af_pop_key("nfe"), "AF_joint_nfe");
        assert_eq!(names.ac_pop_key("nfe"), "AC_joint_nfe");
        assert_eq!(names.an_pop_key("nfe"), "AN_joint_nfe");
        assert_eq!(names.nhomalt_pop_key("nfe"), "nhomalt_joint_nfe");
    }

    #[test]
    fn test_parse_info_id_standard() {
        let line = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">";
        assert_eq!(parse_info_id(line), Some("AF"));
    }

    #[test]
    fn test_parse_info_id_trailing_angle_bracket() {
        // ID is the last attribute — the closing '>' must not be captured.
        let line = "##INFO=<Number=A,Type=Float,ID=AF>";
        assert_eq!(parse_info_id(line), Some("AF"));
    }

    #[test]
    fn test_parse_info_id_reordered_with_quoted_comma() {
        // Description quoted string contains commas — must not split inside it.
        let line = "##INFO=<Number=A,Type=Float,Description=\"AF, joint, multi-pop\",ID=AF_joint>";
        assert_eq!(parse_info_id(line), Some("AF_joint"));
    }

    #[test]
    fn test_parse_info_id_quoted_id_value() {
        let line = "##INFO=<ID=\"weird_id\",Number=A,Type=Float,Description=\"x\">";
        assert_eq!(parse_info_id(line), Some("weird_id"));
    }

    #[test]
    fn test_parse_info_id_non_info_line() {
        assert_eq!(parse_info_id("##fileformat=VCFv4.2"), None);
        assert_eq!(parse_info_id("#CHROM\tPOS\tID"), None);
    }

    #[test]
    fn test_parse_info_id_malformed_no_closing_bracket() {
        // Missing trailing '>' — refuse to parse rather than guess.
        let line = "##INFO=<ID=AF,Number=A";
        assert_eq!(parse_info_id(line), None);
    }

    // ---- v2 (.osa2) encoder ----

    fn chr1_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".to_string(), 0u16);
        m
    }

    #[test]
    fn test_osa2_fields_order_matches_value_layout() {
        // The iterator indexes `fields[0..4]` for the global stats and
        // `fields[4 + pi * COLUMNS_PER_POPULATION + ..]` per population, so
        // this ordering is load-bearing.
        let fields = gnomad_osa2_fields();
        assert_eq!(
            fields.len(),
            4 + POPULATIONS.len() * COLUMNS_PER_POPULATION + extended_fields().len()
        );
        assert_eq!(fields[0].alias, "allAf");
        assert_eq!(fields[1].alias, "allAn");
        assert_eq!(fields[2].alias, "allAc");
        assert_eq!(fields[3].alias, "allHc");
        for (pi, pop) in POPULATIONS.iter().enumerate() {
            let base = 4 + pi * COLUMNS_PER_POPULATION;
            assert_eq!(fields[base].alias, format!("{pop}Af"));
            for (ci, stat) in POP_COUNTS.iter().enumerate() {
                assert_eq!(
                    fields[base + 1 + ci].alias,
                    format!("{pop}{}", stat.alias_suffix())
                );
            }
        }
    }

    #[test]
    fn test_iter_gnomad_osa2_encodes_values() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele number\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_afr=0.002;AF_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF=0.01,0.005;AN=140000;AC=1400,700;nhomalt=10,3;AF_eas=0.02,0.01
";
        let map = chr1_map();
        let recs: Vec<Osa2Record> = iter_gnomad_osa2(vcf.as_bytes(), &map)
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(recs.len(), 3); // 1 SNV + 2 from the multi-allelic site

        let fields = gnomad_osa2_fields();
        // First record: AF=0.001 -> 0.001 * 2_000_000 = 2000.
        assert_eq!(recs[0].position, 10001);
        assert_eq!(recs[0].values[0], 2000);
        assert_eq!(recs[0].values[1], 150000); // AN
        assert_eq!(recs[0].values[2], 150); // AC
        assert_eq!(recs[0].values[3], 2); // nhomalt
                                          // afr AF present, sas AF missing.
        let pop_af_idx = |pop: &str| {
            4 + POPULATIONS.iter().position(|p| *p == pop).unwrap() * COLUMNS_PER_POPULATION
        };
        let afr_idx = pop_af_idx("afr");
        let sas_idx = pop_af_idx("sas");
        assert_eq!(recs[0].values[afr_idx], fields[afr_idx].encode_float(0.002));
        assert_eq!(recs[0].values[sas_idx], u32::MAX); // missing sentinel

        // Multi-allelic second alt picks the second per-allele values.
        assert_eq!(recs[2].position, 20000);
        assert_eq!(recs[2].alt_allele, b"A");
        assert_eq!(recs[2].values[2], 700); // AC second alt
    }

    #[test]
    fn test_iter_gnomad_osa2_handles_joint_release() {
        // v4.1 joint naming must be detected here too (regression parity with
        // the v1 parser).
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF_joint,Number=A,Type=Float,Description=\"Joint AF\">
##INFO=<ID=AN_joint,Number=1,Type=Integer,Description=\"Joint AN\">
##INFO=<ID=AC_joint,Number=A,Type=Integer,Description=\"Joint AC\">
##INFO=<ID=AF_joint_nfe,Number=A,Type=Float,Description=\"Joint AF NFE\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF_joint=0.001;AN_joint=150000;AC_joint=150;AF_joint_nfe=0.0005
";
        let map = chr1_map();
        let recs: Vec<Osa2Record> = iter_gnomad_osa2(vcf.as_bytes(), &map)
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].values[0], 2000); // AF_joint 0.001 encoded
        assert_eq!(recs[0].values[1], 150000); // AN_joint
        let nfe_idx =
            4 + POPULATIONS.iter().position(|p| *p == "nfe").unwrap() * COLUMNS_PER_POPULATION;
        assert_ne!(recs[0].values[nfe_idx], u32::MAX); // nfe populated from joint
    }

    // ---- extended QC / stratified columns ----

    /// A chrX non-PAR record with the full set of extended fields populated:
    /// failed VQSR, inside a segdup, hemizygous counts, and both FAFs.
    const CHRX_VCF: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"AN\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t500\t.\tA\tG\t.\tAS_VQSR;InbreedingCoeff\tAF=0.2;AN=1000;AC=200;nhomalt=20;segdup;lcr
chr1\t600\t.\tA\tG\t.\tPASS\tAF=0.2;AN=1000;AC=200;nhomalt=20;non_par;AC_XY=37;AN_XY=400;nhomalt_XY=5;faf95=1.7e-01;fafmax_faf95_max=2.4e-01
";

    #[test]
    fn test_extended_fields_match_extended_values() {
        // The v1 JSON builder and the v2 u32 encoder both walk these two lists
        // positionally. If they ever fall out of order, every extended column
        // in a v2 database silently holds another column's value.
        let fields = extended_fields();
        let values = extended_values(&HashMap::new(), "PASS", 0);
        assert_eq!(fields.len(), values.len());
        for (f, (alias, _)) in fields.iter().zip(values.iter()) {
            assert_eq!(&f.alias, alias, "extended schema and extraction diverged");
        }
        // And the offset the v2 encoder uses must land on the first of them.
        let all = gnomad_osa2_fields();
        assert_eq!(all[extended_offset()].alias, fields[0].alias);
        assert_eq!(all.len(), extended_offset() + fields.len());
    }

    #[test]
    fn test_v1_json_carries_filter_and_region_flags() {
        let map = chr1_map();
        let records = parse_gnomad_vcf(CHRX_VCF.as_bytes(), &map).unwrap();
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(v.get("filterVqsr").and_then(|x| x.as_bool()), Some(true));
        assert_eq!(
            v.get("filterInbreeding").and_then(|x| x.as_bool()),
            Some(true)
        );
        assert_eq!(v.get("segdup").and_then(|x| x.as_bool()), Some(true));
        assert_eq!(v.get("lcr").and_then(|x| x.as_bool()), Some(true));
        // Unset flags are omitted rather than written as false.
        assert!(v.get("filterAc0").is_none());
        assert!(v.get("nonPar").is_none());
    }

    #[test]
    fn test_v1_json_carries_xy_counts_and_faf() {
        let map = chr1_map();
        let records = parse_gnomad_vcf(CHRX_VCF.as_bytes(), &map).unwrap();
        let v: serde_json::Value = serde_json::from_str(&records[1].json).unwrap();
        assert_eq!(v.get("allAcXY").and_then(|x| x.as_i64()), Some(37));
        assert_eq!(v.get("allAnXY").and_then(|x| x.as_i64()), Some(400));
        assert_eq!(v.get("nonPar").and_then(|x| x.as_bool()), Some(true));
        let faf_max = v.get("faf95Max").and_then(|x| x.as_f64()).unwrap();
        assert!((faf_max - 0.24).abs() < 1e-6, "faf95Max was {faf_max}");
    }

    #[test]
    fn test_osa2_encodes_extended_columns() {
        let map = chr1_map();
        let recs: Vec<Osa2Record> = iter_gnomad_osa2(CHRX_VCF.as_bytes(), &map)
            .collect::<Result<_>>()
            .unwrap();
        let fields = gnomad_osa2_fields();
        let idx = |alias: &str| fields.iter().position(|f| f.alias == alias).unwrap();

        // Flags are 0/1, never the u32::MAX missing sentinel, so an unset flag
        // decodes as `false` rather than `null`.
        assert_eq!(recs[0].values[idx("filterVqsr")], 1);
        assert_eq!(recs[0].values[idx("filterAc0")], 0);
        assert_eq!(recs[0].values[idx("segdup")], 1);

        // Absent integers keep the missing sentinel so they decode as null.
        assert_eq!(recs[0].values[idx("allAcXY")], u32::MAX);

        assert_eq!(recs[1].values[idx("allAcXY")], 37);
        assert_eq!(recs[1].values[idx("allAnXY")], 400);
        assert_eq!(recs[1].values[idx("nonPar")], 1);
        assert_eq!(
            recs[1].values[idx("faf95Max")],
            fields[idx("faf95Max")].encode_float(0.24)
        );
    }

    // ---- per-population counts and the grpmax FAF group ----

    /// A record carrying the per-population counts and both grpmax FAFs.
    const POP_COUNT_VCF: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"AN\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_nfe=0.0005;AN_nfe=60000;AC_nfe=30;nhomalt_nfe=1;fafmax_faf95_max=0.0004;fafmax_faf99_max=0.0003;fafmax_faf95_max_gen_anc=nfe
";

    #[test]
    fn test_v1_json_carries_per_population_counts() {
        // A per-population AF says how common the variant is in that group and
        // nothing about how many alleles the estimate rests on. Without the
        // counts the two are indistinguishable.
        let records = parse_gnomad_vcf(POP_COUNT_VCF.as_bytes(), &chr1_map()).unwrap();
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(v.get("nfeAc").and_then(|x| x.as_i64()), Some(30));
        assert_eq!(v.get("nfeAn").and_then(|x| x.as_i64()), Some(60000));
        assert_eq!(v.get("nfeHc").and_then(|x| x.as_i64()), Some(1));
        assert!(v.get("nfeAf").is_some());
        // A population the record says nothing about carries no keys at all.
        assert!(v.get("sasAc").is_none());
        assert!(v.get("sasAf").is_none());
    }

    #[test]
    fn test_v1_json_carries_both_grpmax_fafs_and_their_group() {
        let records = parse_gnomad_vcf(POP_COUNT_VCF.as_bytes(), &chr1_map()).unwrap();
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        let faf95 = v.get("faf95Max").and_then(|x| x.as_f64()).unwrap();
        let faf99 = v.get("faf99Max").and_then(|x| x.as_f64()).unwrap();
        assert!((faf95 - 0.0004).abs() < 1e-9, "faf95Max was {faf95}");
        assert!((faf99 - 0.0003).abs() < 1e-9, "faf99Max was {faf99}");
        assert_eq!(
            v.get("faf95MaxGenAnc").and_then(|x| x.as_str()),
            Some("nfe")
        );
    }

    #[test]
    fn test_osa2_carries_per_population_counts_and_gen_anc() {
        // The v1 JSON and the v2 columns are what the same consumer reads, so
        // a field present in one and absent from the other is a field that
        // disappears when the build format changes.
        let recs: Vec<Osa2Record> = iter_gnomad_osa2(POP_COUNT_VCF.as_bytes(), &chr1_map())
            .collect::<Result<_>>()
            .unwrap();
        let fields = gnomad_osa2_fields();
        let idx = |alias: &str| fields.iter().position(|f| f.alias == alias).unwrap();

        assert_eq!(recs[0].values[idx("nfeAc")], 30);
        assert_eq!(recs[0].values[idx("nfeAn")], 60000);
        assert_eq!(recs[0].values[idx("nfeHc")], 1);
        assert_eq!(recs[0].values[idx("sasAc")], u32::MAX);
        assert_eq!(
            recs[0].values[idx("faf99Max")],
            fields[idx("faf99Max")].encode_float(0.0003)
        );
        // The group is stored as an index into POPULATIONS, which the archive
        // carries as the field's string table.
        let nfe = POPULATIONS.iter().position(|p| *p == "nfe").unwrap() as u32;
        assert_eq!(recs[0].values[idx("faf95MaxGenAnc")], nfe);
    }

    #[test]
    fn test_gnomad_string_tables_cover_every_categorical_field() {
        // A categorical column whose table is never written decodes to the
        // stored index rather than the string.
        let fields = gnomad_osa2_fields();
        let tables = gnomad_string_tables();
        let categorical: Vec<usize> = fields
            .iter()
            .enumerate()
            .filter(|(_, f)| f.ftype == FieldType::Categorical)
            .map(|(i, _)| i)
            .collect();
        assert!(!categorical.is_empty());
        assert_eq!(
            tables.iter().map(|(i, _)| *i).collect::<Vec<_>>(),
            categorical
        );
        for (_, table) in &tables {
            assert_eq!(table.len(), POPULATIONS.len());
            assert_eq!(table[0], POPULATIONS[0]);
        }
    }

    #[test]
    fn test_gen_anc_outside_the_vocabulary_is_absent_not_misattributed() {
        // The column is an index into POPULATIONS, so a code that is not one
        // has no slot. Encoding it as some other group's index would attribute
        // the frequency to the wrong ancestry.
        let vcf = POP_COUNT_VCF.replace("gen_anc=nfe", "gen_anc=unheard_of");
        let records = parse_gnomad_vcf(vcf.as_bytes(), &chr1_map()).unwrap();
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert!(v.get("faf95MaxGenAnc").is_none());
        assert!(v.get("faf95Max").is_some(), "the frequency itself survives");
    }

    // ---- FILTER vocabulary ----

    #[test]
    fn test_validate_filter_column_accepts_the_declared_vocabulary() {
        for value in ["PASS", "AC0", "AS_VQSR", "InbreedingCoeff", "AC0;AS_VQSR"] {
            assert!(
                validate_filter_column(value).is_ok(),
                "{value} is declared by the gnomAD v4.1 header"
            );
        }
    }

    #[test]
    fn test_validate_filter_column_rejects_a_missing_verdict() {
        // Three unset flags read as PASS, and PASS is what permits benign
        // frequency evidence, so a record with no verdict must not be encoded
        // as one that passed.
        for value in [".", ""] {
            let err = validate_filter_column(value).unwrap_err().to_string();
            assert!(err.contains("no QC verdict"), "{value:?}: {err}");
        }
    }

    #[test]
    fn test_validate_filter_column_rejects_an_unmodelled_filter_name() {
        // The gnomAD v2 vocabulary. A release using it would encode every
        // record as PASS.
        for value in ["RF", "LCR", "SEGDUP", "PASS;RF"] {
            let err = validate_filter_column(value).unwrap_err().to_string();
            assert!(err.contains("Unrecognized"), "{value:?}: {err}");
        }
    }

    #[test]
    fn test_both_encoders_stop_on_an_unmodelled_filter_value() {
        // The guard belongs to the parse, so it holds whichever format the
        // build writes.
        let vcf = POP_COUNT_VCF.replace("\tPASS\t", "\tSEGDUP\t");
        let v1 = parse_gnomad_vcf(vcf.as_bytes(), &chr1_map());
        assert!(v1.is_err(), "v1 parse accepted SEGDUP");
        let v2: Result<Vec<Osa2Record>> = iter_gnomad_osa2(vcf.as_bytes(), &chr1_map()).collect();
        assert!(v2.is_err(), "v2 parse accepted SEGDUP");
    }

    #[test]
    fn test_filter_has_matches_whole_entries_only() {
        // A substring match would make "AS_VQSR" also match a hypothetical
        // "AS_VQSR_2", and "PASS" must never look like a failed filter.
        assert!(filter_has("AC0;AS_VQSR", "AC0"));
        assert!(filter_has("AC0;AS_VQSR", "AS_VQSR"));
        assert!(!filter_has("AC0;AS_VQSR", "InbreedingCoeff"));
        assert!(!filter_has("PASS", "AC0"));
        assert!(!filter_has(".", "AC0"));
        assert!(!filter_has("AS_VQSR_2", "AS_VQSR"));
    }

    #[test]
    fn test_parse_info_records_bare_flags_without_disturbing_valued_keys() {
        let info = parse_info("AF=0.1;segdup;AN=100;.;lcr");
        assert_eq!(info.get("AF").map(|s| s.as_str()), Some("0.1"));
        assert_eq!(info.get("AN").map(|s| s.as_str()), Some("100"));
        assert!(info.contains_key("segdup"));
        assert!(info.contains_key("lcr"));
        assert!(!info.contains_key("non_par"));
        // The "." placeholder is not a flag.
        assert!(!info.contains_key("."));
    }

    #[test]
    fn test_iter_gnomad_osa2_skips_unknown_contigs() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr2\t100\t.\tA\tG\t.\tPASS\tAF=0.5
chr1\t200\t.\tA\tG\t.\tPASS\tAF=0.5
";
        let map = chr1_map(); // only chr1 is valid
        let recs: Vec<Osa2Record> = iter_gnomad_osa2(vcf.as_bytes(), &map)
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].chrom, "chr1");
        assert_eq!(recs[0].position, 200);
    }
}
