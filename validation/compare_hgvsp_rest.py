#!/usr/bin/env python3
"""Compare fastVEP HGVSp against Ensembl VEP for every variant in a VCF.

Unlike ``compare_vep.py`` (which diffs two annotated VCFs and needs the Ensembl
VEP Docker image plus a full GRCh38 GFF3/FASTA), this queries the Ensembl REST
VEP endpoint directly. That makes it runnable with no local reference data, at
the cost of needing network access and honouring Ensembl's rate limits.

It exists to guard one specific failure class: in-frame indels whose protein
description is emitted in a substitution shape. Those render as ``p.Glu117???``
(deletion marker with no three-letter code), ``p.Pro189=`` (surviving residue
compares equal) or ``p.Glu251Gly`` (surviving residue differs) — the last being
the dangerous one, since a wrong missense is indistinguishable from a real call.

Usage:
    # Seed expectations for a new fixture: report what Ensembl returns.
    python3 validation/compare_hgvsp_rest.py --rest-only <input.vcf>

    # Compare a fastVEP run against Ensembl.
    fastvep annotate -i <input.vcf> -o out.json --output-format json --hgvs \\
        --gff3 <gff3> --fasta <fasta>
    python3 validation/compare_hgvsp_rest.py <input.vcf> --fastvep-json out.json

Exit status is non-zero when any variant is malformed (a substitution shape
where Ensembl reports an indel) or missing (Ensembl describes the protein change
and fastVEP does not).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request

REST_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
BATCH = 200  # Ensembl's documented maximum per POST
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
MAX_ATTEMPTS = 5
MIN_RETRY_DELAY_SECONDS = 1.0
MAX_RETRY_DELAY_SECONDS = 60.0

SUBSTITUTION_SHAPE = re.compile(r"^p\.[A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|=|\?{3})$")
INDEL_SHAPE = re.compile(r"(del|dup|ins|fs|ext)")


def variant_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    """Identify a variant by coordinate rather than by VCF ID.

    IDs are not guaranteed unique in a fixture, and Ensembl omits the id
    entirely for some records — keying on either silently merges results.
    """
    return f"{chrom}:{pos}:{ref}:{alt}"


def read_vcf(path: str) -> list[dict[str, str]]:
    variants: list[dict[str, str]] = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, identifier, ref, alt = fields[:5]
            variants.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "id": identifier,
                    "ref": ref,
                    "alt": alt,
                    "key": variant_key(chrom, pos, ref, alt),
                }
            )
    return variants


def retry_delay(response_headers: object, attempt: int) -> float:
    """Seconds to wait before the next attempt, honouring Retry-After."""
    raw = None
    if response_headers is not None:
        raw = getattr(response_headers, "get", lambda _name, _default=None: None)(
            "Retry-After"
        )
    try:
        delay = float(raw)
    except (TypeError, ValueError):
        delay = float(2**attempt)
    return min(max(delay, MIN_RETRY_DELAY_SECONDS), MAX_RETRY_DELAY_SECONDS)


def post(payload: dict) -> list:
    request = urllib.request.Request(
        REST_URL,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    last_error: Exception | None = None
    for attempt in range(MAX_ATTEMPTS):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            last_error = error
            if error.code not in RETRY_STATUSES or attempt == MAX_ATTEMPTS - 1:
                raise
            time.sleep(retry_delay(error.headers, attempt))
        except (urllib.error.URLError, TimeoutError) as error:
            # A transient network fault mid-corpus should not lose the whole run.
            last_error = error
            if attempt == MAX_ATTEMPTS - 1:
                raise
            time.sleep(retry_delay(None, attempt))
    raise RuntimeError(f"Ensembl REST retries exhausted: {last_error}")


def pick_consequence(transcript_consequences: list[dict]) -> dict:
    """Prefer the MANE Select or canonical transcript that carries an HGVSp."""
    preferred = next(
        (
            consequence
            for consequence in transcript_consequences
            if (consequence.get("mane_select") or consequence.get("canonical"))
            and consequence.get("hgvsp")
        ),
        None,
    )
    if preferred is not None:
        return preferred
    return next(
        (
            consequence
            for consequence in transcript_consequences
            if consequence.get("hgvsp")
        ),
        {},
    )


def bare_hgvsp(consequence: dict) -> str:
    """The protein change without its accession prefix; "" when absent or null."""
    return (consequence.get("hgvsp") or "").split(":")[-1]


def query_ensembl(variants: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Return {variant_key: {hgvsp, consequence, transcript}} from Ensembl VEP."""
    results: dict[str, dict[str, str]] = {}
    by_key = {variant["key"]: variant for variant in variants}
    for start in range(0, len(variants), BATCH):
        chunk = variants[start : start + BATCH]
        payload = {
            # Send the coordinate key as the VCF ID so responses can be matched
            # back unambiguously even when the fixture repeats an ID.
            "variants": [
                f"{v['chrom']} {v['pos']} {v['key']} {v['ref']} {v['alt']} . . ."
                for v in chunk
            ],
            "hgvs": 1,
            "canonical": 1,
            "mane": 1,
        }
        for record in post(payload):
            key = record.get("id")
            if key not in by_key:
                continue
            consequence = pick_consequence(record.get("transcript_consequences") or [])
            results[key] = {
                "hgvsp": bare_hgvsp(consequence),
                "consequence": ",".join(consequence.get("consequence_terms") or [])
                or (record.get("most_severe_consequence") or ""),
                "transcript": consequence.get("transcript_id") or "",
            }
        sys.stderr.write(
            f"  queried {min(start + BATCH, len(variants))}/{len(variants)}\n"
        )
    return results


def read_fastvep_hgvsp(path: str) -> dict[str, str]:
    """Pull HGVSp per variant coordinate from fastVEP JSON output."""
    with open(path) as handle:
        data = json.load(handle)
    records = data if isinstance(data, list) else data.get("annotations", data)
    results: dict[str, str] = {}
    for record in records:
        # fastVEP echoes the input allele on each record; fall back to the id
        # when a build predates that field.
        chrom = record.get("seq_region_name")
        start = record.get("start")
        allele = record.get("allele_string") or ""
        if chrom is not None and start is not None and "/" in allele:
            ref, _, alt = allele.partition("/")
            key = variant_key(str(chrom), str(start), ref, alt)
        else:
            key = record.get("id")
        consequence = pick_consequence(record.get("transcript_consequences") or [])
        results[key] = bare_hgvsp(consequence)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="input VCF, the same one fastVEP annotated")
    parser.add_argument(
        "--fastvep-json", help="fastVEP --output-format json result to compare"
    )
    parser.add_argument(
        "--rest-only",
        action="store_true",
        help="report Ensembl output only; do not compare",
    )
    args = parser.parse_args()

    if not args.rest_only and not args.fastvep_json:
        parser.error("--fastvep-json is required unless --rest-only is given")

    variants = read_vcf(args.vcf)
    sys.stderr.write(f"{len(variants)} variants from {args.vcf}\n")
    ensembl = query_ensembl(variants)
    fastvep = read_fastvep_hgvsp(args.fastvep_json) if args.fastvep_json else {}

    malformed: list[str] = []
    missing: list[str] = []
    differing: list[str] = []
    agreed = 0

    print(f"{'variant':<34} {'Ensembl HGVSp':<30} {'fastVEP HGVSp':<30} verdict")
    print("-" * 116)
    for variant in variants:
        key = variant["key"]
        expected = ensembl.get(key, {}).get("hgvsp", "")
        actual = fastvep.get(key, "")
        label = variant["id"] if variant["id"] not in ("", ".") else key

        if args.rest_only:
            terms = ensembl.get(key, {}).get("consequence", "")
            print(f"{label[:33]:<34} {expected:<30} {'-':<30} {terms}")
            continue

        if actual and SUBSTITUTION_SHAPE.match(actual) and INDEL_SHAPE.search(expected):
            verdict = "MALFORMED"
            malformed.append(f"{label}: {actual} (Ensembl {expected})")
        elif expected and not actual:
            # Ensembl describes the protein change and fastVEP emits nothing.
            # This is the regression the terminal-insertion fallback guards.
            verdict = "MISSING"
            missing.append(f"{label}: Ensembl {expected}, fastVEP none")
        elif actual and expected and actual != expected:
            verdict = "differs"
            differing.append(f"{label}: {actual} vs Ensembl {expected}")
        elif actual and actual == expected:
            verdict = "ok"
            agreed += 1
        else:
            verdict = "no-hgvsp"
        print(f"{label[:33]:<34} {expected:<30} {actual:<30} {verdict}")

    if args.rest_only:
        return 0

    print("-" * 116)
    print(f"agreed exactly : {agreed}")
    print(f"differs        : {len(differing)}")
    print(f"MISSING        : {len(missing)}")
    print(f"MALFORMED      : {len(malformed)}")
    for line in malformed + missing:
        print(f"  {line}")

    # `differs` alone does not fail the run: fastVEP falls back to an unshifted
    # description for transcripts with no usable peptide, which is valid HGVS at
    # a position Ensembl may report differently. Malformed and missing do fail —
    # those are the defect classes this script exists to catch.
    return 1 if malformed or missing else 0


if __name__ == "__main__":
    sys.exit(main())
