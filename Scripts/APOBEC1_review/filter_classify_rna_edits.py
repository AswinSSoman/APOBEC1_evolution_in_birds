#!/usr/bin/env python3
"""
Filter and classify candidate RNA editing sites from REDItools output.

Inputs:
  * REDItools table
  * RepeatMasker output, BED, or GFF-like repeat intervals
  * Genome FASTA, indexed with samtools faidx
  * RNA BAM, indexed
  * DNA BAM, indexed
  * GTF annotation

Main filters:
  * Repeat overlap
  * Homopolymer runs in the reference FASTA
  * Low RNA or DNA coverage and base quality
  * Possible genomic variants seen in the DNA counts
  * Poor mapping quality in RNA or DNA BAM pileups
  * Splice-site proximity from GTF exons
  * Indels in flanking regions from RNA or DNA BAM pileups
  * Absence from annotated gene models

Classification:
  * Uses the strand of overlapping gene models.
  * Reports the exact transcript-oriented substitution type, for example A>G, C>T, G>A.
  * Colors all possible ADAR-compatible changes (A>G / T>C) in blue,
    all APOBEC-compatible changes (C>T / G>A) in orange,
    and all remaining change types in grey.

Example:
  python filter_classify_rna_edits.py \
    --reditools sample.redi.tsv \
    --repeatmasker genome.fa.out \
    --genome genome.fa \
    --rna-bam rna.sorted.bam \
    --dna-bam dna.sorted.bam \
    --gtf annotation.gtf \
    --out-prefix results/sample

Dependencies:
  Python 3.6+
  pip install pysam matplotlib
"""

# Compatibility: tested to parse with Python 3.6 grammar; no future annotations import.
import argparse
import csv
import gzip
import logging
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Sequence, Tuple

try:
    import pysam
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency: pysam. Install it with: pip install pysam"
    ) from exc


BASES = ("A", "C", "G", "T")
BASE_TO_I = {b: i for i, b in enumerate(BASES)}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

TRANSCRIPT_SUBSTITUTION_ORDER = [
    "A>G",  # canonical ADAR signal
    "C>T",  # canonical APOBEC signal
    "A>C",
    "A>T",
    "C>A",
    "C>G",
    "G>A",
    "G>C",
    "G>T",
    "T>A",
    "T>C",
    "T>G",
]


REDITOOLS_COLUMNS = [
    "Region",
    "Position",
    "Reference",
    "Strand",
    "Coverage-q30",
    "MeanQ",
    "BaseCount[A,C,G,T]",
    "AllSubs",
    "Frequency",
    "gCoverage-q30",
    "gMeanQ",
    "gBaseCount[A,C,G,T]",
    "gAllSubs",
    "gFrequency",
]

NUMBER_TOKEN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+|nan|NaN|NA|-|\.)"
FREQ_TOKEN = r"(?:-|\.|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[,;][+-]?(?:\d+(?:\.\d*)?|\.\d+))*)"
COUNT_TOKEN = r"(?:-|\[[^\]]+\])"

REDITOOLS_SPACE_RE = re.compile(
    r"^(\S+)\s+"                  # Region
    r"(\d+)\s+"                  # Position
    r"([ACGTNacgtn])\s+"         # Reference
    r"(\S+)\s+"                  # Strand
    r"(\d+|-)\s+"               # Coverage-q30
    r"(" + NUMBER_TOKEN + r")\s+"  # MeanQ
    r"(" + COUNT_TOKEN + r")\s+"   # BaseCount
    r"(.+?)\s+"                  # AllSubs
    r"(" + FREQ_TOKEN + r")\s+"    # Frequency
    r"(\d+|-)\s+"               # gCoverage-q30
    r"(" + NUMBER_TOKEN + r")\s+"  # gMeanQ
    r"(" + COUNT_TOKEN + r")\s+"   # gBaseCount
    r"(.+?)\s+"                  # gAllSubs
    r"(" + FREQ_TOKEN + r")\s*$"   # gFrequency
)


class EditPipelineError(RuntimeError):
    """Raised for user-fixable pipeline errors."""


class Interval:
    def __init__(self, start, end, data):
        self.start = start
        self.end = end
        self.data = data


class BinnedIntervalIndex:
    """A small point-query interval index based on fixed-size genomic bins."""

    def __init__(self, bin_size: int = 10000) -> None:
        if bin_size <= 0:
            raise ValueError("bin_size must be positive")
        self.bin_size = bin_size
        self._bins: Dict[str, Dict[int, List[Interval]]] = defaultdict(lambda: defaultdict(list))
        self.interval_count = 0

    def add(self, chrom: str, start: int, end: int, data: Optional[dict] = None) -> None:
        if start < 0:
            start = 0
        if end <= start:
            return
        interval = Interval(start=start, end=end, data=data or {})
        first_bin = start // self.bin_size
        last_bin = (end - 1) // self.bin_size
        for bin_id in range(first_bin, last_bin + 1):
            self._bins[chrom][bin_id].append(interval)
        self.interval_count += 1

    def query_point(self, chrom: str, pos0: int) -> List[Interval]:
        if pos0 < 0:
            return []
        bin_id = pos0 // self.bin_size
        hits = []
        for iv in self._bins.get(chrom, {}).get(bin_id, []):
            if iv.start <= pos0 < iv.end:
                hits.append(iv)
        return hits

    def __len__(self) -> int:
        return self.interval_count


class ContigResolver:
    """Resolve common contig naming differences such as chr1 versus 1."""

    def __init__(self, names: Sequence[str]) -> None:
        self.names = set(names)
        self.alias: Dict[str, str] = {}
        for name in names:
            self.alias[name] = name
            if name.startswith("chr"):
                self.alias.setdefault(name[3:], name)
            else:
                self.alias.setdefault("chr" + name, name)
        if "MT" in self.names:
            self.alias.setdefault("M", "MT")
            self.alias.setdefault("chrM", "MT")
        if "M" in self.names:
            self.alias.setdefault("MT", "M")
            self.alias.setdefault("chrM", "M")
        if "chrM" in self.names:
            self.alias.setdefault("M", "chrM")
            self.alias.setdefault("MT", "chrM")

    def resolve(self, chrom: str) -> Optional[str]:
        return self.alias.get(chrom)


class BamContext:
    def __init__(self, label, bam, resolver):
        self.label = label
        self.bam = bam
        self.resolver = resolver


class GtfData:
    def __init__(self, genes, splice_sites, loaded_gene_intervals, loaded_splice_intervals):
        self.genes = genes
        self.splice_sites = splice_sites
        self.loaded_gene_intervals = loaded_gene_intervals
        self.loaded_splice_intervals = loaded_splice_intervals


class Evaluation:
    def __init__(self, passed, reasons, record, event_class=None, editing_level=None, event_group=None):
        self.passed = passed
        self.reasons = reasons
        self.record = record
        self.event_class = event_class
        self.editing_level = editing_level
        self.event_group = event_group


def open_text(path: str):
    """Open plain text, gzip text, or stdin."""
    if path == "-":
        return sys.stdin
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def is_int_string(value: str) -> bool:
    try:
        int(value)
        return True
    except ValueError:
        return False


def is_number_string(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


MISSING_FIELD_VALUES = {"", "-", ".", "NA", "N/A", "nan", "NaN"}


def is_missing_field(value: str) -> bool:
    return value is None or str(value).strip() in MISSING_FIELD_VALUES


def parse_int(value: str, field: str, line_no: int, missing_as=None) -> int:
    if is_missing_field(value):
        if missing_as is not None:
            return missing_as
        raise ValueError(f"line {line_no}: missing integer field {field}={value!r}")
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"line {line_no}: could not parse integer field {field}={value!r}") from exc


def parse_float(value: str, field: str, line_no: int) -> float:
    value = str(value).strip()
    if is_missing_field(value):
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"line {line_no}: could not parse float field {field}={value!r}") from exc


def parse_base_count(value: str, line_no: int, missing_as_zero: bool = False) -> List[int]:
    value = str(value).strip()
    if is_missing_field(value):
        if missing_as_zero:
            return [0, 0, 0, 0]
        raise ValueError(f"line {line_no}: missing base-count field {value!r}")
    if not (value.startswith("[") and value.endswith("]")):
        raise ValueError(f"line {line_no}: malformed base-count field {value!r}")
    parts = [x.strip() for x in value[1:-1].split(",")]
    if len(parts) != 4:
        raise ValueError(f"line {line_no}: expected 4 base counts but found {len(parts)}")
    try:
        counts = [int(x) for x in parts]
    except ValueError as exc:
        raise ValueError(f"line {line_no}: base counts must be integers: {value!r}") from exc
    if any(x < 0 for x in counts):
        raise ValueError(f"line {line_no}: base counts cannot be negative: {value!r}")
    return counts


def canonical_reditools_name(name):
    """Normalize REDItools column names while keeping the public output names stable."""
    key = str(name).strip().lstrip("#").strip().lower()
    key = re.sub(r"\s+", "", key)
    key = key.replace("_", "-")
    aliases = {
        "region": "Region",
        "chrom": "Region",
        "chromosome": "Region",
        "seqid": "Region",
        "position": "Position",
        "pos": "Position",
        "reference": "Reference",
        "ref": "Reference",
        "strand": "Strand",
        "coverage-q30": "Coverage-q30",
        "coverage": "Coverage-q30",
        "meanq": "MeanQ",
        "basecount[a,c,g,t]": "BaseCount[A,C,G,T]",
        "basecount": "BaseCount[A,C,G,T]",
        "allsubs": "AllSubs",
        "frequency": "Frequency",
        "gcov": "gCoverage-q30",
        "gcoverage": "gCoverage-q30",
        "gcoverage-q30": "gCoverage-q30",
        "gmeanq": "gMeanQ",
        "gbasecount[a,c,g,t]": "gBaseCount[A,C,G,T]",
        "gbasecount": "gBaseCount[A,C,G,T]",
        "gallsubs": "gAllSubs",
        "gfrequency": "gFrequency",
    }
    return aliases.get(key, str(name).strip())


def looks_like_reditools_header(line):
    low = line.lower()
    return ("region" in low and "position" in low and "reference" in low and
            ("basecount" in low or "base count" in low))


def split_delimited_reditools_line(line, delimiter):
    stripped = line.rstrip("\n\r")
    if delimiter == "\t":
        return [x.strip() for x in stripped.split("\t")]
    if delimiter == ",":
        reader = csv.reader([stripped])
        return [x.strip() for x in next(reader)]
    return None


def split_reditools_line(line, delimiter=None):
    """
    Split a REDItools row.

    REDItools tables are usually tab-delimited, but rendered examples are often
    pasted as whitespace-delimited text. BaseCount fields contain spaces inside
    brackets, so a naive line.split() corrupts the columns. This function first
    respects an explicit delimiter from the header, then falls back to a regex
    that treats bracketed base-count vectors as single fields.
    """
    stripped = line.rstrip("\n\r")
    if delimiter in ("\t", ","):
        parts = split_delimited_reditools_line(stripped, delimiter)
        if parts and len(parts) >= 5:
            return parts
    if "\t" in stripped:
        parts = [x.strip() for x in stripped.split("\t")]
        if len(parts) >= 5:
            return parts
    if looks_like_reditools_header(stripped):
        parts = re.split(r"\s+", stripped.strip())
        if len(parts) >= 5:
            return parts
    match = REDITOOLS_SPACE_RE.match(stripped.strip())
    if match:
        return list(match.groups())
    return None


def build_reditools_header_map(header_parts):
    """Return canonical column positions for a REDItools header line."""
    mapping = {}
    for i, name in enumerate(header_parts):
        canonical = canonical_reditools_name(name)
        if canonical in REDITOOLS_COLUMNS and canonical not in mapping:
            mapping[canonical] = i
    required = ["Region", "Position", "Reference", "Coverage-q30", "MeanQ", "BaseCount[A,C,G,T]"]
    missing = [c for c in required if c not in mapping]
    if missing:
        raise EditPipelineError(
            "REDItools header is missing required column(s): " + ", ".join(missing)
        )
    return mapping


def row_from_parts(parts, header_map=None):
    """Build a canonical REDItools row dict from split fields."""
    if header_map is not None:
        row = {}
        for col in REDITOOLS_COLUMNS:
            idx = header_map.get(col)
            row[col] = parts[idx].strip() if idx is not None and idx < len(parts) else "-"
        return row
    # No usable header was seen. Fall back to the classic 14-column REDItools order.
    if len(parts) < len(REDITOOLS_COLUMNS):
        return None
    return dict(zip(REDITOOLS_COLUMNS, [x.strip() for x in parts[:len(REDITOOLS_COLUMNS)]]))


def validate_reditools_row_fields(row, line_no):
    """Parse numeric/count fields in place and annotate missing DNA evidence."""
    row["Position"] = parse_int(row["Position"], "Position", line_no)
    row["Coverage-q30"] = parse_int(row["Coverage-q30"], "Coverage-q30", line_no, missing_as=0)
    row["MeanQ"] = parse_float(row["MeanQ"], "MeanQ", line_no)
    row["BaseCount[A,C,G,T]"] = parse_base_count(row["BaseCount[A,C,G,T]"], line_no)
    # REDItools sometimes writes '-' for genomic DNA fields. Treat those as missing
    # instead of malformed: coverage/counts become zero and quality becomes NA.
    row["_missing_gCoverage_q30"] = is_missing_field(row["gCoverage-q30"])
    row["_missing_gMeanQ"] = is_missing_field(row["gMeanQ"])
    row["_missing_gBaseCount"] = is_missing_field(row["gBaseCount[A,C,G,T]"])
    row["gCoverage-q30"] = parse_int(row["gCoverage-q30"], "gCoverage-q30", line_no, missing_as=0)
    row["gMeanQ"] = parse_float(row["gMeanQ"], "gMeanQ", line_no)
    row["gBaseCount[A,C,G,T]"] = parse_base_count(row["gBaseCount[A,C,G,T]"], line_no, missing_as_zero=True)
    row["Reference"] = str(row["Reference"]).strip().upper()
    row["Region"] = str(row["Region"]).strip()
    row["Strand"] = str(row.get("Strand", ".")).strip()
    row["AllSubs"] = str(row.get("AllSubs", "-")).strip()
    row["Frequency"] = str(row.get("Frequency", "-")).strip()
    row["gAllSubs"] = str(row.get("gAllSubs", "-")).strip()
    row["gFrequency"] = str(row.get("gFrequency", "-")).strip()
    return row


def parse_reditools(path):
    """
    Yield parsed REDItools rows as dictionaries.

    This parser is header-aware. The attached example uses the canonical order:
    Region, Position, Reference, Strand, Coverage-q30, MeanQ,
    BaseCount[A,C,G,T], AllSubs, Frequency, gCoverage-q30, gMeanQ,
    gBaseCount[A,C,G,T], gAllSubs, gFrequency.
    """
    header_map = None
    delimiter = None
    warned_no_header = False
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if looks_like_reditools_header(stripped):
                if "\t" in line:
                    delimiter = "\t"
                elif "," in line and "BaseCount[A,C,G,T]" not in line:
                    delimiter = ","
                else:
                    delimiter = None
                header_parts = split_reditools_line(line, delimiter=delimiter)
                if header_parts is None:
                    raise EditPipelineError("Could not parse REDItools header at line %s" % line_no)
                header_map = build_reditools_header_map(header_parts)
                logging.info(
                    "Detected REDItools columns from header at line %s: %s",
                    line_no,
                    ", ".join([c for c in REDITOOLS_COLUMNS if c in header_map])
                )
                continue
            parts = split_reditools_line(line, delimiter=delimiter)
            if parts is None:
                logging.warning("Skipping malformed REDItools line %s", line_no)
                yield line_no, {"_parse_error": "malformed_reditools_line", "_raw": stripped}
                continue
            if header_map is None and not warned_no_header:
                warned_no_header = True
                logging.warning(
                    "REDItools header was not detected before data line %s; assuming classic 14-column order.",
                    line_no
                )
            row = row_from_parts(parts, header_map=header_map)
            if row is None:
                logging.warning("Skipping malformed REDItools line %s: expected at least %s columns, found %s", line_no, len(REDITOOLS_COLUMNS), len(parts))
                yield line_no, {"_parse_error": "malformed_reditools_line", "_raw": stripped}
                continue
            try:
                validate_reditools_row_fields(row, line_no)
            except ValueError as exc:
                logging.warning("Skipping REDItools line %s: %s", line_no, exc)
                yield line_no, {"_parse_error": str(exc), "_raw": stripped}
                continue
            yield line_no, row


def parse_gtf_attributes(attr_text: str) -> dict:
    attrs = {}
    for field in attr_text.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field and " " not in field.split("=", 1)[0]:
            key, value = field.split("=", 1)
        else:
            parts = field.split(None, 1)
            if len(parts) == 1:
                key, value = parts[0], ""
            else:
                key, value = parts
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def load_gtf(path: str, resolver: ContigResolver, overlap_feature: str, splice_flank: int) -> GtfData:
    genes = BinnedIntervalIndex()
    splice_sites = BinnedIntervalIndex()
    exon_by_transcript: Dict[str, List[Tuple[str, int, int, str, str, str]]] = defaultdict(list)
    exon_by_gene: Dict[str, List[Tuple[str, int, int, str, str]]] = defaultdict(list)
    skipped_contigs = Counter()
    gene_interval_count = 0

    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                logging.warning("Skipping malformed GTF line %s", line_no)
                continue
            seqid, source, feature, start_s, end_s, score, strand, frame, attrs_s = fields
            chrom = resolver.resolve(seqid)
            if chrom is None:
                skipped_contigs[seqid] += 1
                continue
            try:
                start0 = int(start_s) - 1
                end = int(end_s)
            except ValueError:
                logging.warning("Skipping GTF line %s with invalid coordinates", line_no)
                continue
            if start0 < 0 or end <= start0:
                logging.warning("Skipping GTF line %s with invalid interval", line_no)
                continue
            attrs = parse_gtf_attributes(attrs_s)
            gene_id = (
                attrs.get("gene_id")
                or attrs.get("gene")
                or attrs.get("ID")
                or attrs.get("Parent")
                or f"unknown_gene_line_{line_no}"
            )
            gene_name = attrs.get("gene_name") or attrs.get("Name") or gene_id
            transcript_id = attrs.get("transcript_id") or attrs.get("transcript") or attrs.get("Parent")

            if feature == overlap_feature:
                genes.add(
                    chrom,
                    start0,
                    end,
                    {
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "strand": strand,
                        "feature": feature,
                    },
                )
                gene_interval_count += 1

            if feature == "exon":
                exon_by_gene[gene_id].append((chrom, start0, end, strand, gene_name))
                if transcript_id:
                    exon_by_transcript[transcript_id].append((chrom, start0, end, strand, gene_id, gene_name))

    if gene_interval_count == 0:
        logging.warning(
            "No GTF feature=%r intervals found. Falling back to gene spans built from exon records.",
            overlap_feature,
        )
        for gene_id, exons in exon_by_gene.items():
            by_chrom_strand: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)
            for chrom, start0, end, strand, gene_name in exons:
                by_chrom_strand[(chrom, strand, gene_name)].append((start0, end))
            for (chrom, strand, gene_name), intervals in by_chrom_strand.items():
                start0 = min(x[0] for x in intervals)
                end = max(x[1] for x in intervals)
                genes.add(
                    chrom,
                    start0,
                    end,
                    {
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "strand": strand,
                        "feature": "exon_span_fallback",
                    },
                )
                gene_interval_count += 1

    if splice_flank > 0:
        for transcript_id, exons in exon_by_transcript.items():
            if len(exons) < 2:
                continue
            exons_sorted = sorted(exons, key=lambda x: (x[0], x[1], x[2]))
            for left, right in zip(exons_sorted, exons_sorted[1:]):
                left_chrom, left_start, left_end, strand, gene_id, gene_name = left
                right_chrom, right_start, right_end, _, _, _ = right
                if left_chrom != right_chrom:
                    continue
                for boundary in (left_end, right_start):
                    splice_sites.add(
                        left_chrom,
                        max(0, boundary - splice_flank),
                        boundary + splice_flank,
                        {
                            "gene_id": gene_id,
                            "gene_name": gene_name,
                            "transcript_id": transcript_id,
                            "strand": strand,
                        },
                    )

    if skipped_contigs:
        top = ", ".join(f"{k}:{v}" for k, v in skipped_contigs.most_common(5))
        logging.warning("Skipped GTF records on contigs absent from FASTA: %s", top)

    if len(genes) == 0:
        raise EditPipelineError("No usable gene intervals were loaded from the GTF file.")

    return GtfData(
        genes=genes,
        splice_sites=splice_sites,
        loaded_gene_intervals=len(genes),
        loaded_splice_intervals=len(splice_sites),
    )


def parse_repeatmasker_interval(parts: List[str], fmt: str) -> Optional[Tuple[str, int, int, str]]:
    """Return chrom, start0, end, repeat_name or None.

    Supported repeat inputs:
      * BED: chrom start0 end [name ...]
      * GFF/GTF-like: chrom source feature start1 end1 ... attributes
      * RepeatMasker .out:
          SW perc perc perc query begin end (left) strand repeat class/family ... ID

    RepeatMasker coordinates are 1-based closed, so begin is converted to
    0-based half-open for internal overlap tests.
    """
    requested_rmout = fmt in {"auto", "rmout", "repeatmasker"}

    if fmt in {"auto", "bed"} and len(parts) >= 3 and is_int_string(parts[1]) and is_int_string(parts[2]):
        chrom = parts[0]
        start0 = int(parts[1])
        end = int(parts[2])
        name = parts[3] if len(parts) >= 4 else "repeat"
        return chrom, start0, end, name

    if fmt in {"auto", "gff"} and len(parts) >= 9 and is_int_string(parts[3]) and is_int_string(parts[4]):
        chrom = parts[0]
        start0 = int(parts[3]) - 1
        end = int(parts[4])
        attrs = parse_gtf_attributes(parts[8])
        name = attrs.get("Name") or attrs.get("Target") or parts[2]
        return chrom, start0, end, name

    # Standard RepeatMasker .out row, e.g.:
    # 1435 15.5 11.8 2.8 NC_046332.1 2 365 (117953686) C DR0009509 LTR/ERVL (535) 397 2 1
    #   18 14.7  0.0 0.0 NC_046332.1 2065 2094 (117951957) + (T)n Simple_repeat 1 30 (0) 2
    # A trailing asterisk can be present in some .out files, so only the
    # coordinate/name columns needed here are required.
    if requested_rmout and len(parts) >= 11 and is_number_string(parts[0]):
        chrom = parts[4]
        if not (is_int_string(parts[5]) and is_int_string(parts[6])):
            return None
        start1 = int(parts[5])
        end1 = int(parts[6])
        if start1 <= 0 or end1 <= 0:
            return None
        start0 = min(start1, end1) - 1
        end = max(start1, end1)
        repeat_name = parts[9]
        repeat_class = parts[10] if len(parts) > 10 else "repeat"
        name = f"{repeat_name}|{repeat_class}"
        return chrom, start0, end, name

    return None


def load_repeats(path: str, resolver: ContigResolver, fmt: str) -> BinnedIntervalIndex:
    repeats = BinnedIntervalIndex()
    skipped = Counter()
    malformed = 0
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.lower().startswith(("sw ", "score ", "there were", "repeatmasker")):
                continue
            parts = stripped.split()
            parsed = parse_repeatmasker_interval(parts, fmt)
            if parsed is None:
                malformed += 1
                if malformed <= 5:
                    logging.warning("Skipping unrecognized repeat line %s", line_no)
                continue
            chrom_raw, start0, end, name = parsed
            chrom = resolver.resolve(chrom_raw)
            if chrom is None:
                skipped[chrom_raw] += 1
                continue
            if end <= start0:
                malformed += 1
                continue
            repeats.add(chrom, start0, end, {"repeat_name": name})
    if skipped:
        top = ", ".join(f"{k}:{v}" for k, v in skipped.most_common(5))
        logging.warning("Skipped repeat records on contigs absent from FASTA: %s", top)
    logging.info("Loaded %s repeat intervals", len(repeats))
    return repeats


def validate_file(path: str, label: str, allow_stdin: bool = False) -> None:
    if allow_stdin and path == "-":
        return
    if not path:
        raise EditPipelineError(f"Missing path for {label}.")
    if not os.path.exists(path):
        raise EditPipelineError(f"{label} does not exist: {path}")
    if not os.path.isfile(path):
        raise EditPipelineError(f"{label} is not a regular file: {path}")


def open_fasta(path: str) -> pysam.FastaFile:
    try:
        fasta = pysam.FastaFile(path)
    except Exception as exc:
        raise EditPipelineError(
            f"Could not open genome FASTA {path!r}. Make sure it is indexed with: samtools faidx {path}"
        ) from exc
    if not fasta.references:
        raise EditPipelineError("The genome FASTA has no references.")
    return fasta


def open_bam(path: str, label: str) -> pysam.AlignmentFile:
    try:
        bam = pysam.AlignmentFile(path, "rb")
    except Exception as exc:
        raise EditPipelineError(f"Could not open {label} BAM {path!r}.") from exc
    try:
        has_index = bam.has_index()
    except Exception:
        has_index = False
    if not has_index:
        raise EditPipelineError(
            f"{label} BAM is not indexed or the index is not readable: {path}. Run: samtools index {path}"
        )
    return bam


def fasta_base(fasta: pysam.FastaFile, chrom: str, pos0: int) -> Optional[str]:
    try:
        if pos0 < 0 or pos0 >= fasta.get_reference_length(chrom):
            return None
        seq = fasta.fetch(chrom, pos0, pos0 + 1).upper()
    except Exception:
        return None
    return seq if len(seq) == 1 else None


def homopolymer_run_length(
    fasta: pysam.FastaFile,
    chrom: str,
    pos0: int,
    threshold: int,
) -> int:
    if threshold <= 1:
        return 1
    contig_len = fasta.get_reference_length(chrom)
    start = max(0, pos0 - threshold + 1)
    end = min(contig_len, pos0 + threshold)
    seq = fasta.fetch(chrom, start, end).upper()
    idx = pos0 - start
    if idx < 0 or idx >= len(seq):
        return 0
    base = seq[idx]
    if base not in BASE_TO_I:
        return 0
    run = 1
    left = idx - 1
    while left >= 0 and seq[left] == base:
        run += 1
        left -= 1
    right = idx + 1
    while right < len(seq) and seq[right] == base:
        run += 1
        right += 1
    return run


def usable_alignment(aln, ignore_duplicates: bool, ignore_qcfail: bool) -> bool:
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        return False
    if ignore_duplicates and aln.is_duplicate:
        return False
    if ignore_qcfail and aln.is_qcfail:
        return False
    return True


def site_mapq_stats(
    bam_ctx: BamContext,
    chrom: str,
    pos0: int,
    min_mapq: int,
    min_mean_mapq: float,
    min_good_frac: float,
    max_depth: int,
    ignore_duplicates: bool,
    ignore_qcfail: bool,
) -> Tuple[bool, dict, str]:
    bam_chrom = bam_ctx.resolver.resolve(chrom)
    if bam_chrom is None:
        return False, {}, "contig_not_found_in_bam"
    total = 0
    good = 0
    mapq_sum = 0
    try:
        pileup_iter = bam_ctx.bam.pileup(
            bam_chrom,
            pos0,
            pos0 + 1,
            truncate=True,
            stepper="all",
            max_depth=max_depth,
            min_base_quality=0,
        )
        for col in pileup_iter:
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                aln = pr.alignment
                if not usable_alignment(aln, ignore_duplicates, ignore_qcfail):
                    continue
                if pr.is_del or pr.is_refskip:
                    continue
                total += 1
                mq = int(aln.mapping_quality)
                mapq_sum += mq
                if mq >= min_mapq:
                    good += 1
    except ValueError as exc:
        return False, {}, f"bam_pileup_failed:{exc}"
    except Exception as exc:
        return False, {}, f"bam_pileup_failed:{type(exc).__name__}"

    if total == 0:
        return False, {"depth": 0, "mean_mapq": math.nan, "good_mapq_fraction": math.nan}, "no_usable_bam_reads"
    mean_mapq = mapq_sum / total
    good_frac = good / total
    stats = {"depth": total, "mean_mapq": mean_mapq, "good_mapq_fraction": good_frac}
    if mean_mapq < min_mean_mapq:
        return False, stats, "low_mean_mapq"
    if good_frac < min_good_frac:
        return False, stats, "low_good_mapq_fraction"
    return True, stats, "pass"


def bam_base_counts_at_site(
    bam_ctx: BamContext,
    chrom: str,
    pos0: int,
    min_mapq: int,
    min_baseq: int,
    max_depth: int,
    ignore_duplicates: bool,
    ignore_qcfail: bool,
) -> Tuple[bool, dict, str]:
    """Count A/C/G/T bases in a BAM at one site for DNA fallback genotyping."""
    bam_chrom = bam_ctx.resolver.resolve(chrom)
    if bam_chrom is None:
        return False, {}, "contig_not_found_in_bam"
    counts = [0, 0, 0, 0]
    qual_sum = 0
    try:
        pileup_iter = bam_ctx.bam.pileup(
            bam_chrom,
            pos0,
            pos0 + 1,
            truncate=True,
            stepper="all",
            max_depth=max_depth,
            min_base_quality=0,
        )
        for col in pileup_iter:
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                aln = pr.alignment
                if not usable_alignment(aln, ignore_duplicates, ignore_qcfail):
                    continue
                if int(aln.mapping_quality) < min_mapq:
                    continue
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                seq = aln.query_sequence
                if not seq or pr.query_position >= len(seq):
                    continue
                base = seq[pr.query_position].upper()
                if base not in BASE_TO_I:
                    continue
                q = 0
                if aln.query_qualities is not None and pr.query_position < len(aln.query_qualities):
                    q = int(aln.query_qualities[pr.query_position])
                if q < min_baseq:
                    continue
                counts[BASE_TO_I[base]] += 1
                qual_sum += q
    except ValueError as exc:
        return False, {}, f"bam_pileup_failed:{exc}"
    except Exception as exc:
        return False, {}, f"bam_pileup_failed:{type(exc).__name__}"

    depth = sum(counts)
    if depth == 0:
        return False, {"counts": counts, "coverage": 0, "mean_baseq": math.nan}, "no_usable_dna_bases"
    return True, {"counts": counts, "coverage": depth, "mean_baseq": qual_sum / depth}, "pass"


def indel_stats_in_window(
    bam_ctx: BamContext,
    chrom: str,
    pos0: int,
    flank: int,
    max_indel_frac: float,
    max_indel_reads: int,
    max_depth: int,
    ignore_duplicates: bool,
    ignore_qcfail: bool,
) -> Tuple[bool, dict, str]:
    if flank <= 0:
        return True, {"max_indel_fraction": 0.0, "max_indel_reads": 0}, "pass"
    bam_chrom = bam_ctx.resolver.resolve(chrom)
    if bam_chrom is None:
        return False, {}, "contig_not_found_in_bam"
    start = max(0, pos0 - flank)
    end = pos0 + flank + 1
    max_frac_seen = 0.0
    max_indel_seen = 0
    columns_seen = 0
    try:
        pileup_iter = bam_ctx.bam.pileup(
            bam_chrom,
            start,
            end,
            truncate=True,
            stepper="all",
            max_depth=max_depth,
            min_base_quality=0,
        )
        for col in pileup_iter:
            if not (start <= col.reference_pos < end):
                continue
            columns_seen += 1
            depth = 0
            indel_reads = 0
            for pr in col.pileups:
                aln = pr.alignment
                if not usable_alignment(aln, ignore_duplicates, ignore_qcfail):
                    continue
                depth += 1
                if pr.indel != 0 or pr.is_del or pr.is_refskip:
                    indel_reads += 1
            if depth > 0:
                frac = indel_reads / depth
                max_frac_seen = max(max_frac_seen, frac)
                max_indel_seen = max(max_indel_seen, indel_reads)
    except ValueError as exc:
        return False, {}, f"bam_pileup_failed:{exc}"
    except Exception as exc:
        return False, {}, f"bam_pileup_failed:{type(exc).__name__}"

    stats = {"max_indel_fraction": max_frac_seen, "max_indel_reads": max_indel_seen, "columns_seen": columns_seen}
    if max_indel_seen > max_indel_reads:
        return False, stats, "too_many_indel_reads"
    if max_frac_seen > max_indel_frac:
        return False, stats, "high_indel_fraction"
    return True, stats, "pass"


def mismatch_stats_in_window(
    bam_ctx: BamContext,
    fasta: pysam.FastaFile,
    chrom: str,
    pos0: int,
    flank: int,
    min_mapq: int,
    min_baseq: int,
    max_window_mismatch_fraction: float,
    max_column_mismatch_fraction: float,
    max_high_mismatch_columns: int,
    max_depth: int,
    ignore_duplicates: bool,
    ignore_qcfail: bool,
    gene_strand: str,
    ignore_event_type: str,
    count_same_event_mismatches: bool,
) -> Tuple[bool, dict, str]:
    """Detect locally mismapped regions by measuring mismatch density around a site.

    The candidate position itself is excluded so a genuine edit does not count
    against the window. By default, adjacent mismatches matching the same
    transcript-oriented substitution are not counted, which keeps true local
    A>G or C>T clusters from looking like mapping artefacts. A site fails when
    adjacent columns show too many other bases that disagree with the reference FASTA.
    """
    if flank <= 0:
        return True, {
            "window_mismatch_fraction": 0.0,
            "max_mismatch_fraction": 0.0,
            "high_mismatch_columns": 0,
            "mismatch_columns_seen": 0,
            "mismatch_bases": 0,
            "window_bases": 0,
        }, "pass"
    bam_chrom = bam_ctx.resolver.resolve(chrom)
    if bam_chrom is None:
        return False, {}, "contig_not_found_in_bam"

    try:
        contig_len = fasta.get_reference_length(chrom)
    except Exception:
        return False, {}, "fasta_reference_length_failed"
    start = max(0, pos0 - flank)
    end = min(contig_len, pos0 + flank + 1)
    if end <= start:
        return False, {}, "empty_mismatch_window"
    try:
        ref_seq = fasta.fetch(chrom, start, end).upper()
    except Exception as exc:
        return False, {}, "fasta_fetch_failed:%s" % type(exc).__name__

    total_bases = 0
    mismatch_bases = 0
    max_frac_seen = 0.0
    high_columns = 0
    columns_seen = 0
    try:
        pileup_iter = bam_ctx.bam.pileup(
            bam_chrom,
            start,
            end,
            truncate=True,
            stepper="all",
            max_depth=max_depth,
            min_base_quality=0,
        )
        for col in pileup_iter:
            ref_pos = col.reference_pos
            if not (start <= ref_pos < end):
                continue
            if ref_pos == pos0:
                continue
            idx = ref_pos - start
            if idx < 0 or idx >= len(ref_seq):
                continue
            ref_base = ref_seq[idx]
            if ref_base not in BASE_TO_I:
                continue
            depth = 0
            mismatches = 0
            for pr in col.pileups:
                aln = pr.alignment
                if not usable_alignment(aln, ignore_duplicates, ignore_qcfail):
                    continue
                if int(aln.mapping_quality) < min_mapq:
                    continue
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                seq = aln.query_sequence
                if not seq or pr.query_position >= len(seq):
                    continue
                base = seq[pr.query_position].upper()
                if base not in BASE_TO_I:
                    continue
                q = 0
                if aln.query_qualities is not None and pr.query_position < len(aln.query_qualities):
                    q = int(aln.query_qualities[pr.query_position])
                if q < min_baseq:
                    continue
                depth += 1
                total_bases += 1
                if base != ref_base:
                    if (not count_same_event_mismatches) and gene_strand in {"+", "-"} and ignore_event_type:
                        try:
                            local_change = transcript_change(ref_base, base, gene_strand)
                        except Exception:
                            local_change = ""
                        if local_change == ignore_event_type:
                            continue
                    mismatches += 1
                    mismatch_bases += 1
            if depth > 0:
                columns_seen += 1
                frac = mismatches / depth
                max_frac_seen = max(max_frac_seen, frac)
                if frac > max_column_mismatch_fraction:
                    high_columns += 1
    except ValueError as exc:
        return False, {}, "bam_pileup_failed:%s" % exc
    except Exception as exc:
        return False, {}, "bam_pileup_failed:%s" % type(exc).__name__

    window_frac = (mismatch_bases / total_bases) if total_bases > 0 else 0.0
    stats = {
        "window_mismatch_fraction": window_frac,
        "max_mismatch_fraction": max_frac_seen,
        "high_mismatch_columns": high_columns,
        "mismatch_columns_seen": columns_seen,
        "mismatch_bases": mismatch_bases,
        "window_bases": total_bases,
    }
    if total_bases == 0:
        return True, stats, "pass"
    if window_frac > max_window_mismatch_fraction:
        return False, stats, "high_window_mismatch_fraction"
    if max_frac_seen > max_column_mismatch_fraction:
        return False, stats, "high_column_mismatch_fraction"
    if high_columns > max_high_mismatch_columns:
        return False, stats, "too_many_high_mismatch_columns"
    return True, stats, "pass"


def major_rna_alt(ref: str, counts: List[int]) -> Tuple[Optional[str], int, float, List[str]]:
    if ref not in BASE_TO_I:
        return None, 0, 0.0, []
    total = sum(counts)
    if total <= 0:
        return None, 0, 0.0, []
    ref_i = BASE_TO_I[ref]
    alt_items = [(base, counts[i]) for i, base in enumerate(BASES) if i != ref_i and counts[i] > 0]
    if not alt_items:
        return None, 0, 0.0, []
    alt_items.sort(key=lambda x: (-x[1], x[0]))
    alt_base, alt_count = alt_items[0]
    secondary = [f"{base}:{count}" for base, count in alt_items[1:]]
    return alt_base, alt_count, alt_count / total, secondary


def transcript_change(ref: str, alt: str, strand: str) -> str:
    if strand == "+":
        return f"{ref}>{alt}"
    if strand == "-":
        tref = ref.translate(COMPLEMENT).upper()
        talt = alt.translate(COMPLEMENT).upper()
        return f"{tref}>{talt}"
    return "NA"


def classify_event(ref: str, alt: str, gene_hits: List[Interval]) -> Tuple[str, str, str, str, str, str]:
    """Classify by transcript-oriented substitution without collapsing Other types.

    Returns:
      event_group: ADAR, APOBEC, or Other
      event_type: transcript-oriented change such as A>G or C>T
      event_label: label used in plots. Usually same as event_type; ambiguous
                   strand cases are labelled ambiguous_REF>ALT.
      gene_ids, gene_names, gene_strand_used
    """
    gene_ids = sorted({iv.data.get("gene_id", "NA") for iv in gene_hits})
    gene_names = sorted({iv.data.get("gene_name", "NA") for iv in gene_hits})
    strands = sorted({iv.data.get("strand", ".") for iv in gene_hits if iv.data.get("strand") in {"+", "-"}})
    if len(strands) != 1:
        genomic_change = f"{ref}>{alt}"
        event_type = "ambiguous_" + genomic_change
        return "Other", event_type, event_type, ",".join(gene_ids), ",".join(gene_names), "ambiguous_or_missing_gene_strand"

    strand = strands[0]
    event_type = transcript_change(ref, alt, strand)
    if event_type == "A>G":
        event_group = "ADAR"
    elif event_type == "C>T":
        event_group = "APOBEC"
    else:
        event_group = "Other"
    return event_group, event_type, event_type, ",".join(gene_ids), ",".join(gene_names), strand

def fmt_float(value: float, digits: int = 6) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.{digits}f}"


def evaluate_site(
    row: dict,
    line_no: int,
    fasta: pysam.FastaFile,
    fasta_resolver: ContigResolver,
    repeats: BinnedIntervalIndex,
    gtf_data: GtfData,
    bam_contexts: List[BamContext],
    args,
) -> Evaluation:
    if "_parse_error" in row:
        return Evaluation(False, ["parse_error"], None)

    reasons: List[str] = []
    chrom_raw = str(row["Region"])
    chrom = fasta_resolver.resolve(chrom_raw)
    if chrom is None:
        return Evaluation(False, ["contig_not_found_in_fasta"], None)

    pos1 = int(row["Position"])
    pos0 = pos1 - 1
    ref = str(row["Reference"]).upper()
    if ref not in BASE_TO_I:
        return Evaluation(False, ["non_acgt_reference"], None)

    genome_ref = fasta_base(fasta, chrom, pos0)
    if genome_ref is None:
        return Evaluation(False, ["position_out_of_range_or_fasta_fetch_failed"], None)
    if genome_ref != ref and not args.allow_ref_mismatch:
        return Evaluation(False, [f"reference_mismatch_fasta:{genome_ref}_reditools:{ref}"], None)

    rna_counts: List[int] = row["BaseCount[A,C,G,T]"]
    gdna_counts: List[int] = list(row["gBaseCount[A,C,G,T]"])
    alt, alt_count, edit_freq, secondary = major_rna_alt(ref, rna_counts)
    if alt is None:
        return Evaluation(False, ["no_rna_substitution"], None)

    rna_count_depth = sum(rna_counts)
    gdna_count_depth = sum(gdna_counts)
    rna_cov = int(row["Coverage-q30"])
    gdna_cov = int(row["gCoverage-q30"])
    rna_meanq = float(row["MeanQ"])
    gdna_meanq = float(row["gMeanQ"])
    missing_reditools_dna = bool(row.get("_missing_gCoverage_q30") or row.get("_missing_gBaseCount"))
    dna_count_source = "reditools"

    if row.get("_missing_gCoverage_q30") and gdna_count_depth > 0:
        gdna_cov = gdna_count_depth

    if missing_reditools_dna:
        dna_count_source = "missing_reditools_dna"
        if args.use_dna_bam_counts_if_reditools_dna_missing:
            dna_ctx = None
            for bam_ctx in bam_contexts:
                if bam_ctx.label == "dna":
                    dna_ctx = bam_ctx
                    break
            if dna_ctx is not None:
                ok, dna_stats, dna_reason = bam_base_counts_at_site(
                    bam_ctx=dna_ctx,
                    chrom=chrom,
                    pos0=pos0,
                    min_mapq=args.min_mapq,
                    min_baseq=args.dna_variant_min_baseq,
                    max_depth=args.max_pileup_depth,
                    ignore_duplicates=not args.include_duplicates,
                    ignore_qcfail=not args.include_qcfail,
                )
                if ok:
                    gdna_counts = list(dna_stats["counts"])
                    gdna_count_depth = sum(gdna_counts)
                    gdna_cov = int(dna_stats["coverage"])
                    gdna_meanq = float(dna_stats["mean_baseq"])
                    dna_count_source = "dna_bam_fallback"
                else:
                    dna_count_source = "dna_bam_fallback_failed:" + dna_reason


    if rna_cov < args.min_rna_coverage:
        reasons.append("low_rna_coverage")
    if gdna_cov < args.min_dna_coverage:
        reasons.append("low_dna_coverage")
    if not math.isnan(rna_meanq) and rna_meanq < args.min_rna_meanq:
        reasons.append("low_rna_mean_base_quality")
    if not math.isnan(gdna_meanq) and gdna_meanq < args.min_dna_meanq:
        reasons.append("low_dna_mean_base_quality")
    if alt_count < args.min_alt_count:
        reasons.append("low_rna_alt_count")
    if edit_freq < args.min_edit_freq:
        reasons.append("low_editing_frequency")

    alt_i = BASE_TO_I[alt]
    gdna_alt_count = gdna_counts[alt_i]
    gdna_alt_freq = (gdna_alt_count / gdna_count_depth) if gdna_count_depth > 0 else 0.0
    if gdna_alt_freq > args.max_dna_alt_freq:
        reasons.append("high_dna_alt_frequency_possible_snv")
    if gdna_alt_count > args.max_dna_alt_count:
        reasons.append("high_dna_alt_count_possible_snv")

    repeat_hits = repeats.query_point(chrom, pos0)
    if repeat_hits:
        reasons.append("repeat_overlap")

    hp_len = homopolymer_run_length(fasta, chrom, pos0, args.homopolymer_length)
    if args.homopolymer_length > 1 and hp_len >= args.homopolymer_length:
        reasons.append("homopolymer_run")

    splice_hits = gtf_data.splice_sites.query_point(chrom, pos0) if args.splice_flank > 0 else []
    if splice_hits:
        reasons.append("splice_site_proximity")

    gene_hits = gtf_data.genes.query_point(chrom, pos0)
    if not gene_hits:
        reasons.append("no_gene_model_overlap")

    if reasons:
        return Evaluation(False, reasons, None)

    event_group, event_type, event_label, gene_ids, gene_names, gene_strand = classify_event(ref, alt, gene_hits)

    mapq_output = {}
    if not args.skip_bam_mapq_filter:
        for bam_ctx in bam_contexts:
            ok, stats, reason = site_mapq_stats(
                bam_ctx=bam_ctx,
                chrom=chrom,
                pos0=pos0,
                min_mapq=args.min_mapq,
                min_mean_mapq=args.min_mean_mapq,
                min_good_frac=args.min_good_mapq_fraction,
                max_depth=args.max_pileup_depth,
                ignore_duplicates=not args.include_duplicates,
                ignore_qcfail=not args.include_qcfail,
            )
            mapq_output[f"{bam_ctx.label}_bam_depth"] = stats.get("depth", "NA")
            mapq_output[f"{bam_ctx.label}_mean_mapq"] = fmt_float(stats.get("mean_mapq", math.nan), 3)
            mapq_output[f"{bam_ctx.label}_good_mapq_fraction"] = fmt_float(stats.get("good_mapq_fraction", math.nan), 3)
            if not ok:
                reasons.append(f"{bam_ctx.label}_{reason}")

    indel_output = {}
    if not args.skip_bam_indel_filter:
        for bam_ctx in bam_contexts:
            ok, stats, reason = indel_stats_in_window(
                bam_ctx=bam_ctx,
                chrom=chrom,
                pos0=pos0,
                flank=args.indel_flank,
                max_indel_frac=args.max_indel_fraction,
                max_indel_reads=args.max_indel_reads,
                max_depth=args.max_pileup_depth,
                ignore_duplicates=not args.include_duplicates,
                ignore_qcfail=not args.include_qcfail,
            )
            indel_output[f"{bam_ctx.label}_max_indel_fraction"] = fmt_float(
                stats.get("max_indel_fraction", math.nan), 3
            )
            indel_output[f"{bam_ctx.label}_max_indel_reads"] = stats.get("max_indel_reads", "NA")
            if not ok:
                reasons.append(f"{bam_ctx.label}_{reason}")

    mismatch_output = {}
    if not args.skip_bam_mismatch_filter:
        for bam_ctx in bam_contexts:
            if args.mismatch_bams != "both" and args.mismatch_bams != bam_ctx.label:
                continue
            ok, stats, reason = mismatch_stats_in_window(
                bam_ctx=bam_ctx,
                fasta=fasta,
                chrom=chrom,
                pos0=pos0,
                flank=args.mismatch_flank,
                min_mapq=args.min_mapq,
                min_baseq=args.mismatch_min_baseq,
                max_window_mismatch_fraction=args.max_window_mismatch_fraction,
                max_column_mismatch_fraction=args.max_column_mismatch_fraction,
                max_high_mismatch_columns=args.max_high_mismatch_columns,
                max_depth=args.max_pileup_depth,
                ignore_duplicates=not args.include_duplicates,
                ignore_qcfail=not args.include_qcfail,
                gene_strand=gene_strand,
                ignore_event_type=event_type,
                count_same_event_mismatches=args.count_same_event_mismatches,
            )
            mismatch_output[f"{bam_ctx.label}_window_mismatch_fraction"] = fmt_float(
                stats.get("window_mismatch_fraction", math.nan), 4
            )
            mismatch_output[f"{bam_ctx.label}_max_mismatch_fraction"] = fmt_float(
                stats.get("max_mismatch_fraction", math.nan), 4
            )
            mismatch_output[f"{bam_ctx.label}_high_mismatch_columns"] = stats.get("high_mismatch_columns", "NA")
            mismatch_output[f"{bam_ctx.label}_mismatch_columns_seen"] = stats.get("mismatch_columns_seen", "NA")
            mismatch_output[f"{bam_ctx.label}_mismatch_bases"] = stats.get("mismatch_bases", "NA")
            mismatch_output[f"{bam_ctx.label}_mismatch_window_bases"] = stats.get("window_bases", "NA")
            if not ok:
                reasons.append(f"{bam_ctx.label}_{reason}")

    if reasons:
        return Evaluation(False, reasons, None)

    record = {
        "chrom": chrom,
        "position_1based": pos1,
        "position_0based": pos0,
        "reference": ref,
        "genome_reference": genome_ref,
        "alt_base": alt,
        "genomic_edit": f"{ref}>{alt}",
        "transcript_edit": event_type,
        "event_type": event_type,
        "event_group": event_group,
        "event_label": event_label,
        "event_class": event_label,
        "gene_ids": gene_ids,
        "gene_names": gene_names,
        "gene_strand_used": gene_strand,
        "rna_coverage_q30": rna_cov,
        "rna_count_depth": rna_count_depth,
        "rna_alt_count": alt_count,
        "editing_level": fmt_float(edit_freq, 6),
        "secondary_rna_alts": ";".join(secondary) if secondary else "-",
        "rna_mean_base_quality": fmt_float(rna_meanq, 3),
        "dna_coverage_q30": gdna_cov,
        "dna_count_depth": gdna_count_depth,
        "dna_alt_count_same_base": gdna_alt_count,
        "dna_alt_frequency_same_base": fmt_float(gdna_alt_freq, 6),
        "dna_mean_base_quality": fmt_float(gdna_meanq, 3),
        "dna_count_source": dna_count_source,
        "homopolymer_length": hp_len,
        "reditools_allsubs": row.get("AllSubs", "NA"),
        "reditools_frequency": row.get("Frequency", "NA"),
        "reditools_gallsubs": row.get("gAllSubs", "NA"),
        "reditools_gfrequency": row.get("gFrequency", "NA"),
    }
    record.update(mapq_output)
    record.update(indel_output)
    record.update(mismatch_output)
    return Evaluation(True, [], record, event_class=event_label, editing_level=edit_freq, event_group=event_group)



def event_label_sort_key(label: str):
    if label in TRANSCRIPT_SUBSTITUTION_ORDER:
        return (0, TRANSCRIPT_SUBSTITUTION_ORDER.index(label), label)
    if str(label).startswith("ambiguous_"):
        return (2, 0, label)
    return (1, 0, label)


def sorted_event_labels(class_counts: Counter) -> List[str]:
    observed = [label for label, count in class_counts.items() if count > 0]
    return sorted(observed, key=event_label_sort_key)


def _extract_change_from_label(label: str) -> str:
    """Return the substitution pattern embedded in a plot/event label.

    Examples:
      A>G -> A>G
      ambiguous_T>C -> T>C
      genomic_G>A -> G>A
    """
    if not label:
        return ""
    if ">" in label:
        candidate = label.split("_")[-1]
        if len(candidate) == 3 and candidate[1] == ">":
            return candidate
    return label


ADAR_CHANGES = {"A>G", "T>C"}
APOBEC_CHANGES = {"C>T", "G>A"}


def event_group_from_label(label: str) -> str:
    change = _extract_change_from_label(label)
    if change in ADAR_CHANGES:
        return "ADAR"
    if change in APOBEC_CHANGES:
        return "APOBEC"
    return "Other"


def event_color(label: str, args) -> str:
    group = event_group_from_label(label)
    if group == "ADAR":
        return args.adar_color
    if group == "APOBEC":
        return args.apobec_color
    return args.other_color

def make_bar_plot(class_counts: Counter, out_path: Path, args) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = sorted_event_labels(class_counts)
    if not labels:
        labels = ["no_passed_sites"]
    values = [class_counts.get(label, 0) for label in labels]
    colors = [event_color(label, args) for label in labels]

    fig_width = max(7.0, 0.55 * len(labels) + 3.0)
    fig, ax = plt.subplots(figsize=(fig_width, 4.8))
    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Number of filtered editing sites")
    ax.set_title("RNA editing event types")
    ax.set_ylim(0, max(values + [1]) * 1.18)
    ax.tick_params(axis="x", rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(value), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=args.plot_dpi)
    plt.close(fig)


def make_histogram_plot(edit_levels: Dict[str, List[float]], out_path: Path, args) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = sorted([label for label, values in edit_levels.items() if values], key=event_label_sort_key)
    if not labels:
        labels = ["no_passed_sites"]
    bins = [i / args.hist_bins for i in range(args.hist_bins + 1)]
    fig_height = max(3.0, 2.45 * len(labels))
    fig, axes = plt.subplots(len(labels), 1, figsize=(7.5, fig_height), sharex=True)
    if len(labels) == 1:
        axes = [axes]
    for ax, label in zip(axes, labels):
        data = edit_levels.get(label, [])
        ax.hist(data, bins=bins, color=event_color(label, args), edgecolor="black", linewidth=0.5)
        ax.set_ylabel(f"{label}\ncount")
        ax.text(0.98, 0.86, f"n={len(data)}", transform=ax.transAxes, ha="right", va="top")
        ax.set_xlim(0, 1)
    axes[-1].set_xlabel("Editing level")
    fig.suptitle("Editing-level distributions by event type", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, dpi=args.plot_dpi)
    plt.close(fig)


def write_summary(path: Path, stats: Counter, class_counts: Counter, group_counts: Counter, filter_counts: Counter, args) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key in [
            "input_rows",
            "parse_error",
            "failed_rows",
            "passed_rows",
        ]:
            writer.writerow([key, stats.get(key, 0)])
        writer.writerow(["event_group_counts", ""])
        for label in ["ADAR", "APOBEC", "Other"]:
            writer.writerow([f"event_group_{label}", group_counts.get(label, 0)])
        writer.writerow(["event_type_counts", ""])
        for label in sorted_event_labels(class_counts):
            writer.writerow([f"event_type_{label}", class_counts.get(label, 0)])
        writer.writerow(["filter_reason_counts", ""])
        for reason, count in filter_counts.most_common():
            writer.writerow([f"filter_{reason}", count])
        writer.writerow(["parameters", ""])
        for key, value in sorted(vars(args).items()):
            writer.writerow([f"param_{key}", value])


def parse_args(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(
        description="Filter REDItools candidate RNA editing sites and classify exact event types with ADAR/APOBEC grouping.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    required = parser.add_argument_group("required inputs")
    required.add_argument("--reditools", required=True, help="REDItools output table. Use '-' for stdin.")
    required.add_argument("--repeatmasker", required=True, help="RepeatMasker .out, BED, or GFF-like repeat intervals.")
    required.add_argument("--genome", required=True, help="Indexed genome FASTA.")
    required.add_argument("--rna-bam", required=True, help="Indexed RNA-seq BAM.")
    required.add_argument("--dna-bam", required=True, help="Indexed DNA-seq BAM.")
    required.add_argument("--gtf", required=True, help="GTF annotation file.")
    required.add_argument("--out-prefix", required=True, help="Output prefix, for example results/sample.")

    parser.add_argument("--repeat-format", choices=["auto", "rmout", "repeatmasker", "bed", "gff"], default="auto", help="Repeat interval format. Use rmout/repeatmasker for standard RepeatMasker .out files.")
    parser.add_argument("--overlap-feature", default="gene", help="GTF feature used for gene-model overlap. Falls back to exon spans if absent.")

    filters = parser.add_argument_group("site filters")
    filters.add_argument("--min-rna-coverage", type=int, default=10)
    filters.add_argument("--min-dna-coverage", type=int, default=10)
    filters.add_argument("--min-rna-meanq", type=float, default=25.0, help="Minimum REDItools RNA MeanQ.")
    filters.add_argument("--min-dna-meanq", type=float, default=25.0, help="Minimum REDItools DNA gMeanQ.")
    filters.add_argument("--use-dna-bam-counts-if-reditools-dna-missing", action="store_true", help="When REDItools genomic DNA fields are '-', count DNA bases directly from the DNA BAM instead of treating DNA evidence as zero coverage.")
    filters.add_argument("--dna-variant-min-baseq", type=int, default=20, help="Minimum DNA BAM base quality used by --use-dna-bam-counts-if-reditools-dna-missing.")
    filters.add_argument("--min-alt-count", type=int, default=3)
    filters.add_argument("--min-edit-freq", type=float, default=0.01)
    filters.add_argument("--max-dna-alt-freq", type=float, default=0.05)
    filters.add_argument("--max-dna-alt-count", type=int, default=3)
    filters.add_argument("--homopolymer-length", type=int, default=5, help="Exclude sites inside same-base runs at least this long. Use 0 to disable.")
    filters.add_argument("--splice-flank", type=int, default=2, help="Exclude positions within this many bases of exon junction boundaries. Use 0 to disable.")
    filters.add_argument("--allow-ref-mismatch", action="store_true", help="Do not fail sites when REDItools reference differs from FASTA.")

    bam_filters = parser.add_argument_group("BAM mapping and indel filters")
    bam_filters.add_argument("--min-mapq", type=int, default=30, help="Per-read MAPQ threshold used to compute good MAPQ fraction.")
    bam_filters.add_argument("--min-mean-mapq", type=float, default=30.0, help="Minimum mean MAPQ at the site.")
    bam_filters.add_argument("--min-good-mapq-fraction", type=float, default=0.80, help="Minimum fraction of usable reads with MAPQ >= min-mapq.")
    bam_filters.add_argument("--indel-flank", type=int, default=5, help="Check indels in this many bases on each side of the site.")
    bam_filters.add_argument("--max-indel-fraction", type=float, default=0.02, help="Maximum allowed indel-read fraction at any flanking pileup column.")
    bam_filters.add_argument("--max-indel-reads", type=int, default=2, help="Maximum allowed indel-supporting reads at any flanking pileup column.")
    bam_filters.add_argument("--max-pileup-depth", type=int, default=100000)
    bam_filters.add_argument("--include-duplicates", action="store_true", help="Include duplicate reads in BAM pileup filters.")
    bam_filters.add_argument("--include-qcfail", action="store_true", help="Include QC-fail reads in BAM pileup filters.")
    bam_filters.add_argument("--skip-bam-mapq-filter", action="store_true", help="Skip mapping-quality checks from BAMs.")
    bam_filters.add_argument("--skip-bam-indel-filter", action="store_true", help="Skip flanking-indel checks from BAMs.")
    bam_filters.add_argument("--mismatch-flank", type=int, default=10, help="Check local reference-discordant mismatch density this many bases on each side of the candidate. The candidate base itself is excluded.")
    bam_filters.add_argument("--mismatch-bams", choices=["rna", "dna", "both"], default="rna", help="Which BAMs to use for the local mismatch-density filter. RNA is usually most useful for IGV red-read artefacts.")
    bam_filters.add_argument("--mismatch-min-baseq", type=int, default=20, help="Minimum base quality used by the local mismatch-density filter.")
    bam_filters.add_argument("--max-window-mismatch-fraction", type=float, default=0.08, help="Maximum total mismatch fraction across the flanking window.")
    bam_filters.add_argument("--max-column-mismatch-fraction", type=float, default=0.30, help="Maximum mismatch fraction allowed at any one flanking reference column.")
    bam_filters.add_argument("--max-high-mismatch-columns", type=int, default=2, help="Maximum number of flanking columns whose mismatch fraction is above --max-column-mismatch-fraction.")
    bam_filters.add_argument("--skip-bam-mismatch-filter", action="store_true", help="Skip local mismatch-density checks for possible mismapped reads.")
    bam_filters.add_argument("--count-same-event-mismatches", action="store_true", help="Count nearby mismatches that match the candidate's transcript-oriented edit type in the local mismatch-density filter. By default these same-type changes are ignored to avoid penalizing real edited clusters.")

    output = parser.add_argument_group("output and plots")
    output.add_argument("--plot-format", choices=["png", "pdf", "svg"], default="png")
    output.add_argument("--plot-dpi", type=int, default=180)
    output.add_argument("--hist-bins", type=int, default=20)
    output.add_argument("--adar-color", default="#3B82F6")
    output.add_argument("--apobec-color", default="#F97316")
    output.add_argument("--other-color", default="#9CA3AF")
    output.add_argument("--write-failures", action="store_true", help="Write a TSV of filtered-out sites and their reasons. This can be large.")
    output.add_argument("--max-sites", type=int, default=0, help="Process at most this many REDItools rows. Use 0 for all rows.")
    output.add_argument("--progress-every", type=int, default=100000, help="Log progress every N rows. Use 0 to disable.")
    output.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")

    args = parser.parse_args(argv)
    if args.homopolymer_length < 0:
        parser.error("--homopolymer-length cannot be negative")
    if args.splice_flank < 0:
        parser.error("--splice-flank cannot be negative")
    if args.hist_bins < 1:
        parser.error("--hist-bins must be positive")
    if not (0 <= args.min_edit_freq <= 1):
        parser.error("--min-edit-freq must be between 0 and 1")
    if not (0 <= args.max_dna_alt_freq <= 1):
        parser.error("--max-dna-alt-freq must be between 0 and 1")
    if args.dna_variant_min_baseq < 0:
        parser.error("--dna-variant-min-baseq must be non-negative")
    if not (0 <= args.max_indel_fraction <= 1):
        parser.error("--max-indel-fraction must be between 0 and 1")
    if not (0 <= args.min_good_mapq_fraction <= 1):
        parser.error("--min-good-mapq-fraction must be between 0 and 1")
    if args.mismatch_flank < 0:
        parser.error("--mismatch-flank cannot be negative")
    if args.mismatch_min_baseq < 0:
        parser.error("--mismatch-min-baseq must be non-negative")
    if not (0 <= args.max_window_mismatch_fraction <= 1):
        parser.error("--max-window-mismatch-fraction must be between 0 and 1")
    if not (0 <= args.max_column_mismatch_fraction <= 1):
        parser.error("--max-column-mismatch-fraction must be between 0 and 1")
    if args.max_high_mismatch_columns < 0:
        parser.error("--max-high-mismatch-columns cannot be negative")
    return args


def run(args) -> None:
    validate_file(args.reditools, "REDItools input", allow_stdin=True)
    validate_file(args.repeatmasker, "RepeatMasker input")
    validate_file(args.genome, "genome FASTA")
    validate_file(args.rna_bam, "RNA BAM")
    validate_file(args.dna_bam, "DNA BAM")
    validate_file(args.gtf, "GTF annotation")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fasta = open_fasta(args.genome)
    fasta_resolver = ContigResolver(fasta.references)
    rna_bam = open_bam(args.rna_bam, "RNA")
    dna_bam = open_bam(args.dna_bam, "DNA")
    bam_contexts = [
        BamContext("rna", rna_bam, ContigResolver(rna_bam.references)),
        BamContext("dna", dna_bam, ContigResolver(dna_bam.references)),
    ]

    logging.info("Loading GTF annotation")
    gtf_data = load_gtf(args.gtf, fasta_resolver, args.overlap_feature, args.splice_flank)
    logging.info(
        "Loaded %s gene intervals and %s splice-site intervals",
        gtf_data.loaded_gene_intervals,
        gtf_data.loaded_splice_intervals,
    )
    logging.info("Loading repeats")
    repeats = load_repeats(args.repeatmasker, fasta_resolver, args.repeat_format)

    passed_path = out_prefix.with_suffix(".filtered_edits.tsv")
    summary_path = out_prefix.with_suffix(".summary.tsv")
    failures_path = out_prefix.with_suffix(".filtered_out_sites.tsv")
    bar_path = out_prefix.with_suffix(f".event_class_counts.{args.plot_format}")
    hist_path = out_prefix.with_suffix(f".editing_level_histograms.{args.plot_format}")

    stats = Counter()
    filter_counts = Counter()
    class_counts = Counter()
    group_counts = Counter()
    edit_levels: Dict[str, List[float]] = defaultdict(list)

    output_fields = [
        "chrom",
        "position_1based",
        "position_0based",
        "reference",
        "genome_reference",
        "alt_base",
        "genomic_edit",
        "transcript_edit",
        "event_type",
        "event_group",
        "event_label",
        "event_class",
        "gene_ids",
        "gene_names",
        "gene_strand_used",
        "rna_coverage_q30",
        "rna_count_depth",
        "rna_alt_count",
        "editing_level",
        "secondary_rna_alts",
        "rna_mean_base_quality",
        "dna_coverage_q30",
        "dna_count_depth",
        "dna_alt_count_same_base",
        "dna_alt_frequency_same_base",
        "dna_mean_base_quality",
        "dna_count_source",
        "homopolymer_length",
        "rna_bam_depth",
        "rna_mean_mapq",
        "rna_good_mapq_fraction",
        "dna_bam_depth",
        "dna_mean_mapq",
        "dna_good_mapq_fraction",
        "rna_max_indel_fraction",
        "rna_max_indel_reads",
        "dna_max_indel_fraction",
        "dna_max_indel_reads",
        "rna_window_mismatch_fraction",
        "rna_max_mismatch_fraction",
        "rna_high_mismatch_columns",
        "rna_mismatch_columns_seen",
        "rna_mismatch_bases",
        "rna_mismatch_window_bases",
        "dna_window_mismatch_fraction",
        "dna_max_mismatch_fraction",
        "dna_high_mismatch_columns",
        "dna_mismatch_columns_seen",
        "dna_mismatch_bases",
        "dna_mismatch_window_bases",
        "reditools_allsubs",
        "reditools_frequency",
        "reditools_gallsubs",
        "reditools_gfrequency",
    ]

    failure_fields = ["line_no", "chrom", "position_1based", "reference", "filter_reasons"]

    with open(passed_path, "w", newline="") as passed_handle:
        passed_writer = csv.DictWriter(passed_handle, delimiter="\t", fieldnames=output_fields, extrasaction="ignore")
        passed_writer.writeheader()

        failure_handle = None
        failure_writer = None
        if args.write_failures:
            failure_handle = open(failures_path, "w", newline="")
            failure_writer = csv.DictWriter(failure_handle, delimiter="\t", fieldnames=failure_fields, extrasaction="ignore")
            failure_writer.writeheader()

        try:
            for line_no, row in parse_reditools(args.reditools):
                stats["input_rows"] += 1
                if args.max_sites and stats["input_rows"] > args.max_sites:
                    break
                if args.progress_every and stats["input_rows"] % args.progress_every == 0:
                    logging.info(
                        "Processed %s rows, passed %s",
                        stats["input_rows"],
                        stats["passed_rows"],
                    )
                evaluation = evaluate_site(
                    row=row,
                    line_no=line_no,
                    fasta=fasta,
                    fasta_resolver=fasta_resolver,
                    repeats=repeats,
                    gtf_data=gtf_data,
                    bam_contexts=bam_contexts,
                    args=args,
                )
                if evaluation.passed and evaluation.record:
                    passed_writer.writerow(evaluation.record)
                    stats["passed_rows"] += 1
                    class_counts[evaluation.event_class or "Other"] += 1
                    group_counts[evaluation.event_group or event_group_from_label(evaluation.event_class or "Other")] += 1
                    if evaluation.editing_level is not None:
                        edit_levels[evaluation.event_class or "Other"].append(evaluation.editing_level)
                else:
                    stats["failed_rows"] += 1
                    if "parse_error" in evaluation.reasons:
                        stats["parse_error"] += 1
                    filter_counts.update(evaluation.reasons)
                    if failure_writer is not None:
                        failure_writer.writerow(
                            {
                                "line_no": line_no,
                                "chrom": row.get("Region", "NA"),
                                "position_1based": row.get("Position", "NA"),
                                "reference": row.get("Reference", "NA"),
                                "filter_reasons": ";".join(evaluation.reasons),
                            }
                        )
        finally:
            if failure_handle is not None:
                failure_handle.close()
            fasta.close()
            rna_bam.close()
            dna_bam.close()

    write_summary(summary_path, stats, class_counts, group_counts, filter_counts, args)
    make_bar_plot(class_counts, bar_path, args)
    make_histogram_plot(edit_levels, hist_path, args)

    logging.info("Wrote filtered edits: %s", passed_path)
    logging.info("Wrote summary: %s", summary_path)
    if args.write_failures:
        logging.info("Wrote filtered-out sites: %s", failures_path)
    logging.info("Wrote plot: %s", bar_path)
    logging.info("Wrote plot: %s", hist_path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        run(args)
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        return 130
    except EditPipelineError as exc:
        logging.error("%s", exc)
        return 2
    except BrokenPipeError:
        return 1
    except Exception as exc:
        logging.exception("Unexpected error: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
