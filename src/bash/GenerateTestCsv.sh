#!/usr/bin/env bash
set -euo pipefail

readonly DEFAULT_REFERENCE='/mnt/experimental/ont/source_data/fc_reference.fa'
readonly DEFAULT_ZONES='/mnt/experimental/ont/source_data/DMSZones.txt'

usage() {
  cat >&2 <<'EOF'
Usage:
  GenerateTestCsv.sh <template_csv> <source_dir> <test_input> <output_csv> [options]

Description:
  Populate a test CSV template while preserving its exact header and row order.
  Each row name is resolved as either:
    - a custom Dorado model under <source_dir>
    - the built-in Dorado "sup" model
    - an existing FASTQ under <source_dir>

Required positional arguments:
  <template_csv>   Input CSV template with "name" as the first column
  <source_dir>     Root directory to search recursively for models / FASTQs
  <test_input>     POD5 file or directory used for model/sup rows
  <output_csv>     Output CSV path

Options:
  --reference PATH          Reference FASTA/MMI for minimap2.
                            Default: /mnt/experimental/ont/source_data/fc_reference.fa
  --zones PATH              DMS zones file. Default:
                            /mnt/experimental/ont/source_data/DMSZones.txt
  --work-dir PATH           Directory for cached row artifacts. Default:
                            <output_csv without .csv>.work
  --kit-name NAME           Dorado kit name. Default: SQK-NBD114-24
  --sup-model MODEL         Built-in Dorado model spec for the "sup" row.
                            Default: sup
  --threads N               Worker threads for fastplong/minimap2/samtools.
                            Default: 16
  --min-length N            Minimum read length. Default: 600
  --max-length N            Maximum read length. Default: 800
  --max-indel-events N      DualSiteDMSFilter max indel events. Default: 3
  --max-indel-bases N       DualSiteDMSFilter max indel bases. Default: 3
  --mut-codon-library LIB   DualSiteDMSFilter mutation codon library.
                            Default: NNK
  --filter-bin PATH         Explicit path to DualSiteDMSFilter
  --force                   Re-run rows even if cached artifacts exist
  -h, --help                Show this help message
EOF
}

die() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

require_command() {
  local command_name=$1
  command -v "$command_name" >/dev/null 2>&1 || die "Required command not found: $command_name"
}

require_file() {
  local path=$1
  local label=${2:-File}
  [[ -f "$path" ]] || die "$label does not exist: $path"
}

require_directory() {
  local path=$1
  local label=${2:-Directory}
  [[ -d "$path" ]] || die "$label does not exist: $path"
}

require_path() {
  local path=$1
  local label=${2:-Path}
  [[ -e "$path" ]] || die "$label does not exist: $path"
}

sanitize_filename_component() {
  local value=$1
  value=$(printf '%s' "$value" | sed 's/[^A-Za-z0-9._-]/_/g')
  [[ -n "$value" ]] || value='row'
  printf '%s\n' "$value"
}

dataset_stem_from_name() {
  local dataset_name=$1
  local stem=${dataset_name##*/}
  stem=${stem%.fastq.gz}
  stem=${stem%.fq.gz}
  stem=${stem%.fastq}
  stem=${stem%.fq}
  stem=${stem%.pod5}
  stem=${stem%.bam}
  printf '%s\n' "$stem"
}

fastq_marker_path() {
  local fastq_path=$1
  printf '%s.ok\n' "$fastq_path"
}

mark_fastq_valid() {
  local fastq_path=$1
  : >"$(fastq_marker_path "$fastq_path")"
}

clear_fastq_state() {
  local fastq_path=$1
  rm -f -- "$fastq_path" "$(fastq_marker_path "$fastq_path")"
}

fastq_read_count() {
  local fastq_path=$1

  python3 - "$fastq_path" <<'PY'
import gzip
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    raise SystemExit(f"FASTQ file does not exist: {path}")


def open_fastq(file_path: Path):
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, "rt", encoding="utf-8", newline="")
    return open(file_path, "rt", encoding="utf-8", newline="")


with open_fastq(path) as handle:
    record_count = 0
    while True:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record_count += 1

        if not sequence or not plus or not quality:
            raise SystemExit(f"Incomplete FASTQ record in {path} at record {record_count}")
        if not header.startswith("@"):
            raise SystemExit(f"FASTQ header does not start with @ in {path} at record {record_count}")
        if not plus.startswith("+"):
            raise SystemExit(f"FASTQ separator does not start with + in {path} at record {record_count}")

        sequence = sequence.rstrip("\r\n")
        quality = quality.rstrip("\r\n")
        if len(sequence) != len(quality):
            raise SystemExit(
                f"Sequence/quality length mismatch in {path} at record {record_count}: "
                f"{len(sequence)} != {len(quality)}"
            )

print(record_count)
PY
}

validate_fastq_file() {
  local fastq_path=$1
  fastq_read_count "$fastq_path" >/dev/null
}

bam_record_count() {
  local bam_path=$1
  samtools view -c "$bam_path"
}

alignment_fastq_matches_bam() {
  local bam_path=$1
  local fastq_path=$2
  local bam_records
  local fastq_records

  bam_records=$(bam_record_count "$bam_path")
  fastq_records=$(fastq_read_count "$fastq_path")

  if [[ "$bam_records" != "$fastq_records" ]]; then
    log "Detected BAM/FASTQ inconsistency: $bam_path has $bam_records records, but $fastq_path has $fastq_records reads"
    return 1
  fi

  return 0
}

reuse_fastq_if_valid() {
  local fastq_path=$1
  local force=$2
  local label=$3
  local marker_path
  marker_path=$(fastq_marker_path "$fastq_path")

  if [[ -f "$fastq_path" && "$force" -eq 0 ]]; then
    if [[ -f "$marker_path" ]] || validate_fastq_file "$fastq_path"; then
      mark_fastq_valid "$fastq_path"
      log "Reusing existing $label: $fastq_path"
      return 0
    fi
    log "Existing $label is invalid and will be regenerated: $fastq_path"
  fi

  return 1
}

validate_bam_file() {
  local bam_path=$1
  samtools quickcheck "$bam_path"
}

read_zone_string() {
  local zone_file=$1
  awk '
    /^[[:space:]]*#/ { next }
    /^[[:space:]]*$/ { next }
    { print; exit }
  ' "$zone_file"
}

resolve_filter_binary() {
  local requested_path=$1
  local script_directory
  local project_root

  script_directory=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  project_root=$(cd "$script_directory/../.." && pwd)

  if [[ -n "$requested_path" ]]; then
    [[ -x "$requested_path" ]] || die "DualSiteDMSFilter is not executable: $requested_path"
    printf '%s\n' "$requested_path"
    return 0
  fi

  if command -v DualSiteDMSFilter >/dev/null 2>&1; then
    command -v DualSiteDMSFilter
    return 0
  fi

  if [[ -x "$project_root/build/DualSiteDMSFilter" ]]; then
    printf '%s\n' "$project_root/build/DualSiteDMSFilter"
    return 0
  fi

  die "DualSiteDMSFilter was not found on PATH or under $project_root/build"
}

json_zero_fastq_metrics() {
  local prefix=$1
  python3 - "$prefix" <<'PY'
import json
import sys

prefix = sys.argv[1]
payload = {
    f"{prefix}reads": "0",
    f"{prefix}bases": "0",
    f"{prefix}minimum_length": "0",
    f"{prefix}maximum_length": "0",
    f"{prefix}median_length": "0",
    f"{prefix}mean_length": "0.000000",
    f"{prefix}n50_length": "0",
    f"{prefix}gc_content": "0.000000",
    f"{prefix}q5_bases": "0",
    f"{prefix}q7_bases": "0",
    f"{prefix}q10_bases": "0",
    f"{prefix}q15_bases": "0",
    f"{prefix}q20_bases": "0",
    f"{prefix}q30_bases": "0",
    f"{prefix}q40_bases": "0",
}
json.dump(payload, sys.stdout, sort_keys=True)
PY
}

compute_fastq_metrics_json() {
  local fastq_path=$1
  local prefix=$2

  python3 - "$fastq_path" "$prefix" <<'PY'
import gzip
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
prefix = sys.argv[2]
thresholds = (5, 7, 10, 15, 20, 30, 40)


def open_fastq(file_path: Path):
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, "rt", encoding="utf-8", newline="")
    return open(file_path, "rt", encoding="utf-8", newline="")


def fmt_decimal(value: float) -> str:
    return f"{value:.6f}"


def fmt_number(value) -> str:
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return str(value)


def zero_payload() -> dict[str, str]:
    return {
        f"{prefix}reads": "0",
        f"{prefix}bases": "0",
        f"{prefix}minimum_length": "0",
        f"{prefix}maximum_length": "0",
        f"{prefix}median_length": "0",
        f"{prefix}mean_length": "0.000000",
        f"{prefix}n50_length": "0",
        f"{prefix}gc_content": "0.000000",
        f"{prefix}q5_bases": "0",
        f"{prefix}q7_bases": "0",
        f"{prefix}q10_bases": "0",
        f"{prefix}q15_bases": "0",
        f"{prefix}q20_bases": "0",
        f"{prefix}q30_bases": "0",
        f"{prefix}q40_bases": "0",
    }


if not path.exists():
    raise SystemExit(f"FASTQ file does not exist: {path}")

lengths = []
total_reads = 0
total_bases = 0
gc_bases = 0
quality_counts = {threshold: 0 for threshold in thresholds}

with open_fastq(path) as handle:
    while True:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()

        if not sequence or not plus or not quality:
            raise SystemExit(f"Incomplete FASTQ record in {path}")

        sequence = sequence.rstrip("\r\n")
        quality = quality.rstrip("\r\n")
        if len(sequence) != len(quality):
            raise SystemExit(
                f"Sequence/quality length mismatch in {path}: "
                f"{len(sequence)} != {len(quality)}"
            )

        sequence_upper = sequence.upper()
        read_length = len(sequence_upper)
        lengths.append(read_length)
        total_reads += 1
        total_bases += read_length
        gc_bases += sequence_upper.count("G") + sequence_upper.count("C")

        for character in quality:
            score = ord(character) - 33
            if score >= 40:
                quality_counts[40] += 1
                quality_counts[30] += 1
                quality_counts[20] += 1
                quality_counts[15] += 1
                quality_counts[10] += 1
                quality_counts[7] += 1
                quality_counts[5] += 1
            elif score >= 30:
                quality_counts[30] += 1
                quality_counts[20] += 1
                quality_counts[15] += 1
                quality_counts[10] += 1
                quality_counts[7] += 1
                quality_counts[5] += 1
            elif score >= 20:
                quality_counts[20] += 1
                quality_counts[15] += 1
                quality_counts[10] += 1
                quality_counts[7] += 1
                quality_counts[5] += 1
            elif score >= 15:
                quality_counts[15] += 1
                quality_counts[10] += 1
                quality_counts[7] += 1
                quality_counts[5] += 1
            elif score >= 10:
                quality_counts[10] += 1
                quality_counts[7] += 1
                quality_counts[5] += 1
            elif score >= 7:
                quality_counts[7] += 1
                quality_counts[5] += 1
            elif score >= 5:
                quality_counts[5] += 1

if total_reads == 0 or total_bases == 0:
    json.dump(zero_payload(), sys.stdout, sort_keys=True)
    raise SystemExit(0)

lengths.sort()
minimum_length = lengths[0]
maximum_length = lengths[-1]
middle = total_reads // 2
if total_reads % 2 == 1:
    median_length = float(lengths[middle])
else:
    median_length = (lengths[middle - 1] + lengths[middle]) / 2.0

target_bases = total_bases / 2.0
cumulative_bases = 0
n50_length = 0
for read_length in reversed(lengths):
    cumulative_bases += read_length
    if cumulative_bases >= target_bases:
        n50_length = read_length
        break

mean_length = total_bases / total_reads
gc_content = gc_bases / total_bases

payload = {
    f"{prefix}reads": str(total_reads),
    f"{prefix}bases": str(total_bases),
    f"{prefix}minimum_length": str(minimum_length),
    f"{prefix}maximum_length": str(maximum_length),
    f"{prefix}median_length": fmt_number(median_length),
    f"{prefix}mean_length": fmt_decimal(mean_length),
    f"{prefix}n50_length": str(n50_length),
    f"{prefix}gc_content": fmt_decimal(gc_content),
    f"{prefix}q5_bases": str(quality_counts[5]),
    f"{prefix}q7_bases": str(quality_counts[7]),
    f"{prefix}q10_bases": str(quality_counts[10]),
    f"{prefix}q15_bases": str(quality_counts[15]),
    f"{prefix}q20_bases": str(quality_counts[20]),
    f"{prefix}q30_bases": str(quality_counts[30]),
    f"{prefix}q40_bases": str(quality_counts[40]),
}

json.dump(payload, sys.stdout, sort_keys=True)
PY
}

json_zero_alignment_metrics() {
  local prefix=$1
  python3 - "$prefix" <<'PY'
import json
import sys

prefix = sys.argv[1]
payload = {
    f"{prefix}inserted_bases": "0",
    f"{prefix}deleted_bases": "0",
    f"{prefix}mismatches": "0",
    f"{prefix}bases_mapped": "0",
    f"{prefix}error_rate": "0.00000000",
    f"{prefix}overall_error_rate": "0.00000000",
}
json.dump(payload, sys.stdout, sort_keys=True)
PY
}

compute_alignment_metrics_json() {
  local bam_path=$1
  local prefix=$2
  local filter_primary=$3

  python3 - "$bam_path" "$prefix" "$filter_primary" <<'PY'
import json
import subprocess
import sys
import tempfile
from pathlib import Path

bam_path = Path(sys.argv[1])
prefix = sys.argv[2]
filter_primary = sys.argv[3] == "1"


def zero_payload() -> dict[str, str]:
    return {
        f"{prefix}inserted_bases": "0",
        f"{prefix}deleted_bases": "0",
        f"{prefix}mismatches": "0",
        f"{prefix}bases_mapped": "0",
        f"{prefix}error_rate": "0.00000000",
        f"{prefix}overall_error_rate": "0.00000000",
    }


def extract_sn_metric(stats_path: Path, metric_name: str) -> str | None:
    with stats_path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] == "SN" and parts[1] == metric_name:
                return parts[2]
    return None


def count_indel_bases(path: Path) -> tuple[int, int]:
    view = subprocess.Popen(
        ["samtools", "view", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    inserted = 0
    deleted = 0
    try:
        assert view.stdout is not None
        for line in view.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            cigar = fields[5]
            number = ""
            for ch in cigar:
                if ch.isdigit():
                    number += ch
                    continue
                if not number:
                    continue
                length = int(number)
                if ch == "I":
                    inserted += length
                elif ch == "D":
                    deleted += length
                number = ""
    finally:
        if view.stdout is not None:
            view.stdout.close()
    stderr = view.stderr.read() if view.stderr is not None else ""
    if view.stderr is not None:
        view.stderr.close()
    exit_code = view.wait()
    if exit_code != 0:
        raise SystemExit(stderr.strip() or f"samtools view failed for {path}")
    return inserted, deleted


if not bam_path.exists():
    raise SystemExit(f"BAM file does not exist: {bam_path}")

working_bam = bam_path
temporary_dir = None

try:
    if filter_primary:
        temporary_dir = tempfile.TemporaryDirectory()
        working_bam = Path(temporary_dir.name) / "primary_only.bam"
        subprocess.run(
            ["samtools", "view", "-F", "0x904", "-b", "-o", str(working_bam), str(bam_path)],
            check=True,
            stdout=subprocess.DEVNULL,
        )

    stats_path = Path(tempfile.mkstemp(prefix="alignment_stats_", suffix=".txt")[1])
    try:
        with stats_path.open("w") as handle:
            subprocess.run(["samtools", "stats", str(working_bam)], check=True, stdout=handle)

        bases_mapped = extract_sn_metric(stats_path, "bases mapped (cigar):") or "0"
        mismatches = extract_sn_metric(stats_path, "mismatches:") or "0"
        inserted_bases, deleted_bases = count_indel_bases(working_bam)

        bases_mapped_value = int(float(bases_mapped))
        mismatches_value = int(float(mismatches))
        if bases_mapped_value == 0:
            error_rate = "0.00000000"
            overall_error_rate = "0.00000000"
        else:
            error_rate = f"{mismatches_value / bases_mapped_value:.8f}"
            overall_error_rate = (
                f"{(mismatches_value + inserted_bases + deleted_bases) / bases_mapped_value:.8f}"
            )

        payload = {
            f"{prefix}inserted_bases": str(inserted_bases),
            f"{prefix}deleted_bases": str(deleted_bases),
            f"{prefix}mismatches": str(mismatches_value),
            f"{prefix}bases_mapped": str(bases_mapped_value),
            f"{prefix}error_rate": error_rate,
            f"{prefix}overall_error_rate": overall_error_rate,
        }
        json.dump(payload, sys.stdout, sort_keys=True)
    finally:
        stats_path.unlink(missing_ok=True)
finally:
    if temporary_dir is not None:
        temporary_dir.cleanup()
PY
}

json_zero_dms_metrics() {
  python3 - <<'PY'
import json
import sys

payload = {
    "failed_mapq": "0",
    "failed_indel": "0",
    "failed_coverage": "0",
    "failed_dual_site": "0",
    "failed_mut_codon_library": "0",
    "pass_reads_0_codon_mutations": "0",
    "pass_reads_1_codon_mutations": "0",
    "pass_reads_2_codon_mutations": "0",
    "pass_rate": "0",
    "unique_variants": "0",
    "avg_dms_nucleotide_mismatches_per_read": "0",
    "avg_framework_nucleotide_mismatches_per_read": "0",
    "avg_net_nucleotide_mismatches_per_read": "0",
    "avg_dms_codon_mismatches_per_read": "0",
    "avg_framework_codon_mismatches_per_read": "0",
    "avg_net_codon_mismatches_per_read": "0",
    "dms_bases_compared": "0",
    "framework_bases_compared": "0",
    "dms_codons_compared": "0",
    "framework_codons_compared": "0",
    "dms_nt_mismatch_rate_per_base": "0",
    "framework_nt_mismatch_rate_per_base": "0",
    "net_nt_mismatch_rate_per_base": "0",
    "dms_codon_mismatch_rate_per_codon": "0",
    "framework_codon_mismatch_rate_per_codon": "0",
    "net_codon_mismatch_rate_per_codon": "0",
}
json.dump(payload, sys.stdout, sort_keys=True)
PY
}

parse_dms_filter_log_json() {
  local log_path=$1

  python3 - "$log_path" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    raise SystemExit(f"DMS log does not exist: {path}")

key_map = {
    "fail_mapq": "failed_mapq",
    "fail_indel": "failed_indel",
    "fail_coverage": "failed_coverage",
    "fail_dual_site": "failed_dual_site",
    "fail_mut_codon_library": "failed_mut_codon_library",
    "pass_reads_0_codon_mutations": "pass_reads_0_codon_mutations",
    "pass_reads_1_codon_mutations": "pass_reads_1_codon_mutations",
    "pass_reads_2_codon_mutations": "pass_reads_2_codon_mutations",
    "pass_rate": "pass_rate",
    "unique_variants": "unique_variants",
    "avg_dms_nucleotide_mismatches_per_read": "avg_dms_nucleotide_mismatches_per_read",
    "avg_framework_nucleotide_mismatches_per_read": "avg_framework_nucleotide_mismatches_per_read",
    "avg_net_nucleotide_mismatches_per_read": "avg_net_nucleotide_mismatches_per_read",
    "avg_dms_codon_mismatches_per_read": "avg_dms_codon_mismatches_per_read",
    "avg_framework_codon_mismatches_per_read": "avg_framework_codon_mismatches_per_read",
    "avg_net_codon_mismatches_per_read": "avg_net_codon_mismatches_per_read",
    "dms_bases_compared": "dms_bases_compared",
    "framework_bases_compared": "framework_bases_compared",
    "dms_codons_compared": "dms_codons_compared",
    "framework_codons_compared": "framework_codons_compared",
    "dms_nt_mismatch_rate_per_base": "dms_nt_mismatch_rate_per_base",
    "framework_nt_mismatch_rate_per_base": "framework_nt_mismatch_rate_per_base",
    "net_nt_mismatch_rate_per_base": "net_nt_mismatch_rate_per_base",
    "dms_codon_mismatch_rate_per_codon": "dms_codon_mismatch_rate_per_codon",
    "framework_codon_mismatch_rate_per_codon": "framework_codon_mismatch_rate_per_codon",
    "net_codon_mismatch_rate_per_codon": "net_codon_mismatch_rate_per_codon",
}

payload = {
    "failed_mapq": "0",
    "failed_indel": "0",
    "failed_coverage": "0",
    "failed_dual_site": "0",
    "failed_mut_codon_library": "0",
    "pass_reads_0_codon_mutations": "0",
    "pass_reads_1_codon_mutations": "0",
    "pass_reads_2_codon_mutations": "0",
    "pass_rate": "0",
    "unique_variants": "0",
    "avg_dms_nucleotide_mismatches_per_read": "0",
    "avg_framework_nucleotide_mismatches_per_read": "0",
    "avg_net_nucleotide_mismatches_per_read": "0",
    "avg_dms_codon_mismatches_per_read": "0",
    "avg_framework_codon_mismatches_per_read": "0",
    "avg_net_codon_mismatches_per_read": "0",
    "dms_bases_compared": "0",
    "framework_bases_compared": "0",
    "dms_codons_compared": "0",
    "framework_codons_compared": "0",
    "dms_nt_mismatch_rate_per_base": "0",
    "framework_nt_mismatch_rate_per_base": "0",
    "net_nt_mismatch_rate_per_base": "0",
    "dms_codon_mismatch_rate_per_codon": "0",
    "framework_codon_mismatch_rate_per_codon": "0",
    "net_codon_mismatch_rate_per_codon": "0",
}

in_block = False
found_block = False

with path.open() as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line == "[DmsFilterCount]":
            in_block = True
            found_block = True
            continue
        if not in_block:
            continue
        if "=" not in line:
            if line.startswith("[") and line.endswith("]"):
                break
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key in key_map:
            payload[key_map[key]] = value

if not found_block:
    raise SystemExit(f"DmsFilterCount block not found in {path}")

json.dump(payload, sys.stdout, sort_keys=True)
PY
}

merge_json_files() {
  python3 - "$@" <<'PY'
import json
import sys

merged = {}
for path in sys.argv[1:]:
    with open(path) as handle:
        payload = json.load(handle)
    merged.update(payload)

pre_length_filter_reads = int(merged.get("pre_length_filter_reads", "0") or 0)
post_dms_filter_reads = int(merged.get("post_dms_filter_reads", "0") or 0)
if pre_length_filter_reads == 0:
    merged["pipeline_pass_rate"] = "0.0000000000"
else:
    merged["pipeline_pass_rate"] = f"{post_dms_filter_reads / pre_length_filter_reads:.10f}"

json.dump(merged, sys.stdout, sort_keys=True)
PY
}

emit_template_row_metadata() {
  local template_csv=$1
  python3 - "$template_csv" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open(newline="") as handle:
    reader = csv.reader(handle)
    try:
        header = next(reader)
    except StopIteration:
        raise SystemExit(f"Template CSV is empty: {path}")

    if not header or header[0] != "name":
        raise SystemExit(f"Template CSV must have 'name' as the first column: {path}")

    expected_fields = len(header)
    for index, row in enumerate(reader, start=1):
        if len(row) < expected_fields:
            row = row + [""] * (expected_fields - len(row))
        name = row[0] if row else ""
        complete = (
            bool(name)
            and len(row) >= expected_fields
            and all(cell != "" for cell in row[1:expected_fields])
        )
        print(f"{index}\t{name}\t{1 if complete else 0}")
PY
}

render_output_csv() {
  local template_csv=$1
  local results_dir=$2
  local output_csv=$3

  python3 - "$template_csv" "$results_dir" "$output_csv" <<'PY'
import csv
import json
import sys
from pathlib import Path

template_csv = Path(sys.argv[1])
results_dir = Path(sys.argv[2])
output_csv = Path(sys.argv[3])

with template_csv.open(newline="") as in_handle, output_csv.open("w", newline="") as out_handle:
    reader = csv.reader(in_handle)
    writer = csv.writer(out_handle, lineterminator="\n")

    try:
        header = next(reader)
    except StopIteration:
        raise SystemExit(f"Template CSV is empty: {template_csv}")

    template_header = list(header)
    if "pipeline_pass_rate" not in header:
        header = header + ["pipeline_pass_rate"]

    writer.writerow(header)

    for index, row in enumerate(reader, start=1):
        if len(row) < len(template_header):
            row = row + [""] * (len(template_header) - len(row))
        elif len(row) > len(template_header):
            row = row[: len(template_header)]

        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))

        metrics_path = results_dir / f"{index}.json"
        if metrics_path.exists():
            with metrics_path.open() as handle:
                metrics = json.load(handle)
            for column_index, column_name in enumerate(header):
                if column_name in metrics:
                    value = metrics[column_name]
                    row[column_index] = "" if value is None else str(value)
        writer.writerow(row)
PY
}

resolve_model_dir_from_wrapper() {
  local wrapper_dir=$1
  local row_name=$2
  local preferred_child="$wrapper_dir/${row_name}_v2"
  local exact_child="$wrapper_dir/$row_name"
  local candidates=()
  local candidate

  if [[ -f "$wrapper_dir/config.toml" ]]; then
    printf '%s\n' "$wrapper_dir"
    return 0
  fi

  if [[ -d "$preferred_child" && -f "$preferred_child/config.toml" ]]; then
    printf '%s\n' "$preferred_child"
    return 0
  fi

  if [[ -d "$exact_child" && -f "$exact_child/config.toml" ]]; then
    printf '%s\n' "$exact_child"
    return 0
  fi

  while IFS= read -r candidate; do
    candidates+=("$candidate")
  done < <(find "$wrapper_dir" -mindepth 1 -maxdepth 1 -type d -name '*' -exec test -f '{}/config.toml' ';' -print)

  if (( ${#candidates[@]} == 1 )); then
    printf '%s\n' "${candidates[0]}"
    return 0
  fi

  if (( ${#candidates[@]} == 0 )); then
    return 1
  fi

  printf 'Error: multiple Dorado-ready child directories found for %s:\n' "$row_name" >&2
  printf '  %s\n' "${candidates[@]}" >&2
  return 2
}

resolve_row_source() {
  local row_name=$1
  local source_dir=$2
  local row_stem
  local matches=()
  local exact_fastq_matches=()
  local exact_model_matches=()
  local fastq_matches=()
  local match_path
  local candidate
  local resolved_model
  local resolve_status

  if [[ "$row_name" == "sup" ]]; then
    printf 'sup\t%s\n' "$SUP_MODEL"
    return 0
  fi

  while IFS= read -r match_path; do
    matches+=("$match_path")
  done < <(find "$source_dir" \( -type f -o -type d \) -name "$row_name" -print)

  if (( ${#matches[@]} > 0 )); then
    for match_path in "${matches[@]}"; do
      case "$match_path" in
        *.fastq|*.fq|*.fastq.gz|*.fq.gz)
          exact_fastq_matches+=("$match_path")
          ;;
        *)
          if [[ -d "$match_path" ]]; then
            if resolved_model=$(resolve_model_dir_from_wrapper "$match_path" "$row_name"); then
              exact_model_matches+=("$resolved_model")
            else
              resolve_status=$?
              if (( resolve_status > 1 )); then
                return "$resolve_status"
              fi
            fi
          fi
          ;;
      esac
    done
  fi

  if (( ${#exact_model_matches[@]} == 1 && ${#exact_fastq_matches[@]} == 0 )); then
    printf 'model\t%s\n' "${exact_model_matches[0]}"
    return 0
  fi

  if (( ${#exact_fastq_matches[@]} == 1 && ${#exact_model_matches[@]} == 0 )); then
    printf 'fastq\t%s\n' "${exact_fastq_matches[0]}"
    return 0
  fi

  if (( ${#exact_model_matches[@]} + ${#exact_fastq_matches[@]} > 1 )); then
    printf 'Error: multiple exact resolvable matches found for %s under %s:\n' "$row_name" "$source_dir" >&2
    printf '  %s\n' "${exact_model_matches[@]}" "${exact_fastq_matches[@]}" >&2
    return 1
  fi

  if (( ${#matches[@]} == 1 && ${#exact_model_matches[@]} == 0 && ${#exact_fastq_matches[@]} == 0 )); then
    match_path=${matches[0]}
    case "$match_path" in
      *)
        if [[ -d "$match_path" ]]; then
          die "Exact match exists for row '$row_name' but is not a Dorado-ready model directory: $match_path"
        fi
        die "Unsupported exact match for row '$row_name': $match_path"
        ;;
    esac
  fi

  if (( ${#matches[@]} > 1 && ${#exact_model_matches[@]} == 0 && ${#exact_fastq_matches[@]} == 0 )); then
    printf 'Error: multiple exact matches found for %s under %s:\n' "$row_name" "$source_dir" >&2
    printf '  %s\n' "${matches[@]}" >&2
    return 1
  fi

  row_stem=$(dataset_stem_from_name "$row_name")

  while IFS= read -r candidate; do
    if [[ "$(dataset_stem_from_name "${candidate##*/}")" == "$row_stem" ]]; then
      fastq_matches+=("$candidate")
    fi
  done < <(find "$source_dir" -type f \( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) -print)

  if (( ${#fastq_matches[@]} == 1 )); then
    printf 'fastq\t%s\n' "${fastq_matches[0]}"
    return 0
  fi

  if (( ${#fastq_matches[@]} > 1 )); then
    printf 'Error: multiple FASTQ stem matches found for %s under %s:\n' "$row_name" "$source_dir" >&2
    printf '  %s\n' "${fastq_matches[@]}" >&2
    return 1
  fi

  die "Could not resolve row '$row_name' under $source_dir"
}

run_basecalling_if_needed() {
  local model_spec=$1
  local test_input=$2
  local output_fastq=$3
  local kit_name=$4
  local force=$5

  if reuse_fastq_if_valid "$output_fastq" "$force" "basecalled FASTQ"; then
    return 0
  fi

  mkdir -p "$(dirname "$output_fastq")"

  local dorado_args=(
    basecaller
    --emit-fastq
    --trim all
  )
  if [[ -n "$kit_name" ]]; then
    dorado_args+=(--kit-name "$kit_name")
  fi
  if [[ -d "$test_input" ]]; then
    dorado_args+=(--recursive)
  fi
  dorado_args+=("$model_spec" "$test_input")

  local attempt
  for attempt in 1 2; do
    clear_fastq_state "$output_fastq"
    log "Basecalling $test_input with $model_spec -> $output_fastq (attempt $attempt)"
    if dorado "${dorado_args[@]}" >"$output_fastq"; then
      if validate_fastq_file "$output_fastq"; then
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: dorado produced an invalid FASTQ on attempt $attempt: $output_fastq"
    else
      log "Warning: dorado failed on attempt $attempt"
    fi
  done

  return 1
}

run_length_filter_if_needed() {
  local input_fastq=$1
  local output_fastq=$2
  local report_html=$3
  local report_json=$4
  local min_length=$5
  local max_length=$6
  local threads=$7
  local force=$8

  if reuse_fastq_if_valid "$output_fastq" "$force" "length-filtered FASTQ"; then
    return 0
  fi

  mkdir -p "$(dirname "$output_fastq")"

  local attempt
  for attempt in 1 2; do
    clear_fastq_state "$output_fastq"
    rm -f -- "$report_html" "$report_json"
    log "Length filtering $input_fastq -> $output_fastq (attempt $attempt)"
    if fastplong \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --length_required "$min_length" \
      --length_limit "$max_length" \
      --thread "$threads" \
      --in "$input_fastq" \
      --out "$output_fastq" \
      --html "$report_html" \
      --json "$report_json"
    then
      if validate_fastq_file "$output_fastq"; then
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: length filtering produced an invalid FASTQ on attempt $attempt: $output_fastq"
    else
      log "Warning: length filtering failed on attempt $attempt"
    fi
  done

  return 1
}

run_primary_alignment_if_needed() {
  local input_fastq=$1
  local reference_path=$2
  local aligned_bam=$3
  local primary_bam=$4
  local aligned_fastq=$5
  local threads=$6
  local force=$7

  if [[ -f "$primary_bam" && -f "$aligned_bam" ]] && reuse_fastq_if_valid "$aligned_fastq" "$force" "post-alignment FASTQ"; then
    if validate_bam_file "$primary_bam" && validate_bam_file "$aligned_bam"; then
      if alignment_fastq_matches_bam "$primary_bam" "$aligned_fastq"; then
        log "Reusing existing primary alignment BAMs: $primary_bam"
        return 0
      fi
      log "Existing post-alignment FASTQ is inconsistent and will be regenerated: $aligned_fastq"
    fi
  fi

  if [[ ! -s "$input_fastq" ]]; then
    clear_fastq_state "$aligned_fastq"
    rm -f -- "$aligned_bam" "$primary_bam"
    : >"$aligned_fastq"
    mark_fastq_valid "$aligned_fastq"
    log "Skipping alignment because the input FASTQ is empty: $input_fastq"
    return 0
  fi

  mkdir -p "$(dirname "$aligned_fastq")"
  clear_fastq_state "$aligned_fastq"
  rm -f -- "$aligned_bam" "$primary_bam"

  log "Aligning $input_fastq -> $primary_bam"
  minimap2 -t "$threads" -ax map-ont "$reference_path" "$input_fastq" \
    | samtools view -@ "$threads" -u -o "$aligned_bam" -

  samtools view -@ "$threads" -F 0x904 -u "$aligned_bam" \
    | samtools sort -@ "$threads" -n -o "$primary_bam" -

  samtools fastq "$primary_bam" >"$aligned_fastq"
  validate_bam_file "$aligned_bam"
  validate_bam_file "$primary_bam"
  validate_fastq_file "$aligned_fastq"
  alignment_fastq_matches_bam "$primary_bam" "$aligned_fastq" \
    || die "Generated post-alignment FASTQ does not match primary BAM: $aligned_fastq"
  mark_fastq_valid "$aligned_fastq"
}

run_dms_filter_if_needed() {
  local primary_bam=$1
  local reference_path=$2
  local polished_fastq=$3
  local variant_tsv=$4
  local dms_log=$5
  local zone_string=$6
  local max_indel_events=$7
  local max_indel_bases=$8
  local mut_codon_library=$9
  local filter_binary=${10}
  local force=${11}

  if [[ -f "$variant_tsv" && -f "$dms_log" ]] && reuse_fastq_if_valid "$polished_fastq" "$force" "DMS-filtered FASTQ"; then
    log "Reusing existing DMS artifacts: $polished_fastq"
    return 0
  fi

  mkdir -p "$(dirname "$polished_fastq")"
  clear_fastq_state "$polished_fastq"
  rm -f -- "$variant_tsv" "$dms_log"

  log "Running DualSiteDMSFilter on $primary_bam"
  if "$filter_binary" \
    "$primary_bam" \
    "$reference_path" \
    "$polished_fastq" \
    "$zone_string" \
    --max-indel-events "$max_indel_events" \
    --max-indel-bases "$max_indel_bases" \
    --mut-codon-library "$mut_codon_library" \
    --out-counts "$variant_tsv" \
    2>&1 | tee "$dms_log"
  then
    validate_fastq_file "$polished_fastq"
    mark_fastq_valid "$polished_fastq"
    return 0
  fi

  return 1
}

run_polished_alignment_if_needed() {
  local polished_fastq=$1
  local reference_path=$2
  local polished_bam=$3
  local threads=$4
  local force=$5

  if [[ -f "$polished_bam" && "$force" -eq 0 ]]; then
    if validate_bam_file "$polished_bam"; then
      log "Reusing existing polished alignment BAM: $polished_bam"
      return 0
    fi
  fi

  rm -f -- "$polished_bam"

  if [[ ! -s "$polished_fastq" ]]; then
    log "Skipping polished alignment because the polished FASTQ is empty: $polished_fastq"
    return 0
  fi

  log "Aligning polished FASTQ $polished_fastq -> $polished_bam"
  minimap2 -t "$threads" -ax map-ont "$reference_path" "$polished_fastq" \
    | samtools view -@ "$threads" -u - \
    | samtools sort -@ "$threads" -o "$polished_bam" -

  validate_bam_file "$polished_bam"
}

build_row_metrics_json() {
  local row_name=$1
  local base_fastq=$2
  local length_filtered_fastq=$3
  local aligned_fastq=$4
  local primary_bam=$5
  local polished_fastq=$6
  local polished_bam=$7
  local dms_log=$8
  local row_dir=$9
  local output_json=${10}

  local pre_json="$row_dir/pre_length.json"
  local post_length_json="$row_dir/post_length.json"
  local post_alignment_json="$row_dir/post_alignment.json"
  local primary_json="$row_dir/primary_alignment.json"
  local post_dms_json="$row_dir/post_dms.json"
  local dms_json="$row_dir/dms.json"
  local polished_json="$row_dir/polished_alignment.json"

  if [[ -f "$primary_bam" && -f "$aligned_fastq" ]]; then
    alignment_fastq_matches_bam "$primary_bam" "$aligned_fastq" \
      || die "Post-alignment FASTQ is inconsistent with primary BAM: $aligned_fastq"
  fi

  compute_fastq_metrics_json "$base_fastq" "pre_length_filter_" >"$pre_json"
  compute_fastq_metrics_json "$length_filtered_fastq" "post_length_filter_" >"$post_length_json"

  if [[ -f "$aligned_fastq" ]]; then
    compute_fastq_metrics_json "$aligned_fastq" "post_alignment_filter_" >"$post_alignment_json"
  else
    json_zero_fastq_metrics "post_alignment_filter_" >"$post_alignment_json"
  fi

  if [[ -f "$primary_bam" ]]; then
    compute_alignment_metrics_json "$primary_bam" "primary_" 0 >"$primary_json"
  else
    json_zero_alignment_metrics "primary_" >"$primary_json"
  fi

  if [[ -f "$polished_fastq" ]]; then
    compute_fastq_metrics_json "$polished_fastq" "post_dms_filter_" >"$post_dms_json"
  else
    json_zero_fastq_metrics "post_dms_filter_" >"$post_dms_json"
  fi

  if [[ -f "$dms_log" ]]; then
    parse_dms_filter_log_json "$dms_log" >"$dms_json"
  else
    json_zero_dms_metrics >"$dms_json"
  fi

  if [[ -f "$polished_bam" ]]; then
    compute_alignment_metrics_json "$polished_bam" "polished_" 1 >"$polished_json"
  else
    json_zero_alignment_metrics "polished_" >"$polished_json"
  fi

  merge_json_files \
    "$pre_json" \
    "$post_length_json" \
    "$post_alignment_json" \
    "$primary_json" \
    "$post_dms_json" \
    "$dms_json" \
    "$polished_json" \
    >"$output_json"

  python3 - "$output_json" "$row_name" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
row_name = sys.argv[2]
payload = json.loads(path.read_text())
payload["name"] = row_name
path.write_text(json.dumps(payload, sort_keys=True))
PY
}

process_row() {
  local row_index=$1
  local row_name=$2
  local metrics_json_path=$3
  local row_id
  local row_dir
  local resolution_file
  local source_type
  local resolved_path
  local base_fastq
  local length_filtered_fastq
  local aligned_bam
  local primary_bam
  local aligned_fastq
  local polished_fastq
  local variant_tsv
  local dms_log
  local polished_bam
  local merged_metrics_tmp

  if [[ -s "$metrics_json_path" && "$FORCE" -eq 0 ]]; then
    log "Reusing cached metrics for row $row_name"
    return 0
  fi

  row_id=$(printf '%03d_%s' "$row_index" "$(sanitize_filename_component "$row_name")")
  row_dir="$WORK_ROOT/rows/$row_id"
  resolution_file="$row_dir/source_resolution.tsv"
  mkdir -p "$row_dir"

  if ! resolve_row_source "$row_name" "$SOURCE_DIR" >"$resolution_file"; then
    return 1
  fi

  IFS=$'\t' read -r source_type resolved_path <"$resolution_file"
  log "Resolved row $row_name as $source_type: $resolved_path"

  case "$source_type" in
    model|sup)
      base_fastq="$row_dir/${row_id}.basecalled.fq"
      run_basecalling_if_needed "$resolved_path" "$TEST_INPUT" "$base_fastq" "$KIT_NAME" "$FORCE" || return 1
      ;;
    fastq)
      base_fastq="$resolved_path"
      ;;
    *)
      die "Unsupported resolved row type for $row_name: $source_type"
      ;;
  esac

  length_filtered_fastq="$row_dir/${row_id}.len${MIN_LENGTH}_${MAX_LENGTH}.fq"
  aligned_bam="$row_dir/${row_id}.aligned.bam"
  primary_bam="$row_dir/${row_id}.primary.sorted.bam"
  aligned_fastq="$row_dir/${row_id}.aligned.fq"
  polished_fastq="$row_dir/${row_id}.polished.fq"
  variant_tsv="$row_dir/${row_id}.variants.tsv"
  dms_log="$row_dir/${row_id}.dms.log"
  polished_bam="$row_dir/${row_id}.polished.sorted.bam"

  run_length_filter_if_needed \
    "$base_fastq" \
    "$length_filtered_fastq" \
    "$row_dir/length_filter.html" \
    "$row_dir/length_filter.json" \
    "$MIN_LENGTH" \
    "$MAX_LENGTH" \
    "$THREADS" \
    "$FORCE" \
    || return 1

  if [[ -s "$length_filtered_fastq" ]]; then
    run_primary_alignment_if_needed \
      "$length_filtered_fastq" \
      "$REFERENCE_PATH" \
      "$aligned_bam" \
      "$primary_bam" \
      "$aligned_fastq" \
      "$THREADS" \
      "$FORCE" \
      || return 1
  else
    clear_fastq_state "$aligned_fastq"
    rm -f -- "$aligned_bam" "$primary_bam" "$polished_fastq" "$variant_tsv" "$dms_log" "$polished_bam"
    : >"$aligned_fastq"
    mark_fastq_valid "$aligned_fastq"
  fi

  if [[ -f "$primary_bam" && -s "$aligned_fastq" ]]; then
    run_dms_filter_if_needed \
      "$primary_bam" \
      "$REFERENCE_PATH" \
      "$polished_fastq" \
      "$variant_tsv" \
      "$dms_log" \
      "$ZONE_STRING" \
      "$MAX_INDEL_EVENTS" \
      "$MAX_INDEL_BASES" \
      "$MUT_CODON_LIBRARY" \
      "$FILTER_BINARY" \
      "$FORCE" \
      || return 1
  else
    clear_fastq_state "$polished_fastq"
    rm -f -- "$variant_tsv" "$dms_log" "$polished_bam"
    : >"$polished_fastq"
    mark_fastq_valid "$polished_fastq"
  fi

  if [[ -s "$polished_fastq" ]]; then
    run_polished_alignment_if_needed \
      "$polished_fastq" \
      "$REFERENCE_PATH" \
      "$polished_bam" \
      "$THREADS" \
      "$FORCE" \
      || return 1
  else
    rm -f -- "$polished_bam"
  fi

  merged_metrics_tmp=$(mktemp "$row_dir/metrics.XXXXXX.json")
  build_row_metrics_json \
    "$row_name" \
    "$base_fastq" \
    "$length_filtered_fastq" \
    "$aligned_fastq" \
    "$primary_bam" \
    "$polished_fastq" \
    "$polished_bam" \
    "$dms_log" \
    "$row_dir" \
    "$merged_metrics_tmp" \
    || {
      rm -f -- "$merged_metrics_tmp"
      return 1
    }

  mv -f -- "$merged_metrics_tmp" "$metrics_json_path"
}

main() {
  local template_csv=
  local source_dir=
  local test_input=
  local output_csv=
  local reference_path=$DEFAULT_REFERENCE
  local zones_path=$DEFAULT_ZONES
  local work_root=
  local kit_name='SQK-NBD114-24'
  local sup_model='sup'
  local threads=16
  local min_length=600
  local max_length=800
  local max_indel_events=3
  local max_indel_bases=3
  local mut_codon_library='NNK'
  local filter_binary_override=
  local force=0

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -h|--help)
        usage
        exit 0
        ;;
      --reference)
        [[ $# -ge 2 ]] || die "--reference requires a value"
        reference_path=$2
        shift 2
        ;;
      --zones)
        [[ $# -ge 2 ]] || die "--zones requires a value"
        zones_path=$2
        shift 2
        ;;
      --work-dir)
        [[ $# -ge 2 ]] || die "--work-dir requires a value"
        work_root=$2
        shift 2
        ;;
      --kit-name)
        [[ $# -ge 2 ]] || die "--kit-name requires a value"
        kit_name=$2
        shift 2
        ;;
      --sup-model)
        [[ $# -ge 2 ]] || die "--sup-model requires a value"
        sup_model=$2
        shift 2
        ;;
      --threads)
        [[ $# -ge 2 ]] || die "--threads requires a value"
        threads=$2
        shift 2
        ;;
      --min-length)
        [[ $# -ge 2 ]] || die "--min-length requires a value"
        min_length=$2
        shift 2
        ;;
      --max-length)
        [[ $# -ge 2 ]] || die "--max-length requires a value"
        max_length=$2
        shift 2
        ;;
      --max-indel-events)
        [[ $# -ge 2 ]] || die "--max-indel-events requires a value"
        max_indel_events=$2
        shift 2
        ;;
      --max-indel-bases)
        [[ $# -ge 2 ]] || die "--max-indel-bases requires a value"
        max_indel_bases=$2
        shift 2
        ;;
      --mut-codon-library)
        [[ $# -ge 2 ]] || die "--mut-codon-library requires a value"
        mut_codon_library=$2
        shift 2
        ;;
      --filter-bin)
        [[ $# -ge 2 ]] || die "--filter-bin requires a value"
        filter_binary_override=$2
        shift 2
        ;;
      --force)
        force=1
        shift
        ;;
      --*)
        die "Unknown option: $1"
        ;;
      *)
        if [[ -z "$template_csv" ]]; then
          template_csv=$1
        elif [[ -z "$source_dir" ]]; then
          source_dir=$1
        elif [[ -z "$test_input" ]]; then
          test_input=$1
        elif [[ -z "$output_csv" ]]; then
          output_csv=$1
        else
          die "Unexpected extra argument: $1"
        fi
        shift
        ;;
    esac
  done

  if [[ -z "$template_csv" || -z "$source_dir" || -z "$test_input" || -z "$output_csv" ]]; then
    usage
    exit 1
  fi

  require_file "$template_csv" "Template CSV"
  require_directory "$source_dir" "Source directory"
  require_path "$test_input" "Test input"
  require_file "$reference_path" "Reference"
  require_file "$zones_path" "Zones file"

  require_command dorado
  require_command fastplong
  require_command minimap2
  require_command python3
  require_command samtools
  require_command tee
  require_command awk

  [[ "$threads" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer"
  [[ "$min_length" =~ ^[0-9]+$ ]] || die "--min-length must be a non-negative integer"
  [[ "$max_length" =~ ^[0-9]+$ ]] || die "--max-length must be a non-negative integer"
  [[ "$max_indel_events" =~ ^[0-9]+$ ]] || die "--max-indel-events must be a non-negative integer"
  [[ "$max_indel_bases" =~ ^[0-9]+$ ]] || die "--max-indel-bases must be a non-negative integer"
  (( threads > 0 )) || die "--threads must be greater than zero"
  (( max_length >= min_length )) || die "--max-length must be >= --min-length"

  if [[ -z "$work_root" ]]; then
    work_root=${output_csv%.csv}
    if [[ "$work_root" == "$output_csv" ]]; then
      work_root="${output_csv}.work"
    else
      work_root="${work_root}.work"
    fi
  fi

  mkdir -p "$work_root" "$(dirname "$output_csv")"

  TEMPLATE_CSV=$template_csv
  SOURCE_DIR=$source_dir
  TEST_INPUT=$test_input
  OUTPUT_CSV=$output_csv
  REFERENCE_PATH=$reference_path
  ZONES_PATH=$zones_path
  WORK_ROOT=$work_root
  KIT_NAME=$kit_name
  SUP_MODEL=$sup_model
  THREADS=$threads
  MIN_LENGTH=$min_length
  MAX_LENGTH=$max_length
  MAX_INDEL_EVENTS=$max_indel_events
  MAX_INDEL_BASES=$max_indel_bases
  MUT_CODON_LIBRARY=$mut_codon_library
  FORCE=$force
  FILTER_BINARY=$(resolve_filter_binary "$filter_binary_override")
  ZONE_STRING=$(read_zone_string "$ZONES_PATH")
  [[ -n "$ZONE_STRING" ]] || die "Zone file does not contain a usable zone string: $ZONES_PATH"

  local results_dir="$WORK_ROOT/completed_rows"
  local failed_rows_file="${OUTPUT_CSV}.failed.txt"
  local output_tmp
  local row_index
  local row_name
  local row_complete
  local failure_count=0

  mkdir -p "$results_dir"
  : >"$failed_rows_file"

  log "Template CSV      = $TEMPLATE_CSV"
  log "Source directory  = $SOURCE_DIR"
  log "Test input        = $TEST_INPUT"
  log "Output CSV        = $OUTPUT_CSV"
  log "Work directory    = $WORK_ROOT"
  log "Reference         = $REFERENCE_PATH"
  log "Zones file        = $ZONES_PATH"
  log "Zone string       = $ZONE_STRING"
  log "Filter binary     = $FILTER_BINARY"

  while IFS=$'\t' read -r row_index row_name row_complete; do
    if [[ -z "$row_name" ]]; then
      rm -f -- "$results_dir/${row_index}.json"
      log "Skipping blank row at template index $row_index"
      continue
    fi

    if [[ "$row_complete" == "1" ]]; then
      rm -f -- "$results_dir/${row_index}.json"
      log "Leaving already-complete template row unchanged: $row_name"
      continue
    fi

    if process_row "$row_index" "$row_name" "$results_dir/${row_index}.json"; then
      log "Completed row: $row_name"
    else
      rm -f -- "$results_dir/${row_index}.json"
      failure_count=$((failure_count + 1))
      printf '%s\n' "$row_name" >>"$failed_rows_file"
      log "Warning: row failed and will be left unchanged: $row_name"
    fi
  done < <(emit_template_row_metadata "$TEMPLATE_CSV")

  output_tmp=$(mktemp "$WORK_ROOT/output.XXXXXX.csv")
  render_output_csv "$TEMPLATE_CSV" "$results_dir" "$output_tmp"
  mv -f -- "$output_tmp" "$OUTPUT_CSV"

  log "Completed CSV written to $OUTPUT_CSV"
  if (( failure_count > 0 )); then
    log "Warning: $failure_count row(s) failed. See $failed_rows_file"
  fi
}

main "$@"
