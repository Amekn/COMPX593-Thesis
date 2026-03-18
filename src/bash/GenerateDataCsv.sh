#!/usr/bin/env bash
set -euo pipefail

readonly CSV_HEADER='name,pre_length_filter_reads,pre_length_filter_bases,pre_length_filter_minimum_length,pre_length_filter_maximum_length,pre_length_filter_median_length,pre_length_filter_mean_length,pre_length_filter_n50_length,pre_length_filter_gc_content,pre_length_filter_q5_bases,pre_length_filter_q7_bases,pre_length_filter_q10_bases,pre_length_filter_q15_bases,pre_length_filter_q20_bases,pre_length_filter_q30_bases,pre_length_filter_q40_bases,post_length_filter_reads,post_length_filter_bases,post_length_filter_minimum_length,post_length_filter_maximum_length,post_length_filter_median_length,post_length_filter_mean_length,post_length_filter_n50_length,post_length_filter_gc_content,post_length_filter_q5_bases,post_length_filter_q7_bases,post_length_filter_q10_bases,post_length_filter_q15_bases,post_length_filter_q20_bases,post_length_filter_q30_bases,post_length_filter_q40_bases,post_quality_filter_reads,post_quality_filter_bases,post_quality_filter_minimum_length,post_quality_filter_maximum_length,post_quality_filter_median_length,post_quality_filter_mean_length,post_quality_filter_n50_length,post_quality_filter_gc_content,post_quality_filter_q5_bases,post_quality_filter_q7_bases,post_quality_filter_q10_bases,post_quality_filter_q15_bases,post_quality_filter_q20_bases,post_quality_filter_q30_bases,post_quality_filter_q40_bases,post_alignment_filter_reads,post_alignment_filter_bases,post_alignment_filter_minimum_length,post_alignment_filter_maximum_length,post_alignment_filter_median_length,post_alignment_filter_mean_length,post_alignment_filter_n50_length,post_alignment_filter_gc_content,post_alignment_filter_q5_bases,post_alignment_filter_q7_bases,post_alignment_filter_q10_bases,post_alignment_filter_q15_bases,post_alignment_filter_q20_bases,post_alignment_filter_q30_bases,post_alignment_filter_q40_bases,inserted_bases,deleted_bases,mismatches,bases_mapped,error_rate,overall_error_rate'

usage() {
  cat >&2 <<'EOF'
Usage:
  GenerateDataCsv.sh <source_dir> <names_file> <output_csv> [options]

Description:
  Run the documented basecalling/filtering/alignment workflow for the datasets
  listed in <names_file> and write a completed CSV matching table/data.csv.

  <names_file> can be either:
    - a plain text file with one dataset name per line, or
    - a CSV whose first column is "name" (for example table/data.csv)

Options:
  --reference PATH       Reference FASTA or minimap2 index. Defaults to
                         <source_dir>/fc_reference.mmi or <source_dir>/fc_reference.fa.
  --work-dir PATH        Directory for intermediate files. Defaults to
                         <source_dir>/.table_pipeline
  --kit-name NAME        Dorado kit name for POD5 inputs. Default:
                         SQK-NBD114-24
  --dorado-model MODEL   Dorado basecaller model. Default: sup@v5.2.0
  --threads N            Thread count for fastplong/minimap2/samtools.
                         Default: 16
  --min-length N         Minimum read length for length filtering. Default: 600
  --max-length N         Maximum read length for length filtering. Default: 800
  --min-mean-q N         Minimum mean quality for quality filtering. Default: 15
  --force                Re-run stages even if intermediate files already exist
  -h, --help             Show this help message
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

trim_cr() {
  local value=$1
  printf '%s' "${value%$'\r'}"
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

default_reference_path() {
  local source_dir=$1
  local candidate
  for candidate in \
    "$source_dir/fc_reference.mmi" \
    "$source_dir/fc_reference.fa"
  do
    if [[ -f "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  local found_reference
  found_reference=$(find "$source_dir" -type f \( -name 'fc_reference.mmi' -o -name 'fc_reference.fa' \) -print | head -n 1 || true)
  [[ -n "$found_reference" ]] || return 1
  printf '%s\n' "$found_reference"
}

extract_dataset_names() {
  local names_file=$1
  local first_line
  if ! IFS= read -r first_line < "$names_file"; then
    return 0
  fi
  first_line=$(trim_cr "$first_line")

  if [[ "$first_line" == name || "$first_line" == name,* ]]; then
    awk -F',' '
      NR == 1 { next }
      {
        sub(/\r$/, "", $1)
        if (length($1) > 0) {
          print $1
        }
      }
    ' "$names_file"
  else
    awk '
      {
        sub(/\r$/, "", $0)
        if (length($0) > 0) {
          print $0
        }
      }
    ' "$names_file"
  fi
}

resolve_source_path() {
  local source_dir=$1
  local dataset_name=$2

  if [[ -e "$dataset_name" ]]; then
    printf '%s\n' "$dataset_name"
    return 0
  fi

  if [[ -e "$source_dir/$dataset_name" ]]; then
    printf '%s\n' "$source_dir/$dataset_name"
    return 0
  fi

  local matches=()
  while IFS= read -r match_path; do
    matches+=("$match_path")
  done < <(find "$source_dir" \( -type f -o -type d \) -name "$dataset_name" -print)

  if (( ${#matches[@]} == 1 )); then
    printf '%s\n' "${matches[0]}"
    return 0
  fi
  if (( ${#matches[@]} == 0 )); then
    return 1
  fi

  printf 'Error: multiple matches found for %s:\n' "$dataset_name" >&2
  printf '  %s\n' "${matches[@]}" >&2
  exit 1
}

compute_fastq_metrics_csv() {
  local fastq_path=$1

  python3 - "$fastq_path" <<'PY'
import gzip
import sys
from pathlib import Path

path = Path(sys.argv[1])
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
    zero_values = [
        "0", "0", "0", "0", "0", "0.000000", "0", "0.000000",
        "0", "0", "0", "0", "0", "0", "0",
    ]
    print(",".join(zero_values))
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

values = [
    total_reads,
    total_bases,
    minimum_length,
    maximum_length,
    fmt_number(median_length),
    fmt_decimal(mean_length),
    n50_length,
    fmt_decimal(gc_content),
    quality_counts[5],
    quality_counts[7],
    quality_counts[10],
    quality_counts[15],
    quality_counts[20],
    quality_counts[30],
    quality_counts[40],
]

print(",".join(str(value) for value in values))
PY
}

extract_samtools_sn_metric() {
  local stats_path=$1
  local metric_name=$2
  awk -F '\t' -v key="$metric_name" '$1 == "SN" && $2 == key { print $3; exit }' "$stats_path"
}

count_primary_bam_indel_bases() {
  local bam_path=$1

  samtools view -F 0x904 "$bam_path" \
    | awk '
        {
          cigar = $6
          while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
            token = substr(cigar, RSTART, RLENGTH)
            length_value = substr(token, 1, length(token) - 1) + 0
            operation = substr(token, length(token), 1)
            if (operation == "I") {
              inserted_bases += length_value
            } else if (operation == "D") {
              deleted_bases += length_value
            }
            cigar = substr(cigar, RSTART + RLENGTH)
          }
        }
        END {
          print inserted_bases + 0 "," deleted_bases + 0
        }
      '
}

compute_alignment_metrics_csv() {
  local bam_path=$1
  local temporary_stats
  local bases_mapped
  local mismatches
  local inserted_bases
  local deleted_bases
  local error_rate
  local overall_error_rate

  require_file "$bam_path" "Primary alignment BAM"

  temporary_stats=$(mktemp)
  samtools stats "$bam_path" >"$temporary_stats"

  bases_mapped=$(extract_samtools_sn_metric "$temporary_stats" "bases mapped (cigar):")
  mismatches=$(extract_samtools_sn_metric "$temporary_stats" "mismatches:")
  IFS=',' read -r inserted_bases deleted_bases < <(count_primary_bam_indel_bases "$bam_path")
  rm -f -- "$temporary_stats"

  bases_mapped=${bases_mapped:-0}
  mismatches=${mismatches:-0}
  inserted_bases=${inserted_bases:-0}
  deleted_bases=${deleted_bases:-0}

  read -r error_rate overall_error_rate < <(
    awk \
      -v bases_mapped="$bases_mapped" \
      -v mismatches="$mismatches" \
      -v inserted_bases="$inserted_bases" \
      -v deleted_bases="$deleted_bases" '
        BEGIN {
          if (bases_mapped == 0) {
            error_rate = 0
            overall_error_rate = 0
          } else {
            error_rate = mismatches / bases_mapped
            overall_error_rate = (mismatches + inserted_bases + deleted_bases) / bases_mapped
          }
          printf "%.8f %.8f\n", error_rate, overall_error_rate
        }
      '
  )

  printf '%s,%s,%s,%s,%s,%s\n' \
    "$inserted_bases" \
    "$deleted_bases" \
    "$mismatches" \
    "$bases_mapped" \
    "$error_rate" \
    "$overall_error_rate"
}

validate_fastq_file() {
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
PY
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

reuse_stage_fastq_if_valid() {
  local input_fastq=$1
  local output_fastq=$2
  local force=$3
  local label=$4
  shift 4

  if [[ -f "$output_fastq" && "$force" -eq 0 ]]; then
    if [[ "$input_fastq" -nt "$output_fastq" ]]; then
      log "Existing $label is older than its input and will be regenerated: $output_fastq"
      return 1
    fi

    local required_artifact
    for required_artifact in "$@"; do
      if [[ ! -e "$required_artifact" ]]; then
        log "Existing $label is missing a required artifact and will be regenerated: $required_artifact"
        return 1
      fi
      if [[ "$input_fastq" -nt "$required_artifact" || "$output_fastq" -nt "$required_artifact" ]]; then
        log "Existing $label has an outdated dependent artifact and will be regenerated: $required_artifact"
        return 1
      fi
    done
  fi

  reuse_fastq_if_valid "$output_fastq" "$force" "$label"
}

prepare_output_csv_for_resume() {
  local output_csv=$1
  local existing_header

  if [[ -s "$output_csv" ]]; then
    IFS= read -r existing_header < "$output_csv"
    existing_header=$(trim_cr "$existing_header")
    [[ "$existing_header" == "$CSV_HEADER" ]] || die "Existing output CSV header does not match the current schema: $output_csv"
  else
    printf '%s\n' "$CSV_HEADER" >"$output_csv"
  fi
}

extract_completed_output_names() {
  local output_csv=$1

  [[ -s "$output_csv" ]] || return 0

  awk -F',' '
    NR == 1 {
      expected_fields = NF
      next
    }
    {
      sub(/\r$/, "", $0)
      if ($1 == "") {
        next
      }
      if (NF != expected_fields) {
        next
      }
      complete = 1
      for (i = 2; i <= NF; ++i) {
        if ($i == "") {
          complete = 0
          break
        }
      }
      if (complete) {
        print $1
      }
    }
  ' "$output_csv"
}

run_basecalling_if_needed() {
  local source_path=$1
  local output_fastq=$2
  local kit_name=$3
  local dorado_model=$4
  local force=$5

  if reuse_fastq_if_valid "$output_fastq" "$force" "basecalled FASTQ"; then
    return 0
  fi

  require_command dorado
  mkdir -p "$(dirname "$output_fastq")"

  local dorado_args=(
    basecaller
    --emit-fastq
    --trim all
  )
  if [[ -n "$kit_name" ]]; then
    dorado_args+=(--kit-name "$kit_name")
  fi
  if [[ -d "$source_path" ]]; then
    dorado_args+=(--recursive)
  fi
  dorado_args+=("$dorado_model" "$source_path")

  local attempt exit_code
  for attempt in 1 2; do
    clear_fastq_state "$output_fastq"
    log "Basecalling $source_path -> $output_fastq (attempt $attempt)"

    if dorado "${dorado_args[@]}" >"$output_fastq"; then
      if validate_fastq_file "$output_fastq"; then
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: dorado produced an invalid FASTQ on attempt $attempt: $output_fastq"
    else
      exit_code=$?
      if [[ -f "$output_fastq" ]] && validate_fastq_file "$output_fastq"; then
        log "Warning: dorado exited with status $exit_code but produced a valid FASTQ; continuing."
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: dorado failed with status $exit_code on attempt $attempt"
    fi
  done

  return 1
}

run_bam_to_fastq_if_needed() {
  local source_bam=$1
  local output_fastq=$2
  local force=$3

  if reuse_fastq_if_valid "$output_fastq" "$force" "BAM-to-FASTQ conversion"; then
    return 0
  fi

  require_command samtools
  mkdir -p "$(dirname "$output_fastq")"

  local attempt exit_code
  for attempt in 1 2; do
    clear_fastq_state "$output_fastq"
    log "Converting BAM to FASTQ: $source_bam -> $output_fastq (attempt $attempt)"

    if samtools fastq "$source_bam" >"$output_fastq"; then
      if validate_fastq_file "$output_fastq"; then
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: BAM-to-FASTQ produced an invalid FASTQ on attempt $attempt: $output_fastq"
    else
      exit_code=$?
      log "Warning: BAM-to-FASTQ failed with status $exit_code on attempt $attempt"
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

  if reuse_stage_fastq_if_valid \
    "$input_fastq" \
    "$output_fastq" \
    "$force" \
    "length-filtered FASTQ" \
    "$report_html" \
    "$report_json"
  then
    return 0
  fi

  local attempt exit_code
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
      exit_code=$?
      if [[ -f "$output_fastq" ]] && validate_fastq_file "$output_fastq"; then
        log "Warning: length filtering exited with status $exit_code but produced a valid FASTQ; continuing."
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: length filtering failed with status $exit_code on attempt $attempt"
    fi
  done

  return 1
}

run_quality_filter_if_needed() {
  local input_fastq=$1
  local output_fastq=$2
  local report_html=$3
  local report_json=$4
  local min_mean_q=$5
  local threads=$6
  local force=$7

  if reuse_stage_fastq_if_valid \
    "$input_fastq" \
    "$output_fastq" \
    "$force" \
    "quality-filtered FASTQ" \
    "$report_html" \
    "$report_json"
  then
    return 0
  fi

  local attempt exit_code
  for attempt in 1 2; do
    clear_fastq_state "$output_fastq"
    rm -f -- "$report_html" "$report_json"
    log "Quality filtering $input_fastq -> $output_fastq (attempt $attempt)"

    if fastplong \
      --disable_adapter_trimming \
      --disable_length_filtering \
      --mean_qual "$min_mean_q" \
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
      log "Warning: quality filtering produced an invalid FASTQ on attempt $attempt: $output_fastq"
    else
      exit_code=$?
      if [[ -f "$output_fastq" ]] && validate_fastq_file "$output_fastq"; then
        log "Warning: quality filtering exited with status $exit_code but produced a valid FASTQ; continuing."
        mark_fastq_valid "$output_fastq"
        return 0
      fi
      log "Warning: quality filtering failed with status $exit_code on attempt $attempt"
    fi
  done

  return 1
}

run_alignment_if_needed() {
  local input_fastq=$1
  local reference_path=$2
  local aligned_bam=$3
  local primary_bam=$4
  local aligned_fastq=$5
  local report_html=$6
  local report_json=$7
  local threads=$8
  local force=$9

  if [[ -f "$primary_bam" ]] && reuse_stage_fastq_if_valid \
    "$input_fastq" \
    "$aligned_fastq" \
    "$force" \
    "post-alignment FASTQ" \
    "$aligned_bam" \
    "$primary_bam"
  then
    return 0
  fi

  if [[ -f "$aligned_fastq" && "$force" -eq 0 && ! -f "$primary_bam" ]]; then
    log "Existing post-alignment FASTQ found, but the primary BAM is missing; regenerating alignment artifacts: $primary_bam"
  fi

  if [[ ! -s "$input_fastq" ]]; then
    log "Skipping alignment because the quality-filtered FASTQ is empty: $input_fastq"
    clear_fastq_state "$aligned_fastq"
    : >"$aligned_fastq"
    mark_fastq_valid "$aligned_fastq"
    rm -f -- "$aligned_bam" "$primary_bam" "$report_html" "$report_json"
    return 0
  fi

  local attempt exit_code
  for attempt in 1 2; do
    clear_fastq_state "$aligned_fastq"
    rm -f -- "$aligned_bam" "$primary_bam" "$report_html" "$report_json"
    log "Aligning $input_fastq -> $aligned_fastq (attempt $attempt)"

    if minimap2 -t "$threads" -ax map-ont "$reference_path" "$input_fastq" \
      | samtools view -@ "$threads" -u -o "$aligned_bam" -
    then
      :
    else
      exit_code=$?
      log "Warning: minimap2/samtools view failed with status $exit_code on attempt $attempt"
      continue
    fi

    if samtools view -@ "$threads" -F 0x904 -u "$aligned_bam" \
      | samtools sort -@ "$threads" -n -o "$primary_bam" -
    then
      :
    else
      exit_code=$?
      log "Warning: primary-alignment extraction failed with status $exit_code on attempt $attempt"
      continue
    fi

    if samtools fastq "$primary_bam" >"$aligned_fastq"; then
      if validate_fastq_file "$aligned_fastq"; then
        mark_fastq_valid "$aligned_fastq"
        break
      fi
      log "Warning: alignment produced an invalid FASTQ on attempt $attempt: $aligned_fastq"
    else
      exit_code=$?
      log "Warning: samtools fastq failed with status $exit_code on attempt $attempt"
    fi
  done

  if [[ ! -f "$(fastq_marker_path "$aligned_fastq")" ]]; then
    return 1
  fi

  if ! fastplong \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --thread "$threads" \
    --in "$aligned_fastq" \
    --out /dev/null \
    --html "$report_html" \
    --json "$report_json"
  then
    log "Warning: alignment report generation failed for $aligned_fastq; continuing without 03_alignment_filter.{html,json}"
    rm -f -- "$report_html" "$report_json"
  fi

  return 0
}

process_dataset() {
  local source_dir=$1
  local dataset_name=$2
  local work_root=$3
  local reference_path=$4
  local kit_name=$5
  local dorado_model=$6
  local threads=$7
  local min_length=$8
  local max_length=$9
  local min_mean_q=${10}
  local force=${11}
  local output_csv=${12}

  log "Processing dataset: $dataset_name"

  local source_path
  if ! source_path=$(resolve_source_path "$source_dir" "$dataset_name"); then
    log "Warning: could not locate dataset under $source_dir: $dataset_name"
    return 1
  fi

  local dataset_stem
  dataset_stem=$(dataset_stem_from_name "$dataset_name")
  local dataset_work_dir="$work_root/$dataset_stem"
  mkdir -p "$dataset_work_dir"

  local base_fastq=
  case "$source_path" in
    *.pod5)
      base_fastq="$dataset_work_dir/${dataset_stem}.basecalled.fq"
      run_basecalling_if_needed "$source_path" "$base_fastq" "$kit_name" "$dorado_model" "$force" || return 1
      ;;
    *.fastq|*.fq|*.fastq.gz|*.fq.gz)
      base_fastq="$source_path"
      ;;
    *.bam)
      base_fastq="$dataset_work_dir/${dataset_stem}.source.fq"
      run_bam_to_fastq_if_needed "$source_path" "$base_fastq" "$force" || return 1
      ;;
    *)
      if [[ -d "$source_path" ]]; then
        base_fastq="$dataset_work_dir/${dataset_stem}.basecalled.fq"
        run_basecalling_if_needed "$source_path" "$base_fastq" "$kit_name" "$dorado_model" "$force" || return 1
      else
        log "Warning: unsupported dataset type for $dataset_name: $source_path"
        return 1
      fi
      ;;
  esac

  local length_filtered_fastq="$dataset_work_dir/${dataset_stem}.len${min_length}_${max_length}.fq"
  local quality_filtered_fastq="$dataset_work_dir/${dataset_stem}.q${min_mean_q}.fq"
  local aligned_bam="$dataset_work_dir/${dataset_stem}.aligned.bam"
  local primary_bam="$dataset_work_dir/${dataset_stem}.primary.sorted.bam"
  local aligned_fastq="$dataset_work_dir/${dataset_stem}.aligned.fq"

  run_length_filter_if_needed \
    "$base_fastq" \
    "$length_filtered_fastq" \
    "$dataset_work_dir/01_length_filter.html" \
    "$dataset_work_dir/01_length_filter.json" \
    "$min_length" \
    "$max_length" \
    "$threads" \
    "$force" \
    || return 1

  run_quality_filter_if_needed \
    "$length_filtered_fastq" \
    "$quality_filtered_fastq" \
    "$dataset_work_dir/02_quality_filter.html" \
    "$dataset_work_dir/02_quality_filter.json" \
    "$min_mean_q" \
    "$threads" \
    "$force" \
    || return 1

  run_alignment_if_needed \
    "$quality_filtered_fastq" \
    "$reference_path" \
    "$aligned_bam" \
    "$primary_bam" \
    "$aligned_fastq" \
    "$dataset_work_dir/03_alignment_filter.html" \
    "$dataset_work_dir/03_alignment_filter.json" \
    "$threads" \
    "$force" \
    || return 1

  local pre_length_metrics
  local post_length_metrics
  local post_quality_metrics
  local post_alignment_metrics
  local alignment_metrics

  pre_length_metrics=$(compute_fastq_metrics_csv "$base_fastq")
  post_length_metrics=$(compute_fastq_metrics_csv "$length_filtered_fastq")
  post_quality_metrics=$(compute_fastq_metrics_csv "$quality_filtered_fastq")
  post_alignment_metrics=$(compute_fastq_metrics_csv "$aligned_fastq")
  alignment_metrics=$(compute_alignment_metrics_csv "$primary_bam")

  printf '%s,%s,%s,%s,%s,%s\n' \
    "$dataset_name" \
    "$pre_length_metrics" \
    "$post_length_metrics" \
    "$post_quality_metrics" \
    "$post_alignment_metrics" \
    "$alignment_metrics" \
    >>"$output_csv"
}

main() {
  local source_dir=
  local names_file=
  local output_csv=
  local reference_path=
  local work_root=
  local kit_name='SQK-NBD114-24'
  local dorado_model='sup@v5.2.0'
  local threads=16
  local min_length=600
  local max_length=800
  local min_mean_q=15
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
      --dorado-model)
        [[ $# -ge 2 ]] || die "--dorado-model requires a value"
        dorado_model=$2
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
      --min-mean-q)
        [[ $# -ge 2 ]] || die "--min-mean-q requires a value"
        min_mean_q=$2
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
        if [[ -z "$source_dir" ]]; then
          source_dir=$1
        elif [[ -z "$names_file" ]]; then
          names_file=$1
        elif [[ -z "$output_csv" ]]; then
          output_csv=$1
        else
          die "Unexpected extra argument: $1"
        fi
        shift
        ;;
    esac
  done

  if [[ -z "$source_dir" || -z "$names_file" || -z "$output_csv" ]]; then
    usage
    exit 1
  fi

  require_directory "$source_dir" "Source directory"
  require_file "$names_file" "Dataset list file"
  require_command awk
  require_command fastplong
  require_command mktemp
  require_command minimap2
  require_command python3
  require_command samtools

  [[ "$threads" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer"
  [[ "$min_length" =~ ^[0-9]+$ ]] || die "--min-length must be a non-negative integer"
  [[ "$max_length" =~ ^[0-9]+$ ]] || die "--max-length must be a non-negative integer"
  [[ "$min_mean_q" =~ ^[0-9]+$ ]] || die "--min-mean-q must be a non-negative integer"
  (( threads > 0 )) || die "--threads must be greater than zero"
  (( max_length == 0 || max_length >= min_length )) || die "--max-length must be zero or >= --min-length"

  if [[ -z "$reference_path" ]]; then
    reference_path=$(default_reference_path "$source_dir" || true)
    [[ -n "$reference_path" ]] || die "Could not find fc_reference.mmi or fc_reference.fa under $source_dir. Use --reference."
  fi
  require_file "$reference_path" "Reference"

  if [[ -z "$work_root" ]]; then
    work_root="$source_dir/.table_pipeline"
  fi
  mkdir -p "$work_root"
  mkdir -p "$(dirname "$output_csv")"

  mapfile -t dataset_names < <(extract_dataset_names "$names_file")
  (( ${#dataset_names[@]} > 0 )) || die "No dataset names were found in $names_file"

  prepare_output_csv_for_resume "$output_csv"

  declare -A completed_dataset_names=()
  local existing_dataset_name
  while IFS= read -r existing_dataset_name; do
    completed_dataset_names["$existing_dataset_name"]=1
  done < <(extract_completed_output_names "$output_csv")

  local failed_datasets_file="${output_csv}.failed.txt"
  : >"$failed_datasets_file"

  local dataset_name
  local failure_count=0
  for dataset_name in "${dataset_names[@]}"; do
    if [[ -n "${completed_dataset_names[$dataset_name]+x}" ]]; then
      log "Skipping already completed dataset: $dataset_name"
      continue
    fi

    if process_dataset \
      "$source_dir" \
      "$dataset_name" \
      "$work_root" \
      "$reference_path" \
      "$kit_name" \
      "$dorado_model" \
      "$threads" \
      "$min_length" \
      "$max_length" \
      "$min_mean_q" \
      "$force" \
      "$output_csv"
    then
      :
    else
      failure_count=$((failure_count + 1))
      printf '%s\n' "$dataset_name" >>"$failed_datasets_file"
      log "Warning: dataset failed and was skipped: $dataset_name"
    fi
  done

  if (( failure_count == 0 )); then
    log "Completed CSV written to $output_csv"
  else
    log "Completed CSV written to $output_csv"
    log "Warning: $failure_count dataset(s) failed. See $failed_datasets_file"
  fi
}

main "$@"
