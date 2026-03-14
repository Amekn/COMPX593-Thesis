#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage:
  CalcStats.sh <alignment.bam>

Description:
  Summarise a BAM file using `samtools stats` and report mismatch, indel, and
  mapping-rate metrics for primary alignments.
EOF
}

die() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

require_command() {
  local command_name=$1
  command -v "$command_name" >/dev/null 2>&1 || die "Required command not found: $command_name"
}

extract_sn_metric() {
  local stats_path=$1
  local metric_name=$2
  awk -F '\t' -v key="$metric_name" '$1 == "SN" && $2 == key { print $3; exit }' "$stats_path"
}

count_primary_indel_bases() {
  local bam_path=$1

  # Primary alignments only:
  #   0x4   unmapped
  #   0x100 secondary
  #   0x800 supplementary
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
          print inserted_bases + 0, deleted_bases + 0
        }
      '
}

main() {
  if [[ $# -ne 1 ]]; then
    usage
    exit 1
  fi
  if [[ ${1:-} == "--help" || ${1:-} == "-h" ]]; then
    usage
    exit 0
  fi

  local bam_path=$1
  [[ -f "$bam_path" ]] || die "Input BAM file does not exist: $bam_path"

  require_command samtools
  require_command awk
  require_command mktemp

  local temporary_stats
  temporary_stats=$(mktemp)
  trap 'rm -f "$temporary_stats"' EXIT

  samtools stats "$bam_path" >"$temporary_stats"

  local bases_mapped
  local mismatch_count
  local reads_mapped
  local raw_total_sequences
  bases_mapped=$(extract_sn_metric "$temporary_stats" "bases mapped (cigar):")
  mismatch_count=$(extract_sn_metric "$temporary_stats" "mismatches:")
  reads_mapped=$(extract_sn_metric "$temporary_stats" "reads mapped:")
  raw_total_sequences=$(extract_sn_metric "$temporary_stats" "raw total sequences:")

  bases_mapped=${bases_mapped:-0}
  mismatch_count=${mismatch_count:-0}
  reads_mapped=${reads_mapped:-0}
  raw_total_sequences=${raw_total_sequences:-0}

  local inserted_bases deleted_bases
  read -r inserted_bases deleted_bases < <(count_primary_indel_bases "$bam_path")

  awk \
    -v raw_total_sequences="$raw_total_sequences" \
    -v reads_mapped="$reads_mapped" \
    -v inserted_bases="$inserted_bases" \
    -v deleted_bases="$deleted_bases" \
    -v mismatch_count="$mismatch_count" \
    -v bases_mapped="$bases_mapped" '
      BEGIN {
        if (raw_total_sequences == 0) {
          percentage_reads_mapped = 0
        } else {
          percentage_reads_mapped = reads_mapped / raw_total_sequences
        }

        if (reads_mapped == 0) {
          average_mismatches_per_mapped_read = 0
        } else {
          average_mismatches_per_mapped_read = mismatch_count / reads_mapped
        }

        if (bases_mapped == 0) {
          mismatch_rate = 0
          overall_error_rate = 0
        } else {
          mismatch_rate = mismatch_count / bases_mapped
          overall_error_rate = (mismatch_count + inserted_bases + deleted_bases) / bases_mapped
        }

        printf "raw_total_sequences: %d\n", raw_total_sequences
        printf "reads_mapped: %d\n", reads_mapped
        printf "inserted_bases: %d\n", inserted_bases
        printf "deleted_bases: %d\n", deleted_bases
        printf "mismatches: %d\n", mismatch_count
        printf "bases_mapped_(cigar): %d\n", bases_mapped
        printf "error_rate: %.8f\n", mismatch_rate
        printf "overall_error_rate: %.8f\n", overall_error_rate
        printf "percentage_reads_mapped: %.8f\n", percentage_reads_mapped
        printf "avg_mismatches_per_pass_read: %.8f\n", average_mismatches_per_mapped_read
      }
    '
}

main "$@"
