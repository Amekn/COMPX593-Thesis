#!/usr/bin/env bash
set -euo pipefail

# Usage: bash CalcStats.sh <bam_file>
bam_file="${1:-}"

if [[ -z "$bam_file" || ! -f "$bam_file" ]]; then
  echo "Usage: $0 <bam_file>" >&2
  exit 1
fi

tmp_stats="$(mktemp)"
trap 'rm -f "$tmp_stats"' EXIT

samtools stats "$bam_file" > "$tmp_stats"

# Pull SN fields from samtools stats
bases=$(awk -F'\t' '$1=="SN" && $2=="bases mapped (cigar):" {print $3}' "$tmp_stats")
mm=$(awk -F'\t' '$1=="SN" && $2=="mismatches:" {print $3}' "$tmp_stats")
reads_mapped=$(awk -F'\t' '$1=="SN" && $2=="reads mapped:" {print $3}' "$tmp_stats")
raw_total_sequences=$(awk -F'\t' '$1=="SN" && $2=="raw total sequences:" {print $3}' "$tmp_stats")

# Default to 0 if missing (prevents empty vars breaking awk math)
bases="${bases:-0}"
mm="${mm:-0}"
reads_mapped="${reads_mapped:-0}"
raw_total_sequences="${raw_total_sequences:-0}"

# Count total inserted/deleted bases across PRIMARY alignments only
# -F 0x904 skips: unmapped(0x4), secondary(0x100), supplementary(0x800)
read ins del < <(
  samtools view -F 0x904 "$bam_file" \
  | awk '{
      cigar=$6;
      while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
        token=substr(cigar, RSTART, RLENGTH);
        len=substr(token, 1, length(token)-1) + 0;
        op=substr(token, length(token), 1);
        if (op=="I") ins+=len;
        else if (op=="D") del+=len;
        cigar=substr(cigar, RSTART+RLENGTH);
      }
    }
    END {print ins+0, del+0}'
)

awk -v mm="$mm" -v ins="$ins" -v del="$del" -v bases="$bases" \
    -v reads_mapped="$reads_mapped" -v raw_total_sequences="$raw_total_sequences" '
BEGIN{
  if (raw_total_sequences == 0) {
    percentage_reads_mapped = 0;
  } else {
    percentage_reads_mapped = reads_mapped / raw_total_sequences;
  }

  if (bases == 0) {
    error_rate = 0;
    overall_error_rate = 0;
  } else {
    error_rate = mm / bases;
    overall_error_rate = (mm + ins + del) / bases;
  }

  printf "raw_total_sequences: %d\nreads_mapped: %d\ninserted_bases: %d\ndeleted_bases: %d\nmismatches: %d\nbases_mapped_(cigar): %d\nerror_rate: %.8f\noverall_error_rate: %.8f\npercentage_reads_mapped: %.8f\n",
         raw_total_sequences, reads_mapped, ins, del, mm, bases, error_rate, overall_error_rate, percentage_reads_mapped
}'
