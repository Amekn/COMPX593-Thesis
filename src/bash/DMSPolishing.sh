#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage:
  Basecall mode:
    DMSPolishing.sh --model <model_path> --test <pod5_path> --ref <reference_fasta> \
      --out <output_fastq> --zones <zone_file> [options]

  Existing BAM mode:
    DMSPolishing.sh --bam <raw_basecall.bam> --ref <reference_fasta> \
      --out <output_fastq> --zones <zone_file> [options]

Options:
  --threads <N>              Number of worker threads to use. Default: 32
  --kit-name <NAME>          Dorado kit name. Default: SQK-NBD114-24
  --min-length <N>           Minimum read length for filtering. Default: 600
  --max-length <N>           Maximum read length for filtering. Default: 800
  --max-indel-events <N>     Maximum indel events for DualSiteDMSFilter. Default: 3
  --max-indel-bases <N>      Maximum indel bases for DualSiteDMSFilter. Default: 3
  --mut-codon-library <LIB>  Library constraint for DualSiteDMSFilter. Default: NNK
  --filter-bin <PATH>        Explicit path to DualSiteDMSFilter.
  -h, --help                 Show this help message.

Notes:
  - Exactly one of `--bam` or (`--model` and `--test`) must be provided.
  - The first non-empty, non-comment line in `--zones` is used as the zone string.
  - Output-derived intermediate files are written next to `--out`.
EOF
}

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

log() {
  printf '[%s] %s\n' "$(timestamp)" "$*"
}

step() {
  printf '\n[%s] ===== %s =====\n' "$(timestamp)" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

require_command() {
  local command_name=$1
  command -v "$command_name" >/dev/null 2>&1 || die "Required command not found: $command_name"
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
  local script_directory=$2
  local project_root
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

require_file() {
  local path=$1
  local label=$2
  [[ -f "$path" ]] || die "$label does not exist: $path"
}

main() {
  local model_path=""
  local pod5_path=""
  local basecall_bam=""
  local reference_fasta=""
  local output_fastq=""
  local zones_file=""
  local threads=32
  local kit_name="SQK-NBD114-24"
  local min_length=600
  local max_length=800
  local max_indel_events=3
  local max_indel_bases=3
  local mutant_codon_library="NNK"
  local filter_binary_override=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --model)
        model_path=${2:-}
        shift 2
        ;;
      --test)
        pod5_path=${2:-}
        shift 2
        ;;
      --bam)
        basecall_bam=${2:-}
        shift 2
        ;;
      --ref)
        reference_fasta=${2:-}
        shift 2
        ;;
      --out)
        output_fastq=${2:-}
        shift 2
        ;;
      --zones)
        zones_file=${2:-}
        shift 2
        ;;
      --threads)
        threads=${2:-}
        shift 2
        ;;
      --kit-name)
        kit_name=${2:-}
        shift 2
        ;;
      --min-length)
        min_length=${2:-}
        shift 2
        ;;
      --max-length)
        max_length=${2:-}
        shift 2
        ;;
      --max-indel-events)
        max_indel_events=${2:-}
        shift 2
        ;;
      --max-indel-bases)
        max_indel_bases=${2:-}
        shift 2
        ;;
      --mut-codon-library)
        mutant_codon_library=${2:-}
        shift 2
        ;;
      --filter-bin)
        filter_binary_override=${2:-}
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        usage
        die "Unknown option: $1"
        ;;
    esac
  done

  [[ -n "$reference_fasta" && -n "$output_fastq" && -n "$zones_file" ]] || {
    usage
    exit 1
  }

  if [[ -n "$basecall_bam" ]]; then
    [[ -z "$model_path" && -z "$pod5_path" ]] || die "Use either --bam or (--model and --test), not both."
  else
    [[ -n "$model_path" && -n "$pod5_path" ]] || die "Basecall mode requires both --model and --test."
  fi

  [[ "$threads" =~ ^[0-9]+$ && "$threads" -gt 0 ]] || die "--threads must be a positive integer."
  [[ "$min_length" =~ ^[0-9]+$ && "$max_length" =~ ^[0-9]+$ && "$min_length" -lt "$max_length" ]] \
    || die "--min-length and --max-length must be positive integers with min < max."
  [[ "$max_indel_events" =~ ^[0-9]+$ ]] || die "--max-indel-events must be a non-negative integer."
  [[ "$max_indel_bases" =~ ^[0-9]+$ ]] || die "--max-indel-bases must be a non-negative integer."

  require_file "$reference_fasta" "Reference FASTA"
  require_file "$zones_file" "Zone file"
  if [[ -n "$basecall_bam" ]]; then
    require_file "$basecall_bam" "Input BAM"
  else
    require_file "$model_path" "Model"
    require_file "$pod5_path" "POD5 input"
    require_command dorado
  fi

  require_command samtools
  require_command minimap2
  require_command fastplong
  require_command awk
  require_command tee

  local script_directory
  script_directory=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  local calc_stats_script="$script_directory/CalcStats.sh"
  require_file "$calc_stats_script" "CalcStats helper"
  local filter_binary
  filter_binary=$(resolve_filter_binary "$filter_binary_override" "$script_directory")

  local zone_string
  zone_string=$(read_zone_string "$zones_file")
  [[ -n "$zone_string" ]] || die "Zone file does not contain a usable zone string: $zones_file"

  local output_prefix=${output_fastq%.fastq}
  output_prefix=${output_prefix%.fq}
  local timeline_log="${output_prefix}.timeline.log"

  mkdir -p "$(dirname "$output_fastq")"
  exec > >(tee -a "$timeline_log") 2>&1

  local raw_bam="${output_prefix}.basecalled.bam"
  local raw_fastq="${output_prefix}.basecalled.fq"
  local filtered_fastq="${output_prefix}.len${min_length}_${max_length}.fq"
  local aligned_sorted_bam="${output_prefix}.aligned.sorted.bam"
  local primary_name_sorted_bam="${output_prefix}.primary.sorted.bam"
  local polished_aligned_bam="${output_prefix}.polished.sorted.bam"

  log "DMS polishing pipeline"
  log "reference_fasta = $reference_fasta"
  log "output_fastq    = $output_fastq"
  log "zone_string     = $zone_string"
  log "threads         = $threads"
  log "timeline_log    = $timeline_log"
  log "filter_binary   = $filter_binary"

  if [[ -n "$basecall_bam" ]]; then
    log "mode            = existing_bam"
    log "input_bam       = $basecall_bam"
    raw_bam=$basecall_bam
  else
    log "mode            = basecall"
    log "model_path      = $model_path"
    log "pod5_path       = $pod5_path"
    log "kit_name        = $kit_name"

    step "Basecalling with dorado"
    dorado basecaller "$model_path" "$pod5_path" \
      --kit-name "$kit_name" \
      --trim all \
      >"$raw_bam"
    log "basecall_bam    = $raw_bam"
  fi

  step "Convert BAM to FASTQ"
  samtools fastq "$raw_bam" >"$raw_fastq"
  log "raw_fastq       = $raw_fastq"

  step "Length filtering"
  fastplong \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --length_required "$min_length" \
    --length_limit "$max_length" \
    --thread "$threads" \
    --in "$raw_fastq" \
    --out "$filtered_fastq"
  log "filtered_fastq  = $filtered_fastq"

  step "Align filtered reads"
  minimap2 -t "$threads" -ax map-ont "$reference_fasta" "$filtered_fastq" \
    | samtools view -u - \
    | samtools sort -@ "$threads" -o "$aligned_sorted_bam" -
  log "aligned_bam     = $aligned_sorted_bam"

  step "Extract primary alignments"
  samtools view -F 0x904 -u "$aligned_sorted_bam" \
    | samtools sort -@ "$threads" -n -o "$primary_name_sorted_bam" -
  log "primary_bam     = $primary_name_sorted_bam"

  step "Summarise primary alignments"
  "$calc_stats_script" "$primary_name_sorted_bam"

  step "Polish reads with DualSiteDMSFilter"
  "$filter_binary" \
    "$primary_name_sorted_bam" \
    "$reference_fasta" \
    "$output_fastq" \
    "$zone_string" \
    --max-indel-events "$max_indel_events" \
    --max-indel-bases "$max_indel_bases" \
    --mut-codon-library "$mutant_codon_library"
  log "polished_fastq  = $output_fastq"

  step "Align polished reads"
  minimap2 -t "$threads" -ax map-ont "$reference_fasta" "$output_fastq" \
    | samtools view -u - \
    | samtools sort -@ "$threads" -o "$polished_aligned_bam" -
  log "polished_bam    = $polished_aligned_bam"

  step "Summarise polished alignments"
  "$calc_stats_script" "$polished_aligned_bam"

  step "Pipeline complete"
  log "primary_bam_pre_polish = $primary_name_sorted_bam"
  log "polished_aligned_bam   = $polished_aligned_bam"
  log "timeline_log           = $timeline_log"
}

main "$@"
