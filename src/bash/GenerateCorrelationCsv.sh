#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage:
  GenerateCorrelationCsv.sh <template_csv> <source_dir> <output_csv> [options]

Description:
  Populate a correlation CSV template by resolving each row's variant-key TSV,
  comparing it against the single ground-truth row with VariantConcordance, and
  writing a completed CSV while preserving the original header and row order.

Required positional arguments:
  <template_csv>   Input CSV with at least a "name" column and one
                   "ground_truth=TRUE" row
  <source_dir>     Root directory to search recursively for variant-key TSVs
  <output_csv>     Output CSV path

Options:
  --concordance-bin PATH   Explicit path to VariantConcordance
  -h, --help               Show this help message

Notes:
  - Each row name is matched only against *.variants.tsv files under
    <source_dir>.
  - A row name matches files shaped like:
      *<name>.variants.tsv
      *<name>.*.variants.tsv
  - "haplotype_threshold" may be a single integer or a "-" separated list,
    for example: 1-2-5-10
  - Output metric columns are suffixed by threshold value, for example:
      source_haplotypes_1
      source_haplotypes_2
  - Missing metric columns are appended automatically when the template does
    not provide enough threshold-specific output columns.
  - If the CSV does not contain "haplotype_threshold", the script uses 1.
    Empty threshold cells also default to 1.
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

trim_whitespace() {
  local value=$1
  value=${value#"${value%%[![:space:]]*}"}
  value=${value%"${value##*[![:space:]]}"}
  printf '%s' "$value"
}

to_upper() {
  local value=$1
  printf '%s' "$value" | tr '[:lower:]' '[:upper:]'
}

parse_threshold_list() {
  local raw_value=$1
  local -n destination=$2
  local normalized_value
  local -a tokens=()
  local token
  local canonical_threshold
  local -A seen_thresholds=()

  destination=()
  raw_value=$(trim_whitespace "$raw_value")

  if [[ -z "$raw_value" ]]; then
    destination=(1)
    return 0
  fi

  normalized_value=${raw_value//[[:space:]]/}
  normalized_value=${normalized_value//;/-}
  normalized_value=${normalized_value//:/-}

  IFS='-' read -r -a tokens <<< "$normalized_value"
  (( ${#tokens[@]} > 0 )) || die "Invalid haplotype_threshold list: $raw_value"

  for token in "${tokens[@]}"; do
    [[ -n "$token" ]] || die "Invalid haplotype_threshold list: $raw_value"
    [[ "$token" =~ ^[0-9]+$ ]] || die "Invalid haplotype_threshold value: $token"

    canonical_threshold=$((10#$token))
    token=$canonical_threshold

    if [[ ! -v "seen_thresholds[$token]" ]]; then
      seen_thresholds[$token]=1
      destination+=("$token")
    fi
  done
}

sort_unique_numeric_values() {
  local -n source_ref=$1
  local -n destination_ref=$2

  destination_ref=()
  (( ${#source_ref[@]} > 0 )) || return 0

  mapfile -t destination_ref < <(printf '%s\n' "${source_ref[@]}" | sort -n -u)
}

resolve_concordance_binary() {
  local requested_path=$1
  local script_directory
  local project_root

  script_directory=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  project_root=$(cd "$script_directory/../.." && pwd)

  if [[ -n "$requested_path" ]]; then
    [[ -x "$requested_path" ]] || die "VariantConcordance is not executable: $requested_path"
    printf '%s\n' "$requested_path"
    return 0
  fi

  if command -v VariantConcordance >/dev/null 2>&1; then
    command -v VariantConcordance
    return 0
  fi

  if [[ -x "$project_root/build/VariantConcordance" ]]; then
    printf '%s\n' "$project_root/build/VariantConcordance"
    return 0
  fi

  die "VariantConcordance was not found on PATH or under $project_root/build"
}

parse_csv_line() {
  local line=$1
  local expected_columns=$2
  local -n destination=$3

  IFS=, read -r -a destination <<< "$(trim_cr "$line"),"

  while (( ${#destination[@]} > expected_columns )); do
    unset 'destination[${#destination[@]}-1]'
  done

  while (( ${#destination[@]} < expected_columns )); do
    destination+=("")
  done
}

emit_csv_line() {
  local -n fields_ref=$1
  local field_count=${#fields_ref[@]}
  local index

  for (( index = 0; index < field_count; ++index )); do
    if (( index > 0 )); then
      printf ','
    fi
    printf '%s' "${fields_ref[index]}"
  done
  printf '\n'
}

resolve_variant_path() {
  local row_name=$1
  local -n variant_files_ref=$2
  local raw_name
  local base_name
  local search_name
  local path
  local file_name
  local matches=()

  raw_name=$(trim_whitespace "$row_name")
  [[ -n "$raw_name" ]] || die "Encountered a row with an empty name"

  base_name=${raw_name##*/}
  search_name=${base_name%.variants.tsv}
  search_name=${search_name%.tsv}

  for path in "${variant_files_ref[@]}"; do
    file_name=${path##*/}
    if [[ "$file_name" == *"${search_name}.variants.tsv" || "$file_name" == *"${search_name}."*.variants.tsv ]]; then
      matches+=("$path")
    fi
  done

  if (( ${#matches[@]} == 1 )); then
    printf '%s\n' "${matches[0]}"
    return 0
  fi
  if (( ${#matches[@]} > 1 )); then
    printf 'Error: multiple %s matches found for row "%s":\n' \
      "\"*${search_name}.variants.tsv\" / \"*${search_name}.*.variants.tsv\"" \
      "$row_name" >&2
    printf '  %s\n' "${matches[@]}" >&2
    exit 1
  fi

  die "No variant-key TSV match found for row \"$row_name\" under $source_dir using *${search_name}.variants.tsv or *${search_name}.*.variants.tsv"
}

main() {
  local concordance_bin_option=''
  local positional_args=()

  if [[ $# -eq 0 ]]; then
    usage
    exit 1
  fi

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -h|--help)
        usage
        exit 0
        ;;
      --concordance-bin)
        [[ $# -ge 2 ]] || die "--concordance-bin requires a value"
        concordance_bin_option=$2
        shift 2
        ;;
      --)
        shift
        while [[ $# -gt 0 ]]; do
          positional_args+=("$1")
          shift
        done
        break
        ;;
      -*)
        die "Unknown option: $1"
        ;;
      *)
        positional_args+=("$1")
        shift
        ;;
    esac
  done

  [[ ${#positional_args[@]} -eq 3 ]] || {
    usage
    exit 1
  }

  local template_csv=${positional_args[0]}
  local source_dir=${positional_args[1]}
  local output_csv=${positional_args[2]}
  local metric_prefixes=(
    source_haplotypes
    groundtruth_haplotypes
    intersection_haplotypes
    union_haplotypes
    exact_overlap_mass
    weighted_jaccard
    jensen_shannon_similarity
    top100_spearman
  )

  require_command awk
  require_command find
  require_command mktemp
  require_command sort

  require_file "$template_csv" "Template CSV"
  require_directory "$source_dir" "Source directory"

  local concordance_bin
  concordance_bin=$(resolve_concordance_binary "$concordance_bin_option")

  local output_dir
  output_dir=$(dirname "$output_csv")
  mkdir -p "$output_dir"

  local temporary_output
  local cleanup_trap
  temporary_output=$(mktemp "$output_dir/.tmp.$(basename "$output_csv").XXXXXX")
  printf -v cleanup_trap 'rm -f -- %q' "$temporary_output"
  trap "$cleanup_trap" EXIT

  local header_line
  if ! IFS= read -r header_line < "$template_csv"; then
    die "Template CSV is empty: $template_csv"
  fi
  header_line=$(trim_cr "$header_line")

  local header_fields=()
  IFS=, read -r -a header_fields <<< "$header_line"
  local header_count=${#header_fields[@]}
  (( header_count > 0 )) || die "Template CSV header is empty: $template_csv"

  local -A column_index=()
  local index
  for index in "${!header_fields[@]}"; do
    header_fields[index]=$(trim_whitespace "$(trim_cr "${header_fields[index]}")")
    column_index["${header_fields[index]}"]=$index
  done

  local required_column
  for required_column in \
    name \
    ground_truth
  do
    [[ -v "column_index[$required_column]" ]] || die "Template CSV is missing required column: $required_column"
  done

  local name_index=${column_index[name]}
  local ground_truth_index=${column_index[ground_truth]}
  local haplotype_threshold_index=-1
  if [[ -v "column_index[haplotype_threshold]" ]]; then
    haplotype_threshold_index=${column_index[haplotype_threshold]}
  fi

  local rows=()
  local parsed_row=()
  local row_line
  local row_name
  local ground_truth_value
  local ground_truth_name=''
  local ground_truth_count=0
  local row_threshold_text
  local row_thresholds=()
  local all_thresholds=()
  local sorted_thresholds=()
  local -A seen_output_thresholds=()
  local threshold
  local metric_prefix
  local metric_column_name

  while IFS= read -r row_line || [[ -n "$row_line" ]]; do
    row_line=$(trim_cr "$row_line")
    [[ -n "$row_line" ]] || continue

    rows+=("$row_line")
    parse_csv_line "$row_line" "$header_count" parsed_row

    row_name=$(trim_whitespace "${parsed_row[name_index]}")
    [[ -n "$row_name" ]] || die "Encountered a row with an empty name"

    row_threshold_text=''
    if (( haplotype_threshold_index >= 0 )); then
      row_threshold_text=${parsed_row[haplotype_threshold_index]}
    fi
    parse_threshold_list "$row_threshold_text" row_thresholds
    for threshold in "${row_thresholds[@]}"; do
      if [[ ! -v "seen_output_thresholds[$threshold]" ]]; then
        seen_output_thresholds[$threshold]=1
        all_thresholds+=("$threshold")
      fi
    done

    ground_truth_value=$(to_upper "$(trim_whitespace "${parsed_row[ground_truth_index]}")")
    if [[ -z "$ground_truth_value" || "$ground_truth_value" == "FALSE" ]]; then
      continue
    fi
    if [[ "$ground_truth_value" != "TRUE" ]]; then
      die "Invalid ground_truth value for row \"$row_name\": ${parsed_row[ground_truth_index]}"
    fi

    (( ++ground_truth_count ))
    ground_truth_name=$row_name
  done < <(tail -n +2 "$template_csv")

  (( ${#rows[@]} > 0 )) || die "Template CSV has no data rows: $template_csv"
  (( ground_truth_count == 1 )) || die "Template CSV must contain exactly one ground_truth=TRUE row (found $ground_truth_count)"

  sort_unique_numeric_values all_thresholds sorted_thresholds
  for threshold in "${sorted_thresholds[@]}"; do
    for metric_prefix in "${metric_prefixes[@]}"; do
      metric_column_name=${metric_prefix}_${threshold}
      if [[ ! -v "column_index[$metric_column_name]" ]]; then
        column_index["$metric_column_name"]=${#header_fields[@]}
        header_fields+=("$metric_column_name")
      fi
    done
  done
  header_count=${#header_fields[@]}

  local variant_files=()
  while IFS= read -r -d '' row_line; do
    variant_files+=("$row_line")
  done < <(find "$source_dir" -type f -name '*.variants.tsv' -print0)
  (( ${#variant_files[@]} > 0 )) || die "No *.variants.tsv files were found under $source_dir"

  local ground_truth_path
  ground_truth_path=$(resolve_variant_path "$ground_truth_name" variant_files)
  log "Ground truth row \"$ground_truth_name\" resolved to $ground_truth_path"

  emit_csv_line header_fields >"$temporary_output"

  local source_path
  local concordance_output
  local metrics_line
  local metrics_fields=()
  local thresholds_summary

  for row_line in "${rows[@]}"; do
    parse_csv_line "$row_line" "$header_count" parsed_row
    row_name=$(trim_whitespace "${parsed_row[name_index]}")

    source_path=$(resolve_variant_path "$row_name" variant_files)

    row_threshold_text=''
    if (( haplotype_threshold_index >= 0 )); then
      row_threshold_text=${parsed_row[haplotype_threshold_index]}
    fi
    parse_threshold_list "$row_threshold_text" row_thresholds
    thresholds_summary=$(IFS='-'; printf '%s' "${row_thresholds[*]}")

    for threshold in "${row_thresholds[@]}"; do
      log "Comparing row \"$row_name\" using $(basename "$source_path") against $(basename "$ground_truth_path") with threshold $threshold (row thresholds: $thresholds_summary)"
      if ! concordance_output=$("$concordance_bin" "$source_path" "$ground_truth_path" "$threshold"); then
        die "VariantConcordance failed for row \"$row_name\" at threshold $threshold"
      fi

      metrics_line=$(printf '%s\n' "$concordance_output" | awk 'NR == 2 { print; exit }')
      [[ -n "$metrics_line" ]] || die "VariantConcordance did not return a metrics row for \"$row_name\" at threshold $threshold"

      IFS=$'\t' read -r -a metrics_fields <<< "$metrics_line"
      (( ${#metrics_fields[@]} == 9 )) || die "Unexpected VariantConcordance output for row \"$row_name\" at threshold $threshold: $metrics_line"

      parsed_row[${column_index[source_haplotypes_$threshold]}]=${metrics_fields[1]}
      parsed_row[${column_index[groundtruth_haplotypes_$threshold]}]=${metrics_fields[2]}
      parsed_row[${column_index[intersection_haplotypes_$threshold]}]=${metrics_fields[3]}
      parsed_row[${column_index[union_haplotypes_$threshold]}]=${metrics_fields[4]}
      parsed_row[${column_index[exact_overlap_mass_$threshold]}]=${metrics_fields[5]}
      parsed_row[${column_index[weighted_jaccard_$threshold]}]=${metrics_fields[6]}
      parsed_row[${column_index[jensen_shannon_similarity_$threshold]}]=${metrics_fields[7]}
      parsed_row[${column_index[top100_spearman_$threshold]}]=${metrics_fields[8]}
    done

    emit_csv_line parsed_row >>"$temporary_output"
  done

  mv -f -- "$temporary_output" "$output_csv"
  log "Wrote completed correlation CSV to $output_csv"
}

main "$@"
