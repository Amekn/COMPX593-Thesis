#!/usr/bin/env bash
set -euo pipefail

# GenerateCorrelationCsv.sh -- populates tables/correlation.csv for the thesis.
# Reads a CSV template that lists one row per (model, filter-stage) along with
# the per-row "variant_threshold" (single int or dash-separated list, e.g.
# "1-2-5-10"). For each row, locates the matching *.variants.tsv produced by
# DualSiteDMSFilter under <source_dir>, invokes VariantConcordance against
# the single ground_truth=TRUE row (MGI UDP0057), and fills in eight metrics
# per threshold (source/groundtruth/intersection/union variant counts plus
# exact_overlap_mass, weighted_jaccard, jensen_shannon_similarity,
# top100_spearman). Header order and row order are preserved; missing
# threshold-suffixed metric columns are appended. Dependencies: awk, find,
# sort, mktemp, and the project's VariantConcordance executable.

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
  - "variant_threshold" may be a single integer or a "-" separated list,
    for example: 1-2-5-10
  - Output metric columns are suffixed by threshold value, for example:
      source_variants_1
      source_variants_2
  - Missing metric columns are appended automatically when the template does
    not provide enough threshold-specific output columns.
  - Legacy "haplotype_*" column names are normalised automatically.
  - If the CSV does not contain "variant_threshold", the script uses 1.
    Empty threshold cells also default to 1.
EOF
}

# Print to stderr and exit 1. $*: message.
die() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

# Timestamped info line to stderr (stdout is reserved for CSV emission).
log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

# Die unless $1 is on PATH.
require_command() {
  local command_name=$1
  command -v "$command_name" >/dev/null 2>&1 || die "Required command not found: $command_name"
}

# Die unless $1 is an existing regular file; $2 is a human label.
require_file() {
  local path=$1
  local label=${2:-File}
  [[ -f "$path" ]] || die "$label does not exist: $path"
}

# Die unless $1 is an existing directory; $2 is a human label.
require_directory() {
  local path=$1
  local label=${2:-Directory}
  [[ -d "$path" ]] || die "$label does not exist: $path"
}

# Strip a trailing CR so CSVs authored on Windows behave identically to
# those authored on Unix. $1: raw string. Stdout: CR-stripped value.
trim_cr() {
  local value=$1
  printf '%s' "${value%$'\r'}"
}

# Strip leading + trailing whitespace from $1; echo the trimmed value.
# Uses parameter expansion only so it is safe on arbitrary bytes.
trim_whitespace() {
  local value=$1
  value=${value#"${value%%[![:space:]]*}"}
  value=${value%"${value##*[![:space:]]}"}
  printf '%s' "$value"
}

# Upper-case $1 for case-insensitive comparison of the ground_truth column.
to_upper() {
  local value=$1
  printf '%s' "$value" | tr '[:lower:]' '[:upper:]'
}

# Parse the per-row variant_threshold cell into a de-duplicated integer
# array. Accepts a single int, dash/semicolon/colon-separated lists, and
# empty (defaults to (1)). $1: raw cell text, $2: nameref to destination
# array. Dies on any non-integer token.
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

  # Collapse spaces and accept `;` or `:` as aliases for `-` so the
  # user-facing spec is forgiving.
  normalized_value=${raw_value//[[:space:]]/}
  normalized_value=${normalized_value//;/-}
  normalized_value=${normalized_value//:/-}

  IFS='-' read -r -a tokens <<< "$normalized_value"
  (( ${#tokens[@]} > 0 )) || die "Invalid variant_threshold list: $raw_value"

  for token in "${tokens[@]}"; do
    [[ -n "$token" ]] || die "Invalid variant_threshold list: $raw_value"
    [[ "$token" =~ ^[0-9]+$ ]] || die "Invalid variant_threshold value: $token"

    # 10# forces decimal interpretation so leading-zero tokens (e.g. "08")
    # do not trip bash's octal parser.
    canonical_threshold=$((10#$token))
    token=$canonical_threshold

    if [[ ! -v "seen_thresholds[$token]" ]]; then
      seen_thresholds[$token]=1
      destination+=("$token")
    fi
  done
}

# Rename "haplotype" column families to "variant" so old templates work
# unchanged against the current VariantConcordance metric set. $1: header
# field; stdout: normalised field.
normalize_legacy_header_field() {
  local field=$1

  case "$field" in
    haplotype_threshold)
      printf 'variant_threshold\n'
      ;;
    source_haplotypes_*)
      printf 'source_variants_%s\n' "${field#source_haplotypes_}"
      ;;
    groundtruth_haplotypes_*)
      printf 'groundtruth_variants_%s\n' "${field#groundtruth_haplotypes_}"
      ;;
    intersection_haplotypes_*)
      printf 'intersection_variants_%s\n' "${field#intersection_haplotypes_}"
      ;;
    union_haplotypes_*)
      printf 'union_variants_%s\n' "${field#union_haplotypes_}"
      ;;
    *)
      printf '%s\n' "$field"
      ;;
  esac
}

# Numeric-sort-unique a nameref'd array into a second nameref'd array.
# $1: source array name, $2: destination array name. Used to produce a
# stable, ascending global threshold order for header backfill.
sort_unique_numeric_values() {
  local -n source_ref=$1
  local -n destination_ref=$2

  destination_ref=()
  (( ${#source_ref[@]} > 0 )) || return 0

  mapfile -t destination_ref < <(printf '%s\n' "${source_ref[@]}" | sort -n -u)
}

# Resolve the VariantConcordance binary: honour --concordance-bin, else
# search PATH, else fall back to ../../build/VariantConcordance relative to
# this script. Mirrors resolve_filter_binary() in DMSPolishing.sh.
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

# Split a CSV line into a nameref'd array of exactly $2 cells, padding with
# empty strings when the template row is short. The trailing "," trick
# protects against `read -a` dropping a final empty field. $1: raw line,
# $2: expected column count, $3: destination array name.
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

# Print a CSV row from a nameref'd array with commas between fields and a
# trailing newline. No quoting (template values are plain).
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

# Find the single *.variants.tsv matching a CSV row name. Accepts either
# `<stem>.variants.tsv` or `<stem>.<suffix>.variants.tsv` so per-threshold
# / per-stage suffixes (e.g. `<stem>.t3.variants.tsv`) still resolve.
# Errors on ambiguity (>1 match) or absence. $1: row name, $2: nameref
# containing the full list of variant files under source_dir.
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
  # Order matches VariantConcordance's TSV output: fields[1..8] are pulled
  # index-by-index into column <prefix>_<threshold>. Keep in sync with the
  # assignments inside the per-threshold loop at the end of main().
  local metric_prefixes=(
    source_variants
    groundtruth_variants
    intersection_variants
    union_variants
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

  # Write into a hidden sibling file of $output_csv then atomic-rename on
  # success. %q is used so special characters in the path are safely
  # embedded into the EXIT trap string.
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

  # Build column_index lookup after normalising any legacy "haplotype_*"
  # column names. This keeps `column_index[source_variants_3]` valid even
  # when the on-disk template still used the old "haplotype_" spelling.
  local -A column_index=()
  local index
  local normalized_header
  for index in "${!header_fields[@]}"; do
    header_fields[index]=$(trim_whitespace "$(trim_cr "${header_fields[index]}")")
    normalized_header=$(normalize_legacy_header_field "${header_fields[index]}")
    [[ ! -v "column_index[$normalized_header]" ]] || die "Template CSV contains duplicate or conflicting column: $normalized_header"
    header_fields[index]=$normalized_header
    column_index["$normalized_header"]=$index
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
  # variant_threshold is optional; -1 sentinel means "assume threshold=1 for
  # every row" per the documented default.
  local variant_threshold_index=-1
  if [[ -v "column_index[variant_threshold]" ]]; then
    variant_threshold_index=${column_index[variant_threshold]}
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

  # First pass over data rows: cache the raw lines, collect every distinct
  # threshold seen anywhere in the template (used to extend the header with
  # <prefix>_<threshold> columns), and enforce exactly one ground_truth=TRUE
  # row.
  while IFS= read -r row_line || [[ -n "$row_line" ]]; do
    row_line=$(trim_cr "$row_line")
    [[ -n "$row_line" ]] || continue

    rows+=("$row_line")
    parse_csv_line "$row_line" "$header_count" parsed_row

    row_name=$(trim_whitespace "${parsed_row[name_index]}")
    [[ -n "$row_name" ]] || die "Encountered a row with an empty name"

    row_threshold_text=''
    if (( variant_threshold_index >= 0 )); then
      row_threshold_text=${parsed_row[variant_threshold_index]}
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

  # Extend header with <prefix>_<threshold> columns for any thresholds not
  # already present. Ascending threshold order yields a stable schema.
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

  # Walk <source_dir> once with null-delimited find -print0 so the
  # file-name list is robust to spaces/newlines; reused for every row.
  local variant_files=()
  while IFS= read -r -d '' row_line; do
    variant_files+=("$row_line")
  done < <(find "$source_dir" -type d -name '*.failed_*' -prune -o -type f -name '*.variants.tsv' -print0)
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
    if (( variant_threshold_index >= 0 )); then
      row_threshold_text=${parsed_row[variant_threshold_index]}
    fi
    parse_threshold_list "$row_threshold_text" row_thresholds
    thresholds_summary=$(IFS='-'; printf '%s' "${row_thresholds[*]}")

    for threshold in "${row_thresholds[@]}"; do
      log "Comparing row \"$row_name\" using $(basename "$source_path") against $(basename "$ground_truth_path") with threshold $threshold (row thresholds: $thresholds_summary)"
      if ! concordance_output=$("$concordance_bin" "$source_path" "$ground_truth_path" "$threshold"); then
        die "VariantConcordance failed for row \"$row_name\" at threshold $threshold"
      fi

      # VariantConcordance writes a header line then a single tab-separated
      # metrics line; pick NR==2. Expect 9 fields: [0] is threshold echo,
      # [1..8] are the metrics in metric_prefixes order.
      metrics_line=$(printf '%s\n' "$concordance_output" | awk 'NR == 2 { print; exit }')
      [[ -n "$metrics_line" ]] || die "VariantConcordance did not return a metrics row for \"$row_name\" at threshold $threshold"

      IFS=$'\t' read -r -a metrics_fields <<< "$metrics_line"
      (( ${#metrics_fields[@]} == 9 )) || die "Unexpected VariantConcordance output for row \"$row_name\" at threshold $threshold: $metrics_line"

      parsed_row[${column_index[source_variants_$threshold]}]=${metrics_fields[1]}
      parsed_row[${column_index[groundtruth_variants_$threshold]}]=${metrics_fields[2]}
      parsed_row[${column_index[intersection_variants_$threshold]}]=${metrics_fields[3]}
      parsed_row[${column_index[union_variants_$threshold]}]=${metrics_fields[4]}
      parsed_row[${column_index[exact_overlap_mass_$threshold]}]=${metrics_fields[5]}
      parsed_row[${column_index[weighted_jaccard_$threshold]}]=${metrics_fields[6]}
      parsed_row[${column_index[jensen_shannon_similarity_$threshold]}]=${metrics_fields[7]}
      parsed_row[${column_index[top100_spearman_$threshold]}]=${metrics_fields[8]}
    done

    emit_csv_line parsed_row >>"$temporary_output"
  done

  # Atomic rename; EXIT trap cleans up temp if we die before this point.
  mv -f -- "$temporary_output" "$output_csv"
  log "Wrote completed correlation CSV to $output_csv"
}

main "$@"
