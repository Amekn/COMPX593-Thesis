// UmiExtractor: dual-UMI extraction stage of the thesis pipeline.
//
// Produces the `UmiExtractor` executable. For each nanopore FASTQ record the
// tool locates the 1st-round PCR primer pair (forward + reverse), reads off the
// degenerate (IUPAC ambiguous) UMI bases embedded in each primer, appends a
// ":UMI_<fwd>_<rev>" tag to the header, and trims the read down to the insert
// between the two UMI stretches. Reads are tried first in their original
// orientation, then reverse-complemented; orientation-specific UMIs are
// canonicalised back into the forward frame so downstream tools (Benchmarker,
// DualSiteDMSFilter, VariantConcordance) always see one key per template.
//
// Usage: UmiExtractor <input.fastq> <output.fastq> <fwd_primer> <rev_primer>
//                    <max_mismatch>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Sequence.hpp"

namespace {

// Immutable per-run description of both 1st-round primers plus their reverse
// complements, and the 0-based offsets of every ambiguous (UMI) base within
// each of those four strings. Precomputing the RC strings and UMI offsets
// once lets the per-read hot loop avoid recomputing IUPAC bookkeeping.
struct PrimerLayout {
    std::string forward_primer;
    std::string reverse_primer;
    std::string forward_primer_reverse_complement;
    std::string reverse_primer_reverse_complement;
    std::vector<int> forward_umi_positions;
    std::vector<int> reverse_umi_positions;
    std::vector<int> forward_reverse_complement_umi_positions;
    std::vector<int> reverse_reverse_complement_umi_positions;
};

// Emit the command-line synopsis when arg count / parsing fails.
void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name
        << " <input.fastq> <output.fastq> <forward_primer> <reverse_primer> <max_mismatch>\n\n"
        << "Extract dual UMIs from ambiguous primer positions, append them to the FASTQ header,\n"
        << "and trim the read to the insert region bounded by the primer pair.\n";
}

// Locate a primer pair within `uppercase_sequence` using the leftmost match
// for `left_primer` and the rightmost match for `right_primer`, each allowing
// up to `max_mismatches` Hamming-style errors. Writes primer start offsets
// into *left_start / *right_start. Returns false if either primer is absent
// or the primers overlap (would leave no insert).
bool TryFindPrimerSpan(
    const std::string& uppercase_sequence,
    const std::string& left_primer,
    const std::string& right_primer,
    const int max_mismatches,
    std::size_t* left_start,
    std::size_t* right_start) {
    const int left_index = ont::seq::FindPrimerFromLeft(uppercase_sequence, left_primer, max_mismatches);
    const int right_index = ont::seq::FindPrimerFromRight(uppercase_sequence, right_primer, max_mismatches);
    if (left_index < 0 || right_index < 0) {
        return false;
    }

    const std::size_t left_start_index = static_cast<std::size_t>(left_index);
    const std::size_t right_start_index = static_cast<std::size_t>(right_index);
    if (left_start_index + left_primer.size() > right_start_index) {
        return false;
    }

    *left_start = left_start_index;
    *right_start = right_start_index;
    return true;
}

// Forward-orientation path: read goes forward-primer ... insert ... RC(reverse-primer).
// On success, rewrites `record` in place: header gets the dual UMI tag, and
// sequence/quality are trimmed to the insert (first base after the forward
// UMI block through the last base before the reverse UMI block).
bool ExtractForwardRead(
    ont::fastq::Record* record,
    const std::string& uppercase_sequence,
    const PrimerLayout& layout,
    const int max_mismatches) {
    std::size_t left_start = 0;
    std::size_t right_start = 0;
    if (!TryFindPrimerSpan(
            uppercase_sequence,
            layout.forward_primer,
            layout.reverse_primer_reverse_complement,
            max_mismatches,
            &left_start,
            &right_start)) {
        return false;
    }

    const std::string left_primer_region =
        uppercase_sequence.substr(left_start, layout.forward_primer.size());
    const std::string right_primer_region =
        uppercase_sequence.substr(right_start, layout.reverse_primer_reverse_complement.size());

    const std::string forward_umi =
        ont::seq::ExtractIndexedBases(left_primer_region, layout.forward_umi_positions);
    const std::string reverse_umi = ont::seq::ExtractIndexedBases(
        right_primer_region, layout.reverse_reverse_complement_umi_positions);

    // Trim boundaries: drop everything up to and including the last UMI base of
    // the forward primer; stop just before the first UMI base of the RC'd
    // reverse primer. `.back()+1U` and `.front()` give half-open bounds on the
    // insert. A false return means the primers landed in an impossible order.
    const std::size_t trimmed_start =
        left_start + static_cast<std::size_t>(layout.forward_umi_positions.back()) + 1U;
    const std::size_t trimmed_end =
        right_start + static_cast<std::size_t>(layout.reverse_reverse_complement_umi_positions.front());
    if (trimmed_end < trimmed_start) {
        return false;
    }

    record->header += ":UMI_" + forward_umi + "_" + reverse_umi;
    record->sequence = record->sequence.substr(trimmed_start, trimmed_end - trimmed_start);
    record->quality = record->quality.substr(trimmed_start, trimmed_end - trimmed_start);
    return true;
}

// Reverse-orientation path: read goes reverse-primer ... insert ... RC(forward-primer).
// The UMIs are extracted in read space, then reverse-complemented so the
// ":UMI_<fwd>_<rev>" tag is always written in the canonical forward frame.
// The payload is also reverse-complemented (with quality reversed) before
// being trimmed, so downstream consumers see a forward-oriented insert.
bool ExtractReverseRead(
    ont::fastq::Record* record,
    const std::string& uppercase_sequence,
    const PrimerLayout& layout,
    const int max_mismatches) {
    std::size_t left_start = 0;
    std::size_t right_start = 0;
    if (!TryFindPrimerSpan(
            uppercase_sequence,
            layout.reverse_primer,
            layout.forward_primer_reverse_complement,
            max_mismatches,
            &left_start,
            &right_start)) {
        return false;
    }

    const std::string left_primer_region =
        uppercase_sequence.substr(left_start, layout.reverse_primer.size());
    const std::string right_primer_region =
        uppercase_sequence.substr(right_start, layout.forward_primer.size());

    const std::string left_umi =
        ont::seq::ExtractIndexedBases(left_primer_region, layout.reverse_umi_positions);
    const std::string right_umi = ont::seq::ExtractIndexedBases(
        right_primer_region, layout.forward_reverse_complement_umi_positions);

    // Canonicalise to forward orientation: the "right" UMI in read space was
    // actually derived from the forward primer (RC'd), so reverse-complement
    // it to recover the true forward UMI; mirror for the reverse UMI.
    const std::string forward_umi = ont::seq::ReverseComplement(right_umi);
    const std::string reverse_umi = ont::seq::ReverseComplement(left_umi);

    const std::string reverse_complement_sequence = ont::seq::ReverseComplement(record->sequence);
    const std::string reverse_quality = ont::seq::ReverseQuality(record->quality);

    // Map read-space primer coordinates into RC-space coordinates so the same
    // insert-trimming offsets used in the forward path apply here. The length
    // subtraction flips a [start, start+primer_size) range around the sequence.
    const std::size_t reverse_left_start =
        reverse_complement_sequence.size() - (right_start + layout.forward_primer.size());
    const std::size_t reverse_right_start =
        reverse_complement_sequence.size() - (left_start + layout.reverse_primer.size());

    const std::size_t trimmed_start =
        reverse_left_start + static_cast<std::size_t>(layout.forward_umi_positions.back()) + 1U;
    const std::size_t trimmed_end =
        reverse_right_start + static_cast<std::size_t>(layout.reverse_reverse_complement_umi_positions.front());
    if (trimmed_end < trimmed_start) {
        return false;
    }

    record->header += ":UMI_" + forward_umi + "_" + reverse_umi;
    record->sequence = reverse_complement_sequence.substr(trimmed_start, trimmed_end - trimmed_start);
    record->quality = reverse_quality.substr(trimmed_start, trimmed_end - trimmed_start);
    return true;
}

// Try forward orientation first, then reverse. Returns false if neither path
// can locate a valid primer pair within the allowed mismatch budget, in which
// case the caller drops the read.
bool ExtractAndAnnotateRead(
    ont::fastq::Record* record,
    const PrimerLayout& layout,
    const int max_mismatches) {
    const std::string uppercase_sequence = ont::seq::ToUpperCopy(record->sequence);
    if (ExtractForwardRead(record, uppercase_sequence, layout, max_mismatches)) {
        return true;
    }
    return ExtractReverseRead(record, uppercase_sequence, layout, max_mismatches);
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc != 6) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string input_fastq = argv[1];
    const std::string output_fastq = argv[2];

    // Precompute orientations and UMI index lists. AmbiguousPositions returns
    // the 0-based offsets of every non-ACGT base (N, K, etc.); those positions
    // are where Illumina-style degenerate bases were incorporated during PCR
    // and therefore carry the per-molecule UMI content.
    PrimerLayout layout;
    layout.forward_primer = ont::seq::ToUpperCopy(argv[3]);
    layout.reverse_primer = ont::seq::ToUpperCopy(argv[4]);
    layout.forward_primer_reverse_complement = ont::seq::ReverseComplement(layout.forward_primer);
    layout.reverse_primer_reverse_complement = ont::seq::ReverseComplement(layout.reverse_primer);
    layout.forward_umi_positions = ont::seq::AmbiguousPositions(layout.forward_primer);
    layout.reverse_umi_positions = ont::seq::AmbiguousPositions(layout.reverse_primer);
    layout.forward_reverse_complement_umi_positions =
        ont::seq::AmbiguousPositions(layout.forward_primer_reverse_complement);
    layout.reverse_reverse_complement_umi_positions =
        ont::seq::AmbiguousPositions(layout.reverse_primer_reverse_complement);

    // A primer with zero ambiguous bases carries no UMI content; there is
    // nothing to extract and downstream tools expect a populated dual UMI.
    if (layout.forward_umi_positions.empty() ||
        layout.reverse_umi_positions.empty() ||
        layout.forward_reverse_complement_umi_positions.empty() ||
        layout.reverse_reverse_complement_umi_positions.empty()) {
        std::cerr << "Error: both primer sequences must include at least one ambiguous UMI position.\n";
        return 1;
    }

    int max_mismatches = 0;
    try {
        max_mismatches = std::stoi(argv[5]);
    } catch (const std::exception&) {
        std::cerr << "Error: max_mismatch must be an integer.\n";
        return 1;
    }
    if (max_mismatches < 0) {
        std::cerr << "Error: max_mismatch must be non-negative.\n";
        return 1;
    }

    std::uint64_t total_reads = 0;
    std::uint64_t kept_reads = 0;
    std::uint64_t dropped_reads = 0;

    // Stream the input FASTQ record by record: avoid mutating the source
    // record so the ExtractAndAnnotateRead failure path leaves caller state
    // untouched; only annotated copies for successful reads reach the writer.
    try {
        ont::fastq::Reader reader(input_fastq);
        ont::fastq::Writer writer(output_fastq);

        ont::fastq::Record record;
        while (reader.Next(record)) {
            ++total_reads;
            ont::fastq::Record annotated_record = record;
            if (ExtractAndAnnotateRead(&annotated_record, layout, max_mismatches)) {
                writer.Write(annotated_record);
                ++kept_reads;
            } else {
                ++dropped_reads;
            }
        }
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    // Single-line stderr summary. Consumed by batch wrapper scripts to report
    // the retention rate of the UMI extraction step.
    std::cerr << "[UmiExtractor] total=" << total_reads
              << " kept=" << kept_reads
              << " dropped=" << dropped_reads
              << " max_mismatch=" << max_mismatches << "\n";
    return 0;
}
