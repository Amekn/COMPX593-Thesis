// PrimerTrimmer.cpp
// Produces executable: PrimerTrimmer
// Usage: PrimerTrimmer <input.fastq> <output.fastq> <forward_primer> <reverse_primer> <max_mismatch>
//
// Pipeline role: locate the forward primer (or its RC) at the 5' end and
// the reverse primer (or its RC) at the 3' end of each read, then emit the
// bracketed span. Reads whose orientation is reverse get reverse-complemented
// (and their quality string reversed) so the output is all in forward-read
// orientation relative to the IgA Fc reference before downstream alignment.
#include <cstdint>
#include <iostream>
#include <string>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Sequence.hpp"

namespace {

// Print CLI contract to stderr when args are invalid.
void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name
        << " <input.fastq> <output.fastq> <forward_primer> <reverse_primer> <max_mismatch>\n\n"
        << "Trim reads to the span bracketed by the primer pair.\n"
        << "Forward-oriented reads are emitted unchanged; reverse-oriented reads are reverse-complemented.\n";
}

// Anchor a primer pair against an uppercased sequence and return their
// start offsets via out-params. Returns true only when both primers hit
// and the left match actually ends before the right match begins (so the
// bracketed span is non-overlapping and ordered).
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
    // Reject overlap: the left primer must finish strictly before the right
    // primer starts, otherwise the "span" is a chance double-match.
    if (left_start_index + left_primer.size() > right_start_index) {
        return false;
    }

    *left_start = left_start_index;
    *right_start = right_start_index;
    return true;
}

// Trim a record in place to the span bracketed by a matching primer pair.
// Tries the forward orientation first (fwd...rc(rev)); on miss, falls back
// to the reverse orientation (rev...rc(fwd)) and reverse-complements the
// kept span so output is always in forward orientation. Returns false if
// neither orientation yields a valid span.
bool TrimRead(
    ont::fastq::Record* record,
    const std::string& forward_primer,
    const std::string& reverse_primer,
    const std::string& forward_primer_reverse_complement,
    const std::string& reverse_primer_reverse_complement,
    const int max_mismatches) {
    // Uppercase once: primer search uses IUPAC masks that are
    // case-insensitive after this copy.
    const std::string uppercase_sequence = ont::seq::ToUpperCopy(record->sequence);

    std::size_t left_start = 0;
    std::size_t right_start = 0;

    // Forward orientation: fwd primer at 5', reverse-complement of reverse
    // primer at 3'. Span length includes both primer bases so the trimmed
    // read still carries primer flanks for downstream alignment.
    if (TryFindPrimerSpan(
            uppercase_sequence,
            forward_primer,
            reverse_primer_reverse_complement,
            max_mismatches,
            &left_start,
            &right_start)) {
        const std::size_t span_length =
            right_start + reverse_primer_reverse_complement.size() - left_start;
        record->sequence = record->sequence.substr(left_start, span_length);
        record->quality = record->quality.substr(left_start, span_length);
        return true;
    }

    // Reverse orientation: the read was sequenced from the opposite strand,
    // so we reverse-complement the kept span and reverse the quality string
    // to keep sequence and quality indices aligned.
    if (TryFindPrimerSpan(
            uppercase_sequence,
            reverse_primer,
            forward_primer_reverse_complement,
            max_mismatches,
            &left_start,
            &right_start)) {
        const std::size_t span_length =
            right_start + forward_primer_reverse_complement.size() - left_start;
        record->sequence =
            ont::seq::ReverseComplement(record->sequence.substr(left_start, span_length));
        record->quality = ont::seq::ReverseQuality(record->quality.substr(left_start, span_length));
        return true;
    }

    return false;
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
    std::string forward_primer = ont::seq::ToUpperCopy(argv[3]);
    std::string reverse_primer = ont::seq::ToUpperCopy(argv[4]);

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

    // Pre-compute the primer reverse complements once: they are checked
    // against every read, so recomputing per-read would be wasteful.
    const std::string forward_primer_reverse_complement =
        ont::seq::ReverseComplement(forward_primer);
    const std::string reverse_primer_reverse_complement =
        ont::seq::ReverseComplement(reverse_primer);

    std::uint64_t total_reads = 0;
    std::uint64_t kept_reads = 0;
    std::uint64_t dropped_reads = 0;

    try {
        ont::fastq::Reader reader(input_fastq);
        ont::fastq::Writer writer(output_fastq);

        ont::fastq::Record record;
        while (reader.Next(record)) {
            ++total_reads;
            ont::fastq::Record trimmed_record = record;
            if (TrimRead(
                    &trimmed_record,
                    forward_primer,
                    reverse_primer,
                    forward_primer_reverse_complement,
                    reverse_primer_reverse_complement,
                    max_mismatches)) {
                writer.Write(trimmed_record);
                ++kept_reads;
            } else {
                ++dropped_reads;
            }
        }
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    std::cerr << "[PrimerTrimmer] total=" << total_reads
              << " kept=" << kept_reads
              << " dropped=" << dropped_reads << "\n";
    return 0;
}
