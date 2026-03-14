#include <cstdint>
#include <iostream>
#include <string>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Sequence.hpp"

namespace {

void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name
        << " <input.fastq> <output.fastq> <forward_primer> <reverse_primer> <max_mismatch>\n\n"
        << "Trim reads to the span bracketed by the primer pair.\n"
        << "Forward-oriented reads are emitted unchanged; reverse-oriented reads are reverse-complemented.\n";
}

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

bool TrimRead(
    ont::fastq::Record* record,
    const std::string& forward_primer,
    const std::string& reverse_primer,
    const std::string& forward_primer_reverse_complement,
    const std::string& reverse_primer_reverse_complement,
    const int max_mismatches) {
    const std::string uppercase_sequence = ont::seq::ToUpperCopy(record->sequence);

    std::size_t left_start = 0;
    std::size_t right_start = 0;

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
