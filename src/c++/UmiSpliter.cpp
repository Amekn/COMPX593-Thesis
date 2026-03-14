#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "ont_tools/Fastq.hpp"

namespace {

struct ParsedDualHeader {
    std::string prefix;
    std::string first_umi;
    std::string second_umi;
    std::string suffix;
};

struct UmiSplitStatistics {
    std::uint64_t total_reads = 0;
    std::uint64_t parsed_reads = 0;
    std::uint64_t dropped_reads = 0;

    std::unordered_set<std::string> distinct_dual_keys;
    std::unordered_set<std::string> duplicated_dual_keys;
    std::unordered_set<std::string> distinct_first_umis;
    std::unordered_set<std::string> duplicated_first_umis;
    std::unordered_set<std::string> distinct_second_umis;
    std::unordered_set<std::string> duplicated_second_umis;
};

void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " <dual_umi.fastq> <umi1.fastq> <umi2.fastq>\n\n"
        << "Split dual-UMI FASTQ records into two single-UMI FASTQ files.\n";
}

bool ParseDualUmiHeader(const std::string& header, ParsedDualHeader* parsed_header) {
    constexpr char kTag[] = ":UMI_";
    const std::size_t tag_position = header.find(kTag);
    if (tag_position == std::string::npos) {
        return false;
    }

    const std::size_t first_umi_begin = tag_position + sizeof(kTag) - 1U;
    const std::size_t separator = header.find('_', first_umi_begin);
    if (separator == std::string::npos) {
        return false;
    }

    const std::size_t second_umi_begin = separator + 1U;
    std::size_t second_umi_end = header.find_first_of(" \t\r\n", second_umi_begin);
    if (second_umi_end == std::string::npos) {
        second_umi_end = header.size();
        parsed_header->suffix.clear();
    } else {
        parsed_header->suffix = header.substr(second_umi_end);
    }

    parsed_header->prefix = header.substr(0, tag_position);
    parsed_header->first_umi = header.substr(first_umi_begin, separator - first_umi_begin);
    parsed_header->second_umi = header.substr(second_umi_begin, second_umi_end - second_umi_begin);
    return !parsed_header->first_umi.empty() && !parsed_header->second_umi.empty();
}

void TrackUniqueness(UmiSplitStatistics* statistics, const ParsedDualHeader& parsed_header) {
    const std::string dual_key = parsed_header.first_umi + "_" + parsed_header.second_umi;
    if (!statistics->distinct_dual_keys.insert(dual_key).second) {
        statistics->duplicated_dual_keys.insert(dual_key);
    }

    if (!statistics->distinct_first_umis.insert(parsed_header.first_umi).second) {
        statistics->duplicated_first_umis.insert(parsed_header.first_umi);
    }

    if (!statistics->distinct_second_umis.insert(parsed_header.second_umi).second) {
        statistics->duplicated_second_umis.insert(parsed_header.second_umi);
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc != 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    UmiSplitStatistics statistics;
    statistics.distinct_dual_keys.reserve(1U << 22);
    statistics.duplicated_dual_keys.reserve(1U << 22);
    statistics.distinct_first_umis.reserve(1U << 22);
    statistics.duplicated_first_umis.reserve(1U << 22);
    statistics.distinct_second_umis.reserve(1U << 22);
    statistics.duplicated_second_umis.reserve(1U << 22);

    try {
        ont::fastq::Reader reader(argv[1]);
        ont::fastq::Writer first_writer(argv[2]);
        ont::fastq::Writer second_writer(argv[3]);

        ont::fastq::Record input_record;
        while (reader.Next(input_record)) {
            ++statistics.total_reads;

            ParsedDualHeader parsed_header;
            if (!ParseDualUmiHeader(input_record.header, &parsed_header)) {
                ++statistics.dropped_reads;
                continue;
            }

            TrackUniqueness(&statistics, parsed_header);

            ont::fastq::Record first_record = input_record;
            ont::fastq::Record second_record = input_record;
            first_record.header =
                parsed_header.prefix + ":UMI_" + parsed_header.first_umi + parsed_header.suffix;
            second_record.header =
                parsed_header.prefix + ":UMI_" + parsed_header.second_umi + parsed_header.suffix;

            first_writer.Write(first_record);
            second_writer.Write(second_record);
            ++statistics.parsed_reads;
        }
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    const std::uint64_t singleton_dual_keys =
        statistics.distinct_dual_keys.size() - statistics.duplicated_dual_keys.size();
    const std::uint64_t singleton_first_umis =
        statistics.distinct_first_umis.size() - statistics.duplicated_first_umis.size();
    const std::uint64_t singleton_second_umis =
        statistics.distinct_second_umis.size() - statistics.duplicated_second_umis.size();

    std::cerr << "[UmiSpliter]\n"
              << "input_reads=" << statistics.total_reads << "\n"
              << "parsed_reads=" << statistics.parsed_reads << "\n"
              << "dropped_reads=" << statistics.dropped_reads << "\n"
              << "distinct_dual_umi=" << statistics.distinct_dual_keys.size() << "\n"
              << "duplicate_dual_umi=" << statistics.duplicated_dual_keys.size() << "\n"
              << "singleton_dual_umi=" << singleton_dual_keys << "\n"
              << "distinct_umi1=" << statistics.distinct_first_umis.size() << "\n"
              << "duplicate_umi1=" << statistics.duplicated_first_umis.size() << "\n"
              << "singleton_umi1=" << singleton_first_umis << "\n"
              << "distinct_umi2=" << statistics.distinct_second_umis.size() << "\n"
              << "duplicate_umi2=" << statistics.duplicated_second_umis.size() << "\n"
              << "singleton_umi2=" << singleton_second_umis << "\n";
    return 0;
}
