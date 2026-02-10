/*
Usage: UmiOverlap <fastq1> <fastq2>
Description:
    Compute the number of overlapping UMIs between two FASTQ files.
    Output:
        Number of reads in FASTQ1
        Number of reads in FASTQ2
        Number of overlapping UMIs
*/
/*
Usage: Benchmarker <ground_truth_fastq> <test_fastq1> <test_fastq2> ...

Description:
    Compares multiple FASTQ files containing processed nanopore reads against a ground truth illumina
    reads FASTQ file. Outputs precision, recall, and F1-score for each test FASTQ file.

    Each read in the ground_truth_fastq must contain dual UMI tags in the header
    Each read in the test FASTQ file must contain dual UMI tags in the header (same format as ground truth)
    A map is created from the ground truth FASTQ file using the dual UMI tags as keys.

    Each read in the test FASTQ file is looked up in the map using the dual UMI tags.
    If a match is found, the read sequences are compared to determine statistics of correctness.

Important:
    If a ground-truth quality character is '#', that reference position is IGNORED in match/sub/del totals.

Output (TSV) columns for each test FASTQ:
    Test FASTQ filename
    Number of processed reads
    Number of nucleotides compared
    Number of matched nucleotides
    Number of deleted nucleotides
    Number of inserted nucleotides
    Number of substituted nucleotides
    Normalized Matched Nucleotides (matched / total)
    Normalized Deletions (deletions / total)
    Normalized Insertions (insertions / total)
    Normalized Substitutions (substitutions / total)
    Precision
    Recall
    F1-score

Precision/Recall definition (nucleotide-level, alignment-based):
    TP = matched
    FP = substitutions + insertions
    FN = substitutions + deletions
    precision = TP / (TP + FP) = matched / (matched + subs + ins)
    recall    = TP / (TP + FN) = matched / (matched + subs + del)
*/
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;

static inline char up(char c) {
    return static_cast<char>(toupper(static_cast<unsigned char>(c)));
}

static bool read_fastq_record(istream& in, string& header, string& seq, string& plus, string& qual) {
    header.clear(); seq.clear(); plus.clear(); qual.clear();
    if (!std::getline(in, header)) return false;
    if (!std::getline(in, seq))    throw runtime_error("Malformed FASTQ: missing sequence line");
    if (!std::getline(in, plus))   throw runtime_error("Malformed FASTQ: missing '+' line");
    if (!std::getline(in, qual))   throw runtime_error("Malformed FASTQ: missing quality line");
    if (header.empty() || header[0] != '@')
        throw runtime_error("Malformed FASTQ: header does not start with '@': " + header);
    if (plus.empty() || plus[0] != '+')
        throw runtime_error("Malformed FASTQ: third line does not start with '+': " + plus);
    if (seq.size() != qual.size())
        throw runtime_error("Malformed FASTQ: sequence and quality lengths differ: " + to_string(seq.size()) + " vs " + to_string(qual.size()));
    return true;
}

// Parse ":UMI_<FWD>_<REV>" from header (either with or without the leading '@').
static string parse_umi(const string& header_line) {
    size_t start = header_line.find_first_of(":UMI_");
    if (start == string::npos) {
        throw runtime_error("Header missing ':UMI_' tag: " + header_line);
    }
    // The UMI segment ends at whitespace
    size_t end = header_line.find_first_of(" \t\r\n", start);
    if (end == string::npos) end = header_line.size();
    string umi_part = header_line.substr(start, end - start);
    return umi_part;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        if (argc < 3) {
            cerr << "Usage: " << argv[0] << " <fastq1> <fastq2>\n";
            return 1;
        }

        string fastq1 = argv[1];
        ifstream in(fastq1);

        unordered_set<string> fastq1_umis;
        fastq1_umis.reserve(3000000);
        fastq1_umis.max_load_factor(0.70f);

        string h, s, p, q;
        while (read_fastq_record(in, h, s, p, q)){
            string key = parse_umi(h);
            fastq1_umis.insert(key);
        }

        string fastq2 = argv[2];
        ifstream in2(fastq2);
        uint64_t fastq1_umi_count = fastq1_umis.size();
        uint64_t fastq2_umi_count = 0;
        uint64_t overlap_count = 0;

        while (read_fastq_record(in2, h, s, p, q)){
            string key = parse_umi(h);
            fastq2_umi_count++;
            if (fastq1_umis.find(key) != fastq1_umis.end()){
                overlap_count++;
            }
        }

        cout << "FASTQ1_reads=" << fastq1_umi_count
                << " FASTQ2_reads=" << fastq2_umi_count
                << " Overlapping_UMIs=" << overlap_count
                << "\n";
        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}
