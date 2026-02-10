#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>

using namespace std;

static inline void trim_cr(string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

static bool read_fastq_record(istream& in, string& header, string& seq, string& plus, string& qual) {
    header.clear(); seq.clear(); plus.clear(); qual.clear();

    if (!std::getline(in, header)) return false;
    if (!std::getline(in, seq))    throw runtime_error("Malformed FASTQ: missing sequence line");
    if (!std::getline(in, plus))   throw runtime_error("Malformed FASTQ: missing '+' line");
    if (!std::getline(in, qual))   throw runtime_error("Malformed FASTQ: missing quality line");

    trim_cr(header); trim_cr(seq); trim_cr(plus); trim_cr(qual);

    if (header.empty() || header[0] != '@')
        throw runtime_error("Malformed FASTQ: header does not start with '@': " + header);
    if (plus.empty() || plus[0] != '+')
        throw runtime_error("Malformed FASTQ: third line does not start with '+': " + plus);
    if (seq.size() != qual.size())
        throw runtime_error("Malformed FASTQ: sequence and quality lengths differ: " +
                            to_string(seq.size()) + " vs " + to_string(qual.size()));
    return true;
}

// Extract the UMI code after ":UMI_" up to whitespace.
// Returns UMI code only (not including ":UMI_"), normalized to uppercase.
static string parse_umi_code(const string& header_line) {
    static constexpr const char* TAG = ":UMI_";
    const size_t tag_pos = header_line.find(TAG);
    if (tag_pos == string::npos) {
        throw runtime_error("Header missing ':UMI_' tag: " + header_line);
    }

    const size_t start = tag_pos + std::char_traits<char>::length(TAG);
    size_t end = header_line.find_first_of(" \t\r\n", start);
    if (end == string::npos) end = header_line.size();

    if (end <= start) {
        throw runtime_error("Empty UMI after ':UMI_' tag in header: " + header_line);
    }

    string umi = header_line.substr(start, end - start);
    for (char& c : umi) c = static_cast<char>(toupper(static_cast<unsigned char>(c)));
    return umi;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        if (argc != 3) {
            cerr << "Usage: " << argv[0] << " <fastq1> <fastq2>\n";
            return 1;
        }

        const string fastq1 = argv[1];
        const string fastq2 = argv[2];

        ifstream in1(fastq1);
        if (!in1) throw runtime_error("Failed to open FASTQ1: " + fastq1);

        unordered_set<string> umis1;
        umis1.reserve(3'000'000);
        umis1.max_load_factor(0.70f);

        string h, s, p, q;
        uint64_t fastq1_reads = 0;

        while (read_fastq_record(in1, h, s, p, q)) {
            ++fastq1_reads;
            umis1.insert(parse_umi_code(h));
        }

        ifstream in2(fastq2);
        if (!in2) throw runtime_error("Failed to open FASTQ2: " + fastq2);

        unordered_set<string> seen2; // to count unique UMIs in FASTQ2
        seen2.reserve(3'000'000);
        seen2.max_load_factor(0.70f);

        uint64_t fastq2_reads = 0;
        uint64_t overlap_reads_in_fastq2 = 0;

        uint64_t fastq2_unique_umis_only = 0; // unique UMIs present only in FASTQ2
        uint64_t overlapping_unique_umis = 0; // unique UMIs in intersection

        while (read_fastq_record(in2, h, s, p, q)) {
            ++fastq2_reads;
            const string umi = parse_umi_code(h);

            const bool in_fastq1 = (umis1.find(umi) != umis1.end());
            if (in_fastq1) ++overlap_reads_in_fastq2;

            // Count unique UMIs in FASTQ2 exactly once
            if (seen2.insert(umi).second) {
                if (in_fastq1) ++overlapping_unique_umis;
                else ++fastq2_unique_umis_only;
            }
        }

        const uint64_t fastq1_unique_umis = umis1.size();
        const uint64_t fastq1_unique_umis_only =
            (fastq1_unique_umis >= overlapping_unique_umis)
                ? (fastq1_unique_umis - overlapping_unique_umis)
                : 0; // defensive, should never happen now

        cout
            << "FASTQ1_Reads=" << fastq1_reads
            << " FASTQ2_Reads=" << fastq2_reads
            << " FASTQ1_UniqueUMIs=" << fastq1_unique_umis
            << " FASTQ2_UniqueUMIsOnly=" << fastq2_unique_umis_only
            << " Overlapping_UniqueUMIs=" << overlapping_unique_umis
            << " FASTQ2_OverlapReads=" << overlap_reads_in_fastq2
            << "\n";

        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}
