/*
Usage: UmiSpliter <DualUMI.fastq/fq> <UMI1.fastq/fq> <UMI2.fastq/fq>
Description: Split dual-UMI encoded reads into two single-UMI datasets.
Build: g++ -O3 -std=c++17 -o UmiSpliter UmiSpliter.cpp
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <unordered_set>

using namespace std;

// ──────────────────────────────────────────────────────────────────────────────
// Processing Tracker
// Note: these count UNIQUE UMI KEYS, not reads.
// distinct_*_set.size() == number of unique keys seen
// duplicate_*_set.size() == number of unique keys seen >= 2 times
// ──────────────────────────────────────────────────────────────────────────────
static unordered_set<string> distinct_dual_umi_set;
static unordered_set<string> duplicate_dual_umi_set;
static unordered_set<string> distinct_umi1_set;
static unordered_set<string> duplicate_umi1_set;
static unordered_set<string> distinct_umi2_set;
static unordered_set<string> duplicate_umi2_set;

// ──────────────────────────────────────────────────────────────────────────────
// FASTQ I/O
// ──────────────────────────────────────────────────────────────────────────────
struct FastqRecord {
    string id;
    string seq;
    string plus;
    string qual;
};

static inline void chomp_cr(string &s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

class FastqReader {
public:
    explicit FastqReader(const string &path) : in_(path) {
        if (!in_.is_open()) throw runtime_error("Failed to open input: " + path);
    }

    bool next(FastqRecord &rec) {
        if (!getline(in_, rec.id)) return false;
        if (!getline(in_, rec.seq) || !getline(in_, rec.plus) || !getline(in_, rec.qual)) {
            throw runtime_error("Truncated FASTQ record (EOF mid-record).");
        }

        chomp_cr(rec.id);
        chomp_cr(rec.seq);
        chomp_cr(rec.plus);
        chomp_cr(rec.qual);

        if (rec.id.empty() || rec.id[0] != '@')
            throw runtime_error("FASTQ header is empty or does not start with '@'.");

        if (rec.plus.empty() || rec.plus[0] != '+')
            throw runtime_error("FASTQ '+' line is missing or does not start with '+'.");

        if (rec.seq.size() != rec.qual.size())
            throw runtime_error("Sequence length does not match quality length.");

        return true;
    }

private:
    ifstream in_;
};

class FastqWriter {
public:
    explicit FastqWriter(const string &path) : out_(path) {
        if (!out_.is_open()) throw runtime_error("Failed to open output: " + path);
    }

    void write(const FastqRecord &rec) {
        out_ << rec.id << '\n'
             << rec.seq << '\n'
             << rec.plus << '\n'
             << rec.qual << '\n';
    }

private:
    ofstream out_;
};

// ──────────────────────────────────────────────────────────────────────────────
// Header parsing
// Expected: <prefix>:UMI_<umi1>_<umi2><suffix>
// suffix = optional whitespace + rest of header
// ──────────────────────────────────────────────────────────────────────────────
static bool parse_dual_umi_header(
    const string &header,
    string &prefix,
    string &umi1,
    string &umi2,
    string &suffix
) {
    const string tag = ":UMI_";
    size_t tag_pos = header.find(tag);
    if (tag_pos == string::npos) return false;

    prefix = header.substr(0, tag_pos); // includes leading '@'

    size_t start = tag_pos + tag.size();
    size_t sep = header.find('_', start);
    if (sep == string::npos) return false;

    umi1 = header.substr(start, sep - start);
    if (umi1.empty()) return false;

    size_t umi2_start = sep + 1;
    size_t end_umi2 = header.find_first_of(" \t\r\n", umi2_start);
    if (end_umi2 == string::npos) {
        umi2 = header.substr(umi2_start);
        suffix.clear();
    } else {
        umi2 = header.substr(umi2_start, end_umi2 - umi2_start);
        suffix = header.substr(end_umi2); // keep the whitespace + rest
    }

    return !umi2.empty();
}

static bool process_record(const FastqRecord &dual_rec, FastqRecord &umi1_rec, FastqRecord &umi2_rec) {
    string prefix, umi1, umi2, suffix;
    if (!parse_dual_umi_header(dual_rec.id, prefix, umi1, umi2, suffix)) {
        return false;
    }

    // Build output records (copy seq/plus/qual, replace header)
    umi1_rec = dual_rec;
    umi2_rec = dual_rec;

    umi1_rec.id = prefix + ":UMI_" + umi1 + suffix;
    umi2_rec.id = prefix + ":UMI_" + umi2 + suffix;

    // Stats keys
    const string dual_key = umi1 + "_" + umi2;

    const bool dual_first = distinct_dual_umi_set.insert(dual_key).second;
    if (!dual_first) duplicate_dual_umi_set.insert(dual_key);

    const bool u1_first = distinct_umi1_set.insert(umi1).second;
    if (!u1_first) duplicate_umi1_set.insert(umi1);

    const bool u2_first = distinct_umi2_set.insert(umi2).second;
    if (!u2_first) duplicate_umi2_set.insert(umi2);

    return true;
}

// ──────────────────────────────────────────────────────────────────────────────
// main
// ──────────────────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 4) {
        cerr << "Usage: " << argv[0]
             << " <DualUMI.fastq/fq> <UMI1.fastq/fq> <UMI2.fastq/fq>\n";
        return 1;
    }

    const string dual_umi_fastq = argv[1];
    const string umi1_fastq = argv[2];
    const string umi2_fastq = argv[3];

    uint64_t n_dual_umi_read = 0, n_parsed = 0, n_dropped = 0;

    // Reserve if you expect big files. Adjust as needed.
    distinct_dual_umi_set.reserve(1 << 22);
    duplicate_dual_umi_set.reserve(1 << 22);
    distinct_umi1_set.reserve(1 << 22);
    duplicate_umi1_set.reserve(1 << 22);
    distinct_umi2_set.reserve(1 << 22);
    duplicate_umi2_set.reserve(1 << 22);

    try {
        FastqReader reader(dual_umi_fastq);
        FastqWriter umi1_writer(umi1_fastq);
        FastqWriter umi2_writer(umi2_fastq);

        FastqRecord dual_rec, umi1_rec, umi2_rec;

        while (reader.next(dual_rec)) {
            ++n_dual_umi_read;

            if (process_record(dual_rec, umi1_rec, umi2_rec)) {
                umi1_writer.write(umi1_rec);
                umi2_writer.write(umi2_rec);
                ++n_parsed;
            } else {
                ++n_dropped;
            }
        }

        uint64_t distinct_dual_umi = distinct_dual_umi_set.size();
        uint64_t duplicate_dual_umi = duplicate_dual_umi_set.size();
        uint64_t singleton_dual_umi = distinct_dual_umi - duplicate_dual_umi;

        uint64_t distinct_umi1 = distinct_umi1_set.size();
        uint64_t duplicate_umi1 = duplicate_umi1_set.size();
        uint64_t singleton_umi1 = distinct_umi1 - duplicate_umi1;

        uint64_t distinct_umi2 = distinct_umi2_set.size();
        uint64_t duplicate_umi2 = duplicate_umi2_set.size();
        uint64_t singleton_umi2 = distinct_umi2 - duplicate_umi2;

        cerr << "[UmiSpliter]\n"
             << "input_reads=" << n_dual_umi_read << "\n"
             << "parsed=" << n_parsed << "\n"
             << "dropped=" << n_dropped << "\n"
             << "distinct_dual_umi=" << distinct_dual_umi << "\n"
             << "duplicate_dual_umi=" << duplicate_dual_umi << "\n"
             << "singleton_dual_umi=" << singleton_dual_umi << "\n"
             << "distinct_umi1=" << distinct_umi1 << "\n"
             << "duplicate_umi1=" << duplicate_umi1 << "\n"
             << "singleton_umi1=" << singleton_umi1 << "\n"
             << "distinct_umi2=" << distinct_umi2 << "\n"
             << "duplicate_umi2=" << duplicate_umi2 << "\n"
             << "singleton_umi2=" << singleton_umi2 << endl;
    }
    catch (const exception &e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    return 0;
}
