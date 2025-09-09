#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <unordered_set>

using namespace std;

static inline char up(char c) {
    return static_cast<char>(toupper(static_cast<unsigned char>(c)));
}

static string extract_umi(const string& header) {
    const string tag = "|UMI:";
    size_t pos = header.find(tag);
    if (pos == string::npos) return string();
    size_t start = pos + tag.size();
    size_t end_bar = header.find('|', start);
    string umi = (end_bar == string::npos) ? header.substr(start) : header.substr(start, end_bar - start);
    for (char& c : umi) c = up(c);
    return umi;
}

struct FastqRecord {
    string h, s, p, q; // header, seq, plus, qual
};

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " <input.fastq> <output.fastq>\n";
        return 1;
    }
    const string input_fastq  = argv[1];
    const string output_fastq = argv[2];

    ifstream in(input_fastq);
    if (!in) {
        cerr << "ERROR: cannot open input: " << input_fastq << "\n";
        return 2;
    }
    ofstream out(output_fastq);
    if (!out) {
        cerr << "ERROR: cannot open output: " << output_fastq << "\n";
        return 2;
    }

    // Reserve a big-ish set to reduce rehashing if you expect many UMIs.
    unordered_set<string> seen;
    seen.reserve(1 << 17); // ~100k buckets; adjust to your scale
    int total = 0, kept = 0, dropped=0;

    FastqRecord r;
    r.h.reserve(256); r.s.reserve(1024); r.p.reserve(4); r.q.reserve(1024);

    while (true) {
        if (!getline(in, r.h)) break;            // EOF ok here
        if (!getline(in, r.s) || !getline(in, r.p) || !getline(in, r.q)) {
            cerr << "ERROR: truncated FASTQ near record " << (total + 1) << "\n";
            return 3;
        }
        ++total;

        // Minimal FASTQ sanity check
        if (r.h.empty() || r.h[0] != '@' || r.p.empty() || r.p[0] != '+') {
            cerr << "WARN: malformed FASTQ record at " << total << "\n";
        }
        if (r.s.size() != r.q.size()) {
            continue;
        }

        string umi = extract_umi(r.h);
        if (umi.empty()) {
            ++dropped;
            continue;
        }

        auto [it, inserted] = seen.insert(move(umi));
        if (inserted) {
            // First time seeing this UMI => keep
            out << r.h << '\n' << r.s << '\n' << r.p << '\n' << r.q << '\n';
            ++kept;
        } else {
            // Duplicate UMI => drop
            ++dropped;
        }
    }

    cerr << "[umi_dedupe] total=" << total
         << " kept=" << kept
         << " dropped_dupe=" << dropped
         << " kept_no_umi=" << kept << "\n";
    return 0;
}