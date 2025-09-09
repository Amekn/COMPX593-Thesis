#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cctype>
#include <cstdint>

using namespace std;

// ──────────────────────────────────────────────────────────────────────────────
// Utility helpers
// ──────────────────────────────────────────────────────────────────────────────
static uint8_t IUPAC_bits[256];
static char IUPAC_comp[256];

static inline char up(char c)
{
    return static_cast<char>(toupper(static_cast<unsigned char>(c)));
}


static inline uint8_t base_bit(char b)
{
    switch (b)
    {
    case 'A':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 4;
    case 'T':
        return 8;
    default:
        return 0; // non-ACGT in read -> mismatch unless primer allows none
    }
}

static string to_rc(const string &s)
{
    const size_t n = s.size();
    string out;
    out.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        char c = up(s[n - 1 - i]);
        char cc = IUPAC_comp[static_cast<unsigned char>(c)];
        out[i] = cc;
    }
    return out;
}

// ──────────────────────────────────────────────────────────────────────────────
// FASTQ I/O
// ──────────────────────────────────────────────────────────────────────────────
struct FastqRecord
{
    string id;
    string seq;
    string plus;
    string qual;
};

class FastqReader
{
public:
    explicit FastqReader(const string &path) : in_(path)
    {
        if (!in_.is_open())
            throw std::runtime_error("Failed to open input: " + path);
    }
    bool next(FastqRecord &rec)
    {
        string l1, l2, l3, l4;
        if (!std::getline(in_, l1))
            return false;
        if (!std::getline(in_, l2))
            return false;
        if (!std::getline(in_, l3))
            return false;
        if (!std::getline(in_, l4))
            return false;
        if (l1.empty() || l1[0] != '@')
            return false;
        rec.id = l1;
        rec.seq = l2;
        rec.plus = l3;
        rec.qual = l4;
        return true;
    }

private:
    std::ifstream in_;
};

class FastqWriter
{
public:
    explicit FastqWriter(const string &path) : out_(path)
    {
        if (!out_.is_open())
            throw std::runtime_error("Failed to open output: " + path);
    }
    void write(const FastqRecord &rec)
    {
        out_ << rec.id << '\n'
             << rec.seq << '\n'
             << rec.plus << '\n'
             << rec.qual << '\n';
    }

private:
    std::ofstream out_;
};

// ──────────────────────────────────────────────────────────────────────────────
// Processors
// ──────────────────────────────────────────────────────────────────────────────
// Initiate IUPAC bitmasks and complements
static void init_iupac()
{
    auto set_code = [](char c, const char *allowed)
    {
        uint8_t bits = 0;
        for (const char *p = allowed; *p; ++p)
            bits |= base_bit(*p);
        IUPAC_bits[static_cast<unsigned char>(c)] = bits;
    };

    // Zero-init then set
    for (int i = 0; i < 256; ++i)
    {
        IUPAC_bits[i] = 0;
        IUPAC_comp[i] = '?';
    }

    set_code('A', "A");
    set_code('C', "C");
    set_code('G', "G");
    set_code('T', "T");
    set_code('R', "AG");
    set_code('Y', "CT");
    set_code('S', "GC");
    set_code('W', "AT");
    set_code('K', "GT");
    set_code('M', "AC");
    set_code('B', "CGT");
    set_code('D', "AGT");
    set_code('H', "ACT");
    set_code('V', "ACG");
    set_code('N', "ACGT");

    IUPAC_comp[(unsigned char)'A'] = 'T';
    IUPAC_comp[(unsigned char)'C'] = 'G';
    IUPAC_comp[(unsigned char)'G'] = 'C';
    IUPAC_comp[(unsigned char)'T'] = 'A';
    IUPAC_comp[(unsigned char)'R'] = 'Y';
    IUPAC_comp[(unsigned char)'Y'] = 'R';
    IUPAC_comp[(unsigned char)'S'] = 'S';
    IUPAC_comp[(unsigned char)'W'] = 'W';
    IUPAC_comp[(unsigned char)'K'] = 'M';
    IUPAC_comp[(unsigned char)'M'] = 'K';
    IUPAC_comp[(unsigned char)'B'] = 'V';
    IUPAC_comp[(unsigned char)'D'] = 'H';
    IUPAC_comp[(unsigned char)'H'] = 'D';
    IUPAC_comp[(unsigned char)'V'] = 'B';
    IUPAC_comp[(unsigned char)'N'] = 'N';
}



// Any non-ACGT in primer marks a UMI position
static vector<int> umi_mask(const string &primer)
{
    vector<int> pos;
    for (int i = 0; i < (int)primer.size(); ++i)
    {
        char c = up(primer[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            pos.push_back(i);
    }
    return pos;
}

// IUPAC-aware match with <= max_mismatch
static inline bool iupac_match_leq_mismatches(const string &read_seq, int start,
                                              const string &primer, int max_mismatch)
{
    int mm = 0;
    const int L = (int)primer.size();
    for (int i = 0; i < L; ++i)
    {
        char r = up(read_seq[start + i]);
        char p = up(primer[i]);
        uint8_t allowed = IUPAC_bits[static_cast<unsigned char>(p)];
        uint8_t rb = base_bit(r);
        if ((rb & allowed) == 0)
        {
            if (++mm > max_mismatch)
                return false;
        }
    }
    return true;
}

// Slide left->right: first match index, else -1
static int find_from_start(const string &read_seq, const string &primer, int max_mismatch)
{
    int L = (int)primer.size();
    int R = (int)read_seq.size();
    if (L > R)
        return -1;
    for (int i = 0; i <= R - L; ++i)
    {
        if (iupac_match_leq_mismatches(read_seq, i, primer, max_mismatch))
            return i;
    }
    return -1;
}

// Slide right->left: last match index, else -1
static int find_from_end(const string &read_seq, const string &primer, int max_mismatch)
{
    int L = (int)primer.size();
    int R = (int)read_seq.size();
    if (L > R)
        return -1;
    for (int i = R - L; i >= 0; --i)
    {
        if (iupac_match_leq_mismatches(read_seq, i, primer, max_mismatch))
            return i;
    }
    return -1;
}

// Extract UMI bases from region sequence given UMI positions
static string extract_umi_from_region(const string &region_seq, const vector<int> &umi_pos)
{
    string out;
    out.reserve(umi_pos.size());
    for (int pos : umi_pos)
    {
        if (pos >= 0 && pos < (int)region_seq.size())
            out.push_back(region_seq[pos]);
    }
    return out;
}

// Core read processing function
static bool process_record(
    FastqRecord &rec,
    const string &fwd, const string &rev,
    const string &fwd_rc, const string &rev_rc,
    int max_mismatch,
    const vector<int> &f_umi_pos,
    const vector<int> &r_umi_pos,
    const vector<int> &frc_umi_pos,
    const vector<int> &rrc_umi_pos)
{
    string seq = rec.seq;
    for (char &c : seq)
        c = up(c);

    const int Lf = (int)fwd.size();
    const int Lr = (int)rev.size();

    // Orientation A: F ... rc(R)
    {
        int left_idx = find_from_start(seq, fwd, max_mismatch);
        int right_idx = find_from_end(seq, rev_rc, max_mismatch);
        if (left_idx >= 0 && right_idx >= 0 && left_idx + Lf <= right_idx)
        {
            string start_region = seq.substr(left_idx, Lf);
            string end_region = seq.substr(right_idx, Lr);
            string umi1 = extract_umi_from_region(start_region, f_umi_pos);
            string umi2 = extract_umi_from_region(end_region, rrc_umi_pos);
            string umi = umi1 + "-" + umi2;
            rec.id += "|UMI:" + umi + "|STR:FWD|ORIENT:F-rcR";
            int left_end = left_idx + f_umi_pos.front();
            int extend_length = right_idx - left_end + rrc_umi_pos.back() + 1;
            rec.seq = rec.seq.substr(left_end, extend_length);
            rec.qual = rec.qual.substr(left_end, extend_length);
            return true;
        }
    }
    // Orientation B: R ... rc(F)
    {
        int left_idx = find_from_start(seq, rev, max_mismatch);
        int right_idx = find_from_end(seq, fwd_rc, max_mismatch);
        if (left_idx >= 0 && right_idx >= 0 && left_idx + Lr <= right_idx)
        {
            string start_region = seq.substr(left_idx, Lr);
            string end_region = seq.substr(right_idx, Lf);
            string umi1 = extract_umi_from_region(start_region, r_umi_pos);
            string umi2 = extract_umi_from_region(end_region, frc_umi_pos);
            string umi = umi1 + "-" + umi2;
            rec.id += "|UMI:" + umi + "|STR:REV|ORIENT:R-rcF";
            int left_end = left_idx + r_umi_pos.front();
            int extend_length = right_idx - left_end + frc_umi_pos.back() + 1;
            rec.seq = rec.seq.substr(left_end, extend_length);
            rec.qual = rec.qual.substr(left_end, extend_length);
            return true;
        }
    }
    return false;
}

// ──────────────────────────────────────────────────────────────────────────────
// main
// ──────────────────────────────────────────────────────────────────────────────
int main(int argc, char **argv)
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Parse CLI args
    if (argc < 5)
    {
        cerr << "Usage: " << argv[0] << " <input.fastq> <output.fastq> <forward_primer> <reverse_primer> [--max-mismatch N]\n";
        exit(1);
    }
    const string in_fastq = argv[1];
    const string out_fastq = argv[2];
    string fwd = argv[3];
    string rev = argv[4];
    int max_mismatch = 0; // default
    if (argc > 5)
    {
        if (argv[5] == string("--max-mismatch") && argc == 7)
        {
            max_mismatch = stoi(argv[6]);
        }
        else
        {
            cerr << "Usage: " << argv[0] << " <input.fastq> <output.fastq> <forward_primer> <reverse_primer> [--max-mismatch N]\n";
            exit(1);
        }
    }

    // Uppercase primers
    transform(fwd.begin(), fwd.end(), fwd.begin(), up);
    transform(rev.begin(), rev.end(), rev.begin(), up);

    // Build IUPAC bit masks.
    init_iupac();

    // Precompute reverse-complements and UMI positions
    const string fwd_rc = to_rc(fwd);
    const string rev_rc = to_rc(rev);

    const vector<int> f_umi_pos = umi_mask(fwd);
    const vector<int> r_umi_pos = umi_mask(rev);
    const vector<int> frc_umi_pos = umi_mask(fwd_rc);
    const vector<int> rrc_umi_pos = umi_mask(rev_rc);

    uint64_t n_total = 0, n_kept = 0, n_dropped = 0;

    try
    {
        FastqReader reader(in_fastq);
        FastqWriter writer(out_fastq);

        FastqRecord rec;
        while (reader.next(rec))
        {
            ++n_total;
            if (rec.seq.size() != rec.qual.size())
            {
                ++n_dropped;
                continue;
            }
            FastqRecord out = rec;
            if (process_record(out, fwd, rev, fwd_rc, rev_rc, max_mismatch,
                               f_umi_pos, r_umi_pos, frc_umi_pos, rrc_umi_pos))
            {
                writer.write(out);
                ++n_kept;
            }
            else
            {
                ++n_dropped;
            }
        }

        std::cerr << "[UMI-filter] total=" << n_total
                  << " kept=" << n_kept
                  << " dropped=" << n_dropped
                  << " (mismatches per primer <= " << max_mismatch << ")\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
