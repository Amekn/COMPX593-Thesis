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
// main
// ──────────────────────────────────────────────────────────────────────────────
int main(int argc, char **argv)
{
    ios::sync_with_stdio(false);

    // Parse CLI args
    if (argc < 5)
    {
        cerr << "Usage: " << argv[0] << " <input.fastq> <output.fastq> <minimum_length> <maximum_length>\n";
        exit(1);
    }
    const string in_fastq = argv[1];
    const string out_fastq = argv[2];
    int minimum_length = stoi(argv[3]);
    int maximum_length = stoi(argv[4]);

    int n_total = 0;
    int n_dropped = 0;
    int n_kept = 0;

    try
    {
        FastqReader reader(in_fastq);
        FastqWriter writer(out_fastq);

        FastqRecord rec;
        while (reader.next(rec))
        {
            ++n_total;
            if(rec.seq.size() > minimum_length && rec.seq.size() < maximum_length){
                writer.write(rec);
                ++n_kept;
            }
            else{
                ++n_dropped;
            }
        }
        std::cerr << "[UMI-filter] total=" << n_total
                  << " kept=" << n_kept
                  << " dropped=" << n_dropped << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
