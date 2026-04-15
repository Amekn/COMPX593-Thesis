// MQAAR.cpp
// Produces executable: MQAAR
// Usage: MQAAR <report.tsv> [--convert]
//
// Pipeline role: read a Dorado/Bonito sequencing_summary TSV and compute the
// arithmetic mean of 'mean_qscore_template'. With --convert the mean Phred
// score is reported as accuracy percentage via (1 - 10^(-Q/10)) * 100, which
// is the form used in the thesis basecaller-comparison tables. MQAAR stands
// for "Mean Q-score to Approximate Accuracy Reporter".
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

// Print CLI contract to stderr.
static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <report.tsv> [--convert]\n\n"
         <<
R"(Usage:
Description:
  Compute the mean of the mean_qscore_template column from a TSV report.

Options:
  --convert    Convert the mean Phred Q-score into an accuracy percentage:
               accuracy = (1 - 10^(-Q/10)) * 100
)";
}

// Resolve the 0-based index of a named TSV column by scanning the header
// row. Returns -1 if no column matches; the caller treats that as a fatal
// configuration error rather than silently picking the wrong column.
static int find_column_index(const string& header, const string& wanted) {
    string token;
    stringstream ss(header);
    int idx = 0;
    while (getline(ss, token, '\t')) {
        if (token == wanted) return idx;
        ++idx;
    }
    return -1;
}

// Extract the value at column_index from a tab-separated line into out.
// Returns false if the line has fewer columns than expected (short row),
// in which case the caller skips it rather than aborting the scan.
static bool read_column_value(const string& line, int column_index, string& out) {
    out.clear();
    string token;
    stringstream ss(line);
    for (int i = 0; getline(ss, token, '\t'); ++i) {
        if (i == column_index) {
            out = token;
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    bool convert = false;
    string report_path;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--convert" || arg == "-c") {
            convert = true;
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else if (report_path.empty()) {
            report_path = std::move(arg);
        } else {
            cerr << "Unexpected argument: " << arg << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    if (report_path.empty()) {
        usage(argv[0]);
        return 1;
    }

    ifstream in(report_path);
    if (!in) {
        cerr << "Error: cannot open '" << report_path << "'\n";
        return 2;
    }

    string header;
    if (!getline(in, header)) {
        cerr << "Error: input file is empty: '" << report_path << "'\n";
        return 2;
    }

    const string column_name = "mean_qscore_template";
    int column_index = find_column_index(header, column_name);
    if (column_index < 0) {
        cerr << "Missing required column: '" << column_name << "'\n";
        return 2;
    }

    // Accumulate in long double to keep precision stable over millions of
    // reads; stod-failed or partially-parsed fields are skipped silently so
    // a handful of malformed rows do not abort the pass.
    long double sum = 0.0;
    uint64_t count = 0;
    string line;
    string field;

    while (getline(in, line)) {
        if (!read_column_value(line, column_index, field)) continue;
        try {
            size_t used = 0;
            double value = stod(field, &used);
            if (used != field.size()) continue;
            sum += value;
            ++count;
        } catch (const exception&) {
            continue;
        }
    }

    if (count == 0) {
        cerr << "No valid numeric entries in '" << column_name << "'\n";
        return 2;
    }

    double mean_q = static_cast<double>(sum / static_cast<long double>(count));
    cout << fixed << setprecision(2);
    if (convert) {
        // Phred-to-accuracy conversion; note this is the naive per-base
        // accuracy derived from the mean Q, not the geometric mean.
        double mean_pct = (1.0 - std::pow(10.0, -mean_q / 10.0)) * 100.0;
        cout << "Mean accuracy percentage based on mean Q-score (" << mean_q << "): "
             << mean_pct << "%\n";
    } else {
        cout << "Mean Q-score of '" << column_name << "': " << mean_q << "\n";
    }

    return 0;
}
