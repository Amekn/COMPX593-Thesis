#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

struct Row {
    string metric;
    string value;
    string comment;
};

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " [stats_file]\n\n"
         <<
R"(Usage:
Description:
  Parse samtools-style .stats input and print the SN lines in a formatted table.
  If no file is provided, input is read from stdin.
)";
}

static vector<Row> parse_rows(istream& in) {
    vector<Row> rows;
    string line;
    while (getline(in, line)) {
        if (line.rfind("SN\t", 0) != 0) continue;

        vector<string> fields;
        string token;
        size_t start = 0;
        while (true) {
            size_t pos = line.find('\t', start);
            if (pos == string::npos) {
                fields.push_back(line.substr(start));
                break;
            }
            fields.push_back(line.substr(start, pos - start));
            start = pos + 1;
        }

        if (fields.size() < 3) continue;
        Row row;
        row.metric = fields[1];
        row.value = fields[2];
        if (fields.size() > 3) {
            row.comment = fields[3];
            if (!row.comment.empty() && row.comment[0] == '#') {
                row.comment.erase(0, 1);
                while (!row.comment.empty() && row.comment[0] == ' ') row.comment.erase(0, 1);
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc > 2) {
        usage(argv[0]);
        return 1;
    }
    if (argc == 2) {
        string arg = argv[1];
        if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        }
    }

    istream* in = &cin;
    ifstream file;
    if (argc == 2) {
        file.open(argv[1]);
        if (!file) {
            cerr << "Error: cannot open '" << argv[1] << "'\n";
            return 2;
        }
        in = &file;
    }

    vector<Row> rows = parse_rows(*in);
    if (rows.empty()) {
        cout << "No valid 'SN' entries found in the input.\n";
        return 0;
    }

    size_t metric_w = string("Metric").size();
    size_t value_w = string("Value").size();
    size_t comment_w = string("Comment").size();
    for (const auto& row : rows) {
        metric_w = max(metric_w, row.metric.size());
        value_w = max(value_w, row.value.size());
        comment_w = max(comment_w, row.comment.size());
    }

    cout << "\nAlignment Statistics:\n\n";
    cout << left << setw(static_cast<int>(metric_w)) << "Metric" << "  "
         << right << setw(static_cast<int>(value_w)) << "Value" << "  "
         << left << setw(static_cast<int>(comment_w)) << "Comment" << "\n";
    cout << string(metric_w, '-') << "  "
         << string(value_w, '-') << "  "
         << string(comment_w, '-') << "\n";

    for (const auto& row : rows) {
        cout << left << setw(static_cast<int>(metric_w)) << row.metric << "  "
             << right << setw(static_cast<int>(value_w)) << row.value << "  "
             << left << setw(static_cast<int>(comment_w)) << row.comment << "\n";
    }

    return 0;
}
