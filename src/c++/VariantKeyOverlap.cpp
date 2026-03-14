#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

struct Dataset {
    string id;
    string label;
    string path;
    unordered_set<string> variants;
};

static string basename_of(const string& path) {
    size_t pos = path.find_last_of("/\\");
    return (pos == string::npos) ? path : path.substr(pos + 1);
}

static string trim_copy(const string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == string::npos) return string();
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " [--exclude-wt] <counts1.tsv> <counts2.tsv> [counts3.tsv ...]\n\n"
         <<
R"(Usage:
Description:
  Compare DualSiteDMSFilter --out-counts TSV files across any number of datasets.
  The program uses the variant_key column as a unique set per dataset and reports:
    - unique variants per dataset
    - pairwise unique-variant intersections
    - overlap counts by number of datasets containing a variant
    - exact dataset-sharing patterns

Options:
  --exclude-wt    Ignore the WT variant_key during comparison.
)";
}

static unordered_set<string> load_variant_keys(const string& path, bool exclude_wt) {
    ifstream in(path);
    if (!in) {
        throw runtime_error("cannot open input file: " + path);
    }

    unordered_set<string> variants;
    variants.reserve(1 << 15);

    string line;
    uint64_t line_no = 0;
    while (getline(in, line)) {
        ++line_no;
        line = trim_copy(line);
        if (line.empty()) continue;

        size_t tab = line.find('\t');
        string key = trim_copy((tab == string::npos) ? line : line.substr(0, tab));
        if (key.empty()) continue;
        if (line_no == 1 && key == "variant_key") continue;
        if (exclude_wt && key == "WT") continue;
        variants.insert(std::move(key));
    }

    return variants;
}

static uint64_t intersection_size(const unordered_set<string>& a, const unordered_set<string>& b) {
    const unordered_set<string>* small = &a;
    const unordered_set<string>* large = &b;
    if (a.size() > b.size()) {
        small = &b;
        large = &a;
    }

    uint64_t count = 0;
    for (const auto& key : *small) {
        if (large->find(key) != large->end()) ++count;
    }
    return count;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        bool exclude_wt = false;
        vector<string> paths;
        paths.reserve(static_cast<size_t>(max(0, argc - 1)));

        for (int i = 1; i < argc; ++i) {
            string arg = argv[i];
            if (arg == "--exclude-wt") {
                exclude_wt = true;
            } else if (arg == "--help" || arg == "-h") {
                usage(argv[0]);
                return 0;
            } else {
                paths.push_back(std::move(arg));
            }
        }

        if (paths.size() < 2) {
            usage(argv[0]);
            return 1;
        }

        vector<Dataset> datasets;
        datasets.reserve(paths.size());
        for (size_t i = 0; i < paths.size(); ++i) {
            Dataset ds;
            ds.id = "D" + to_string(i + 1);
            ds.path = paths[i];
            ds.label = basename_of(paths[i]);
            ds.variants = load_variant_keys(paths[i], exclude_wt);
            datasets.push_back(std::move(ds));
        }

        size_t n = datasets.size();
        uint64_t total_unique_entries = 0;
        for (const auto& ds : datasets) total_unique_entries += ds.variants.size();

        unordered_map<string, vector<uint8_t>> membership;
        membership.reserve(static_cast<size_t>(total_unique_entries * 1.3) + 1);

        for (size_t i = 0; i < n; ++i) {
            for (const auto& key : datasets[i].variants) {
                auto [it, inserted] = membership.emplace(key, vector<uint8_t>(n, 0));
                it->second[i] = 1;
            }
        }

        vector<uint64_t> multiplicity_counts(n + 1, 0);
        unordered_map<string, uint64_t> pattern_counts;
        pattern_counts.reserve(membership.size());
        uint64_t shared_in_all = 0;
        uint64_t shared_in_at_least_two = 0;

        for (const auto& kv : membership) {
            const vector<uint8_t>& present = kv.second;
            size_t count = 0;
            string pattern;
            for (size_t i = 0; i < n; ++i) {
                if (!present[i]) continue;
                ++count;
                if (!pattern.empty()) pattern += "&";
                pattern += datasets[i].id;
            }

            multiplicity_counts[count]++;
            pattern_counts[pattern]++;
            if (count == n) ++shared_in_all;
            if (count >= 2) ++shared_in_at_least_two;
        }

        vector<pair<string, uint64_t>> pattern_items(pattern_counts.begin(), pattern_counts.end());
        sort(pattern_items.begin(), pattern_items.end(), [](const auto& x, const auto& y) {
            if (x.second != y.second) return x.second > y.second;
            return x.first < y.first;
        });

        cout << fixed << setprecision(6);
        cout << "# VariantKeyOverlap\n";
        cout << "exclude_wt=" << (exclude_wt ? "true" : "false") << "\n";
        cout << "dataset_count=" << n << "\n";
        cout << "union_unique_variants=" << membership.size() << "\n";
        cout << "shared_in_at_least_2_datasets=" << shared_in_at_least_two << "\n";
        cout << "shared_in_all_datasets=" << shared_in_all << "\n";

        cout << "\n# Datasets\n";
        cout << "dataset_id\tlabel\tunique_variants\tpath\n";
        for (const auto& ds : datasets) {
            cout << ds.id << '\t' << ds.label << '\t' << ds.variants.size() << '\t' << ds.path << '\n';
        }

        cout << "\n# PairwiseIntersections\n";
        cout << "dataset_a\tdataset_b\tshared_unique_variants\tjaccard\n";
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                uint64_t shared = intersection_size(datasets[i].variants, datasets[j].variants);
                uint64_t union_size = datasets[i].variants.size() + datasets[j].variants.size() - shared;
                double jaccard = union_size ? static_cast<double>(shared) / static_cast<double>(union_size) : 0.0;
                cout << datasets[i].id << '\t' << datasets[j].id << '\t' << shared << '\t' << jaccard << '\n';
            }
        }

        cout << "\n# OverlapByMultiplicity\n";
        cout << "datasets_present\tvariant_count\n";
        for (size_t k = 1; k <= n; ++k) {
            cout << k << '\t' << multiplicity_counts[k] << '\n';
        }

        cout << "\n# PatternCounts\n";
        cout << "dataset_pattern\tvariant_count\n";
        for (const auto& item : pattern_items) {
            cout << item.first << '\t' << item.second << '\n';
        }

        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}
