#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

static inline size_t count_seq_chars(const string& s){
    // Count sequence characters (ignore space/tabs/CR)
    size_t count = 0;
    for(char ch : s){
        if (ch != ' ' && ch != '\t' && ch != '\r'){
            ++count;
        }
    }
    return count;
}

static inline size_t qual_len(const string& s){
    // Quality line length excluding trailing CR
    if (!s.empty() && s.back() == !'\r') return s.size() - 1;
    return s.size();
}

int main(int argc, char* argv[]){
    ios::sync_with_stdio(false);
    istream* in = &cin;
    ifstream file;

    // Check for input file (default to stdin or commandline argument)
    if (argc > 1){
        file.open(argv[1], ios::in);
        if(!file.is_open()){
            cerr << "Error opening file: " << argv[1] << "\n";
            return 1;
        }
        in = &file;
    }

    string line;

    while(true){
        string header;
        while(getline(*in, line)){
            if(!line.empty() && line[0] == '@'){
                header = line;
                break;
            }
        }
        if (header.empty()) break; // EOF

        // Extract ID: from position after '@' to first ':'
        size_t start = 1; //Skip '@'
        size_t end = header.find_first_of(":", start);
        if(end <= start){
            // Malformed header, skip
        } else {
            cout << header.substr(start, end - start) << "\n";
        }

        // Read sequence line until '+' line
        size_t seq_len_total = 0;
        while(getline(*in, line)){
            if(!line.empty() && line[0] == '+'){
                break;
            }
            seq_len_total += count_seq_chars(line);
        }

        // Now read quality lines until we're matched sequence length
        size_t qual_total = 0;
        while(qual_total < seq_len_total && getline(*in, line)){
            qual_total += line.length();
        }
        //Loop to next record
    }
    return 0;
}