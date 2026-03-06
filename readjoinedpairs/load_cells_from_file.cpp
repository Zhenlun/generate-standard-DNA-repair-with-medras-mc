#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <any>
#include <vector>
#include <queue>
#include <filesystem>
#include <regex>

using namespace std;

// Function to load cells from a .joinedpairs file
vector<pair<string, string>> load_cells_from_file(const string& filename) {
    vector<pair<string, string>> cells;

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return cells;
    }

    vector<string> lines;
    string line;
    while (getline(file, line)) {
        lines.push_back(line);
    }
    file.close();

    // Find delimiter lines (lines consisting only of -)
    vector<size_t> delimiters;
    regex delim_regex(R"(^-+$)");
    for (size_t i = 0; i < lines.size(); ++i) {
        if (regex_match(lines[i], delim_regex)) {
            delimiters.push_back(i);
        }
    }

    if (delimiters.empty()) {
        cerr << "No delimiters found in file." << endl;
        return cells;
    }

    // The sections are between delimiters
    // First section: 0 to delimiters[0]
    // Second: delimiters[0]+1 to delimiters[1]
    // ...
    // Last section: delimiters.back()+1 to end, but wait, no: the lists start after the last delimiter.

    // Actually, the sections are from start to first delim, then after first delim to second, etc.
    // But the last section ends at the last delim, and then lists after.

    // Number of cells = number of sections = delimiters.size()

    size_t num_cells = delimiters.size();

    // Now, collect cell data
    vector<string> cell_data_texts;
    size_t start = 0;
    for (size_t i = 0; i < num_cells; ++i) {
        size_t end = delimiters[i];
        string cell_text;
        for (size_t j = start; j < end; ++j) {
            if (!cell_text.empty()) cell_text += "\n";
            cell_text += lines[j];
        }
        cell_data_texts.push_back(cell_text);
        start = end + 1;
    }

    // Now, the repair data starts after the last delimiter
    size_t repair_start = delimiters.back() + 1;
    vector<string> repair_texts;
    for (size_t i = repair_start; i < lines.size(); ++i) {
        repair_texts.push_back(lines[i]);
    }

    if (repair_texts.size() != num_cells) {
        cerr << "Number of repair data lines (" << repair_texts.size() << ") does not match number of cells (" << num_cells << ")" << endl;
        return cells;
    }

    // Now, pair them
    for (size_t i = 0; i < num_cells; ++i) {
        cells.emplace_back(cell_data_texts[i], repair_texts[i]);
    }

    return cells;
}

/*
int main() {
    // try to open the sample file relative to cwd; tests in output/ may need to go up one level
    filesystem::path sample = "examplesdd.joinedpairs";
    if (!filesystem::exists(sample)) {
        sample = filesystem::current_path() / ".." / "examplesdd.joinedpairs";
    }

    cout << "Looking for file at " << sample.string() << endl;
    auto cells = load_cells_from_file(sample.string());
    if (!cells.empty()) {
        cout << "First cell data:" << endl;
        cout << cells[0].first << endl;
        cout << "First cell repair data:" << endl;
        cout << cells[0].second << endl;
        cout << "Last cell data:" << endl;
        cout << cells.back().first << endl;
        cout << "Last cell repair data:" << endl;
        cout << cells.back().second << endl;
    } else {
        cout << "No cells loaded." << endl;
    }
    return 0;
}
*/