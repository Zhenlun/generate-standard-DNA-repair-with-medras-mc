#include "load_cells_from_file.cpp"
#include <sstream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>


using namespace std;

// Structure to represent a double strand break end
struct dsb_end {
    int break_index;
    int is_complex;
    int chromosome_id;
    int dna_structure;
    int chromosome_set;
    int chromatid_strand;
    int chromosome_arm;
    double position_in_chromosome;
    int upstream_downstream;
};

// Structure to represent a DNA repair pair
struct dna_repair_pair {
    int inter_chrom;  // 0 or 1
    dsb_end end1;
    dsb_end end2;
};

// Segment representation if you want to keep track of pieces
struct segment {
    int old_chromosome_id;
    int old_chromatid_strand;
    double old_chromatid_start_position;
    double old_chromatid_end_position;
    int chromosome_arm;
    int upstream_downstream; // from original ends
};


// Class to represent repaired chromosome geometry
class repaired_chromosome {
public:

    // imported information from medras
    std::string cell_id;
    int dna_dsb_count;
    int unrepaired_dsb;
    int misrepairs;
    int large_misrepairs;
    int inter_chrom_misrepair;
    int dicentrics;
    int rings;
    int excess_linear;
    int total_aberrations;
    int viability;

    // after reconstruction analysis to make sure its the same as imported from medras
    int after_reconstruction_dicentrics;
    int after_reconstruction_rings;
    int after_reconstruction_excess_linear;
    int after_reconstruction_multicentrics;
    int after_reconstruction_total_aberration;

    std::vector<dna_repair_pair> repair_pairs;

    repaired_chromosome(const std::pair<std::string, std::string>& cell);
    void parse_cell_data(const std::string& cell_data);

    // New: analyze geometry of repaired chromosomes
    void analyze_repair_geometry();

    void visualize_graph();

private:
    // Helper: build graph and classify components
    void build_graph_and_classify();
};



void repaired_chromosome::analyze_repair_geometry() {
    // Reset counts in case this is called multiple times
    after_reconstruction_dicentrics = 0;
    after_reconstruction_rings = 0;
    after_reconstruction_excess_linear = 0;
    after_reconstruction_total_aberration = 0;
    after_reconstruction_multicentrics = 0;

    build_graph_and_classify();
}




// Core logic
void repaired_chromosome::build_graph_and_classify() {
    using namespace std;

    // Local counters (fully initialized)
    int dicentrics = 0;
    int rings = 0;
    int excess_linear = 0;
    int multicentrics = 0;

    // Map each break_index to its dsb_end
    unordered_map<int, dsb_end> ends;
    // Adjacency list: break_index -> neighbors
    unordered_map<int, vector<int>> adj;

    // Build nodes and edges from repair_pairs
    for (const auto& rp : repair_pairs) {
        const dsb_end& e1 = rp.end1;
        const dsb_end& e2 = rp.end2;

        ends[e1.break_index] = e1;
        ends[e2.break_index] = e2;

        adj[e1.break_index].push_back(e2.break_index);
        adj[e2.break_index].push_back(e1.break_index);
    }

    unordered_set<int> visited;

    // Traverse each connected component
    for (const auto& kv : ends) {
        int start = kv.first;
        if (visited.count(start)) continue;

        vector<int> component;
        queue<int> q;
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            component.push_back(u);
            for (int v : adj[u]) {
                if (!visited.count(v)) {
                    visited.insert(v);
                    q.push(v);
                }
            }
        }

        if (component.empty()) continue;

        // Degree analysis
        bool is_ring = true;

        for (int node : component) {
            int deg = (int)adj[node].size();
            if (deg != 2) is_ring = false;
        }

        // Centromere counting:
        // For each chromosome_id in this component, check if we have both upstream and downstream ends.
        unordered_map<int, pair<bool,bool>> chrom_up_down; // chrom -> (has_upstream, has_downstream)

        for (int node : component) {
            const dsb_end& e = ends[node];
            auto& flags = chrom_up_down[e.chromosome_id];
            if (e.upstream_downstream == 0) {
                flags.first = true;   // upstream
            } else {
                flags.second = true;  // downstream
            }
        }

        int centromere_count = 0;
        for (const auto& ckv : chrom_up_down) {
            const auto& flags = ckv.second;
            if (flags.first && flags.second) {
                centromere_count++;
            }
        }

        // Classify this component
        if (is_ring && !component.empty()) {
            rings++;
            if (centromere_count == 2) {
                dicentrics++; // dicentric ring
            } else if (centromere_count > 2) {
                multicentrics++;
            }
        } else {
            if (centromere_count == 0) {
                excess_linear++; // acentric linear fragment
            } else if (centromere_count == 2) {
                dicentrics++;
            } else if (centromere_count > 2) {
                multicentrics++;
            }
        }
    }

    // Commit to members (no +=, no reuse of old values)
    after_reconstruction_dicentrics = dicentrics;
    after_reconstruction_rings = rings;
    after_reconstruction_excess_linear = excess_linear;
    after_reconstruction_multicentrics = multicentrics;
    after_reconstruction_total_aberration =
        dicentrics + rings + excess_linear + multicentrics;
}


void repaired_chromosome::visualize_graph() {
    using namespace std;

    unordered_map<int, vector<int>> adj;
    unordered_map<int, dsb_end> ends;

    // Build adjacency
    for (const auto& rp : repair_pairs) {
        ends[rp.end1.break_index] = rp.end1;
        ends[rp.end2.break_index] = rp.end2;

        adj[rp.end1.break_index].push_back(rp.end2.break_index);
        adj[rp.end2.break_index].push_back(rp.end1.break_index);
    }

    unordered_set<int> visited;
    int comp_id = 1;

    cout << "\n=== Graph Visualization for Cell " << cell_id << " ===\n";

    for (auto& kv : ends) {
        int start = kv.first;
        if (visited.count(start)) continue;

        cout << "\nComponent " << comp_id++ << ":\n";

        // BFS to list nodes
        queue<int> q;
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            int u = q.front(); q.pop();

            const dsb_end& e = ends[u];
            cout << "  Node " << u
                 << "  (Chr " << e.chromosome_id
                 << ", Arm " << e.chromosome_arm
                 << ", Pos " << e.position_in_chromosome
                 << ", Up/Down " << e.upstream_downstream
                 << ")  -->  ";

            for (int v : adj[u]) cout << v << " ";
            cout << "\n";

            for (int v : adj[u]) {
                if (!visited.count(v)) {
                    visited.insert(v);
                    q.push(v);
                }
            }
        }
    }

    cout << "=============================================\n";
}




// Helper prototype -- defined below

dsb_end parse_dsb_end(const std::string& end_str);

// -----------------------------------------------------------------------------
// Implementation of helper and class methods
// -----------------------------------------------------------------------------

dsb_end parse_dsb_end(const std::string& end_str) {
    dsb_end end = {};
    
    // Manual parsing: find key values in the string
    // Format: [break_index, array([x, y, z]), is_complex, [dna_struct, chrom_set, chromatid, arm], position, upstream/downstream, ...]
    
    std::string clean = end_str;
    // Remove outer brackets
    if (clean.front() == '[') clean.erase(0, 1);
    if (clean.back() == ']') clean.pop_back();
    
    // Split by comma, but be careful with nested brackets
    std::vector<std::string> tokens;
    std::string token;
    int bracket_depth = 0;
    
    for (size_t i = 0; i < clean.size(); ++i) {
        char c = clean[i];
        
        if (c == '[' || c == '(') {
            bracket_depth++;
            token += c;
        } else if (c == ']' || c == ')') {
            bracket_depth--;
            token += c;
        } else if (c == ',' && bracket_depth == 0) {
            tokens.push_back(token);
            token.clear();
        } else {
            token += c;
        }
    }
    if (!token.empty()) tokens.push_back(token);
    
    // Trim whitespace from tokens
    for (auto& t : tokens) {
        t.erase(t.begin(), std::find_if(t.begin(), t.end(), [](unsigned char ch) { return !std::isspace(ch); }));
        t.erase(std::find_if(t.rbegin(), t.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), t.end());
    }
    
    // Parse tokens
    if (tokens.size() >= 6) {
        try {
            end.break_index = std::stoi(tokens[0]);
            // tokens[1] is array(...), skip
            end.is_complex = std::stoi(tokens[2]);
            
            // Parse chromosome info from tokens[3]: [dna_struct, chrom_set, chromatid, arm]
            std::string chrom_str = tokens[3];
            if (chrom_str.front() == '[') chrom_str.erase(0, 1);
            if (chrom_str.back() == ']') chrom_str.pop_back();
            
            std::vector<int> chrom_vals;
            std::stringstream ss(chrom_str);
            std::string val;
            while (std::getline(ss, val, ',')) {
                val.erase(val.begin(), std::find_if(val.begin(), val.end(), [](unsigned char ch) { return !std::isspace(ch); }));
                val.erase(std::find_if(val.rbegin(), val.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), val.end());
                chrom_vals.push_back(std::stoi(val));
            }
            
            if (chrom_vals.size() >= 4) {
                end.dna_structure = chrom_vals[0];
                end.chromosome_set = chrom_vals[1];
                end.chromatid_strand = chrom_vals[2];
                end.chromosome_arm = chrom_vals[3];
            }
            
            end.position_in_chromosome = std::stod(tokens[4]);
            end.upstream_downstream = std::stoi(tokens[5]);
        } catch (...) {
            // Set defaults on parse error
        }
    }
    return end;
}

// constructor
repaired_chromosome::repaired_chromosome(const std::pair<std::string, std::string>& cell) {
    // Parse repair data
    std::string repair_str = cell.second;

    // Remove leading [ and trailing ]
    if (!repair_str.empty() && repair_str.front() == '[' && repair_str.back() == ']') {
        repair_str = repair_str.substr(1, repair_str.size() - 2);
    }

    // Split by ", "
    std::vector<std::string> parts;
    std::stringstream ss(repair_str);
    std::string token;
    while (std::getline(ss, token, ',')) {
        // Trim spaces and single quotes
        token.erase(token.begin(), std::find_if(token.begin(), token.end(), [](unsigned char ch) { return !std::isspace(ch) && ch != '\''; }));
        token.erase(std::find_if(token.rbegin(), token.rend(), [](unsigned char ch) { return !std::isspace(ch) && ch != '\''; }).base(), token.end());
        parts.push_back(token);
    }

    if (parts.size() == 11) {
        // Normal case
        cell_id = parts[0];
        dna_dsb_count = std::stoi(parts[1]);
        unrepaired_dsb = std::stoi(parts[2]);
        misrepairs = std::stoi(parts[3]);
        large_misrepairs = std::stoi(parts[4]);
        inter_chrom_misrepair = std::stoi(parts[5]);
        dicentrics = std::stoi(parts[6]);
        rings = std::stoi(parts[7]);
        excess_linear = std::stoi(parts[8]);
        total_aberrations = std::stoi(parts[9]);
        viability = std::stoi(parts[10]);
    } else if (parts.size() == 6 && parts[4] == "No" && parts[5] == "misrepair!") {
        // Special case: no misrepair
        cell_id = parts[0];
        dna_dsb_count = std::stoi(parts[1]);
        unrepaired_dsb = std::stoi(parts[2]);
        misrepairs = std::stoi(parts[3]);
        large_misrepairs = 0;
        inter_chrom_misrepair = 0;
        dicentrics = 0;
        rings = 0;
        excess_linear = 0;
        total_aberrations = 0;
        viability = 1;  // Assume alive since no misrepair
    } else {
        // Set defaults on parse error
        cell_id = "error";
        dna_dsb_count = 0;
        unrepaired_dsb = 0;
        misrepairs = 0;
        large_misrepairs = 0;
        inter_chrom_misrepair = 0;
        dicentrics = 0;
        rings = 0;
        excess_linear = 0;
        total_aberrations = 0;
        viability = 0;
    }

    // Parse cell data
    parse_cell_data(cell.first);

    // Run reconstruction
    analyze_repair_geometry();

}

void repaired_chromosome::parse_cell_data(const std::string& cell_data) {
    // Split cell data by newlines
    std::vector<std::string> lines;
    std::stringstream ss(cell_data);
    std::string line;
    while (std::getline(ss, line)) {
        if (!line.empty()) {
            lines.push_back(line);
        }
    }

    // Process every 3 lines as a pair
    for (size_t i = 0; i + 2 < lines.size(); i += 3) {
        dna_repair_pair pair;

        // Parse interchrom line
        std::string interchrom_line = lines[i];
        if (interchrom_line.find("interChrom = 0") != std::string::npos) {
            pair.inter_chrom = 0;
        } else if (interchrom_line.find("interChrom = 1") != std::string::npos) {
            pair.inter_chrom = 1;
        }

        // Parse end1 and end2
        pair.end1 = parse_dsb_end(lines[i + 1]);
        pair.end2 = parse_dsb_end(lines[i + 2]);

        repair_pairs.push_back(pair);
    }
}



int main() {
    // Load cells from the sample file
    filesystem::path sample = "examplesdd.joinedpairs";
    if (!filesystem::exists(sample)) {
        sample = filesystem::current_path() / ".." / "examplesdd.joinedpairs";
    }

    cout << "Loading cells from " << sample.string() << endl;
    auto cells = load_cells_from_file(sample.string());
    if (cells.empty()) {
        cout << "No cells loaded." << endl;
        return 1;
    }

    // Process and print data for each cell
    for (size_t i = 0; i < cells.size(); ++i) {
        repaired_chromosome geom(cells[i]);
        cout << "Cell " << (i) << ":" << endl;
        cout << "Cell ID: " << geom.cell_id << endl;
        cout << "DNA DSB Count: " << geom.dna_dsb_count << endl;
        cout << "Unrepaired DSB: " << geom.unrepaired_dsb << endl;
        cout << "Misrepairs: " << geom.misrepairs << endl;
        cout << "Large Misrepairs: " << geom.large_misrepairs << endl;
        cout << "Inter-Chromosome Misrepair: " << geom.inter_chrom_misrepair << endl;
        cout << "Dicentrics: " << geom.dicentrics << endl;
        cout << "Rings: " << geom.rings << endl;
        cout << "Excess Linear Fragments: " << geom.excess_linear << endl;
        cout << "Total Aberrations: " << geom.total_aberrations << endl;
        cout << "Viability: " << geom.viability << endl;
        cout << "Repair Pairs: " << geom.repair_pairs.size() << endl;
        cout << "After Reconstruction - Dicentrics: " << geom.after_reconstruction_dicentrics << endl;
        cout << "After Reconstruction - Rings: " << geom.after_reconstruction_rings << endl;
        cout << "After Reconstruction - Excess Linear: " << geom.after_reconstruction_excess_linear << endl;
        cout << "After Reconstruction - Multicentrics: " << geom.after_reconstruction_multicentrics << endl;
        cout << "After Reconstruction - Total Aberrations: " << geom.after_reconstruction_total_aberration << endl;

        geom.visualize_graph();

        // Print repair pairs
        /*
        for (size_t j = 0; j < geom.repair_pairs.size(); ++j) {
            cout << "  Pair " << (j + 1) << ":" << endl;
            cout << "    Inter-Chromosome: " << geom.repair_pairs[j].inter_chrom << endl;
            cout << "    End1:" << endl;
            cout << "      Break Index: " << geom.repair_pairs[j].end1.break_index << endl;
            cout << "      Is Complex: " << geom.repair_pairs[j].end1.is_complex << endl;
            cout << "      DNA Structure: " << geom.repair_pairs[j].end1.dna_structure << endl;
            cout << "      Chromosome Set: " << geom.repair_pairs[j].end1.chromosome_set << endl;
            cout << "      Chromatid Strand: " << geom.repair_pairs[j].end1.chromatid_strand << endl;
            cout << "      Chromosome Arm: " << geom.repair_pairs[j].end1.chromosome_arm << endl;
            cout << "      Position in Chromosome: " << geom.repair_pairs[j].end1.position_in_chromosome << endl;
            cout << "      Upstream/Downstream: " << geom.repair_pairs[j].end1.upstream_downstream << endl;
            cout << "    End2:" << endl;
            cout << "      Break Index: " << geom.repair_pairs[j].end2.break_index << endl;
            cout << "      Is Complex: " << geom.repair_pairs[j].end2.is_complex << endl;
            cout << "      DNA Structure: " << geom.repair_pairs[j].end2.dna_structure << endl;
            cout << "      Chromosome Set: " << geom.repair_pairs[j].end2.chromosome_set << endl;
            cout << "      Chromatid Strand: " << geom.repair_pairs[j].end2.chromatid_strand << endl;
            cout << "      Chromosome Arm: " << geom.repair_pairs[j].end2.chromosome_arm << endl;
            cout << "      Position in Chromosome: " << geom.repair_pairs[j].end2.position_in_chromosome << endl;
            cout << "      Upstream/Downstream: " << geom.repair_pairs[j].end2.upstream_downstream << endl;
        }
            */
        cout << "------------------------------" << endl;
    }

    return 0;
}
