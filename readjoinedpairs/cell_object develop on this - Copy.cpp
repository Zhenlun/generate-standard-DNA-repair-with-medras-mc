#include "load_cells_from_file.cpp"
#include <sstream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>


using namespace std;

// Structure to represent a double strand break end
struct dsb_end {
    int break_index;
    int is_complex;
    //int chromosome_id;
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
    int id; // unique ID across all segments
    int old_chromosome_id;
    int old_chromatid_strand;
    double old_chromatid_start_position;
    double old_chromatid_end_position;
    int has_centromere;
};

enum SegmentSide {
    LEFT = 0,   // at old_chromatid_start_position
    RIGHT = 1   // at old_chromatid_end_position
};

struct SegmentSideNode {
    int segment_id;
    SegmentSide side;

    bool operator==(const SegmentSideNode& other) const {
        return segment_id == other.segment_id && side == other.side;
    }
};

struct SegmentSideNodeHash {
    std::size_t operator()(const SegmentSideNode& n) const {
        std::size_t h1 = std::hash<int>()(n.segment_id);
        std::size_t h2 = std::hash<int>()(static_cast<int>(n.side));
        return (h1 ^ (h2 << 1));
    }
};

struct StrandSegment {
    int segment_id;
    int orientation; // +1: start→end, -1: end→start
};

struct ChimericStrand {
    std::vector<StrandSegment> segments;
    bool is_ring;
};


struct ChromChromatidKey {
    int chromosome_set;
    int chromatid_strand;

    bool operator<(const ChromChromatidKey& other) const {
        if (chromosome_set != other.chromosome_set)
            return chromosome_set < other.chromosome_set;
        return chromatid_strand < other.chromatid_strand;
    }
};

// Structure to store affected chromosome and chromatid strand combinations
struct affected_chrom_chromatid {
    int chromosome_set;
    int chromatid_strand;
    
    affected_chrom_chromatid(int chrom, int strand) 
        : chromosome_set(chrom), chromatid_strand(strand) {}
    
    bool operator<(const affected_chrom_chromatid& other) const {
        if (chromosome_set != other.chromosome_set) {
            return chromosome_set < other.chromosome_set;
        }
        return chromatid_strand < other.chromatid_strand;
    }
    
    bool operator==(const affected_chrom_chromatid& other) const {
        return chromosome_set == other.chromosome_set && 
               chromatid_strand == other.chromatid_strand;
    }
};

// Function to get all affected chromosomes and chromatid strands
std::set<affected_chrom_chromatid> get_affected_chromosomes(
    const std::vector<dna_repair_pair>& repair_pairs) {
    
    std::set<affected_chrom_chromatid> affected;
    
    // Iterate through all repair pairs
    for (const auto& pair : repair_pairs) {
        // Add chromosome-chromatid combination for end1
        affected.insert(affected_chrom_chromatid(pair.end1.chromosome_set, pair.end1.chromatid_strand));
        
        // Add chromosome-chromatid combination for end2
        affected.insert(affected_chrom_chromatid(pair.end2.chromosome_set, pair.end2.chromatid_strand));
    }
    
    return affected;
}

// Function to get all broken segments for a given chromosome-chromatid combination
std::vector<segment> get_broken_segments(
    const affected_chrom_chromatid& chr_chromatid,
    const std::vector<dna_repair_pair>& repair_pairs) {
    
    std::vector<segment> segments;
    std::vector<double> break_positions;
    
    // Collect all break positions for this chromosome-chromatid combination
    for (const auto& pair : repair_pairs) {
        // Check end1
        if (pair.end1.chromosome_set == chr_chromatid.chromosome_set &&
            pair.end1.chromatid_strand == chr_chromatid.chromatid_strand) {
            break_positions.push_back(pair.end1.position_in_chromosome);
        }
        
        // Check end2
        if (pair.end2.chromosome_set == chr_chromatid.chromosome_set &&
            pair.end2.chromatid_strand == chr_chromatid.chromatid_strand) {
            break_positions.push_back(pair.end2.position_in_chromosome);
        }
    }
    
    // Sort break positions
    std::sort(break_positions.begin(), break_positions.end());
    
    // Remove duplicates
    break_positions.erase(std::unique(break_positions.begin(), break_positions.end()), break_positions.end());
    
    // Centromere location
    const double centromere_position = 0.5;
    
    // DNA strand spans from 0 to 1
    double start_pos = 0.0;
    const double chromosome_end = 1.0;
    
    // Create segments between breaks
    for (size_t i = 0; i < break_positions.size(); ++i) {
        segment seg;
        seg.old_chromosome_id = chr_chromatid.chromosome_set;
        seg.old_chromatid_strand = chr_chromatid.chromatid_strand;
        seg.old_chromatid_start_position = start_pos;
        seg.old_chromatid_end_position = break_positions[i];
        
        // Check if centromere is in this segment
        seg.has_centromere = (start_pos <= centromere_position && centromere_position < break_positions[i]) ? 1 : 0;
        
        segments.push_back(seg);
        start_pos = break_positions[i];
    }
    
    // Add final segment from last break to end of chromosome (1.0)
    if (!break_positions.empty()) {
        segment final_seg;
        final_seg.old_chromosome_id = chr_chromatid.chromosome_set;
        final_seg.old_chromatid_strand = chr_chromatid.chromatid_strand;
        final_seg.old_chromatid_start_position = break_positions.back();
        final_seg.old_chromatid_end_position = chromosome_end;
        
        // Check if centromere is in this final segment
        final_seg.has_centromere = (break_positions.back() <= centromere_position && centromere_position < chromosome_end) ? 1 : 0;
        
        segments.push_back(final_seg);
    }
    
    return segments;
}

std::map<ChromChromatidKey, std::vector<segment>>
build_all_segments(const std::vector<dna_repair_pair>& repair_pairs) {
    std::map<ChromChromatidKey, std::vector<segment>> segments_by_cc;

    auto affected = get_affected_chromosomes(repair_pairs);
    int next_id = 0;

    for (const auto& chr_chromatid : affected) {
        ChromChromatidKey key{chr_chromatid.chromosome_set,
                              chr_chromatid.chromatid_strand};
        auto segs = get_broken_segments(chr_chromatid, repair_pairs);
        for (auto& s : segs) {
            s.id = next_id++;
            s.old_chromosome_id = chr_chromatid.chromosome_set;
            s.old_chromatid_strand = chr_chromatid.chromatid_strand;
        }
        segments_by_cc[key] = std::move(segs);
    }

    return segments_by_cc;
}

struct SegmentSideKey {
    int chromosome_set;
    int chromatid_strand;
    double position;
    SegmentSide side;

    bool operator==(const SegmentSideKey& other) const {
        return chromosome_set == other.chromosome_set &&
               chromatid_strand == other.chromatid_strand &&
               position == other.position &&
               side == other.side;
    }
};

struct SegmentSideKeyHash {
    std::size_t operator()(const SegmentSideKey& k) const {
        std::size_t h1 = std::hash<int>()(k.chromosome_set);
        std::size_t h2 = std::hash<int>()(k.chromatid_strand);
        std::size_t h3 = std::hash<double>()(k.position);
        std::size_t h4 = std::hash<int>()(static_cast<int>(k.side));
        return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1)) ^ (h4 << 1);
    }
};

std::unordered_map<SegmentSideKey, SegmentSideNode, SegmentSideKeyHash>
build_segment_side_map(const std::map<ChromChromatidKey, std::vector<segment>>& segments_by_cc) {
    std::unordered_map<SegmentSideKey, SegmentSideNode, SegmentSideKeyHash> map;

    for (const auto& kv : segments_by_cc) {
        const auto& key = kv.first;
        const auto& segs = kv.second;

        for (const auto& s : segs) {
            // internal break at start (not telomere)
            if (s.old_chromatid_start_position > 0.0) {
                SegmentSideKey sk{key.chromosome_set, key.chromatid_strand,
                                  s.old_chromatid_start_position, LEFT};
                map[sk] = SegmentSideNode{s.id, LEFT};
            }
            // internal break at end (not telomere)
            if (s.old_chromatid_end_position < 1.0) {
                SegmentSideKey sk{key.chromosome_set, key.chromatid_strand,
                                  s.old_chromatid_end_position, RIGHT};
                map[sk] = SegmentSideNode{s.id, RIGHT};
            }
        }
    }

    return map;
}

bool map_dsb_end_to_segment_side(
    const dsb_end& end,
    const std::unordered_map<SegmentSideKey, SegmentSideNode, SegmentSideKeyHash>& side_map,
    SegmentSideNode& out_node)
{
    // upstream_downstream: -1 = upstream, +1 = downstream
    SegmentSide side = (end.upstream_downstream == -1) ? RIGHT : LEFT;

    SegmentSideKey key{end.chromosome_set,
                       end.chromatid_strand,
                       end.position_in_chromosome,
                       side};

    auto it = side_map.find(key);
    if (it == side_map.end()) return false;
    out_node = it->second;
    return true;
}

using AdjacencyMap = std::unordered_map<SegmentSideNode, SegmentSideNode, SegmentSideNodeHash>;

AdjacencyMap build_adjacency(
    const std::vector<dna_repair_pair>& repair_pairs,
    const std::unordered_map<SegmentSideKey, SegmentSideNode, SegmentSideKeyHash>& side_map)
{
    AdjacencyMap adj;

    for (const auto& pair : repair_pairs) {
        SegmentSideNode n1, n2;
        if (!map_dsb_end_to_segment_side(pair.end1, side_map, n1)) continue;
        if (!map_dsb_end_to_segment_side(pair.end2, side_map, n2)) continue;

        adj[n1] = n2;
        adj[n2] = n1;
    }

    return adj;
}

SegmentSide opposite_side(SegmentSide s) {
    return (s == LEFT) ? RIGHT : LEFT;
}

std::vector<ChimericStrand> build_chimeric_strands(
    const std::map<ChromChromatidKey, std::vector<segment>>& segments_by_cc,
    const AdjacencyMap& adj)
{
    std::vector<ChimericStrand> strands;

    // segment_id -> segment
    std::unordered_map<int, segment> seg_by_id;
    for (const auto& kv : segments_by_cc)
        for (const auto& s : kv.second)
            seg_by_id[s.id] = s;

    std::unordered_set<int> used_segments;

    auto walk_from = [&](SegmentSideNode start_node,
                         std::vector<StrandSegment>& path,
                         bool& hit_ring)
    {
        std::unordered_set<SegmentSideNode, SegmentSideNodeHash> visited_nodes;
        SegmentSideNode current = start_node;
        visited_nodes.insert(current);

        while (true) {
            SegmentSide through_side = current.side;
            SegmentSide exit_side = opposite_side(through_side);

            int orientation = (through_side == LEFT && exit_side == RIGHT) ? +1 : -1;

            if (path.empty() || path.back().segment_id != current.segment_id) {
                path.push_back(StrandSegment{current.segment_id, orientation});
            }

            SegmentSideNode exit_node{current.segment_id, exit_side};

            auto it = adj.find(exit_node);
            if (it == adj.end()) {
                // telomere: no further connection
                break;
            }

            SegmentSideNode next_node = it->second;
            if (visited_nodes.count(next_node)) {
                hit_ring = true;
                break;
            }

            visited_nodes.insert(next_node);
            current = next_node;
        }
    };

    for (const auto& kv : seg_by_id) {
        int start_id = kv.first;
        if (used_segments.count(start_id)) continue;

        SegmentSideNode left_start{start_id, LEFT};
        SegmentSideNode right_start{start_id, RIGHT};

        std::vector<StrandSegment> left_path, right_path;
        bool left_ring = false, right_ring = false;

        walk_from(left_start, left_path, left_ring);
        walk_from(right_start, right_path, right_ring);

        std::reverse(left_path.begin(), left_path.end());

        ChimericStrand strand;
        strand.is_ring = left_ring || right_ring;

        for (const auto& seg : left_path)
            strand.segments.push_back(seg);

        for (size_t i = 0; i < right_path.size(); ++i) {
            if (!strand.segments.empty() &&
                strand.segments.back().segment_id == right_path[i].segment_id)
                continue;
            strand.segments.push_back(right_path[i]);
        }

        for (const auto& seg : strand.segments)
            used_segments.insert(seg.segment_id);

        strands.push_back(std::move(strand));
    }

    return strands;
}



// Class to represent repaired chromosome geometry
class repaired_chromosome {
public:
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
    std::vector<dna_repair_pair> repair_pairs;

    repaired_chromosome(const std::pair<std::string, std::string>& cell);
    void parse_cell_data(const std::string& cell_data);
};


// Structure to store DSB end index mapping and pairing information
struct dsb_end_pairing {
    int dsb_end_index;       // Global index for this DSB end
    int paired_with_index;   // Index of the DSB end it's paired with
    int inter_chrom;         // 0 or 1: whether pairing is inter-chromosomal
    dsb_end_pairing(int idx, int paired, int inter) 
        : dsb_end_index(idx), paired_with_index(paired), inter_chrom(inter) {}
};



// Function to index all DSB ends and return their pairings
std::vector<dsb_end_pairing> index_dsb_ends_and_pairings(
    const std::vector<dna_repair_pair>& repair_pairs) {
    
    std::vector<dsb_end_pairing> pairings;
    int global_index = 0;
    
    // Iterate through all repair pairs
    for (size_t i = 0; i < repair_pairs.size(); ++i) {
        const auto& pair = repair_pairs[i];
        
        // Index for end1
        int end1_index = global_index++;
        // Index for end2
        int end2_index = global_index++;
        
        // Record the pairing: end1 paired with end2
        pairings.push_back(dsb_end_pairing(end1_index, end2_index, pair.inter_chrom));
        pairings.push_back(dsb_end_pairing(end2_index, end1_index, pair.inter_chrom));
    }
    
    return pairings;
}





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



int main_notused() {
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

        // Index DSB ends and get pairing information
        auto dsb_pairings = index_dsb_ends_and_pairings(geom.repair_pairs);
        cout << "\nDSB End Pairings:" << endl;
        for (const auto& pairing : dsb_pairings) {
            // Determine which repair pair and which end this index corresponds to
            int pair_index = pairing.dsb_end_index / 2;
            int end_in_pair = pairing.dsb_end_index % 2;  // 0 for end1, 1 for end2
            
            const dsb_end& current_end = (end_in_pair == 0) 
                ? geom.repair_pairs[pair_index].end1 
                : geom.repair_pairs[pair_index].end2;
            
            // Determine the paired end
            int paired_pair_index = pairing.paired_with_index / 2;
            int paired_end_in_pair = pairing.paired_with_index % 2;
            
            const dsb_end& paired_end = (paired_end_in_pair == 0) 
                ? geom.repair_pairs[paired_pair_index].end1 
                : geom.repair_pairs[paired_pair_index].end2;
            
            cout << "  DSB End Index " << pairing.dsb_end_index 
                 << " (Chr: " << current_end.chromosome_set 
                 << ", Pos: " << current_end.position_in_chromosome
                 << ", Up/Down: " << current_end.upstream_downstream << ")"
                 << " is paired with Index " << pairing.paired_with_index
                 << " (Chr: " << paired_end.chromosome_set
                 << ", Pos: " << paired_end.position_in_chromosome
                 << ", Up/Down: " << paired_end.upstream_downstream << ")"
                 << " (Inter-chrom: " << pairing.inter_chrom << ")" << endl;
        }

        // Get affected chromosomes and chromatids
        auto affected = get_affected_chromosomes(geom.repair_pairs);
        cout << "\nAffected Chromosomes and Chromatid Strands:" << endl;
        for (const auto& chr_chromatid : affected) {
            cout << "  Chromosome " << chr_chromatid.chromosome_set 
                 << ", Chromatid Strand " << chr_chromatid.chromatid_strand << endl;
            
            // Get broken segments for this chromosome-chromatid combination
            auto broken_segs = get_broken_segments(chr_chromatid, geom.repair_pairs);
            cout << "    Broken Segments:" << endl;
            for (size_t k = 0; k < broken_segs.size(); ++k) {
                cout << "      Segment " << (k + 1) << ": [" 
                     << broken_segs[k].old_chromatid_start_position << ", "
                     << broken_segs[k].old_chromatid_end_position << "] "
                     << "(Centromere: " << broken_segs[k].has_centromere << ")" << endl;
            }
        }

        auto segments_by_cc = build_all_segments(geom.repair_pairs);
        auto side_map       = build_segment_side_map(segments_by_cc);
        auto adj            = build_adjacency(geom.repair_pairs, side_map);
        auto strands        = build_chimeric_strands(segments_by_cc, adj);

        cout << "\nChimeric strands:\n";

        for (size_t s = 0; s < strands.size(); ++s) {
            const auto& strand = strands[s];
            cout << "  Strand " << s + 1 << " (" << (strand.is_ring ? "ring" : "linear") << "):\n";

            // rebuild seg_by_id for printing
            std::unordered_map<int, segment> seg_by_id;
            for (const auto& kv : segments_by_cc)
                for (const auto& seg : kv.second)
                    seg_by_id[seg.id] = seg;

            for (const auto& seg_inst : strand.segments) {
                const auto& seg = seg_by_id[seg_inst.segment_id];
                cout << "    SegID " << seg.id
                    << " Chr " << seg.old_chromosome_id
                    << " Strand " << seg.old_chromatid_strand
                    << " [" << seg.old_chromatid_start_position
                    << ", " << seg.old_chromatid_end_position << "]"
                    << " Centromere " << seg.has_centromere
                    << " Orientation " << (seg_inst.orientation > 0 ? "+" : "-")
                    << "\n";
            }
        }

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
        }*/

        cout << "------------------------------" << endl;
    }

    return 0;
}




int main() {
    cout << "=== Synthetic Multi‑Segment Ring Test ===\n";

    // Helper to build DSB ends
    auto make_end = [](int break_index,
                       int chrom,
                       double pos,
                       int updown) -> dsb_end {
        dsb_end e{};
        e.break_index = break_index;
        e.is_complex = 0;
        e.dna_structure = 0;
        e.chromosome_set = chrom;
        e.chromatid_strand = 1;
        e.chromosome_arm = 0;
        e.position_in_chromosome = pos;
        e.upstream_downstream = updown;
        return e;
    };

    // Breakpoints for 6‑segment ring
    std::vector<double> bp = {0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 0.95};

    // For each breakpoint, create two ends:
    // upstream = -1, downstream = +1
    // We will pair downstream of i with upstream of i+1 (mod N)
    std::vector<dna_repair_pair> repair_pairs;
    int N = bp.size();

    for (int i = 0; i < N; i++) {
        int j = (i + 1) % N;  // next breakpoint in ring

        dna_repair_pair p;
        p.inter_chrom = 0;

        // downstream end of breakpoint i
        p.end1 = make_end(2*i,   1, bp[i], +1);

        // upstream end of breakpoint j
        p.end2 = make_end(2*j+1, 1, bp[j], -1);

        repair_pairs.push_back(p);
    }

    // Print DSB end pairings
    auto dsb_pairings = index_dsb_ends_and_pairings(repair_pairs);
    cout << "\nDSB End Pairings:\n";
    for (const auto& pairing : dsb_pairings) {
        int pair_index = pairing.dsb_end_index / 2;
        int end_in_pair = pairing.dsb_end_index % 2;

        const dsb_end& current_end =
            (end_in_pair == 0) ? repair_pairs[pair_index].end1
                               : repair_pairs[pair_index].end2;

        int paired_pair_index = pairing.paired_with_index / 2;
        int paired_end_in_pair = pairing.paired_with_index % 2;

        const dsb_end& paired_end =
            (paired_end_in_pair == 0) ? repair_pairs[paired_pair_index].end1
                                      : repair_pairs[paired_pair_index].end2;

        cout << "  DSB End Index " << pairing.dsb_end_index
             << " (Pos: " << current_end.position_in_chromosome
             << ", Up/Down: " << current_end.upstream_downstream << ")"
             << " paired with " << pairing.paired_with_index
             << " (Pos: " << paired_end.position_in_chromosome
             << ", Up/Down: " << paired_end.upstream_downstream << ")\n";
    }

    // Build segments
    auto affected = get_affected_chromosomes(repair_pairs);
    cout << "\nAffected Chromosomes:\n";
    for (auto& cc : affected) {
        cout << "  Chr " << cc.chromosome_set
             << ", Strand " << cc.chromatid_strand << "\n";

        auto segs = get_broken_segments(cc, repair_pairs);
        cout << "    Segments:\n";
        for (auto& s : segs) {
            cout << "      [" << s.old_chromatid_start_position
                 << ", " << s.old_chromatid_end_position << "]\n";
        }
    }

    // Build chimeric strands
    auto segments_by_cc = build_all_segments(repair_pairs);
    auto side_map       = build_segment_side_map(segments_by_cc);
    auto adj            = build_adjacency(repair_pairs, side_map);
    auto strands        = build_chimeric_strands(segments_by_cc, adj);

    cout << "\nChimeric strands:\n";

    std::unordered_map<int, segment> seg_by_id;
    for (const auto& kv : segments_by_cc)
        for (const auto& seg : kv.second)
            seg_by_id[seg.id] = seg;

    for (size_t s = 0; s < strands.size(); ++s) {
        const auto& strand = strands[s];
        cout << "  Strand " << s + 1 << " ("
             << (strand.is_ring ? "ring" : "linear") << "):\n";

        for (const auto& seg_inst : strand.segments) {
            const auto& seg = seg_by_id[seg_inst.segment_id];
            cout << "    SegID " << seg.id
                 << " [" << seg.old_chromatid_start_position
                 << ", " << seg.old_chromatid_end_position << "]"
                 << " Orientation " << (seg_inst.orientation > 0 ? "+" : "-")
                 << "\n";
        }
    }

    return 0;
}