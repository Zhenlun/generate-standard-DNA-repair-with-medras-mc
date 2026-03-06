#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <any>
#include <vector>
#include <queue>
#include <filesystem>

#include <nlohmann/json.hpp>
#include <regex>



std::string readfile() {
    std::string fn;
    return fn;
}

class dsb_end {
    public:
        // build this by looking at SDD file format
        int break_index; //total of break_index * 2 amount of breaks
        int chromosome_ID; //total of 46 in human
        int chromatid_strand; //1 for not dividing, 1 or 2 for G2 phase
        float chromosome_position; //position on the chromosome in percentage (0.0-1.0)
        int upstream_downstream; //-1 for upstream, 1 for downstream

        // Constructor
        dsb_end(std::string s) {
            s = std::regex_replace(s, std::regex(R"(array\(\s*\[)"), "[");
            s = std::regex_replace(s, std::regex(R"(\]\s*\))"), "]");
            s = std::regex_replace(s, std::regex(R"(\bNone\b)"), "null");
            //s = std::regex_replace(s, std::regex(R"(\bTrue\b)"), "true");
            //s = std::regex_replace(s, std::regex(R"(\bFalse\b)"), "false");

            nlohmann::json j = nlohmann::json::parse(s);

            // extract values
            break_index = j[0];

            auto vec = j[3].get<std::vector<int>>();
            if (vec.size() < 4) throw std::runtime_error("bad list1");
            //int v0 = vec[0], v1 = vec[1], v2 = vec[2], v3 = vec[3];
            chromosome_ID = vec[1];
            chromatid_strand = vec[2];

            chromosome_position = j[4].get<double>();
            upstream_downstream = j[5].get<int>();
        }

};

class rdf_file {
    public:
        // create variable to store all lines in the file
        std::queue<std::string> lines;

        // create variables for header and data
        std::map<std::string, std::any> header;
        std::vector<dsb_end> unrepaired_ends;
        std::vector<std::any> repaired_pairs;

        // Constructor
        rdf_file(const std::string& filename) {
            // read the file line by line and store in lines vector
            std::ifstream file(filename);
            if (!file.is_open()) {throw std::runtime_error("Could not open file");}

            std::string line;
            while (std::getline(file, line)) {
                lines.push(line);
            }
            file.close();
        }

        //get the header and data from lines vector
        int read_lines() {

            for (std::string l = lines.front(); !lines.empty(); lines.pop()) {
                std::string line = lines.front();
                // process each line
                // std::cout << line << std::endl;
            
                // Implementation for defining header

                // Implementation for defining unrepaired_ends
                unrepaired_ends.push_back(dsb_end(line));

                // Implementation for defining repaired_pairs

            }
            return 0;
        }

};



int main() {
    // Use current working directory to construct file path
    std::filesystem::path cwd = std::filesystem::current_path();

    //actual line to be used is commented out to for testing
    //std::string filename = (cwd.parent_path() / "n5MeV_inner_alpha_3_SDDOutput.txtn").string();
    std::string filename = (cwd.parent_path() / "test.txtn").string();

    rdf_file rdf(filename);

    rdf.read_lines();
    
    // print all lines of data
    for (const auto& line : rdf.unrepaired_ends) {
        // print each line
        std::cout << line.break_index << ", ";
        std::cout << line.chromosome_ID << ", ";
        std::cout << line.chromatid_strand << ", ";
        std::cout << line.chromosome_position << ", ";
        std::cout << line.upstream_downstream << std::endl;
    }
    

    return 0;
}