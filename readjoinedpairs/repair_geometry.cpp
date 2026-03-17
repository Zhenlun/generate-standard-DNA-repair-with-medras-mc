#include "cell_object.cpp"
#include "load_cells_from_file.cpp"
#include <filesystem>
#include <iostream>
#include <string>

using namespace std;



//list all segments first?


class geometry_file_format { // use the repaired_chromosome class to output each chromosome's geometry
public:
    // header section capturing metadata about the file/generation run
    struct Header {
        std::string author;
        std::string sdd_file_name;
        std::string simulation_detail;
        std::string source;
        std::string incident_particles;
        double mean_particle_energy;
        std::string energy_distribution;
        double dose;
        double dose_rate;
        std::string irradiation_target;
        int old_chromatid_count;
        std::vector<double> old_chromatid_sizes;   // floats
        std::vector<int> intact_strands_id;        // integers
        int damaged_chromatid_count;
        std::vector<double> new_chromatid_sizes;   // floats
        std::string cell_cycle_phase;
        int total_dsb_count;
        int number_of_misrepairs;
        int number_of_interchromosome_misrepairs;
    } header;


    // additional members for chromosome geometry could go here
};
