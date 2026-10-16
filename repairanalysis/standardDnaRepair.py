"""
This file is the extension to Medras-MC to output a file containing information on the post-repair DNA.

The post-repair DNA should be a map that maps intact DNA's fragments to create new, aberrant DNA strands.

This extension should create a new "Standard for DNA Repair" (SDR) file per SDD file, in the same directory of the SDD file.

The SDR file should end in ".sdr" for the file extension.

The SDR file is inspired by the SDD file format.
"""


from dataclasses import field
from email import header

import numpy as np
import os

keep_intact_strand = False
print_SDR_to_console = False
save_SDR_to_cwd = True

class Master_header:
    def __init__(self):
        self.author = ""
        self.associated_sdd_file = ""
        self.cell_cycle_phase = ""
        self.old_chromosome_sizes = []

    def add_sdd_file_name(self, sddFileName):
        self.associated_sdd_file = sddFileName

    def add_old_chromosome_sizes(self, old_chromosome_sizes):
        self.old_chromosome_sizes = old_chromosome_sizes

    def add_author(self, author):
        self.author = author

    def add_cell_cycle_phase(self, phase):
        if type(phase) == float or type(phase) == int:
            self.cell_cycle_phase = f"{str(int(phase))}"
        elif type(phase) == list:
            ph, progress = phase
            self.cell_cycle_phase = f"{str(int(ph))}, {str(float(progress))}"

    def list_of_header_lines(self):
        # return a list of strings, each string is a line in the header section of the SDR file.
        header_lines = []
        header_lines.append("SDR Version, SDRv2.2;")
        header_lines.append(f"Author, {self.author};")
        header_lines.append(f"Associated SDD File, {self.associated_sdd_file};")
        header_lines.append(f"Cell Cycle Phase, {self.cell_cycle_phase};")
        header_lines.append(f"Intact Chromosome Sizes, {self.old_chromosome_sizes};")
        header_lines.append("***end of master header***")
        return header_lines
    
    def __str__(self):
        # print each header line on a new line
        return "\n".join(self.list_of_header_lines())


class Repaired_cell:
    '''
    contains the per-cell header including information on:
    - size of post-repair chromosomes
    - ID of intact strands
    - total DSB count
    - total misrepair count

    also contains all the DNA strands in the cell, which are represented as Strand objects.
    '''
    def __init__(self, cell_id):

        #per cell header fields:
        self.cell_id = cell_id
        self.new_chromosome_sizes = {}
        self.intact_strands_ID = []
        self.total_dsb_count = 0
        self.total_misrepair_count = 0
        self.medras_log_fields = [] # this will contain the fields from the medras log for this cell, which contains the misrepair spectrum information. we can use this to add more fields to the header if needed.
        
        # DNA strands in the cell:
        self.strands = []
    
    def add_strand(self, strand):
        self.strands = self.strands + [strand]

    def calculate_new_chromosome_sizes(self):
        # calculate the new chromosome sizes
        new_chromosome_sizes = {}
        for strand in self.strands:
            chromID = strand.new_strand_ID
            frag_size = sum(abs(fragment.end - fragment.start) for fragment in strand.fragments)
            if chromID in new_chromosome_sizes:
                new_chromosome_sizes[chromID] += frag_size
            else:
                new_chromosome_sizes[chromID] = frag_size
        self.new_chromosome_sizes = new_chromosome_sizes

        #sanity check
        sanity_check = False
        if sanity_check:
            print(len(self.strands) == len(new_chromosome_sizes))

    # def add_single_intact_strands_ID(self, intact_strands_ID):
    #     self.intact_strands_ID = self.intact_strands_ID + [intact_strands_ID]

    def add_all_intact_strands_ID(self, list_of_intact_strands_ID):
        if len(list_of_intact_strands_ID) > 0:
            self.intact_strands_ID = list_of_intact_strands_ID
        elif len(list_of_intact_strands_ID) == 0:
            self.intact_strands_ID = list(range(1, 47))
    
    def add_total_dsb_count(self, total_dsb_count):
        self.total_dsb_count = total_dsb_count
    
    def add_total_misrepair_count(self, total_misrepair_count):
        self.total_misrepair_count = total_misrepair_count

    def __str__(self):
        start_of_header = f"\n***subheader - cell {self.cell_id}***"
        # header fields, each field is a separate line:
        header_str = f"""Cell ID, {self.cell_id};
Mutated Chromosome Sizes, {self.new_chromosome_sizes};
Intact Strands ID, {len(self.intact_strands_ID)}, {', '.join(map(str, self.intact_strands_ID))};
Total DSB Count, {self.total_dsb_count};
Total Misrepair Count, {self.total_misrepair_count}; \n"""
        header_medras_log_str = f"Medras-MC Log, {', '.join(self.medras_log_fields)};"
        header_str = header_str + header_medras_log_str

        start_of_data = f"\n***data - cell {self.cell_id}***"
        # strand fields, each strand is a separate line:
        strands_str = "\n".join([str(strand) for strand in self.strands])

        return f"{start_of_header}\n{header_str}\n{start_of_data}\n{strands_str}"


    

class Fragment:
    def __init__(self, old_strand_id, start, end, has_centromere):
        # strand ID is the ID of the old, intact strand.
        self.old_strand_id = old_strand_id
        self.start = start
        self.end = end
        self.has_centromere = has_centromere # 1 is true, 0 is false

    def __str__(self):
        # Formats: ID/Start/End/Centromere
        return f"{self.old_strand_id}/{self.start}/{self.end}/{self.has_centromere}"


class Strand:
    # a line in the data section. 
    def __init__(self, cell_ID, new_strand_ID, is_linear):
        self.cell_ID = cell_ID
        self.new_strand_ID = new_strand_ID
        self.fragments = []
        self.is_linear = 1 if is_linear else 0
    
    def add_fragment(self, old_strand_id, start, end, has_centromere):
        new_fragment = Fragment(old_strand_id, start, end, has_centromere)
        self.fragments = self.fragments + [new_fragment]

    def write_data_line(self):
        # 1. Join fragments with a comma
        frag_string = ", ".join([str(f) for f in self.fragments])
        # 2. Join all fields with a semicolon
        # Format: cell; strand; fragments; linear/ring;
        return f"{self.cell_ID}; {self.new_strand_ID}; {frag_string}; {self.is_linear};"
    
    def __str__(self):
        return self.write_data_line()


class outputFromMisrepairSpectrum:
    def __init__(self):
        self.break_set = 0
        self.dsb = 0
        self.remaining_dsb = 0
        self.misrepairs = 0
        self.large_misrepairs = 0
        self.interchromosome_misrepairs = 0
        self.dicentrics = 0
        self.rings = 0
        self.excess_linear_fragments = 0
        self.total_aberrations = 0
        self.viability = 0

    def add_fields(self, fields):
        if len(fields) < 7:
            # theres no misrepair from medras
            self.break_set = fields[0]
            self.dsb = fields[1]
            self.remaining_dsb = fields[2]
            self.misrepairs = fields[3]
        else: 
            self.break_set = fields[0]
            self.dsb = fields[1]
            self.remaining_dsb = fields[2]
            self.misrepairs = fields[3]
            self.large_misrepairs = fields[4]
            self.interchromosome_misrepairs = fields[5]
            self.dicentrics = fields[6]
            self.rings = fields[7]
            self.excess_linear_fragments = fields[8]
            self.total_aberrations = fields[9]
            self.viability = fields[10]

    def __str__(self):
        return f"""
Break Set: {self.break_set}, 
DSB: {self.dsb}, 
Remaining DSB: {self.remaining_dsb}, 
Misrepairs: {self.misrepairs}, 
Large Misrepairs: {self.large_misrepairs}, 
Interchromosome Misrepairs: {self.interchromosome_misrepairs}, 
Dicentrics: {self.dicentrics}, 
Rings: {self.rings}, 
Excess Linear Fragments: {self.excess_linear_fragments}, 
Total Aberrations: {self.total_aberrations}, 
Viability: {self.viability}"""

class Import_cell:
    def __init__(self, cell_id):
        self.cell_id = cell_id
        self.strands = []
        self.intact_strand_ids = []


def import_strands(repair_results):
    cells = []

    for cell_id, (chromosomes, linear_chroms, ring_chroms) in enumerate(repair_results):
        cell = Import_cell(cell_id)

        if keep_intact_strand:
            strand_id = 1
        else:
            strand_id = 47

        # --- Process linear chromosomes ---
        temp_chrom_ID = 1
        for chrom in linear_chroms:
            strand = Strand(cell_id, strand_id, is_linear=True)

            
            intact_strand = False
            for seg in chrom:
                chromID, start, end, _, _ = seg
                chrom_size = chromosomes[chromID][2]
                centromere = 0.5 * chrom_size

                # Check if centromere lies inside this fragment
                has_cent = int(start <= centromere <= end or end <= centromere <= start)
                strand.add_fragment(chromID + 1, start, end, has_cent)

                # Detect intact linear chromosome (handle possible reversed coordinates)
                if len(chrom) == 1 and min(start, end) == 0 and max(start, end) == chrom_size:
                    # cell.intact_strand_ids.append(strand_id)
                    cell.intact_strand_ids.append(chromID + 1)
                    temp_chrom_ID = chromID + 1
                    intact_strand = True

            # check if its intact strand out of the seg loop
            if keep_intact_strand and intact_strand:
                strand.new_strand_ID = temp_chrom_ID
                cell.strands.append(strand)
            elif keep_intact_strand and (not intact_strand):
                temp_chrom_ID += 1
                strand.new_strand_ID = temp_chrom_ID
                cell.strands.append(strand)

            elif (not keep_intact_strand) and intact_strand:
                continue
            else:
                temp_chrom_ID += 1
                cell.strands.append(strand)


            strand_id += 1

        # --- Process ring chromosomes ---
        for chrom in ring_chroms:
            temp_chrom_ID += 1
            strand = Strand(cell_id, temp_chrom_ID, is_linear=False)

            for seg in chrom:
                chromID, start, end, _, _ = seg
                chrom_size = chromosomes[chromID][2]
                centromere = 0.5 * chrom_size

                has_cent = int(start <= centromere <= end or end <= centromere <= start)
                strand.add_fragment(chromID + 1, start, end, has_cent)

            cell.strands.append(strand)
            strand_id += 1

        cells.append(cell)

    return cells



# the final function to be called by medras
def medras_bridge(repair_results, imported_header, sddFileName, captured_logs):
    
    # medras default output:
    first_log_fields = captured_logs[0].split()
    log_first_line_cleaned = " ".join(first_log_fields[18:])
    new_captured_logs = [log_first_line_cleaned] + captured_logs[1:]

    old_chromosome_sizes_str = str(imported_header["Chromosomes"][0]) + ", "
    old_chromosome_sizes_list = imported_header["Chromosomes"][1]
    for size in old_chromosome_sizes_list:
        old_chromosome_sizes_str += str(size) + ", "
    #old_chromosome_sizes_str = ", ".join(str(size) for size in imported_header["Chromosomes"])
    old_chromosome_sizes_str = old_chromosome_sizes_str[:-2]

    # declare classes
    sdr_header = Master_header()
    sdr_header.add_author("Zhenlun Dai")
    sdr_header.add_sdd_file_name(sddFileName)
    sdr_header.add_cell_cycle_phase(imported_header["Cell Cycle Phase"])
    sdr_header.add_old_chromosome_sizes(old_chromosome_sizes_str) # this is old chromosome sizes in MBPs. we can also add the individual chromosome sizes if needed.

    if print_SDR_to_console:
        print(sdr_header)

    if save_SDR_to_cwd:
        SDRfilename = sddFileName[:-4] + ".sdr"
        file = open(SDRfilename, "w")
        file.write(str(sdr_header) + "\n")
    
    # add data to sdr
    imported_cells = import_strands(repair_results)

    # add imported_cells to repaired_cells
    repaired_cells = []
    for cell in imported_cells:
        repaired_cell = Repaired_cell(cell.cell_id)
        repaired_cell.add_all_intact_strands_ID(cell.intact_strand_ids)
        for strand in cell.strands:
            repaired_cell.add_strand(strand)

        #

        # add medras log
        current_cell_captured_log = new_captured_logs[cell.cell_id]
        log_fields = current_cell_captured_log.split()
        misrepair_spectrum = outputFromMisrepairSpectrum()
        misrepair_spectrum.add_fields(log_fields)
        repaired_cell.medras_log_fields = log_fields

        #calculate new chromosome size
        repaired_cell.calculate_new_chromosome_sizes()
        repaired_cell.add_total_dsb_count(misrepair_spectrum.dsb)
        repaired_cell.add_total_misrepair_count(misrepair_spectrum.misrepairs)

        repaired_cells.append(repaired_cell)

        if print_SDR_to_console:
            print(repaired_cell)

        if save_SDR_to_cwd:
            file.write(str(repaired_cell) + "\n")
    
    if save_SDR_to_cwd:
        file.close()
