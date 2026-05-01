"""
This file is the extension to Medras-MC to output a file containing information on the post-repair DNA.

The post-repair DNA should be a map that maps intact DNA's fragments to create new, aberrant DNA strands.

This extension should create a new "Standard for DNA Repair" (SDR) file per SDD file, in the same directory of the SDD file.

The SDR file should end in ".sdr" for the file extension.

The SDR file is inspired by the SDD file format.
"""


from email import header

import numpy as np
import os

class Header:
    def __init__(self):
        self.associated_sdd_file = ""
        self.old_chromosome_sizes = []
        self.intact_strands_ID = []
        self.new_chromosome_sizes = []
        self.number_of_misrepairs = 0
        self.number_of_inter_chromosome_misrepairs = 0

    def add_sdd_file_name(self, sddFileName):
        self.associated_sdd_file = sddFileName

    def add_old_chromosome_sizes(self, old_chromosome_sizes):
        self.old_chromosome_sizes = old_chromosome_sizes

    def add_intact_strands_ID(self, intact_strands_ID):
        number_of_intact_strands = len(intact_strands_ID)
        self.intact_strands_ID = [number_of_intact_strands] + intact_strands_ID

    def add_new_chromosome_sizes(self, new_chromosome_sizes):
        self.new_chromosome_sizes = new_chromosome_sizes

    def add_number_of_misrepairs(self, number_of_misrepairs):
        self.number_of_misrepairs = number_of_misrepairs

    def add_number_of_inter_chromosome_misrepairs(self, number_of_inter_chromosome_misrepairs):
        self.number_of_inter_chromosome_misrepairs = number_of_inter_chromosome_misrepairs

    def list_of_header_lines(self):
        # return a list of strings, each string is a line in the header section of the SDR file.
        header_lines = []
        header_lines.append(f"Associated SDD File: {self.associated_sdd_file}")
        header_lines.append(f"Old Chromosome Sizes: {self.old_chromosome_sizes}")
        header_lines.append(f"Intact Strands ID: {self.intact_strands_ID}")
        header_lines.append(f"New Chromosome Sizes: {self.new_chromosome_sizes}")
        header_lines.append(f"Number of Misrepairs: {self.number_of_misrepairs}")
        header_lines.append(f"Number of Inter-Chromosome Misrepairs: {self.number_of_inter_chromosome_misrepairs}")
        return header_lines
    
    def __str__(self):
        return "\n".join(self.list_of_header_lines())


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


class SDRwriter:
    def __init__(self, sddfilename):
        self.filename = sddfilename.replace(".txt", ".sdr")
        self.header = []
        self.data = []
    
    def add_header(self, list_of_header):
        self.header = list_of_header

    def add_entry(self, entry):
        self.data = self.data + [entry]

    def add_all_entries(self, list_of_entries):
        self.data = list_of_entries

    def write_file(self, destination=None): # default is to write where the script is running
        if destination is None:
            destination = self.filename
        if not destination.endswith(".sdr"):
            raise ValueError("Destination filename must end with .sdr")
        
        with open(destination, 'w') as f:

            # Write Header fields
            for line in self.header:
                f.write(line + ";\n")
            
            f.write("***EndOfHeader***;")

            #write data lines
            for line in self.data:
                f.write(line.write_data_line() + "\n")


def import_strands(repair_results):
    strands = []
    cell_id = 0
    intact_strands_ID = []

    for res in repair_results:
        chromosomes, linearchrom, ringchrom = res
        strandID = 0

        for chrom in linearchrom:
            new_strand = Strand(cell_id, strandID, is_linear=True)
            for seg in chrom:
                chromID, start, end, _, _ = seg
                hasCent = 1 if (start <= 0.5*chromosomes[chromID][2] <= end) or (end <= 0.5*chromosomes[chromID][2] <= start) else 0
                new_strand.add_fragment(chromID, start, end, hasCent)

                # check if this is an intact strand
                # if theres only 1 segment in the chromosome and start is 0 and end is the chromosome size:
                if len(chrom) == 1 and start == 0 and end == chromosomes[chromID][2]:
                    intact_strands_ID.append(strandID)

            strandID += 1
            strands.append(new_strand)

        for chrom in ringchrom:
            new_strand = Strand(cell_id, strandID, is_linear=False)
            for seg in chrom:
                chromID, start, end, _, _ = seg
                hasCent = 1 if (start <= 0.5*chromosomes[chromID][2] <= end) or (end <= 0.5*chromosomes[chromID][2] <= start) else 0
                new_strand.add_fragment(chromID, start, end, hasCent)

            strandID += 1
            strands.append(new_strand)

        cell_id += 1
    return strands, intact_strands_ID


# the final function to be called by medras
def medras_bridge_working(repair_results, imported_header, sddFileName, captured_logs):
    
    # medras default output:
    first_log_fields = captured_logs[0].split()
    log_first_line_cleaned = " ".join(first_log_fields[18:])
    new_captured_logs = [log_first_line_cleaned] + captured_logs[1:]

    # declare classes
    sdr_header = Header()
    SDR = SDRwriter(sddFileName)
    
    # add data to sdr
    all_strands_in_file, intact_strands_ID = import_strands(repair_results)
    SDR.add_all_entries(all_strands_in_file)

    #add header fields to header
    sdr_header.add_sdd_file_name(sddFileName)
    sdr_header.add_old_chromosome_sizes(imported_header["Chromosomes"])
    sdr_header.add_intact_strands_ID(intact_strands_ID)

    #WIP
    # these fields are different per cell. think of a way to represent these per cell
    sdr_header.add_new_chromosome_sizes()
    sdr_header.add_number_of_misrepairs()
    sdr_header.add_number_of_inter_chromosome_misrepairs()


# use this one for tests
def medras_bridge(repair_results, imported_header, sddFileName, captured_logs):
    # ---------------------------------------------------------
    # OUTPUT: Pre-repair segments, repair joins, final strands
    # ---------------------------------------------------------
    print_results = False
    print_header = False
    print_filename = False
    print_captured_logs = True

    # 1. Take the first log and split it into a list of words
    first_log_fields = captured_logs[0].split()

    # 2. Keep only the data (everything after the 18th word)
    # Then join it back into a string so it matches the format of the other logs
    log_first_line_cleaned = " ".join(first_log_fields[18:])

    # 3. Create the new list with the cleaned first line
    new_captured_logs = [log_first_line_cleaned] + captured_logs[1:]

    if print_captured_logs:
        for i, log in enumerate(new_captured_logs):
            print(f"--- Captured Log for Cell {i} ---")
            print(log.split())

    if print_filename:
        print(sddFileName)

    if print_results:
        cell_id = 0
        for res in repair_results:
            chromosomes, linearchrom, ringchrom = res

            strandID = 0
            # Linear chromosomes
            for chrom in linearchrom:
                print(f"STRAND {strandID} LINEAR")
                for seg in chrom:
                    chromID, start, end, _, _ = seg
                    hasCent = 1 if (start <= 0.5*chromosomes[chromID][2] <= end) or (end <= 0.5*chromosomes[chromID][2] <= start) else 0
                    print(f"  SEG {chromID} {start} {end} {hasCent}")
                strandID += 1

            # Ring chromosomes
            for chrom in ringchrom:
                print(f"STRAND {strandID} RING")
                for seg in chrom:
                    chromID, start, end, _, _ = seg
                    hasCent = 1 if (start <= 0.5*chromosomes[chromID][2] <= end) or (end <= 0.5*chromosomes[chromID][2] <= start) else 0
                    print(f"  SEG {chromID} {start} {end} {hasCent}")
                strandID += 1
            cell_id += 1
            
    if print_header:
        print(imported_header)
        print(imported_header["Chromosomes"]) # this is old chromosome sizes in MBPs

#----- added -----



# # --- USAGE EXAMPLE ---

# # 1. Create a Strand Entry (Cell 1, Strand 15, Linear)
# row1 = StrandEntry(cell_id=1, strand_id=15, is_linear=True)
# row1.add_fragment(42, 0, 0.231, 0)
# row1.add_fragment(13, 0.55, 0.45, 1)
# row1.add_fragment(6, 0.35, 0, 0)

# # 2. Create another entry (Cell 2, Strand 7, Ring)
# row2 = StrandEntry(cell_id=2, strand_id=7, is_linear=False)
# row2.add_fragment(32, 0.57, 0.76, 0)

# # 3. Save to file
# writer = SDRWriter("examplesdd.txt")
# writer.add_entry(row1)
# writer.add_entry(row2)

# header = ["Format: SDD-Inspired", "Generated by Python Script", "Version: 1.0"]
# writer.write_file(header)