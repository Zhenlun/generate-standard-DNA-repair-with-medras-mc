# generate-post-repair-genome-with-medras-mc

This repo builds on top of Medras-MC to output a Standard DNA Repair (SDR) file.

The SDR file contains information on the repaired genome. For more details, see this proposal: https://www.overleaf.com/read/szzqsthxtjpq#1d4e64

The SDR file is meant as an input to RadiSeq sequencing simulator, aiming to detect post-repair DNA structural variants after sequencing.

Functions modified in Medras-MC are:
  - medrasrepair.misrepairSpectrum_withoutput() as a drop in replacement of misrepairSpectrum()
  - The simulation wrapper medrasrepair.repairSimulation() to add functionality to output the SDR file.

Functions added are:
  - medrasrepair.postRepairDNA()
  - The entire standardDnaRepair.py file.

To output an SDR file, use medrasrepair.repairSimulation() with analysisFunction="postRepairDNA".

SDRv2.1 Change log:
- Changed headers to use comma and semi-colon to look more similar to SDD format.
- Strand index changed to have 1 to 46 to be always intact strands and 47+ are always mutated strands
- The fragment index are also changed to have the +1, just like strand index
- Fixed a bug that when intact strands are displayed, ring chromosomes would always start at 46 despite it not being the 46th chromosome

SDRv2.0 Change log:
- Added master header - SDR version
- Changed master header - OLD chromosome size to INTACT chromosome size
- Changed subheader - NEW chromosome size to MUTATED chromosome size
- Added toggle at beginning of standardDnaRepair.py to toggle keep intact strand, print to console and save to current work directory
- Changed formatting for subheader entries to not use python arrays and have unified format
- Strand index is no longer continuous. 0 to 45 are always intact strands and 46+ are always mutated strands
