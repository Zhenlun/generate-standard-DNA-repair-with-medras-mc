# generate-post-repair-genome-with-medras-mc

This repo builds on top of Medras-MC to output a Standard DNA Repair (SDR) file.

The SDR file contains information on the repaired genome. For more details, see this proposal: https://www.overleaf.com/read/szzqsthxtjpq#1d4e64

The SDR file is meant as an input to RadiSeq sequencing simulator, aiming to detect post-repair DNA structural variants after sequencing.

Functions modified in Medras-MC are:
  - medrasrepair.misrepairSpectrum_withoutput() as a drop in replacement of misrepairSpectrum()
  - The simulation wrapper medrasrepair.repairSimulation() to add functionality to output the SDR file.

Functions added are:
  - medrasrepair.postRepairDNA()
  - The entire standardDnaRepair file.

To output an SDR file, use medrasrepair.repairSimulation() with analysisFunction="postRepairDNA".

status:
  - As of now, SDR file is still in development. Functions and entries will change.
  - SDR file will now print to console with the ipynb notebook.
  - fields arent in their final form, but will only require minimal tweaks.
