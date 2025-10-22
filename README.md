# generate-post-repair-genome-with-medras-mc

aims to output a file containing:
- all the DSB ends before any repair
- DSB ends that werent joined with their original ends and their index
- associated Medras-MC misrepair spectrum information

2025-10-21: now does the first 2 points
2025-10-22: added misrepairSpectrum output into the file. the code now only works with misrepairSpectrum function and 1 cell inside an SDD file.
