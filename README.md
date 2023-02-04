# Aligner_Benchmarking

Additional file 2 was used to list all the FLAG scores in a SAM file. This was used to help write additional file 6.
Additional file 3 was used to make SHRiMP1 output compairable to SAM files.
Additional file 4 was used to make Maq output compairable to SAM files
Additional file 5 was used to compare shrimp1 fixed output files or Maq fixed output files to SAM files. Both files must be sorted with linux sort command first.
Additional file 6 was used to compare SAM files to SAM files. Both files must be sorted with linux sort command first.
parseSam.py was used to count how often 0, 1, 2, 3, 4, 5 or 6 mismatches and/or indels occured in the sequencing reads
countNotAlignedReads.py was used to sort Bowtie2 multiple mapping reads in to categories; primary correct by chance, primary correct by more accurately,
  primary incorrect because of less accurately, secondary correct by chance, secondary correct by more accurately, secondary incorrect because of less 
  accurately, primary and secondary incorrect due to errors, primary and secondary incorrect without errors.
