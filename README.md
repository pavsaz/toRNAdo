toRNAdo is a script for identification of non-coding RNAs (ncRNAs) from bacterial RNA-Seq data. It identifies both antisense and intergenic ncRNAs as well as those that have both antisense and intergenic properties - "mixed" type. In addition, toRNAdo can also be used to identify 5' and 3' untranslated regions (UTRs) as well as regions that belong to operons.

Prior to running toRNAdo, you must generate appropriate input files with Prepare_for_toRNAdo script. It is found at https://github.com/pavsaz/Prepare_for_toRNAdo

Prepare_for_toRNAdo script will generate a directory for each of your sorted bam files. You then need to run toRNAdo on each of these new directories.

Usage: python toRNAdo.py directory

IMPORTANT: When you specify directory, do not put a forward slash at the end.

Output: Within each of the directories there will now be a new directory named "toRNAdo_output". Within this directory there will be multiple files. Explanation below:

1) .wig files are generated for each type of RNA feature identified. These files include toRNAdo-normalised coverage of each nucleotide in a feature.
  a) Separate files are produced for plus and minus strands.
  b) "filtered" .wig files contain those ncRNAs that have been filtered by the presence of expression "peaks" (>5-fold expression difference between the lowest and highest points of a peak). These are the files you would most likely want.
  c) All files can be imported as a "User plot" to Artemis genome browser and thus visualised.

2) .csv file containing a summary of all identified "filtered" antisense, intergenic and mixed ncRNAs. The "Peak height" column shows the toRNAdo-normalised expression value for a nucleotide with the highest read coverage within the ncRNA.
