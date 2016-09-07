# ProcessSomaticMutation
This program is used to process and return a total count of all the somatic mutations in a given patient file.

Goal for this program:
Take commandline arguments for cancer type folder taken from TCGA website
and process the file_manifest.txt, FILE_SAMPLE_MAP.txt, and .maf files

Program will take input of .maf files and create an output for the total count
of somatic mutations per patient per gene, also the total count of mutations
minus silent mutations.

Options to compare old data against newly created data, process entire folders into one file, or run it on a single file.
