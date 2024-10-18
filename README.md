# metabarcoding-pipeline
Code to process metabarcoding results

This contains the code I used to process the raw data from an illumina novaseq 6000 run.
The data was processed with HPC on a linux server.

For step 1, sabre was manually installed on the server and used to cut the adaptors. I used a slurm job submission. This means the entire text file was submitted as a .slurm file.
The output provides split files for the different sample names associated with the tags.

For step 2, anaconda was used to use cutadapt in python. I am unaware of a better solution. For this I used interactive scheduling.

For step 3, I submitted a job that opened R and ran the dada2_pipeline.R script.

The final output is a csv table with ASV numbers, read counts per sample, and the classifications.

Annotation cleanup is used to filter and correct names using the WoRMS database.

Statistical analysis contains R scripts that can be used to analyze the csv file.
