# mpox-apobec3-screening
APOBEC3 Mutation Screening Script
This Python script analyzes annotated variant data (VCF format) against a reference genome and gene annotation (GFF) to identify APOBEC3-related mutational signatures in Mpox virus sequences.

Requirements
Python 3.x
pandas
Biopython
Usage
Edit the file paths in the script or use command-line arguments (future improvement).
Run the script to generate an annotated CSV of variants classified by APOBEC context.

Input files
Reference genome FASTA
Variant VCF file (annotated)
Gene annotation GFF file
Output
CSV file with mutation context and APOBEC classification
Summary statistics printed to console
