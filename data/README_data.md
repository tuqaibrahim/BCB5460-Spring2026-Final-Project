# Data directory

This directory contains metadata and accession information required to reproduce the analysis.

## Files

- `ids.txt`: SRA run accessions.
- `fastq_map.csv`: mapping from SRR runs to biological samples.
- `samples.csv`: paired metadata used for differential expression.
- QC/alignment summary files copied from the analysis.

## External data

Raw FASTQ files should be downloaded from SRA/BioProject PRJNA796350.

The processed lncRNA count matrix should be downloaded from GEO accession GSE193490.

Large raw FASTQ, BAM, and reference files are not included in this repository.
