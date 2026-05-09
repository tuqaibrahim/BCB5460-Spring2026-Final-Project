# Code directory

This folder contains the original numbered scripts used for the reproduction we attempted.

## Script groups

### 1. Data preparation

- `00_validate_fastq_map.py`: checks that the SRR-to-sample mapping is valid.
- `01_merge_fastqs.py`: merges four SRR runs into each biological sample.
- `02_count_fastq_reads.sh`: counts reads in merged FASTQ files.

### 2. QC and trimming

- `03_fastqc_raw_array.slurm`: runs FastQC on raw merged FASTQ files.
- `04_summarize_fastqc_raw.sh`: summarizes raw FastQC output.
- `05_extract_fastqc_raw_details.sh`: extracts detailed FastQC module results.
- `06_cutadapt_and_fastqc_array.slurm`: trims reads with cutadapt and runs post-trim FastQC.
- `07_summarize_cutadapt_logs.sh`: summarizes cutadapt logs.
- `08_summarize_fastqc_trimmed.sh`: summarizes post-trim FastQC.

### 3. STAR branch

- `09_download_gencode_v38.sh`: downloads GENCODE v38 reference.
- `10_build_STAR_index_gencode_v38.slurm`: builds STAR index.
- `11_STAR_align_array.slurm`: aligns trimmed reads with STAR.
- `12_summarize_STAR_logs.sh`: summarizes STAR alignment logs.

### 4. featureCounts and annotation trials

Scripts `13` through `27` test strandedness, GENCODE v38/v37, primary-only counting, multimapping, multi-overlap, and fractional-counting settings.

### 5. Differential expression

- `28_prepare_candidate_counts_for_DE.py`: prepares the best derived count matrix.
- `29_RUVSeq_edgeR_LRT_candidate.R`: runs RUVSeq + edgeR.
- `30_compare_candidate_DE_to_paper.py`: compares derived DE results to paper genes.

### 6. Bowtie1 branch

Scripts `31` through `36` test Bowtie1 alignment and compare it to STAR alignment and GEO-derived counts.

### 7. Audit and figures

Scripts `38` through `41` create audit summaries and paper-style figures.
