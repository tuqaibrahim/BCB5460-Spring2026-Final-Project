# BCB5460 Spring 2026 Final Project: Emerging Role of Exosomal Long Non-coding RNAs in Spaceflight-Associated Risks in Astronauts

## Project goal

In this repository I document mt attempt to reproduce the lncRNA RNA-seq analysis from a published astronaut exosome study. 
The authors reported 27 differentially regulated exosomal lncRNAs comparing three days after return with ten days pre-flight baseline. My goal was to reconstruct the computational work from public GEO/SRA data, compare my derived results to the published results, and document where I find difficulties while reproduction.

## Main findings

1. Downloaded and inspected the raw SRA RNA-seq data.
2. Merged SRR runs into biological samples using GEO/SRA metadata.
3. Reconstructed the RNA-seq workflow using FastQC, cutadapt, STAR, Bowtie1, featureCounts, RUVSeq, and edgeR.
4. Tested multiple plausible pipeline branches:
   - GENCODE v38 versus GENCODE v37
   - STAR versus Bowtie1
   - primary-only counting
   - multimapping reads
   - multi-overlap reads
   - fractional counting4. 
5. The pipeline that gave the closest counts to the original study was using STAR aligner + GENCODE v37 + featureCounts -M -O --fraction.
5. Even the closest raw-derived pipeline did not reproduce the same differentially expressed lncRNAs reported in the paper.
6. The conclusion is that the published methods are incomplete for exact computational reproduction.

## Repository organization

- `code/`: original numbered commented scripts used in the analysis. The code is divided into many short scripts to make the debugging and pipline reconstruction easier.
- `data/`: metadata files SRA accession list, sample mapping, and links to raw/processed data.
- `results/`: tables and generated figures.
- `presentation/`: powerpoint presentation slides.
- `Replication attempt using the published counts file/`: contains scripts used to replicate the paper figure from the GEO publised `GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv.gz` file along with the resulting figures.


## How to reproduce the analysis

As you will see, we exmined more than one option and did many checking steps that were not meant to directly replicate the paper results rather than help us figure out the best pipeline that will allow us to reach the same results as the authors.
In order to replicate our work please run scripts in numerical order from the `code/` directory. The main workflow is:

1. Validate and merge FASTQ files:
   - `00_validate_fastq_map.py`
   - `01_merge_fastqs.py`
   - `02_count_fastq_reads.sh`

2. Quality control and trimming:
   - `03_fastqc_raw_array.slurm`
   - `04_summarize_fastqc_raw.sh`
   - `05_extract_fastqc_raw_details.sh`
   - `06_cutadapt_and_fastqc_array.slurm`
   - `07_summarize_cutadapt_logs.sh`
   - `08_summarize_fastqc_trimmed.sh`

3. Reference setup and STAR alignment:
   - `09_download_gencode_v38.sh`
   - `10_build_STAR_index_gencode_v38.slurm`
   - `11_STAR_align_array.slurm`
   - `12_summarize_STAR_logs.sh`

4. featureCounts and annotation trials:
   - `13_featureCounts_strand_test_array.slurm`
   - `14_summarize_featureCounts_strand_test.sh`
   - `15_make_lncRNA_count_matrix_from_featureCounts.py`
   - `16_compare_raw_counts_to_GEO.py`
   - `17_check_GEO_IDs_in_GTF.py`
   - `18_search_gencode_releases_for_missing_ids.sh`
   - `19_featureCounts_gencode_v37_s0.slurm`
   - `20_make_lncRNA_count_matrix_gencode_v37.py`
   - `21_compare_v37_raw_counts_to_GEO.py`
   - `22_featureCounts_v37_multimap_test.slurm`
   - `23_make_lncRNA_count_matrix_gencode_v37_multimap.py`
   - `24_compare_v37_multimap_modes_to_GEO.py`
   - `25_featureCounts_v37_no_primary_tests.slurm`
   - `26_make_lncRNA_count_matrix_gencode_v37_noPrimary.py`
   - `27_compare_all_featureCounts_modes_to_GEO.py`

5. Differential expression:
   - `28_prepare_candidate_counts_for_DE.py`
   - `29_RUVSeq_edgeR_LRT_candidate.R`
   - `30_compare_candidate_DE_to_paper.py`

6. Bowtie1 branch:
   - `31_build_bowtie1_index.slurm`
   - `32_bowtie1_align_array.slurm`
   - `33_summarize_bowtie1_logs.sh`
   - `34_featureCounts_bowtie1_v37_MO_fraction.slurm`
   - `35_make_lncRNA_count_matrix_bowtie1_v37.py`
   - `36_compare_bowtie1_to_GEO.py`

7. Audit and figures:
   - `38_audit_paper_genes_GEO_vs_best_raw.py`
   - `39_make_pipeline_audit_summary.py`
   - `40_plot_pipeline_audit_summary.R`
   - `41_make_paper_like_figures_our_derived.R`
   - `42_make_paper_comparison_up_down_panel_our_data.R`
   - `43_make_paper_style_figures_from_our_derived_data.R`

## Data availability

Raw data are available from SRA/BioProject PRJNA796350. The processed count matrix was obtained from GEO accession GSE193490.

Large FASTQ, BAM, genome index, and intermediate count files are not included in this repository because of file size. 
The repository includes accession IDs, sample mapping, code, final figures, and audit summaries.

## Important reproducibility note

This study could be approximated but not exactly reproduced from our published methods alone. 
The biggest issues were discrepancies in annotation version, aligner choice, featureCounts parameters, and downstream processing choices.

