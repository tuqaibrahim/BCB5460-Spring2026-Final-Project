library(edgeR)
library(RUVSeq)

# ==========================
# Paths
# ==========================
base_path <- "/lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline"

input_file <- file.path(base_path, "GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv")

out_all_ruv <- file.path(base_path, "DEGs_lncRNA_RUVSeq_edgeR_all_results_NO_FILTER.csv")
out_paper_genes <- file.path(base_path, "paper_genes_in_RUVSeq_edgeR_results.csv")

# ==========================
# Step 1: Read count file
# ==========================
counts_raw <- read.csv(
  input_file,
  check.names = FALSE
)

# ==========================
# Step 2: Prepare count matrix
# Sum duplicated gene IDs/names
# ==========================
count_cols <- !(colnames(counts_raw) %in% c("ID", "Type"))

counts_only <- counts_raw[, count_cols]
gene_ids <- counts_raw$ID

counts_only[] <- lapply(counts_only, as.numeric)

mtx_data <- rowsum(
  as.matrix(counts_only),
  group = gene_ids
)

cat("Original rows:", nrow(counts_raw), "\n")
cat("Rows after summing duplicate genes:", nrow(mtx_data), "\n")

# ==========================
# Step 3: Metadata
# ==========================
metadata <- data.frame(
  sample = c(
    "C1_L-10", "C1_R+3",
    "C2_L-10", "C2_R+3",
    "C3_L-10", "C3_R+3"
  ),
  time_point = c(
    "L-10", "R+3",
    "L-10", "R+3",
    "L-10", "R+3"
  )
)

stopifnot(all(metadata$sample %in% colnames(mtx_data)))

mtx_data <- mtx_data[, metadata$sample]

group <- factor(metadata$time_point, levels = c("L-10", "R+3"))

# ==========================
# Step 4: First edgeR pass
# This is only to choose empirical control genes for RUVSeq
# ==========================
dge0 <- DGEList(counts = mtx_data, group = group)

# No CPM filtering
dge0 <- calcNormFactors(dge0)

design0 <- model.matrix(~ group)

dge0 <- estimateDisp(dge0, design0)

fit0 <- glmFit(dge0, design0)

lrt0 <- glmLRT(fit0, coef = "groupR+3")

res0 <- topTags(lrt0, n = Inf)$table

# Empirical negative controls:
# genes that look least different in first edgeR pass
control_genes <- rownames(res0)[res0$PValue > 0.5]

cat("Number of empirical control genes for RUVSeq:", length(control_genes), "\n")

# If too few controls, relax cutoff
if (length(control_genes) < 100) {
  control_genes <- rownames(res0)[res0$PValue > 0.2]
  cat("Relaxed control genes cutoff to PValue > 0.2\n")
  cat("Number of empirical control genes:", length(control_genes), "\n")
}

# ==========================
# Step 5: RUVSeq correction
# ==========================
set <- newSeqExpressionSet(
  as.matrix(mtx_data),
  phenoData = data.frame(
    group = group,
    row.names = colnames(mtx_data)
  )
)

# Upper-quartile normalization for RUVSeq
set <- betweenLaneNormalization(set, which = "upper")

# Estimate unwanted variation
# k = 1 is safest with only 6 samples
set_ruv <- RUVg(
  set,
  control_genes,
  k = 1
)

# Extract W_1 unwanted variation factor
W <- pData(set_ruv)$W_1

cat("RUVSeq W_1 values:\n")
print(W)

# ==========================
# Step 6: edgeR after RUVSeq
# ==========================
dge <- DGEList(counts = mtx_data, group = group)

# No CPM filtering
dge <- calcNormFactors(dge)

design_ruv <- model.matrix(~ W + group)

dge <- estimateDisp(dge, design_ruv)

fit <- glmFit(dge, design_ruv)

lrt <- glmLRT(fit, coef = "groupR+3")

degs_ruv <- topTags(lrt, n = Inf)$table

# Add gene names as first column
degs_ruv$gene_id <- rownames(degs_ruv)

degs_ruv <- degs_ruv[, c(
  "gene_id",
  setdiff(colnames(degs_ruv), "gene_id")
)]

# Save ALL results, no filtering
write.csv(
  degs_ruv,
  out_all_ruv,
  row.names = FALSE
)
# ==========================
# Step 7: Check paper genes
# ==========================
paper_genes <- c(
  "LINC02049",
  "AL031767.1",
  "AC244131.2",
  "LINC01954",
  "AL732292.2",
  "AC104561.3",
  "AL161935.3",
  "AC005008.2",
  "AC243960.1",
  "AP003486.1",
  "AL031432.4",
  "AC139769.2",
  "AC098654.1",
  "LINC02511",
  "AL354919.1",
  "LINC00891",
  "AL049796.1",
  "BX284668.2",
  "AL592295.4",
  "LINC01783",
  "DLX2-DT",
  "TAB2-AS1",
  "LINC00324",
  "AC026469.1",
  "LINC02022",
  "AC024651.2",
  "AC129926.1"
)

paper_found <- degs_ruv[degs_ruv$gene_id %in% paper_genes, ]

write.csv(
  paper_found,
  out_paper_genes,
  row.names = FALSE
)

cat("\nDone\n")
cat("All RUVSeq + edgeR results saved to:\n", out_all_ruv, "\n")
cat("Paper genes found in RUVSeq + edgeR results saved to:\n", out_paper_genes, "\n")

cat("\nPaper genes in image:", length(paper_genes), "\n")
cat("Paper genes found in RUVSeq + edgeR all results:", nrow(paper_found), "\n")

cat("\nPaper genes found:\n")
print(paper_found[, c("gene_id", "logFC", "PValue", "FDR")])