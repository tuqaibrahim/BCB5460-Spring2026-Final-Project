suppressPackageStartupMessages({
  library(edgeR)
  library(RUVSeq)
})

# ============================================================
# Input/output
# ============================================================
count_file <- "differential_expression_v37_MO_fraction/lncRNA_counts_candidate_rounded.csv"
metadata_file <- "metadata/samples.csv"
outdir <- "differential_expression_v37_MO_fraction"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Read data
# ============================================================
counts_df <- read.csv(count_file, check.names = FALSE)
meta <- read.csv(metadata_file, check.names = FALSE)

sample_cols <- meta$sample

cat("Count matrix file:\n")
cat(count_file, "\n\n")

cat("Input count table dimensions:\n")
print(dim(counts_df))

# Remove duplicated IDs by summing, just in case
counts_df <- aggregate(. ~ ID, data = counts_df[, c("ID", sample_cols)], FUN = sum)

rownames(counts_df) <- counts_df$ID
counts <- as.matrix(counts_df[, sample_cols])
storage.mode(counts) <- "integer"

cat("\nCount matrix dimensions after ID aggregation:\n")
print(dim(counts))

cat("\nMetadata:\n")
print(meta)

# Ensure same order
stopifnot(all(colnames(counts) == meta$sample))

meta$condition <- factor(meta$condition, levels = c("pre", "post"))
meta$astronaut <- factor(meta$astronaut)

cat("\nLibrary sizes:\n")
print(colSums(counts))

# ============================================================
# edgeR object
# ============================================================
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

cat("\nTMM normalization factors:\n")
print(dge$samples$norm.factors)

# ============================================================
# RUVSeq setup
# ============================================================
# RUVSeq expects a SeqExpressionSet
set <- newSeqExpressionSet(
  as.matrix(dge$counts),
  phenoData = data.frame(
    row.names = colnames(dge$counts),
    sample = meta$sample,
    astronaut = meta$astronaut,
    condition = meta$condition,
    timepoint = meta$timepoint
  )
)

# Empirical negative controls:
# Use genes that are not strongly different by a quick initial edgeR test.
design0 <- model.matrix(~ condition, data = meta)
dge0 <- estimateDisp(dge, design0)
fit0 <- glmFit(dge0, design0)
lrt0 <- glmLRT(fit0, coef = "conditionpost")
tab0 <- topTags(lrt0, n = Inf)$table

# Control genes: not strongly associated with condition
controls <- rownames(tab0)[tab0$PValue > 0.5]
cat("\nNumber of empirical control genes for RUVg:\n")
print(length(controls))

if (length(controls) < 100) {
  stop("Too few empirical control genes for RUVg. Need to revisit filtering/control selection.")
}

set_ruv <- RUVg(set, controls, k = 1)

pdat <- pData(set_ruv)

cat("\nRUVSeq phenotype data with W_1:\n")
print(pdat)

# ============================================================
# edgeR with RUV factor
# ============================================================
dge_ruv <- DGEList(counts = counts)
dge_ruv <- calcNormFactors(dge_ruv)

design <- model.matrix(~ W_1 + condition, data = pdat)

cat("\nDesign matrix:\n")
print(design)

dge_ruv <- estimateDisp(dge_ruv, design)

cat("\nCommon dispersion:\n")
print(dge_ruv$common.dispersion)

fit <- glmFit(dge_ruv, design)
lrt <- glmLRT(fit, coef = "conditionpost")

res <- topTags(lrt, n = Inf)$table
res$ID <- rownames(res)

# Put ID first
res <- res[, c("ID", setdiff(colnames(res), "ID"))]

# Save full result
full_out <- file.path(outdir, "RUVSeq_edgeR_LRT_candidate_full_results.csv")
write.csv(res, full_out, row.names = FALSE)

# Paper-like thresholds
sig <- res[
  abs(res$logFC) > 2 &
  res$PValue < 0.001 &
  res$FDR < 0.05,
]

sig <- sig[order(sig$PValue), ]

sig_out <- file.path(outdir, "RUVSeq_edgeR_LRT_candidate_significant_paper_thresholds.csv")
write.csv(sig, sig_out, row.names = FALSE)

cat("\nTop results:\n")
print(head(res, 30))

cat("\nNumber significant using paper-like thresholds:\n")
print(nrow(sig))

cat("\nSignificant genes:\n")
print(sig)

cat("\nSaved:\n")
cat(full_out, "\n")
cat(sig_out, "\n")
