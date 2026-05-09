# ============================================================
# FINAL CONSISTENT paper-style figures from OUR derived data
#
# Panels:
#   B) Table of ALL significant lncRNAs
#   C) Volcano plot, labeling the SAME genes shown in Panel E
#   D) Chromosomal location of ALL significant lncRNAs
#   E) Top 5 up + top 5 down normalized-count barplots
#
# Consistency rule:
#   Panel E genes are selected by logFC:
#      top 5 highest positive logFC
#      top 5 most negative logFC
#   Panel C labels exactly those same 10 genes.
#   Panel B contains all significant genes, including all Panel E genes.
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(grid)
})

# -----------------------------
# File paths
# -----------------------------
de_file <- "differential_expression_v37_MO_fraction/RUVSeq_edgeR_LRT_candidate_full_results.csv"
count_file <- "differential_expression_v37_MO_fraction/lncRNA_counts_candidate_rounded.csv"
gene_info_file <- "counts/from_raw_STAR_featureCounts_gencode_v37_M_O_fraction_noPrimary/lncRNA_gene_info_s0.csv"
metadata_file <- "metadata/samples.csv"

outdir <- "figures_paper_style_our_derived_FINAL"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Thresholds
# -----------------------------
LOGFC_CUT <- 2
PVALUE_CUT <- 0.001
FDR_CUT <- 0.05
N_TOP_UPDOWN <- 5

# -----------------------------
# Load data
# -----------------------------
de <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)
counts_df <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_info <- read.csv(gene_info_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

sample_cols <- meta$sample

# Aggregate duplicate IDs if any
counts_df <- aggregate(. ~ ID, data = counts_df[, c("ID", sample_cols)], FUN = sum)
rownames(counts_df) <- counts_df$ID

count_mat <- as.matrix(counts_df[, sample_cols])
storage.mode(count_mat) <- "integer"

# CPM matrix for Panel E
dge <- DGEList(counts = count_mat)
dge <- calcNormFactors(dge)
cpm_mat <- cpm(dge, log = FALSE)

# -----------------------------
# Significant genes
# -----------------------------
de$significance <- "Not significant"
de$significance[
  de$logFC > LOGFC_CUT &
    de$PValue < PVALUE_CUT &
    de$FDR < FDR_CUT
] <- "Up-regulated"

de$significance[
  de$logFC < -LOGFC_CUT &
    de$PValue < PVALUE_CUT &
    de$FDR < FDR_CUT
] <- "Down-regulated"

sig <- de[de$significance != "Not significant", ]

# Up genes: strongest positive logFC first
up_sig <- sig[sig$significance == "Up-regulated", ]
up_sig <- up_sig[order(-up_sig$logFC, up_sig$PValue), ]

# Down genes: most negative logFC first
down_sig <- sig[sig$significance == "Down-regulated", ]
down_sig <- down_sig[order(down_sig$logFC, down_sig$PValue), ]

# Final ordered table:
# strongest up first, then strongest down first
sig_ordered <- rbind(up_sig, down_sig)

# Panel E genes
top5_up <- head(up_sig, N_TOP_UPDOWN)
top5_down <- head(down_sig, N_TOP_UPDOWN)

panelE_genes <- c(top5_up$ID, top5_down$ID)

cat("Total significant genes:\n")
print(nrow(sig_ordered))

cat("\nUp-regulated significant genes:\n")
print(nrow(up_sig))

cat("\nDown-regulated significant genes:\n")
print(nrow(down_sig))

cat("\nPanel E genes that will also be labeled in Panel C:\n")
print(panelE_genes)

# Save source tables
write.csv(sig_ordered,
          file.path(outdir, "panel_B_all_significant_lncRNAs_ordered.csv"),
          row.names = FALSE)

write.csv(top5_up,
          file.path(outdir, "panel_E_top5_up_genes.csv"),
          row.names = FALSE)

write.csv(top5_down,
          file.path(outdir, "panel_E_top5_down_genes.csv"),
          row.names = FALSE)

# ============================================================
# PANEL B: TABLE OF ALL SIGNIFICANT GENES
# ============================================================

make_table_panel <- function(df, outfile_base) {

  show_df <- data.frame(
    Gene = df$ID,
    Log2FC = sprintf("%.9f", df$logFC),
    PValue = format(df$PValue, digits = 6, scientific = TRUE),
    FDR = format(df$FDR, digits = 6, scientific = TRUE),
    Direction = df$significance,
    stringsAsFactors = FALSE
  )

  nrows <- nrow(show_df)

  draw_table <- function() {
    grid.newpage()

    x_positions <- c(0.12, 0.42, 0.69, 0.88)
    headers <- c("Gene", "Log\u2082 Fold Change\nR+3 vs L-10", "PValue", "FDR")

    # Top/header/bottom lines
    grid.lines(x = unit(c(0.02, 0.98), "npc"),
               y = unit(c(0.98, 0.98), "npc"),
               gp = gpar(lwd = 2))

    grid.lines(x = unit(c(0.02, 0.98), "npc"),
               y = unit(c(0.91, 0.91), "npc"),
               gp = gpar(lwd = 2))

    grid.lines(x = unit(c(0.02, 0.98), "npc"),
               y = unit(c(0.02, 0.02), "npc"),
               gp = gpar(lwd = 2))

    # Headers
    for (i in seq_along(headers)) {
      grid.text(headers[i],
                x = x_positions[i],
                y = 0.945,
                gp = gpar(fontsize = 14, fontface = "bold"))
    }

    # Rows
    y_start <- 0.885
    y_end <- 0.04
    row_step <- (y_start - y_end) / nrows

    for (i in seq_len(nrows)) {
      y <- y_start - (i - 0.5) * row_step

      col_text <- ifelse(
        show_df$Direction[i] == "Up-regulated",
        "darkolivegreen4",
        "firebrick3"
      )

      grid.text(show_df$Gene[i],
                x = x_positions[1],
                y = y,
                gp = gpar(fontsize = 10.5, col = col_text))

      grid.text(show_df$Log2FC[i],
                x = x_positions[2],
                y = y,
                gp = gpar(fontsize = 10.5, col = col_text))

      grid.text(show_df$PValue[i],
                x = x_positions[3],
                y = y,
                gp = gpar(fontsize = 10.5, col = col_text))

      grid.text(show_df$FDR[i],
                x = x_positions[4],
                y = y,
                gp = gpar(fontsize = 10.5, col = col_text))
    }
  }

  png(paste0(outfile_base, ".png"), width = 1700, height = 2200, res = 220)
  draw_table()
  dev.off()

  pdf(paste0(outfile_base, ".pdf"), width = 8.5, height = 11)
  draw_table()
  dev.off()
}

make_table_panel(
  sig_ordered,
  file.path(outdir, "panel_B_all_significant_table_our_derived")
)

# ============================================================
# PANEL C: VOLCANO PLOT
# Labels only Panel E genes
# ============================================================

make_volcano_panel <- function(de, label_genes, outfile_base) {

  de$neglog10P <- -log10(de$PValue)

  point_cols <- rep("gray75", nrow(de))
  point_cols[de$significance == "Up-regulated"] <- "darkolivegreen4"
  point_cols[de$significance == "Down-regulated"] <- "red2"

  label_df <- de[de$ID %in% label_genes, ]

  draw_volcano <- function() {
    par(mar = c(5.2, 5.2, 2.5, 2))

    plot(
      de$logFC,
      de$neglog10P,
      pch = 21,
      bg = point_cols,
      col = "gray25",
      cex = 1,
      xlab = "Fold Change (Log\u2082)",
      ylab = "Significance (-Log\u2081\u2080)",
      main = ""
    )

    abline(v = c(-LOGFC_CUT, LOGFC_CUT), lty = 2, lwd = 1.5)
    abline(h = -log10(PVALUE_CUT), lty = 2, lwd = 1.5)

    # Label the same 10 genes used in Panel E
    text(
      label_df$logFC,
      -log10(label_df$PValue),
      labels = label_df$ID,
      pos = ifelse(label_df$logFC > 0, 4, 2),
      cex = 0.8,
      xpd = TRUE
    )
  }

  png(paste0(outfile_base, ".png"), width = 1500, height = 1100, res = 220)
  draw_volcano()
  dev.off()

  pdf(paste0(outfile_base, ".pdf"), width = 7.5, height = 6)
  draw_volcano()
  dev.off()
}

make_volcano_panel(
  de,
  panelE_genes,
  file.path(outdir, "panel_C_volcano_labeled_by_panel_E_genes")
)

# ============================================================
# PANEL D: CHROMOSOMAL LOCATION OF ALL SIGNIFICANT GENES
# ============================================================

gene_info2 <- gene_info[, c("gene_name", "chrom", "start", "end")]
gene_info2 <- gene_info2[!duplicated(gene_info2$gene_name), ]
colnames(gene_info2)[1] <- "ID"

sig_loc <- merge(sig_ordered, gene_info2, by = "ID", all.x = TRUE)

sig_loc$chrom_clean <- gsub("^chr", "", sig_loc$chrom)

chr_levels <- c(as.character(1:22), "X", "Y", "M", "MT")
sig_loc <- sig_loc[sig_loc$chrom_clean %in% chr_levels, ]
sig_loc$chrom_clean <- factor(sig_loc$chrom_clean, levels = chr_levels)

write.csv(sig_loc,
          file.path(outdir, "panel_D_chromosomal_location_source_data.csv"),
          row.names = FALSE)

make_chromosomal_location_panel <- function(df, outfile_base) {

  point_cols <- ifelse(df$significance == "Up-regulated",
                       "darkolivegreen4", "red2")

  yvals <- df$start / 1e6
  xvals <- as.numeric(df$chrom_clean)

  draw_chr <- function() {
    par(mar = c(6, 5, 2.5, 2))

    plot(
      xvals,
      yvals,
      pch = 19,
      col = point_cols,
      xaxt = "n",
      xlab = "Chromosome",
      ylab = "Genomic start position (Mb)",
      main = ""
    )

    axis(1, at = seq_along(chr_levels), labels = chr_levels, cex.axis = 0.8)
    grid(col = "gray90")

    text(
      xvals,
      yvals,
      labels = df$ID,
      pos = 3,
      cex = 0.62,
      col = point_cols
    )
  }

  png(paste0(outfile_base, ".png"), width = 1600, height = 1100, res = 220)
  draw_chr()
  dev.off()

  pdf(paste0(outfile_base, ".pdf"), width = 8.5, height = 6)
  draw_chr()
  dev.off()
}

make_chromosomal_location_panel(
  sig_loc,
  file.path(outdir, "panel_D_chromosomal_location_all_significant")
)

# ============================================================
# PANEL E: TOP 5 UP + TOP 5 DOWN BY logFC
# ============================================================

pre_samples <- c("C1_L-10", "C2_L-10", "C3_L-10")
post_samples <- c("C1_R+3", "C2_R+3", "C3_R+3")

get_gene_bar_summary <- function(gene_id) {

  if (!(gene_id %in% rownames(cpm_mat))) {
    return(NULL)
  }

  pre_vals <- as.numeric(cpm_mat[gene_id, pre_samples])
  post_vals <- as.numeric(cpm_mat[gene_id, post_samples])

  data.frame(
    ID = gene_id,
    group = c("L-10", "R+3"),
    mean_cpm = c(mean(pre_vals), mean(post_vals)),
    sem = c(sd(pre_vals) / sqrt(length(pre_vals)),
            sd(post_vals) / sqrt(length(post_vals))),
    stringsAsFactors = FALSE
  )
}

make_top5_bar_panel <- function(up_genes, down_genes, outfile_base) {

  plot_one_gene <- function(gene_id, title_col) {

    sub <- get_gene_bar_summary(gene_id)

    if (is.null(sub)) {
      plot.new()
      title(main = paste(gene_id, "(missing)"))
      return()
    }

    ymax <- max(sub$mean_cpm + sub$sem, na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1

    mids <- barplot(
      sub$mean_cpm,
      names.arg = sub$group,
      col = c("#a6cee3", "#4f81bd"),
      border = NA,
      ylim = c(0, ymax * 1.25),
      ylab = "CPM",
      cex.names = 1.2,
      cex.axis = 1.05,
      cex.lab = 1.1
    )

    arrows(
      x0 = mids,
      y0 = pmax(sub$mean_cpm - sub$sem, 0),
      x1 = mids,
      y1 = sub$mean_cpm + sub$sem,
      angle = 90,
      code = 3,
      length = 0.05,
      lwd = 1.5,
      col = "gray60"
    )

    title(main = gene_id, col.main = title_col, font.main = 2, cex.main = 1.3)
  }

  draw_panel <- function() {
    par(mfrow = c(2, 5), mar = c(4.3, 4.1, 3, 1), oma = c(0, 0, 2.5, 0))

    for (g in up_genes) {
      plot_one_gene(g, "darkolivegreen4")
    }

    for (g in down_genes) {
      plot_one_gene(g, "firebrick3")
    }

    mtext(
      "Normalized counts (CPM) for top 5 up- and down-regulated lncRNAs",
      outer = TRUE,
      cex = 1.35,
      font = 2
    )
  }

  png(paste0(outfile_base, ".png"), width = 2300, height = 1350, res = 220)
  draw_panel()
  dev.off()

  pdf(paste0(outfile_base, ".pdf"), width = 14, height = 8)
  draw_panel()
  dev.off()
}

make_top5_bar_panel(
  top5_up$ID,
  top5_down$ID,
  file.path(outdir, "panel_E_top5_up_down_by_logFC")
)

# ============================================================
# Summary
# ============================================================

sink(file.path(outdir, "summary_FINAL_consistent_figures.txt"))

cat("FINAL consistent paper-style figures from our derived data\n")
cat("=========================================================\n\n")

cat("Thresholds:\n")
cat("abs(logFC) >", LOGFC_CUT, "\n")
cat("PValue <", PVALUE_CUT, "\n")
cat("FDR <", FDR_CUT, "\n\n")

cat("Total significant genes:\n")
print(nrow(sig_ordered))

cat("\nUp-regulated genes:\n")
print(nrow(up_sig))

cat("\nDown-regulated genes:\n")
print(nrow(down_sig))

cat("\nPanel E top 5 up genes, also labeled in Panel C:\n")
print(top5_up[, c("ID", "logFC", "PValue", "FDR")])

cat("\nPanel E top 5 down genes, also labeled in Panel C:\n")
print(top5_down[, c("ID", "logFC", "PValue", "FDR")])

cat("\nPanel B full ordered significant table:\n")
print(sig_ordered[, c("ID", "logFC", "PValue", "FDR", "significance")])

sink()

cat("\nDone. Final consistent figures saved in:\n")
cat(outdir, "\n\n")

cat("Files:\n")
print(list.files(outdir, full.names = TRUE))
