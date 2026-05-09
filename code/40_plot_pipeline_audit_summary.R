# No external R packages needed

infile <- "audit_results/pipeline_audit_summary.csv"
outdir <- "audit_results"

df <- read.csv(infile, stringsAsFactors = FALSE)

# Remove GEO reference row
plot_df <- df[df$overlap_with_GEO != "reference", ]

plot_df$overlap_with_GEO <- as.numeric(plot_df$overlap_with_GEO)
plot_df$mean_spearman_vs_GEO <- as.numeric(gsub("~", "", plot_df$mean_spearman_vs_GEO))

# Short labels
plot_df$label <- c(
  "STAR v38\nprimary-only",
  "STAR v37\nprimary-only",
  "STAR v37\n-M primary",
  "STAR v37\n-M -O fraction",
  "Bowtie1 v37\n-M -O fraction"
)

# -------------------------------
# Plot 1: overlap with GEO
# -------------------------------
png(file.path(outdir, "pipeline_audit_overlap_with_GEO.png"),
    width = 1800, height = 1100, res = 200)

par(mar = c(8, 5, 4, 2))

barplot(
  plot_df$overlap_with_GEO,
  names.arg = plot_df$label,
  las = 2,
  ylab = "Number of lncRNA IDs overlapping GEO",
  main = "Pipeline audit: lncRNA ID overlap with GEO processed matrix"
)

dev.off()

pdf(file.path(outdir, "pipeline_audit_overlap_with_GEO.pdf"),
    width = 9, height = 5.5)

par(mar = c(8, 5, 4, 2))

barplot(
  plot_df$overlap_with_GEO,
  names.arg = plot_df$label,
  las = 2,
  ylab = "Number of lncRNA IDs overlapping GEO",
  main = "Pipeline audit: lncRNA ID overlap with GEO processed matrix"
)

dev.off()

# -------------------------------
# Plot 2: Spearman correlation
# -------------------------------
png(file.path(outdir, "pipeline_audit_spearman_vs_GEO.png"),
    width = 1800, height = 1100, res = 200)

par(mar = c(8, 5, 4, 2))

barplot(
  plot_df$mean_spearman_vs_GEO,
  names.arg = plot_df$label,
  las = 2,
  ylim = c(0, 0.6),
  ylab = "Mean Spearman correlation vs GEO",
  main = "Pipeline audit: count similarity to GEO processed matrix"
)

dev.off()

pdf(file.path(outdir, "pipeline_audit_spearman_vs_GEO.pdf"),
    width = 9, height = 5.5)

par(mar = c(8, 5, 4, 2))

barplot(
  plot_df$mean_spearman_vs_GEO,
  names.arg = plot_df$label,
  las = 2,
  ylim = c(0, 0.6),
  ylab = "Mean Spearman correlation vs GEO",
  main = "Pipeline audit: count similarity to GEO processed matrix"
)

dev.off()

cat("Saved plots:\n")
cat(file.path(outdir, "pipeline_audit_overlap_with_GEO.png"), "\n")
cat(file.path(outdir, "pipeline_audit_overlap_with_GEO.pdf"), "\n")
cat(file.path(outdir, "pipeline_audit_spearman_vs_GEO.png"), "\n")
cat(file.path(outdir, "pipeline_audit_spearman_vs_GEO.pdf"), "\n")

cat("\nPlot data:\n")
print(plot_df[, c("label", "overlap_with_GEO", "mean_spearman_vs_GEO")])

