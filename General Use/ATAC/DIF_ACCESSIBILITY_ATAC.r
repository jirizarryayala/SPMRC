###############################################################################
# STANDALONE DIFFERENTIAL ACCESSIBILITY ANALYSIS FOR ATAC_143
# Picks up from Step 8 - Uses EXACT same code as original pipeline
# Run with: Rscript DA_ATAC_143_standalone.r
###############################################################################

## -------------------- LIBRARY PATHS (same as main pipeline) ----------------
lib_paths <- "/lustre/project/crosslin/crosslin_team/Sharifi/LibC"
.libPaths(c(lib_paths, .libPaths()))
cat("===== .libPaths() =====\n"); print(.libPaths()); cat("\n")

## -------------------- LOAD PACKAGES -----------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

## -------------------- SETUP PATHS (identical to main pipeline) -------------
outs_dir <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/ATAC_143"
plotdir <- file.path(outs_dir, "plots")
tabdir  <- file.path(outs_dir, "tables")
rdsdir  <- file.path(outs_dir, "rds")

## -------------------- LOGGING -----------------------------------------------
logfile <- file.path(outs_dir, "DA_analysis.log")
sink(logfile, split = TRUE)
cat("===== STANDALONE DA ANALYSIS - ATAC_143 =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("R.version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")

## -------------------- HELPER FUNCTION (same as main) -----------------------
save_plot <- function(p, stem, w=10, h=6) {
  pdff <- file.path(plotdir, paste0(stem, ".pdf"))
  grDevices::pdf(pdff, width=w, height=h)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  cat("Saved:", pdff, "\n")
}

###############################################################################
# LOAD OBJECT FROM STEP 8
###############################################################################
cat("\n=== Loading ATAC object from Step 8 ===\n")

rds_file <- file.path(rdsdir, "atac_07_labeled.rds")

if (!file.exists(rds_file)) {
  stop("ERROR: Cannot find ", rds_file, "\n",
       "Make sure main pipeline completed through Step 8!\n",
       "Expected file: ", rds_file)
}

pbmc <- readRDS(rds_file)
cat("Successfully loaded object:\n")
print(pbmc)
cat("\n")

###############################################################################
# STEP 9: DIFFERENTIAL ACCESSIBILITY (EXACT CODE FROM ORIGINAL PIPELINE)
###############################################################################
cat("\n=== STEP 9: Differential accessibility ===\n")

DefaultAssay(pbmc) <- "ATAC"

if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  Idents(pbmc) <- pbmc$predicted.id
  cell_counts <- sort(table(Idents(pbmc)), decreasing = TRUE)
  ident_1 <- names(cell_counts)[1]
  ident_2 <- names(cell_counts)[2]
} else {
  Idents(pbmc) <- pbmc$seurat_clusters
  ident_1 <- levels(Idents(pbmc))[1]
  ident_2 <- levels(Idents(pbmc))[2]
}

cat("Comparing:", ident_1, "vs", ident_2, "\n")
cat("Group 1 (", ident_1, "):", sum(Idents(pbmc) == ident_1), "cells\n")
cat("Group 2 (", ident_2, "):", sum(Idents(pbmc) == ident_2), "cells\n")
cat("\nStarting DA analysis...\n")
cat("This may take 30-60 minutes for large datasets.\n")
cat("Progress will be shown below:\n\n")

# EXACT SAME CODE AS ORIGINAL - with added downsampling to prevent choking
da_peaks <- FindMarkers(
  pbmc, 
  ident.1 = ident_1, 
  ident.2 = ident_2,
  test.use = "LR",                   # Same method
  latent.vars = "nCount_ATAC",       # Same covariate correction
  max.cells.per.ident = 2000,        # Downsample to prevent timeout
  min.pct = 0.05,                    # Only test peaks in >5% cells
  min.cells.feature = 3,             # Minimum cells per peak
  min.cells.group = 3,               # Minimum per group
  verbose = TRUE                     # Show progress
)

cat("\n=== DA Analysis Complete ===\n")

da_peaks <- da_peaks[order(da_peaks$p_val_adj), ]
write.csv(da_peaks, file.path(tabdir, "DA_peaks.csv"))

cat("Total peaks tested:", nrow(da_peaks), "\n")
cat("Significant peaks (padj < 0.05):", sum(da_peaks$p_val_adj < 0.05), "\n")
cat("Significant peaks (padj < 0.01):", sum(da_peaks$p_val_adj < 0.01), "\n")
cat("\nResults saved to:", file.path(tabdir, "DA_peaks.csv"), "\n")

# Show top 10 peaks
cat("\nTop 10 differential peaks:\n")
print(head(da_peaks, 10))

# Visualize top peak (EXACT SAME CODE AS ORIGINAL)
cat("\n=== Generating visualization for top peak ===\n")
top_peak <- rownames(da_peaks)[1]
cat("Top peak:", top_peak, "\n")

p_vln <- VlnPlot(pbmc, features = top_peak, pt.size = 0.1, idents = c(ident_1, ident_2))
p_umap <- FeaturePlot(pbmc, features = top_peak, reduction = "umap", pt.size = 0.1)
save_plot(p_vln | p_umap, "top_DA_peak", w=14, h=6)

###############################################################################
# ADDITIONAL ANALYSES (Optional - can comment out if too slow)
###############################################################################
cat("\n=== Additional Analyses ===\n")

# Top 5 peaks visualization
cat("Visualizing top 5 differential peaks...\n")
top5_peaks <- head(rownames(da_peaks), 5)

for (i in 1:length(top5_peaks)) {
  peak <- top5_peaks[i]
  tryCatch({
    p_vln <- VlnPlot(pbmc, features = peak, pt.size = 0.1, idents = c(ident_1, ident_2))
    p_umap <- FeaturePlot(pbmc, features = peak, reduction = "umap", pt.size = 0.1)
    save_plot(p_vln | p_umap, paste0("DA_peak_rank", i), w=14, h=6)
  }, error = function(e) {
    cat("Could not plot peak", peak, "\n")
  })
}

# Summary statistics by chromosome
cat("\nSummary by chromosome:\n")
da_sig <- da_peaks[da_peaks$p_val_adj < 0.05, ]
if (nrow(da_sig) > 0) {
  chr_counts <- table(gsub(":.*", "", rownames(da_sig)))
  write.csv(as.data.frame(chr_counts), file.path(tabdir, "DA_peaks_by_chr.csv"))
  print(head(sort(chr_counts, decreasing = TRUE), 10))
}

# Export for downstream analysis
cat("\nExporting results for downstream tools...\n")

# BED format for visualization
if (nrow(da_sig) > 0) {
  bed_data <- data.frame(
    chr = gsub(":.*", "", rownames(da_sig)),
    start = as.integer(gsub(".*:(.*)-.*", "\\1", rownames(da_sig))),
    end = as.integer(gsub(".*-", "", rownames(da_sig))),
    name = rownames(da_sig),
    score = -log10(da_sig$p_val_adj),
    strand = "."
  )
  bed_data <- bed_data[order(bed_data$chr, bed_data$start), ]
  write.table(bed_data, file.path(tabdir, "DA_peaks_significant.bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  cat("BED file saved:", file.path(tabdir, "DA_peaks_significant.bed"), "\n")
}

###############################################################################
# SUMMARY
###############################################################################
cat("\n\n===== DA ANALYSIS COMPLETE =====\n")
cat("Comparison:", ident_1, "vs", ident_2, "\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("  Group 1:", sum(Idents(pbmc) == ident_1), "cells\n")
cat("  Group 2:", sum(Idents(pbmc) == ident_2), "cells\n")
cat("\nResults:\n")
cat("  Total peaks tested:", nrow(da_peaks), "\n")
cat("  Significant (padj < 0.05):", sum(da_peaks$p_val_adj < 0.05), "\n")
cat("  Significant (padj < 0.01):", sum(da_peaks$p_val_adj < 0.01), "\n")
cat("\nOutput files:\n")
cat("  Full results:", file.path(tabdir, "DA_peaks.csv"), "\n")
cat("  BED file:", file.path(tabdir, "DA_peaks_significant.bed"), "\n")
cat("  Plots:", file.path(plotdir, "top_DA_peak.pdf"), "\n")
cat("  Log file:", logfile, "\n")

cat("\nCompleted:", format(Sys.time()), "\n")
sink()

cat("\n===== DA ANALYSIS FINISHED SUCCESSFULLY =====\n")
cat("Check", logfile, "for full log\n")
###############################################################################
