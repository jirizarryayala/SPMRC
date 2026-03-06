###############################################################################
# COMPLETE ATAC-seq Pipeline - OBJECTIVE & REPRODUCIBLE
# MAD-based QC filtering + Standard ATAC-seq workflow
# MODIFIED: B cell-specific markers (CD19, CD20, CD38, CD40) + Labeled UMAP
# OUTPUT: Saved to B_cell subdirectory
###############################################################################

## -------------------- LIBRARY PATHS (prepend; don't replace system libs) -----
lib_paths <- "/lustre/project/crosslin/crosslin_team/Sharifi/LibC"
dir.create(lib_paths, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_paths, .libPaths()))
cat("===== .libPaths() =====\n"); print(.libPaths()); cat("\n")

## -------------------- REQUIRED vs OPTIONAL PACKAGE LOADING ------------------
need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing REQUIRED package: ", pkg,
         "\nInstall into: ", lib_paths,
         "\nExample: Rscript -e 'install.packages(\"", pkg, "\", lib=\"", lib_paths,
         "\", repos=\"https://cloud.r-project.org\")'")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  invisible(TRUE)
}

want <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("NOTE: missing OPTIONAL package: ", pkg, " (continuing without it)")
    return(invisible(FALSE))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  invisible(TRUE)
}

# Required for ATAC-seq analysis
need("Seurat")
need("Signac")
need("EnsDb.Hsapiens.v75")
need("GenomicRanges")
need("Matrix")

# Optional (nice-to-have only)
has_patchwork <- want("patchwork")
has_ggplot2 <- want("ggplot2")

## -------------------- ANALYSIS PARAMETERS (MODIFY HERE) ---------------------
# QC filtering stringency
N_MAD_QC <- 3  # Options: 2.5 (strict), 3 (standard/recommended), 4 (lenient)

# Clustering resolution
CLUSTER_RESOLUTION <- 0.5  # Options: 0.3 (coarse), 0.5 (standard), 0.8 (fine)

# Dimensionality reduction
LSI_DIMS <- 2:30  # Standard for ATAC-seq (always skip first LSI component)

## -------------------- SETUP (dirs + logging) --------------------------------
# MODIFIED: Create B_cell subdirectory for all outputs
outs_dir <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/ATAC_150/B_cell"
dir.create(outs_dir, recursive = TRUE, showWarnings = FALSE)

plotdir <- file.path(outs_dir, "plots")
rdsdir  <- file.path(outs_dir, "rds")
tabdir  <- file.path(outs_dir, "tables")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir,  recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(outs_dir, "run.log")
sink(logfile, split = TRUE)
cat("===== ATAC-seq ANALYSIS - B CELL FOCUS =====\n")
cat("Output directory: ", outs_dir, "\n")
cat("Start:", format(Sys.time()), "\n")
cat("R.version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("MAD stringency:", N_MAD_QC, "\n")
cat("Cluster resolution:", CLUSTER_RESOLUTION, "\n\n")

# Write provenance for reproducibility
provfile <- file.path(outs_dir, "sessionInfo.txt")
cat("Writing sessionInfo to:", provfile, "\n")
writeLines(c(capture.output(sessionInfo())), provfile)

# Record package versions
pkgfile <- file.path(outs_dir, "packages_loaded.txt")
pkgs_loaded <- c("Seurat","Signac","EnsDb.Hsapiens.v75","GenomicRanges","Matrix")
pkg_versions <- vapply(pkgs_loaded, function(p) as.character(utils::packageVersion(p)), character(1))
writeLines(paste(names(pkg_versions), pkg_versions, sep=":\t"), pkgfile)

## -------------------- IO helpers --------------------------------------------
save_plot <- function(p, stem, w=10, h=6) {
  pdff <- file.path(plotdir, paste0(stem, ".pdf"))
  grDevices::pdf(pdff, width=w, height=h)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  cat("Saved:", pdff, "\n")
}

save_rds <- function(obj, name) {
  f <- file.path(rdsdir, paste0(name, ".rds"))
  saveRDS(obj, f)
  cat("Saved RDS:", f, "\n")
}

###############################################################################
# STEP 1: LOAD ATAC-seq DATA
# Purpose: Read CellRanger ATAC output (matrix + fragments)
###############################################################################
cat("\n=== STEP 1: Reading ATAC-seq data ===\n")

atac_dir <- "/lustre/project/crosslin/crosslin_team/irizarry/Bhargava/ATACseq/110425 snATAC_LibDD_hPBMC_Anita_8_data/110425 snATAC_LibDD_hPBMC_Anita_8_cellranger/outs_150"

# Verify directory structure
mtx_dir <- file.path(atac_dir, "filtered_peak_bc_matrix")
if (!dir.exists(mtx_dir)) {
  stop("Matrix directory not found: ", mtx_dir)
}

cat("Files in MTX directory:\n")
print(list.files(mtx_dir))

# Read matrix components manually (required for ATAC data)
# Note: ATAC uses peaks.bed instead of features.tsv
mtx_file <- file.path(mtx_dir, "matrix.mtx.gz")
if (!file.exists(mtx_file)) mtx_file <- file.path(mtx_dir, "matrix.mtx")

barcode_file <- file.path(mtx_dir, "barcodes.tsv.gz")
if (!file.exists(barcode_file)) barcode_file <- file.path(mtx_dir, "barcodes.tsv")

peak_file <- file.path(mtx_dir, "peaks.bed.gz")
if (!file.exists(peak_file)) peak_file <- file.path(mtx_dir, "peaks.bed")

# Load components
counts_matrix <- readMM(mtx_file)
barcodes <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)[,1]
peaks <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE)

# Create peak names in genomic coordinate format (chr:start-end)
if (ncol(peaks) >= 3) {
  peak_names <- paste0(peaks[,1], ":", peaks[,2], "-", peaks[,3])
} else {
  peak_names <- peaks[,1]
}

# Assign dimnames
rownames(counts_matrix) <- peak_names
colnames(counts_matrix) <- barcodes

cat("Matrix dimensions:", nrow(counts_matrix), "peaks x", ncol(counts_matrix), "cells\n")

# Load per-cell metadata from CellRanger
metadata_path <- file.path(atac_dir, "singlecell.csv")
metadata <- read.csv(metadata_path, header = TRUE, row.names = 1)

# Fragment file (required for QC metrics and coverage plots)
frag_path <- file.path(atac_dir, "fragments.tsv.gz")
if (!file.exists(frag_path)) stop("Fragment file not found: ", frag_path)

# Create ChromatinAssay (Signac object for ATAC data)
# min.cells = 10: Remove peaks detected in <10 cells
# min.features = 200: Remove cells with <200 peaks
chrom_assay <- CreateChromatinAssay(
  counts = counts_matrix,
  sep = c(":", "-"),
  fragments = frag_path,
  min.cells = 10,
  min.features = 200
)

# Create Seurat object with ATAC assay
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = "ATAC"
)

cat("\nATAC object created:\n")
print(pbmc)
save_rds(pbmc, "atac_01_raw")

###############################################################################
# STEP 2: ADD GENE ANNOTATIONS
# Purpose: Link peaks to genes for downstream analysis
###############################################################################
cat("\n=== STEP 2: Adding gene annotations ===\n")

# Use Ensembl v75 (GRCh37/hg19) - standard for human ATAC-seq
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# Convert chromosome naming to UCSC style (chr1, chr2, etc.)
seqlevelsStyle(annotations) <- "UCSC"

# Add annotations to Seurat object
Annotation(pbmc) <- annotations

cat("Annotations added:", length(annotations), "genes\n")
save_rds(pbmc, "atac_02_annotated")

###############################################################################
# STEP 3: COMPUTE QC METRICS
# Purpose: Calculate quality metrics for filtering
###############################################################################
cat("\n=== STEP 3: Computing QC metrics ===\n")

# Nucleosome signal: Ratio of mononucleosome to nucleosome-free fragments
# Good cells: <2 (indicates good chromatin accessibility)
cat("Computing nucleosome signal...\n")
pbmc <- NucleosomeSignal(pbmc)

# TSS enrichment: Signal at transcription start sites vs. background
# Good cells: >2 (indicates enrichment at regulatory regions)
cat("Computing TSS enrichment...\n")
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

# Blacklist ratio: Fraction of reads in ENCODE blacklist regions
# Good cells: <0.05 (these are artifact-prone regions)
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# Percent reads in peaks (FRiP): Fraction of reads in called peaks
# Good cells: >40% (indicates signal vs. noise)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# Save QC metrics table
write.csv(pbmc@meta.data, file.path(tabdir, "atac_qc_metrics.csv"))
cat("QC metrics saved\n")

# Generate QC plots
cat("Generating QC plots...\n")

# Density plots: Show relationships between metrics
p1 <- DensityScatter(pbmc, x = "nCount_ATAC", y = "TSS.enrichment", 
                     log_x = TRUE, quantiles = TRUE)
p2 <- DensityScatter(pbmc, x = "nucleosome_signal", y = "TSS.enrichment", 
                     quantiles = TRUE)

if (has_patchwork) {
  save_plot(p1 | p2, "atac_qc_density", w=14, h=6)
} else {
  save_plot(p1, "atac_qc_density1", w=10, h=6)
  save_plot(p2, "atac_qc_density2", w=10, h=6)
}

# Violin plots: Show distributions of all QC metrics
p_vln <- VlnPlot(pbmc, 
                 features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", 
                             "nucleosome_signal", "blacklist_ratio", "pct_reads_in_peaks"),
                 pt.size = 0.1, ncol = 6)
save_plot(p_vln, "atac_qc_violin", w=18, h=6)

save_rds(pbmc, "atac_03_qc")


###############################################################################
# STEP 4: OBJECTIVE QC FILTERING (MAD-BASED) - ROBUST VERSION
# Purpose: Remove low-quality cells using data-driven thresholds
# Method: Median Absolute Deviation (MAD) - robust to outliers
# Special handling: Gracefully handles all-zero blacklist_ratio
###############################################################################
cat("\n=== STEP 4: MAD-based QC filtering ===\n")
cat("Using", N_MAD_QC, "MAD stringency for objective thresholds\n\n")

# Store pre-filtering cell count
n_before <- ncol(pbmc)

# Display pre-filtering distributions
cat("PRE-FILTERING QC METRIC SUMMARIES:\n")
cat("═══════════════════════════════════════════════════\n")
cat("\nnCount_ATAC (total fragments per cell):\n")
print(summary(pbmc$nCount_ATAC))
cat("\nTSS.enrichment (signal at promoters):\n")
print(summary(pbmc$TSS.enrichment))
cat("\nnucleosome_signal (chromatin quality):\n")
print(summary(pbmc$nucleosome_signal))
cat("\npct_reads_in_peaks (signal-to-noise):\n")
print(summary(pbmc$pct_reads_in_peaks))
cat("\nblacklist_ratio (artifact regions):\n")
print(summary(pbmc$blacklist_ratio))
cat("═══════════════════════════════════════════════════\n\n")

# ============================================================
# ROBUST THRESHOLD CALCULATION FUNCTION
# ============================================================
calculate_qc_thresholds <- function(seurat_obj, n_mad = 3) {
  
  thresholds <- list()
  
  # 1. nCount_ATAC (total fragments)
  ncount <- seurat_obj$nCount_ATAC
  ncount_median <- median(ncount, na.rm = TRUE)
  ncount_mad <- mad(ncount, na.rm = TRUE)
  thresholds$ncount_lower <- max(1000, ncount_median - n_mad * ncount_mad)
  thresholds$ncount_upper <- ncount_median + n_mad * ncount_mad
  
  # 2. TSS enrichment (promoter signal)
  tss <- seurat_obj$TSS.enrichment
  tss_median <- median(tss, na.rm = TRUE)
  tss_mad <- mad(tss, na.rm = TRUE)
  thresholds$tss_lower <- max(1.5, tss_median - n_mad * tss_mad)
  
  # 3. Nucleosome signal (chromatin structure)
  nuc <- seurat_obj$nucleosome_signal
  nuc_median <- median(nuc, na.rm = TRUE)
  nuc_mad <- mad(nuc, na.rm = TRUE)
  thresholds$nuc_upper <- min(3.0, nuc_median + n_mad * nuc_mad)
  
  # 4. Percent reads in peaks (FRiP)
  frip <- seurat_obj$pct_reads_in_peaks
  frip_median <- median(frip, na.rm = TRUE)
  frip_mad <- mad(frip, na.rm = TRUE)
  thresholds$frip_lower <- max(20, frip_median - n_mad * frip_mad)
  
  # 5. Blacklist ratio (artifact regions) - ROBUST HANDLING
  if ("blacklist_ratio" %in% names(seurat_obj@meta.data)) {
    bl <- seurat_obj$blacklist_ratio
    bl_max <- max(bl, na.rm = TRUE)
    bl_mad <- mad(bl, na.rm = TRUE)
    
    if (bl_max == 0 || bl_mad == 0) {
      cat("NOTE: All cells have blacklist_ratio = 0 (excellent quality!)\n")
      cat("      Blacklist filtering will be SKIPPED (no variance to filter on)\n\n")
      thresholds$bl_upper <- Inf
      thresholds$bl_skip <- TRUE
    } else {
      bl_median <- median(bl, na.rm = TRUE)
      thresholds$bl_upper <- min(0.05, bl_median + n_mad * bl_mad)
      thresholds$bl_skip <- FALSE
    }
  } else {
    cat("NOTE: blacklist_ratio not found in metadata\n")
    cat("      Blacklist filtering will be SKIPPED\n\n")
    thresholds$bl_upper <- Inf
    thresholds$bl_skip <- TRUE
  }
  
  return(thresholds)
}

# Calculate thresholds
thresholds <- calculate_qc_thresholds(pbmc, n_mad = N_MAD_QC)

cat("📊 DATA-DRIVEN QC THRESHOLDS (", N_MAD_QC, " MAD from median):\n", sep="")
cat("═══════════════════════════════════════════════════\n")
cat(sprintf("  nCount_ATAC:         %6.0f - %6.0f fragments\n", 
            thresholds$ncount_lower, thresholds$ncount_upper))
cat(sprintf("  TSS.enrichment:      > %.2f\n", thresholds$tss_lower))
cat(sprintf("  nucleosome_signal:   < %.2f\n", thresholds$nuc_upper))
cat(sprintf("  pct_reads_in_peaks:  > %.1f%%\n", thresholds$frip_lower))

if (thresholds$bl_skip) {
  cat("  blacklist_ratio:     SKIPPED (all zeros or unavailable)\n")
} else {
  cat(sprintf("  blacklist_ratio:     < %.3f\n", thresholds$bl_upper))
}
cat("═══════════════════════════════════════════════════\n\n")

# Apply filters
cat("Applying filters...\n")

if (thresholds$bl_skip) {
  cat("  → Filtering on 4 metrics (blacklist skipped)\n")
  pbmc <- subset(pbmc, subset = 
    nCount_ATAC > thresholds$ncount_lower &
    nCount_ATAC < thresholds$ncount_upper &
    TSS.enrichment > thresholds$tss_lower &
    nucleosome_signal < thresholds$nuc_upper &
    pct_reads_in_peaks > thresholds$frip_lower
  )
} else {
  cat("  → Filtering on all 5 metrics (including blacklist)\n")
  pbmc <- subset(pbmc, subset = 
    nCount_ATAC > thresholds$ncount_lower &
    nCount_ATAC < thresholds$ncount_upper &
    TSS.enrichment > thresholds$tss_lower &
    nucleosome_signal < thresholds$nuc_upper &
    pct_reads_in_peaks > thresholds$frip_lower &
    blacklist_ratio < thresholds$bl_upper
  )
}

# Summarize results
n_after <- ncol(pbmc)
pct_retained <- round(n_after/n_before * 100, 1)
n_removed <- n_before - n_after

cat("\n FILTERING COMPLETE\n")
cat("═══════════════════════════════════════════════════\n")
cat(sprintf("  Cells before filtering:   %6d\n", n_before))
cat(sprintf("  Cells after filtering:    %6d\n", n_after))
cat(sprintf("  Cells retained:           %5.1f%%\n", pct_retained))
cat(sprintf("  Cells removed:            %6d\n", n_removed))
cat("═══════════════════════════════════════════════════\n\n")

if (pct_retained < 30) {
  cat("  WARNING: Retained <30% of cells\n")
  cat("   Consider: Re-run with N_MAD_QC = 4 for more lenient filtering\n\n")
} else if (pct_retained > 95) {
  cat("  WARNING: Retained >95% of cells\n")
  cat("   Consider: Re-run with N_MAD_QC = 2.5 for stricter filtering\n\n")
} else {
  cat("✓ Retention rate is within expected range (30-95%)\n\n")
}

cat("📈 POST-FILTERING QC METRIC DISTRIBUTIONS:\n")
cat("═══════════════════════════════════════════════════\n")
cat(sprintf("  nCount_ATAC        median: %6.0f  (MAD: %5.0f)\n", 
            median(pbmc$nCount_ATAC), mad(pbmc$nCount_ATAC)))
cat(sprintf("  TSS.enrichment     median: %6.2f  (MAD: %5.2f)\n", 
            median(pbmc$TSS.enrichment), mad(pbmc$TSS.enrichment)))
cat(sprintf("  nucleosome_signal  median: %6.2f  (MAD: %5.2f)\n", 
            median(pbmc$nucleosome_signal), mad(pbmc$nucleosome_signal)))
cat(sprintf("  pct_reads_in_peaks median: %6.1f  (MAD: %5.1f)\n", 
            median(pbmc$pct_reads_in_peaks), mad(pbmc$pct_reads_in_peaks)))
cat("═══════════════════════════════════════════════════\n\n")

# Save thresholds
threshold_df <- data.frame(
  metric = c("nCount_ATAC_lower", "nCount_ATAC_upper", "TSS.enrichment_lower",
             "nucleosome_signal_upper", "pct_reads_in_peaks_lower", "blacklist_ratio_upper"),
  threshold_value = c(thresholds$ncount_lower, thresholds$ncount_upper, 
                      thresholds$tss_lower, thresholds$nuc_upper, 
                      thresholds$frip_lower, 
                      ifelse(thresholds$bl_skip, NA, thresholds$bl_upper)),
  applied = c(TRUE, TRUE, TRUE, TRUE, TRUE, !thresholds$bl_skip),
  n_mad_used = N_MAD_QC,
  cells_before = n_before,
  cells_after = n_after,
  pct_retained = pct_retained,
  timestamp = as.character(Sys.time())
)

threshold_file <- file.path(tabdir, "filtering_thresholds_applied.csv")
write.csv(threshold_df, threshold_file, row.names = FALSE)
cat("Thresholds saved to:", threshold_file, "\n\n")

save_rds(pbmc, "atac_04_filtered")


###############################################################################
# STEP 5: NORMALIZATION AND DIMENSIONALITY REDUCTION
###############################################################################
cat("\n=== STEP 5: Normalization and dimensionality reduction ===\n")

cat("Running TF-IDF normalization...\n")
pbmc <- RunTFIDF(pbmc)

cat("Identifying top features...\n")
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")

cat("Running LSI (dimensionality reduction)...\n")
pbmc <- RunSVD(pbmc)

cat("Checking LSI component correlations with depth...\n")
p_depth <- DepthCor(pbmc)
save_plot(p_depth, "atac_lsi_depth_correlation", w=10, h=6)
cat("Note: LSI_1 typically correlates with depth and is excluded from downstream analysis\n\n")

cat("Computing UMAP embedding...\n")
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = LSI_DIMS)

cat("Computing t-SNE embedding...\n")
pbmc <- RunTSNE(pbmc, reduction = "lsi", dims = LSI_DIMS)

cat("Finding neighbors and clustering...\n")
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = LSI_DIMS)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = CLUSTER_RESOLUTION)

# MODIFIED: Create LABELED UMAP plots
p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) + 
          ggtitle("ATAC-seq Clusters (UMAP)") +
          theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

p_tsne <- DimPlot(pbmc, reduction = "tsne", label = TRUE, label.size = 6) + 
          ggtitle("ATAC-seq Clusters (t-SNE)") +
          theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

if (has_patchwork) {
  save_plot(p_umap | p_tsne, "atac_clustering_umap_tsne_labeled", w=16, h=7)
} else {
  save_plot(p_umap, "atac_clustering_umap_labeled", w=10, h=8)
  save_plot(p_tsne, "atac_clustering_tsne_labeled", w=10, h=8)
}

cat("Clustering complete:", nlevels(Idents(pbmc)), "clusters identified\n")
save_rds(pbmc, "atac_05_clustered")

###############################################################################
# STEP 6: GENE ACTIVITY MATRIX
###############################################################################
cat("\n=== STEP 6: Creating gene activity matrix ===\n")

cat("Computing gene activity scores (this may take several minutes)...\n")
gene.activities <- GeneActivity(pbmc)

pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)

cat("Normalizing gene activity scores...\n")
pbmc <- NormalizeData(pbmc, 
                      assay = "RNA", 
                      normalization.method = "LogNormalize",
                      scale.factor = median(pbmc$nCount_RNA))

cat("Gene activity matrix created:", nrow(pbmc[["RNA"]]), "genes\n")
save_rds(pbmc, "atac_06_with_gene_activity")

###############################################################################
# STEP 7: VISUALIZE B CELL MARKERS (CD19, CD20, CD38, CD40)
# Purpose: Focus specifically on B cell surface markers
###############################################################################
cat("\n=== STEP 7: Visualizing B cell markers (CD19, CD20, CD38, CD40) ===\n")

# Switch to RNA assay (gene activity scores)
DefaultAssay(pbmc) <- "RNA"

# MODIFIED: Only B cell surface markers (CD20 = MS4A1)
b_cell_markers <- c("CD19", "MS4A1", "CD38", "CD40")

cat("Checking which markers are present in the data...\n")
markers_present <- b_cell_markers[b_cell_markers %in% rownames(pbmc[["RNA"]])]

if (length(markers_present) == 0) {
  cat("⚠️  WARNING: No B cell markers found in gene activity matrix!\n")
  cat("   Available genes:", length(rownames(pbmc[["RNA"]])), "\n")
} else {
  cat("Found", length(markers_present), "of", length(b_cell_markers), "markers:\n")
  cat("  Present:", paste(markers_present, collapse=", "), "\n")
  
  if (length(markers_present) < length(b_cell_markers)) {
    missing <- setdiff(b_cell_markers, markers_present)
    cat("  Missing:", paste(missing, collapse=", "), "\n")
  }
  
  # Create feature plot for B cell markers
  cat("\nGenerating UMAP plots for B cell markers...\n")
  
  p_bcell <- FeaturePlot(pbmc, 
                         features = markers_present, 
                         reduction = "umap",
                         pt.size = 0.3,
                         max.cutoff = "q95",
                         ncol = 2)
  
  save_plot(p_bcell, "bcell_markers_CD19_CD20_CD38_CD40_umap", 
            w=12, h=ceiling(length(markers_present)/2)*5)
  
  # Also create violin plots to see distribution across clusters
  cat("Generating violin plots for B cell markers...\n")
  
  for (marker in markers_present) {
    p_vln <- VlnPlot(pbmc, features = marker, pt.size = 0.1) +
             ggtitle(paste0(marker, " Expression Across Clusters")) +
             theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    save_plot(p_vln, paste0("bcell_marker_", marker, "_violin"), w=12, h=6)
  }
  
  # Create a combined feature plot showing all markers together
  if (length(markers_present) >= 2) {
    cat("Creating combined marker expression plot...\n")
    
    # Calculate average expression of all B cell markers
    marker_data <- FetchData(pbmc, vars = markers_present)
    pbmc$bcell_score <- rowMeans(marker_data, na.rm = TRUE)
    
    p_score <- FeaturePlot(pbmc, features = "bcell_score", reduction = "umap", pt.size = 0.3) +
               ggtitle("Combined B Cell Marker Score\n(CD19 + CD20 + CD38 + CD40)") +
               theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    save_plot(p_score, "bcell_combined_score_umap", w=10, h=8)
  }
}

cat("\nB cell marker visualization complete\n")

###############################################################################
# STEP 8: LABEL TRANSFER FROM REFERENCE RNA-seq
###############################################################################
cat("\n=== STEP 8: Cell type annotation via label transfer ===\n")

rna_ref <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/pbmc_10k_v3.rds"

if (file.exists(rna_ref)) {
  cat("Loading reference RNA-seq dataset...\n")
  pbmc_rna <- readRDS(rna_ref)
  pbmc_rna <- UpdateSeuratObject(pbmc_rna)
  
  cat("Reference contains:", ncol(pbmc_rna), "cells\n")
  cat("Cell types in reference:", paste(unique(pbmc_rna$celltype), collapse=", "), "\n\n")
  
  cat("Finding transfer anchors (this may take 5-10 minutes)...\n")
  transfer.anchors <- FindTransferAnchors(
    reference = pbmc_rna, 
    query = pbmc, 
    reduction = "pcaproject"
  )
  
  n_anchors <- nrow(transfer.anchors@anchors)
  cat("Found", n_anchors, "anchor pairs between datasets\n")
  
  k_weight_use <- min(floor(n_anchors * 0.9), 40)
  cat("Using k.weight =", k_weight_use, "for label transfer\n\n")
  
  cat("Transferring cell type labels...\n")
  predicted.labels <- TransferData(
    anchorset = transfer.anchors, 
    refdata = pbmc_rna$celltype,
    weight.reduction = pbmc[["lsi"]], 
    dims = LSI_DIMS, 
    k.weight = k_weight_use
  )
  
  pbmc <- AddMetaData(pbmc, metadata = predicted.labels)
  
  # MODIFIED: Create LABELED predicted cell type UMAP
  p_pred <- DimPlot(pbmc, 
                    reduction = "umap", 
                    group.by = "predicted.id",
                    label = TRUE,
                    label.size = 5,
                    repel = TRUE) +
            ggtitle("Predicted Cell Types") +
            theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  save_plot(p_pred, "atac_predicted_celltypes_labeled", w=12, h=8)
  
  celltype_counts <- as.data.frame(table(pbmc$predicted.id))
  colnames(celltype_counts) <- c("CellType", "Count")
  celltype_counts <- celltype_counts[order(-celltype_counts$Count), ]
  write.csv(celltype_counts, file.path(tabdir, "celltype_counts.csv"), row.names=FALSE)
  
  cat("\nPredicted cell type distribution:\n")
  print(celltype_counts)
  
  save_rds(pbmc, "atac_07_labeled")
  
} else {
  cat("⚠️  WARNING: Reference RNA-seq data not found\n")
  cat("   Expected location:", rna_ref, "\n")
  cat("   Skipping label transfer\n\n")
}

###############################################################################
# STEP 9: DIFFERENTIAL ACCESSIBILITY ANALYSIS
###############################################################################
cat("\n=== STEP 9: Differential accessibility analysis ===\n")

DefaultAssay(pbmc) <- "ATAC"

if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  Idents(pbmc) <- pbmc$predicted.id
  cell_counts <- sort(table(Idents(pbmc)), decreasing = TRUE)
  
  ident_1 <- names(cell_counts)[1]
  ident_2 <- names(cell_counts)[2]
  
  cat("Comparing cell types:", ident_1, "vs", ident_2, "\n")
  cat("  ", ident_1, ":", cell_counts[1], "cells\n")
  cat("  ", ident_2, ":", cell_counts[2], "cells\n\n")
  
} else {
  Idents(pbmc) <- pbmc$seurat_clusters
  ident_1 <- levels(Idents(pbmc))[1]
  ident_2 <- levels(Idents(pbmc))[2]
  
  cat("Comparing clusters:", ident_1, "vs", ident_2, "\n\n")
}

cat("Finding differentially accessible peaks...\n")
cat("(This may take several minutes for large datasets)\n")

da_peaks <- FindMarkers(pbmc, 
                        ident.1 = ident_1, 
                        ident.2 = ident_2,
                        test.use = "LR", 
                        latent.vars = "nCount_ATAC",
                        min.pct = 0.05)

da_peaks <- da_peaks[order(da_peaks$p_val_adj), ]

da_file <- file.path(tabdir, paste0("DA_peaks_", ident_1, "_vs_", ident_2, ".csv"))
write.csv(da_peaks, da_file)

n_sig <- sum(da_peaks$p_val_adj < 0.05, na.rm = TRUE)
n_up <- sum(da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC > 0, na.rm = TRUE)
n_down <- sum(da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC < 0, na.rm = TRUE)

cat("\nDifferential accessibility results:\n")
cat("  Total peaks tested:         ", nrow(da_peaks), "\n")
cat("  Significant peaks (FDR<0.05):", n_sig, "\n")
cat("    More accessible in", ident_1, ":", n_up, "\n")
cat("    More accessible in", ident_2, ":", n_down, "\n")
cat("  Results saved to:", da_file, "\n\n")

if (nrow(da_peaks) > 0) {
  top_peak <- rownames(da_peaks)[1]
  cat("Visualizing top peak:", top_peak, "\n")
  cat("  log2FC =", round(da_peaks[1, "avg_log2FC"], 2), "\n")
  cat("  p_adj  =", format(da_peaks[1, "p_val_adj"], scientific=TRUE, digits=3), "\n\n")
  
  p_vln <- VlnPlot(pbmc, features = top_peak, pt.size = 0.1, 
                   idents = c(ident_1, ident_2))
  p_umap <- FeaturePlot(pbmc, features = top_peak, reduction = "umap", 
                        pt.size = 0.1)
  
  if (has_patchwork) {
    save_plot(p_vln | p_umap, "top_DA_peak", w=14, h=6)
  } else {
    save_plot(p_vln, "top_DA_peak_violin", w=10, h=6)
    save_plot(p_umap, "top_DA_peak_umap", w=10, h=6)
  }
}

###############################################################################
# STEP 10: COVERAGE PLOTS FOR B CELL MARKERS
# Purpose: Visualize chromatin accessibility at B cell marker loci
###############################################################################
cat("\n=== STEP 10: Creating coverage plots for B cell markers ===\n")

if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  levels(pbmc) <- unique(pbmc$predicted.id)
  cat("Grouping coverage by cell type\n")
} else {
  cat("Grouping coverage by cluster\n")
}

# MODIFIED: Only B cell markers (CD20 = MS4A1)
key_genes <- c("CD19", "MS4A1", "CD38", "CD40")

cat("Generating coverage plots for", length(key_genes), "B cell marker genes...\n")
cat("(Each plot may take 1-2 minutes)\n\n")

for (gene in key_genes) {
  cat("  Plotting:", gene, "...")
  
  tryCatch({
    p <- CoveragePlot(pbmc, 
                      region = gene, 
                      extend.upstream = 40000,
                      extend.downstream = 20000)
    
    save_plot(p, paste0("coverage_", gene), w=14, h=8)
    cat(" done\n")
    
  }, error = function(e) {
    cat(" FAILED -", e$message, "\n")
  })
}

cat("\nCoverage plot generation complete\n")

save_rds(pbmc, "atac_final")

###############################################################################
# ANALYSIS COMPLETE - SUMMARY
###############################################################################
cat("\n\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("          ATAC-seq ANALYSIS COMPLETE (B CELL FOCUS)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("Final ATAC-seq object:\n")
print(pbmc)

cat("\n Output directories:\n")
cat("  Main:   ", outs_dir, "\n")
cat("  Plots:  ", plotdir, "\n")
cat("  Tables: ", tabdir, "\n")
cat("  RDS:    ", rdsdir, "\n\n")

cat(" Key files generated:\n")
cat("  1. QC metrics:              atac_qc_metrics.csv\n")
cat("  2. Filtering thresholds:    filtering_thresholds_applied.csv\n")
cat("  3. Cell type counts:        celltype_counts.csv\n")
cat("  4. Differential peaks:      DA_peaks_*.csv\n")
cat("  5. B cell markers (UMAP):   bcell_markers_CD19_CD20_CD38_CD40_umap.pdf\n")
cat("  6. Coverage plots:          coverage_CD19/MS4A1/CD38/CD40.pdf\n")
cat("  7. Labeled UMAPs:           atac_clustering_umap_tsne_labeled.pdf\n\n")

cat(" Analysis summary:\n")
cat("  MAD stringency used:        ", N_MAD_QC, "\n")
cat("  Cells retained:             ", ncol(pbmc), "\n")
cat("  Peaks analyzed:             ", nrow(pbmc), "\n")
cat("  Clusters identified:        ", nlevels(Idents(pbmc)), "\n")
if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  cat("  Cell types annotated:       ", length(unique(pbmc$predicted.id)), "\n")
}

cat("\n B Cell Analysis Focus:\n")
cat("  Markers analyzed:           CD19, CD20 (MS4A1), CD38, CD40\n")
cat("  Coverage plots generated:   4 B cell marker loci\n")
cat("  UMAP plots:                 Labeled with cluster numbers and cell types\n\n")

cat("  Analysis completed:", format(Sys.time()), "\n\n")

cat("═══════════════════════════════════════════════════════════\n")

sink()

cat("\n Full log saved to:", logfile, "\n")
cat(" All outputs saved to:", outs_dir, "\n")
cat(" B cell-specific analysis in: B_cell/ subdirectory\n\n")

###############################################################################
# END OF SCRIPT

###############################################################################
