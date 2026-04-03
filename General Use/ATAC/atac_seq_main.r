###############################################################################
# COMPLETE ATAC-seq Pipeline - OBJECTIVE & REPRODUCIBLE
# MAD-based QC filtering + Standard ATAC-seq workflow
# a qc based on MAD - median and +_ 3 SD -- this automates qc for outiers cells 
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

## -------------------- ANALYSIS PARAMETERS (MODIFY HERE) ---------------------
# QC filtering stringency
N_MAD_QC <- 3  # Options: 2.5 (strict), 3 (standard/recommended), 4 (lenient)

# Clustering resolution
CLUSTER_RESOLUTION <- 0.5  # Options: 0.3 (coarse), 0.5 (standard), 0.8 (fine)

# Dimensionality reduction
LSI_DIMS <- 2:30  # Standard for ATAC-seq (always skip first LSI component)

## -------------------- SETUP (dirs + logging) --------------------------------
outs_dir <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/ATAC_150"
dir.create(outs_dir, recursive = TRUE, showWarnings = FALSE)

plotdir <- file.path(outs_dir, "plots")
rdsdir  <- file.path(outs_dir, "rds")
tabdir  <- file.path(outs_dir, "tables")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir,  recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(outs_dir, "run.log")
sink(logfile, split = TRUE)
cat("===== ATAC-seq ANALYSIS - OBJECTIVE QC PIPELINE =====\n")
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
cat("---------------------------------------------------\n")
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
cat("---------------------------------------------------\n\n")

# ============================================================
# ROBUST THRESHOLD CALCULATION FUNCTION
# ============================================================
# Uses MAD (Median Absolute Deviation) - more robust than SD
# MAD = median(|X - median(X)|) * 1.4826 (scale factor for normal distribution)
# 
# Why MAD instead of mean ± SD?
# - Robust to outliers (doesn't get inflated by extreme values)
# - Adapts to your specific dataset
# - Standard in single-cell genomics (used by Seurat, Scanpy)
#
# Special feature: Handles edge cases (all-zero blacklist_ratio)
# ============================================================

calculate_qc_thresholds <- function(seurat_obj, n_mad = 3) {
  
  thresholds <- list()
  
  # ----------------------------------------
  # 1. nCount_ATAC (total fragments)
  # ----------------------------------------
  # Strategy: Remove cells with too few or too many fragments
  # Too few = low-quality cells
  # Too many = doublets (two cells in one droplet)
  
  ncount <- seurat_obj$nCount_ATAC
  ncount_median <- median(ncount, na.rm = TRUE)
  ncount_mad <- mad(ncount, na.rm = TRUE)
  
  # Lower threshold: median - n_mad * MAD (but never below 1000)
  thresholds$ncount_lower <- max(1000, ncount_median - n_mad * ncount_mad)
  
  # Upper threshold: median + n_mad * MAD
  thresholds$ncount_upper <- ncount_median + n_mad * ncount_mad
  
  # ----------------------------------------
  # 2. TSS enrichment (promoter signal)
  # ----------------------------------------
  # Strategy: Remove cells with low TSS enrichment (poor chromatin accessibility)
  # Higher = better (only filter lower bound)
  
  tss <- seurat_obj$TSS.enrichment
  tss_median <- median(tss, na.rm = TRUE)
  tss_mad <- mad(tss, na.rm = TRUE)
  
  # Biological minimum for ATAC-seq is ~1.5
  thresholds$tss_lower <- max(1.5, tss_median - n_mad * tss_mad)
  
  # ----------------------------------------
  # 3. Nucleosome signal (chromatin structure)
  # ----------------------------------------
  # Strategy: Remove cells with high nucleosome signal (poor fragmentation)
  # Lower = better (only filter upper bound)
  
  nuc <- seurat_obj$nucleosome_signal
  nuc_median <- median(nuc, na.rm = TRUE)
  nuc_mad <- mad(nuc, na.rm = TRUE)
  
  # Cap at biological maximum of 3.0
  thresholds$nuc_upper <- min(3.0, nuc_median + n_mad * nuc_mad)
  
  # ----------------------------------------
  # 4. Percent reads in peaks (FRiP)
  # ----------------------------------------
  # Strategy: Remove cells with low signal-to-noise ratio
  # Higher = better (only filter lower bound)
  
  frip <- seurat_obj$pct_reads_in_peaks
  frip_median <- median(frip, na.rm = TRUE)
  frip_mad <- mad(frip, na.rm = TRUE)
  
  # Never go below 20% (biological minimum)
  thresholds$frip_lower <- max(20, frip_median - n_mad * frip_mad)
  
  # ----------------------------------------
  # 5. Blacklist ratio (artifact regions) - ROBUST HANDLING
  # ----------------------------------------
  # Strategy: Remove cells with high reads in problematic regions
  # Lower = better (only filter upper bound)
  #
  # Edge case handling: If all cells have blacklist_ratio = 0:
  # - This indicates EXCELLENT quality (no artifact contamination)
  # - MAD would be 0, making threshold calculation unstable
  # - Solution: Skip blacklist filtering (can't distinguish cells anyway)
  
  if ("blacklist_ratio" %in% names(seurat_obj@meta.data)) {
    bl <- seurat_obj$blacklist_ratio
    bl_max <- max(bl, na.rm = TRUE)
    bl_mad <- mad(bl, na.rm = TRUE)
    
    # Check if all values are zero or MAD is zero
    if (bl_max == 0 || bl_mad == 0) {
      cat("NOTE: All cells have blacklist_ratio = 0 (excellent quality!)\n")
      cat("      Blacklist filtering will be SKIPPED (no variance to filter on)\n\n")
      thresholds$bl_upper <- Inf  # Inf means don't filter any cells
      thresholds$bl_skip <- TRUE
      
    } else {
      # Normal case: some variation exists, use MAD-based threshold
      bl_median <- median(bl, na.rm = TRUE)
      thresholds$bl_upper <- min(0.05, bl_median + n_mad * bl_mad)
      thresholds$bl_skip <- FALSE
    }
  } else {
    # Metric not available in metadata
    cat("NOTE: blacklist_ratio not found in metadata\n")
    cat("      Blacklist filtering will be SKIPPED\n\n")
    thresholds$bl_upper <- Inf
    thresholds$bl_skip <- TRUE
  }
  
  return(thresholds)
}

# ============================================================
# CALCULATE THRESHOLDS
# ============================================================
thresholds <- calculate_qc_thresholds(pbmc, n_mad = N_MAD_QC)

cat("?? DATA-DRIVEN QC THRESHOLDS (", N_MAD_QC, " MAD from median):\n", sep="")
cat("---------------------------------------------------\n")
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
cat("---------------------------------------------------\n\n")

# ============================================================
# APPLY FILTERS TO SEURAT OBJECT - CONDITIONAL BLACKLIST
# ============================================================
cat("Applying filters...\n")

# Build filter expression conditionally based on blacklist availability
if (thresholds$bl_skip) {
  # Skip blacklist filtering (all zeros, no variance, or metric unavailable)
  cat("  ? Filtering on 4 metrics (blacklist skipped)\n")
  pbmc <- subset(pbmc, subset = 
    nCount_ATAC > thresholds$ncount_lower &
    nCount_ATAC < thresholds$ncount_upper &
    TSS.enrichment > thresholds$tss_lower &
    nucleosome_signal < thresholds$nuc_upper &
    pct_reads_in_peaks > thresholds$frip_lower
  )
} else {
  # Include blacklist filtering
  cat("  ? Filtering on all 5 metrics (including blacklist)\n")
  pbmc <- subset(pbmc, subset = 
    nCount_ATAC > thresholds$ncount_lower &
    nCount_ATAC < thresholds$ncount_upper &
    TSS.enrichment > thresholds$tss_lower &
    nucleosome_signal < thresholds$nuc_upper &
    pct_reads_in_peaks > thresholds$frip_lower &
    blacklist_ratio < thresholds$bl_upper
  )
}

# ============================================================
# SUMMARIZE FILTERING RESULTS
# ============================================================
n_after <- ncol(pbmc)
pct_retained <- round(n_after/n_before * 100, 1)
n_removed <- n_before - n_after

cat("\n? FILTERING COMPLETE\n")
cat("---------------------------------------------------\n")
cat(sprintf("  Cells before filtering:   %6d\n", n_before))
cat(sprintf("  Cells after filtering:    %6d\n", n_after))
cat(sprintf("  Cells retained:           %5.1f%%\n", pct_retained))
cat(sprintf("  Cells removed:            %6d\n", n_removed))
cat("---------------------------------------------------\n\n")

# ============================================================
# QUALITY CONTROL CHECKS
# ============================================================
# Typical ATAC-seq datasets retain 50-80% of cells
# <30% = thresholds too strict (consider n_mad = 4)
# >95% = thresholds too lenient (consider n_mad = 2.5)

if (pct_retained < 30) {
  cat("??  WARNING: Retained <30% of cells\n")
  cat("   Your data may be lower quality, or thresholds are too strict\n")
  cat("   Consider: Re-run with N_MAD_QC = 4 for more lenient filtering\n\n")
} else if (pct_retained > 95) {
  cat("??  WARNING: Retained >95% of cells\n")
  cat("   Thresholds may be too permissive\n")
  cat("   Consider: Re-run with N_MAD_QC = 2.5 for stricter filtering\n\n")
} else {
  cat("? Retention rate is within expected range (30-95%)\n")
  cat("  This suggests appropriate filtering stringency\n\n")
}

# ============================================================
# POST-FILTERING QC DISTRIBUTIONS
# ============================================================
cat("?? POST-FILTERING QC METRIC DISTRIBUTIONS:\n")
cat("---------------------------------------------------\n")
cat(sprintf("  nCount_ATAC        median: %6.0f  (MAD: %5.0f)\n", 
            median(pbmc$nCount_ATAC), mad(pbmc$nCount_ATAC)))
cat(sprintf("  TSS.enrichment     median: %6.2f  (MAD: %5.2f)\n", 
            median(pbmc$TSS.enrichment), mad(pbmc$TSS.enrichment)))
cat(sprintf("  nucleosome_signal  median: %6.2f  (MAD: %5.2f)\n", 
            median(pbmc$nucleosome_signal), mad(pbmc$nucleosome_signal)))
cat(sprintf("  pct_reads_in_peaks median: %6.1f  (MAD: %5.1f)\n", 
            median(pbmc$pct_reads_in_peaks), mad(pbmc$pct_reads_in_peaks)))
cat("---------------------------------------------------\n\n")

# ============================================================
# SAVE THRESHOLDS FOR REPRODUCIBILITY
# ============================================================
# This CSV documents exactly what thresholds were used
# Critical for methods section and reproducing analysis
# Includes an 'applied' column to track which filters were actually used

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
cat("Thresholds saved to:", threshold_file, "\n")
cat("Include this file in supplementary materials for full reproducibility\n\n")

save_rds(pbmc, "atac_04_filtered")


###############################################################################
# STEP 5: NORMALIZATION AND DIMENSIONALITY REDUCTION
# Purpose: TF-IDF normalization + Latent Semantic Indexing (LSI)
# This is the standard workflow for ATAC-seq (analogous to PCA for RNA-seq)
###############################################################################
cat("\n=== STEP 5: Normalization and dimensionality reduction ===\n")

# TF-IDF: Term Frequency-Inverse Document Frequency
# Normalizes for differences in sequencing depth and peak accessibility
cat("Running TF-IDF normalization...\n")
pbmc <- RunTFIDF(pbmc)

# Find top features (most variable peaks across cells)
# min.cutoff = "q0" retains all features (recommended for ATAC)
cat("Identifying top features...\n")
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")

# LSI: Latent Semantic Indexing (SVD on TF-IDF matrix)
# This is analogous to PCA for RNA-seq
# Note: First LSI component often correlates with sequencing depth (technical)
cat("Running LSI (dimensionality reduction)...\n")
pbmc <- RunSVD(pbmc)

# Check correlation between LSI components and sequencing depth
# Component 1 usually correlates with depth and should be excluded
cat("Checking LSI component correlations with depth...\n")
p_depth <- DepthCor(pbmc)
save_plot(p_depth, "atac_lsi_depth_correlation", w=10, h=6)
cat("Note: LSI_1 typically correlates with depth and is excluded from downstream analysis\n\n")

# UMAP: Non-linear dimensionality reduction for visualization
# dims = 2:30 excludes first component (correlated with depth)
cat("Computing UMAP embedding...\n")
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = LSI_DIMS)

# t-SNE: Alternative visualization (slower but sometimes clearer separation)
cat("Computing t-SNE embedding...\n")
pbmc <- RunTSNE(pbmc, reduction = "lsi", dims = LSI_DIMS)

# Clustering: Find neighbors and cluster cells
cat("Finding neighbors and clustering...\n")
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = LSI_DIMS)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = CLUSTER_RESOLUTION)

# Visualize clustering results
p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
p_tsne <- DimPlot(pbmc, reduction = "tsne", label = TRUE) + NoLegend()

if (has_patchwork) {
  save_plot(p_umap | p_tsne, "atac_clustering_umap_tsne", w=14, h=6)
} else {
  save_plot(p_umap, "atac_clustering_umap", w=10, h=6)
  save_plot(p_tsne, "atac_clustering_tsne", w=10, h=6)
}

cat("Clustering complete:", nlevels(Idents(pbmc)), "clusters identified\n")
save_rds(pbmc, "atac_05_clustered")

###############################################################################
# STEP 6: GENE ACTIVITY MATRIX
# Purpose: Infer gene expression from chromatin accessibility
# Method: Sum accessibility of peaks within gene body + promoter
###############################################################################
cat("\n=== STEP 6: Creating gene activity matrix ===\n")

# Gene Activity: Estimate gene expression from ATAC-seq data
# This allows RNA-based marker genes to be used with ATAC data
# Note: This is an *approximation* - not as accurate as real RNA-seq
cat("Computing gene activity scores (this may take several minutes)...\n")
gene.activities <- GeneActivity(pbmc)

# Add as RNA assay to Seurat object
pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)

# Normalize gene activity scores
# Use median normalization (robust to outliers)
cat("Normalizing gene activity scores...\n")
pbmc <- NormalizeData(pbmc, 
                      assay = "RNA", 
                      normalization.method = "LogNormalize",
                      scale.factor = median(pbmc$nCount_RNA))

cat("Gene activity matrix created:", nrow(pbmc[["RNA"]]), "genes\n")
save_rds(pbmc, "atac_06_with_gene_activity")

###############################################################################
# STEP 7: VISUALIZE MARKER GENES
# Purpose: Identify cell types using known marker genes
###############################################################################
cat("\n=== STEP 7: Visualizing marker genes ===\n")

# Switch to RNA assay (gene activity scores)
DefaultAssay(pbmc) <- "RNA"

# Define marker genes for major cell types and your research focus
markers <- list(
  # B cell markers (your focus: nephritis, glycosylation)
  b_cell = c("MS4A1", "CD19", "CD79A", "PAX5", "CD27"),
  
  # Glycosylation enzymes (SLE-relevant)
  glycosylation = c("MGAT5", "ST6GAL1", "FUT8", "B4GALT1", "MGAT1", "MGAT2"),
  
  # Metabolic genes (B cell activation)
  metabolism = c("IMPDH2", "GLS", "FASN", "SLC2A1", "HK2", "LDHA"),
  
  # SLE/autoimmune markers
  sle = c("IFI44L", "PRDM1", "CD38", "C3", "IRF7", "MX1"),
  
  # T cell markers (context)
  t_cell = c("CD3D", "CD3E", "CD4", "CD8A"),
  
  # Myeloid markers (context)
  myeloid = c("CD14", "FCGR3A", "CD68", "LYZ")
)

# Plot each marker set
for (marker_name in names(markers)) {
  cat("Plotting:", marker_name, "markers...\n")
  
  genes <- markers[[marker_name]]
  genes_present <- genes[genes %in% rownames(pbmc[["RNA"]])]
  
  if (length(genes_present) == 0) {
    cat("  WARNING: No markers found for", marker_name, "\n")
    next
  }
  
  if (length(genes_present) < length(genes)) {
    missing <- setdiff(genes, genes_present)
    cat("  Note: Missing genes:", paste(missing, collapse=", "), "\n")
  }
  
  # Create feature plot (max 2 genes creates better layout)
  if (length(genes_present) >= 2) {
    p <- FeaturePlot(pbmc, 
                     features = genes_present, 
                     reduction = "umap",
                     pt.size = 0.1, 
                     max.cutoff = "q95",  # Cap at 95th percentile for better contrast
                     ncol = 3)
    
    save_plot(p, 
              paste0("markers_", marker_name, "_umap"), 
              w=12, 
              h=ceiling(length(genes_present)/3)*4)
  }
}

cat("\nMarker visualization complete\n")

###############################################################################
# STEP 8: LABEL TRANSFER FROM REFERENCE RNA-seq
# Purpose: Annotate cell types using a reference dataset
# This transfers labels from a well-annotated RNA-seq dataset
###############################################################################
cat("\n=== STEP 8: Cell type annotation via label transfer ===\n")

rna_ref <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/pbmc_10k_v3.rds"

if (file.exists(rna_ref)) {
  cat("Loading reference RNA-seq dataset...\n")
  pbmc_rna <- readRDS(rna_ref)
  pbmc_rna <- UpdateSeuratObject(pbmc_rna)
  
  cat("Reference contains:", ncol(pbmc_rna), "cells\n")
  cat("Cell types in reference:", paste(unique(pbmc_rna$celltype), collapse=", "), "\n\n")
  
  # Find anchors between reference (RNA) and query (ATAC)
  # This identifies corresponding cells across modalities
  cat("Finding transfer anchors (this may take 5-10 minutes)...\n")
  transfer.anchors <- FindTransferAnchors(
    reference = pbmc_rna, 
    query = pbmc, 
    reduction = "pcaproject"
  )
  
  n_anchors <- nrow(transfer.anchors@anchors)
  cat("Found", n_anchors, "anchor pairs between datasets\n")
  
  # Dynamic k.weight based on number of anchors
  # k.weight = number of neighbors used for label transfer
  # Use 90% of anchors, capped at 40
  k_weight_use <- min(floor(n_anchors * 0.9), 40)
  cat("Using k.weight =", k_weight_use, "for label transfer\n\n")
  
  # Transfer cell type labels from reference to query
  cat("Transferring cell type labels...\n")
  predicted.labels <- TransferData(
    anchorset = transfer.anchors, 
    refdata = pbmc_rna$celltype,
    weight.reduction = pbmc[["lsi"]], 
    dims = LSI_DIMS, 
    k.weight = k_weight_use
  )
  
  # Add predictions to ATAC object
  pbmc <- AddMetaData(pbmc, metadata = predicted.labels)
  
  # Visualize predicted cell types
  p_pred <- DimPlot(pbmc, 
                    reduction = "umap", 
                    group.by = "predicted.id",
                    label = TRUE, 
                    repel = TRUE) + 
            NoLegend()
  save_plot(p_pred, "atac_predicted_celltypes", w=10, h=7)
  
  # Save cell type counts
  celltype_counts <- as.data.frame(table(pbmc$predicted.id))
  colnames(celltype_counts) <- c("CellType", "Count")
  celltype_counts <- celltype_counts[order(-celltype_counts$Count), ]
  write.csv(celltype_counts, file.path(tabdir, "celltype_counts.csv"), row.names=FALSE)
  
  cat("\nPredicted cell type distribution:\n")
  print(celltype_counts)
  
  save_rds(pbmc, "atac_07_labeled")
  
} else {
  cat("??  WARNING: Reference RNA-seq data not found\n")
  cat("   Expected location:", rna_ref, "\n")
  cat("   Skipping label transfer - you can manually annotate using marker genes\n\n")
}

###############################################################################
# STEP 9: DIFFERENTIAL ACCESSIBILITY ANALYSIS
# Purpose: Find peaks that differ between cell types or conditions
###############################################################################
cat("\n=== STEP 9: Differential accessibility analysis ===\n")

# Switch to ATAC assay for peak analysis
DefaultAssay(pbmc) <- "ATAC"

# Determine which groups to compare
if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  # If we have cell type labels, compare cell types
  Idents(pbmc) <- pbmc$predicted.id
  cell_counts <- sort(table(Idents(pbmc)), decreasing = TRUE)
  
  # Compare top 2 most abundant cell types
  ident_1 <- names(cell_counts)[1]
  ident_2 <- names(cell_counts)[2]
  
  cat("Comparing cell types:", ident_1, "vs", ident_2, "\n")
  cat("  ", ident_1, ":", cell_counts[1], "cells\n")
  cat("  ", ident_2, ":", cell_counts[2], "cells\n\n")
  
} else {
  # If no cell type labels, compare clusters
  Idents(pbmc) <- pbmc$seurat_clusters
  ident_1 <- levels(Idents(pbmc))[1]
  ident_2 <- levels(Idents(pbmc))[2]
  
  cat("Comparing clusters:", ident_1, "vs", ident_2, "\n\n")
}

# Find differentially accessible peaks
# test.use = "LR": Logistic regression (standard for ATAC)
# latent.vars = "nCount_ATAC": Control for sequencing depth
cat("Finding differentially accessible peaks...\n")
cat("(This may take several minutes for large datasets)\n")

da_peaks <- FindMarkers(pbmc, 
                        ident.1 = ident_1, 
                        ident.2 = ident_2,
                        test.use = "LR", 
                        latent.vars = "nCount_ATAC",
                        min.pct = 0.05)  # Only test peaks in =5% of cells

# Sort by significance
da_peaks <- da_peaks[order(da_peaks$p_val_adj), ]

# Save results
da_file <- file.path(tabdir, paste0("DA_peaks_", ident_1, "_vs_", ident_2, ".csv"))
write.csv(da_peaks, da_file)

# Summarize results
n_sig <- sum(da_peaks$p_val_adj < 0.05, na.rm = TRUE)
n_up <- sum(da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC > 0, na.rm = TRUE)
n_down <- sum(da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC < 0, na.rm = TRUE)

cat("\nDifferential accessibility results:\n")
cat("  Total peaks tested:         ", nrow(da_peaks), "\n")
cat("  Significant peaks (FDR<0.05):", n_sig, "\n")
cat("    More accessible in", ident_1, ":", n_up, "\n")
cat("    More accessible in", ident_2, ":", n_down, "\n")
cat("  Results saved to:", da_file, "\n\n")

# Visualize top differentially accessible peak
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
# STEP 10: CHROMATIN ACCESSIBILITY TRACKS (COVERAGE PLOTS)
# Purpose: Visualize accessibility at specific genes
###############################################################################
cat("\n=== STEP 10: Creating coverage plots for key genes ===\n")

# Set identity classes for coverage plot grouping
if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  levels(pbmc) <- unique(pbmc$predicted.id)
  cat("Grouping coverage by cell type\n")
} else {
  cat("Grouping coverage by cluster\n")
}

# Key genes for your research (B cells, glycosylation, SLE)
key_genes <- c("MS4A1", "CD19", "MGAT5", "ST6GAL1", "FUT8", "B4GALT1", 
               "IMPDH2", "IFI44L", "CD3D", "CD4", "PRDM1", "CD27")

cat("Generating coverage plots for", length(key_genes), "genes...\n")
cat("(Each plot may take 1-2 minutes)\n\n")

for (gene in key_genes) {
  cat("  Plotting:", gene, "...")
  
  tryCatch({
    # Create coverage plot
    # extend.upstream/downstream: Show regulatory regions around gene
    p <- CoveragePlot(pbmc, 
                      region = gene, 
                      extend.upstream = 40000,    # 40kb upstream (promoter)
                      extend.downstream = 20000)  # 20kb downstream
    
    save_plot(p, paste0("coverage_", gene), w=14, h=8)
    cat(" done\n")
    
  }, error = function(e) {
    cat(" FAILED -", e$message, "\n")
  })
}

cat("\nCoverage plot generation complete\n")

# Save final object
save_rds(pbmc, "atac_final")

###############################################################################
# ANALYSIS COMPLETE - SUMMARY
###############################################################################
cat("\n\n")
cat("-----------------------------------------------------------\n")
cat("               ATAC-seq ANALYSIS COMPLETE\n")
cat("-----------------------------------------------------------\n\n")

cat("Final ATAC-seq object:\n")
print(pbmc)

cat("\n?? Output directories:\n")
cat("  Main:   ", outs_dir, "\n")
cat("  Plots:  ", plotdir, "\n")
cat("  Tables: ", tabdir, "\n")
cat("  RDS:    ", rdsdir, "\n\n")

cat("?? Key files generated:\n")
cat("  1. QC metrics:              atac_qc_metrics.csv\n")
cat("  2. Filtering thresholds:    filtering_thresholds_applied.csv\n")
cat("  3. Cell type counts:        celltype_counts.csv\n")
cat("  4. Differential peaks:      DA_peaks_*.csv\n")
cat("  5. Session info:            sessionInfo.txt\n")
cat("  6. Package versions:        packages_loaded.txt\n\n")

cat("?? Analysis summary:\n")
cat("  MAD stringency used:        ", N_MAD_QC, "\n")
cat("  Cells retained:             ", ncol(pbmc), "\n")
cat("  Peaks analyzed:             ", nrow(pbmc), "\n")
cat("  Clusters identified:        ", nlevels(Idents(pbmc)), "\n")
if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  cat("  Cell types annotated:       ", length(unique(pbmc$predicted.id)), "\n")
}

cat("\n?? Next steps for your lupus nephritis project:\n")
cat("  1. ? Basic QC and clustering complete\n")
cat("  2. ? Integrate CITE-seq protein data (see template in tables/)\n")
cat("  3. ? Add clinical metadata (nephritis status, disease activity)\n")
cat("  4. ? Re-run DA analysis comparing nephritis groups:\n")
cat("       pbmc$nephritis <- your_clinical_data\n")
cat("       Idents(pbmc) <- pbmc$nephritis\n")
cat("       FindMarkers(pbmc, ident.1='with', ident.2='without')\n")
cat("  5. ? Focus on B cell glycosylation genes (MGAT5, ST6GAL1, etc.)\n")
cat("  6. ? Correlate accessibility with disease severity\n\n")

cat("?? For Methods section:\n")
cat("  'Single-nucleus ATAC-seq data were processed using Signac v", 
    as.character(packageVersion("Signac")), " and Seurat v", 
    as.character(packageVersion("Seurat")), ". ", sep="")
cat("Quality control filtering employed median absolute deviation (MAD)-based\n")
cat("  thresholds (", N_MAD_QC, " MAD) to objectively identify outlier cells. ", sep="")
cat("Cells were retained\n")
cat("  if they had ", round(thresholds$ncount_lower), "-", 
    round(thresholds$ncount_upper), " fragments, TSS enrichment >", 
    round(thresholds$tss_lower, 2), ",\n", sep="")
cat("  nucleosome signal <", round(thresholds$nuc_upper, 2), 
    ", and percent reads in peaks >", round(thresholds$frip_lower, 1), "%. ", sep="")
cat("This data-driven\n")
cat("  approach ensures reproducibility while adapting to dataset-specific quality\n")
cat("  distributions. Dimensionality reduction used TF-IDF normalization followed by\n")
cat("  LSI (components 2-30), with clustering performed at resolution ", 
    CLUSTER_RESOLUTION, ".'\n\n", sep="")

cat("??  Analysis completed:", format(Sys.time()), "\n")
cat("??  Total runtime:      ", round(difftime(Sys.time(), 
    as.POSIXct(readLines(logfile)[2]), units="mins"), 1), " minutes\n\n", sep="")

cat("-----------------------------------------------------------\n")

# Close log file
sink()

cat("\n? Full log saved to:", logfile, "\n")
cat("? All outputs saved to:", outs_dir, "\n\n")

###############################################################################
# END OF SCRIPT
###############################################################################
