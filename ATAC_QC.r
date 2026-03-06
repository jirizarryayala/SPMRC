###############################################################################
# COMPLETE ATAC + CITE-seq Pipeline - MATRIX APPROACH (no HDF5)
# Reproducible + HPC-safe setup block 
# MODIFIED: STEP 9 (DA) runs separately to avoid session timeout
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

# Required for your pipeline
need("Seurat")
need("Signac")
need("EnsDb.Hsapiens.v75")
need("GenomicRanges")
need("Matrix")

# Optional (nice-to-have only)
has_tidyverse <- want("tidyverse")
has_patchwork <- want("patchwork")

## -------------------- SETUP (dirs + logging) --------------------------------
outs_dir <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/ATAC_143"
dir.create(outs_dir, recursive = TRUE, showWarnings = FALSE)

plotdir <- file.path(outs_dir, "plots")
rdsdir  <- file.path(outs_dir, "rds")
tabdir  <- file.path(outs_dir, "tables")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir,  recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(outs_dir, "run.log")
sink(logfile, split = TRUE)
cat("===== COMPLETE ANALYSIS - YOUR DATA =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("R.version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")

# Write provenance immediately (this is what makes it reproducible)
provfile <- file.path(outs_dir, "sessionInfo.txt")
cat("Writing sessionInfo to:", provfile, "\n")
writeLines(c(capture.output(sessionInfo())), provfile)

# Optional: record package versions explicitly (best-effort)
pkgfile <- file.path(outs_dir, "packages_loaded.txt")
pkgs_loaded <- c("Seurat","Signac","EnsDb.Hsapiens.v75","GenomicRanges","Matrix",
                 if (has_tidyverse) "tidyverse" else NA,
                 if (has_patchwork) "patchwork" else NA)
pkgs_loaded <- pkgs_loaded[!is.na(pkgs_loaded)]
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


cat("\n=== STEP 1: Reading ATAC-seq data (outs_143) ===\n")

library(Seurat)
library(Signac)

atac_dir <- "/lustre/project/crosslin/crosslin_team/irizarry/Bhargava/ATACseq/110425 snATAC_LibDD_hPBMC_Anita_8_data/110425 snATAC_LibDD_hPBMC_Anita_8_cellranger/outs_143"

# Check MTX directory contents
mtx_dir <- file.path(atac_dir, "filtered_peak_bc_matrix")
cat("Files in MTX directory:\n")
print(list.files(mtx_dir))

# Read the matrix components manually for ATAC data
library(Matrix)

# Read matrix
mtx_file <- file.path(mtx_dir, "matrix.mtx.gz")
if (!file.exists(mtx_file)) mtx_file <- file.path(mtx_dir, "matrix.mtx")
counts_matrix <- readMM(mtx_file)

# Read barcodes
barcode_file <- file.path(mtx_dir, "barcodes.tsv.gz")
if (!file.exists(barcode_file)) barcode_file <- file.path(mtx_dir, "barcodes.tsv")
barcodes <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)[,1]

# Read peaks (features) - ATAC uses peaks.bed instead of features.tsv
peak_file <- file.path(mtx_dir, "peaks.bed.gz")
if (!file.exists(peak_file)) peak_file <- file.path(mtx_dir, "peaks.bed")
if (!file.exists(peak_file)) peak_file <- file.path(mtx_dir, "features.tsv.gz")
if (!file.exists(peak_file)) peak_file <- file.path(mtx_dir, "features.tsv")

peaks <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE)

# Create peak names in format chr:start-end
if (ncol(peaks) >= 3) {
  peak_names <- paste0(peaks[,1], ":", peaks[,2], "-", peaks[,3])
} else {
  peak_names <- peaks[,1]
}

# Set dimnames
rownames(counts_matrix) <- peak_names
colnames(counts_matrix) <- barcodes

cat("Matrix dimensions:", nrow(counts_matrix), "peaks x", ncol(counts_matrix), "cells\n")

# Load metadata
metadata_path <- file.path(atac_dir, "singlecell.csv")
metadata <- read.csv(metadata_path, header = TRUE, row.names = 1)

# Fragment file
frag_path <- file.path(atac_dir, "fragments.tsv.gz")
if (!file.exists(frag_path)) stop("Fragment file not found!")

# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts_matrix,
  sep = c(":", "-"),
  fragments = frag_path,
  min.cells = 10,
  min.features = 200
)

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = "ATAC"
)

cat("ATAC object created:\n")
print(pbmc)

# Save
save_rds(pbmc, "atac_01_raw")


###############################################################################
# STEP 2: ADD GENE ANNOTATIONS
###############################################################################
cat("\n=== STEP 2: Adding gene annotations ===\n")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- "UCSC"
Annotation(pbmc) <- annotations

cat("Annotations added:", length(annotations), "genes\n")
save_rds(pbmc, "atac_02_annotated")

###############################################################################
# STEP 3: QC METRICS (with TSS enrichment)
###############################################################################
cat("\n=== STEP 3: Computing QC metrics ===\n")

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

write.csv(pbmc@meta.data, file.path(tabdir, "atac_qc_metrics.csv"))

# QC plots
p1 <- DensityScatter(pbmc, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
p2 <- DensityScatter(pbmc, x = "nucleosome_signal", y = "TSS.enrichment", quantiles = TRUE)
save_plot(p1 | p2, "atac_qc_density", w=14, h=6)

p_vln <- VlnPlot(pbmc, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", 
                                     "nucleosome_signal", "blacklist_ratio", "pct_reads_in_peaks"),
                 pt.size = 0.1, ncol = 6)
save_plot(p_vln, "atac_qc_violin", w=18, h=6)

save_rds(pbmc, "atac_03_qc")





###############################################################################
# STEP 4: FILTER CELLS (FOR HIGH-QUALITY SAMPLE)
###############################################################################
cat("\n=== STEP 4: Filtering cells ===\n")
n_before <- ncol(pbmc)
cat("Before filtering:", n_before, "cells\n")

# Print QC summaries
cat("\nQC metric summaries:\n")
cat("  nCount_ATAC:        ", summary(pbmc$nCount_ATAC), "\n")
cat("  TSS.enrichment:     ", summary(pbmc$TSS.enrichment), "\n")
cat("  nucleosome_signal:  ", summary(pbmc$nucleosome_signal), "\n")
cat("  pct_reads_in_peaks: ", summary(pbmc$pct_reads_in_peaks), "\n")

# Stringent filtering for high-quality data
pbmc <- subset(pbmc, subset = 
  nCount_ATAC > 2000 &           # Higher minimum (remove low tail)
  nCount_ATAC < 20000 &          # Remove extreme outliers/doublets
  pct_reads_in_peaks > 30 &      # Take advantage of high quality
  nucleosome_signal < 2.0 &      # Tight control
  TSS.enrichment > 1.5           # Higher threshold - your data supports it
)

cat("After filtering:", ncol(pbmc), "cells\n")
cat("Retained:", round(ncol(pbmc)/n_before*100, 1), "%\n")

save_rds(pbmc, "atac_04_filtered")

###############################################################################
# STEP 5: NORMALIZATION AND CLUSTERING
###############################################################################
cat("\n=== STEP 5: Normalization and clustering ===\n")

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

p_depth <- DepthCor(pbmc)
save_plot(p_depth, "atac_lsi_depth", w=10, h=6)

pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- RunTSNE(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.5)

p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
p_tsne <- DimPlot(pbmc, reduction = "tsne", label = TRUE) + NoLegend()
save_plot(p_umap | p_tsne, "atac_umap_tsne", w=14, h=6)

save_rds(pbmc, "atac_05_clustered")

###############################################################################
# STEP 6: GENE ACTIVITY MATRIX
###############################################################################
cat("\n=== STEP 6: Creating gene activity matrix ===\n")

gene.activities <- GeneActivity(pbmc)
pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(pbmc, assay = "RNA", normalization.method = "LogNormalize",
                      scale.factor = median(pbmc$nCount_RNA))

save_rds(pbmc, "atac_06_with_gene_activity")

###############################################################################
# STEP 7: VISUALIZE YOUR MARKER GENES
###############################################################################
cat("\n=== STEP 7: Visualizing marker genes ===\n")

DefaultAssay(pbmc) <- "RNA"

# Your research-specific markers
markers <- list(
  b_cell = c("MS4A1", "CD19", "CD79A", "PAX5", "CD27"),
  glycosylation = c("MGAT5", "ST6GAL1", "FUT8", "B4GALT1", "MGAT1", "MGAT2"),
  metabolism = c("IMPDH2", "GLS", "FASN", "SLC2A1", "HK2", "LDHA"),
  sle = c("IFI44L", "PRDM1", "CD38", "C3", "IRF7", "MX1"),
  t_cell = c("CD3D", "CD3E", "CD4", "CD8A"),
  myeloid = c("CD14", "FCGR3A", "CD68", "LYZ")
)

for (marker_name in names(markers)) {
  genes <- markers[[marker_name]]
  genes_present <- genes[genes %in% rownames(pbmc[["RNA"]])]
  
  if (length(genes_present) >= 2) {
    p <- FeaturePlot(pbmc, features = genes_present, reduction = "umap",
                     pt.size = 0.1, max.cutoff = "q95", ncol = 3)
    save_plot(p, paste0("markers_", marker_name, "_umap"), 
              w=12, h=ceiling(length(genes_present)/3)*4)
  }
}




###############################################################################
# STEP 8: LABEL TRANSFER
###############################################################################
cat("\n=== STEP 8: Label transfer ===\n")

rna_ref <- "/lustre/project/crosslin/crosslin_team/Sharifi/notebook/SPMRC/pbmc_10k_v3.rds"

if (file.exists(rna_ref)) {
  pbmc_rna <- readRDS(rna_ref)
  pbmc_rna <- UpdateSeuratObject(pbmc_rna)
  
  transfer.anchors <- FindTransferAnchors(
    reference = pbmc_rna, 
    query = pbmc, 
    reduction = "pcaproject",
    k.anchor = 5
  )
  
  n_anchors <- nrow(transfer.anchors@anchors)
  k_weight_use <- max(5, min(20, floor(n_anchors * 0.2)))
  cat("Found", n_anchors, "anchors. Using k.weight =", k_weight_use, "\n")
  
  predicted.labels <- TransferData(
    anchorset = transfer.anchors, 
    refdata = pbmc_rna$celltype,
    weight.reduction = pbmc[["lsi"]], 
    dims = 2:30, 
    k.weight = k_weight_use
  )
  
  pbmc <- AddMetaData(pbmc, metadata = predicted.labels)
  
  p_pred <- DimPlot(pbmc, reduction = "umap", group.by = "predicted.id",
                    label = TRUE, repel = TRUE) + NoLegend()
  save_plot(p_pred, "atac_predicted_celltypes", w=10, h=7)
  
  write.csv(table(pbmc$predicted.id), file.path(tabdir, "celltype_counts.csv"))
  save_rds(pbmc, "atac_07_labeled")
} else {
  cat("WARNING: Reference RNA-seq data not found. Skipping label transfer.\n")
}


###############################################################################
# STEP 9: DIFFERENTIAL ACCESSIBILITY - SKIPPED (RUN SEPARATELY)
###############################################################################
cat("\n=== STEP 9: Differential accessibility - DEFERRED ===\n")
cat("NOTICE: DA analysis skipped in main pipeline to avoid session timeout.\n")
cat("For large datasets (10K+ cells), DA should be run separately.\n\n")

cat("To run DA analysis:\n")
cat("  Rscript DA_ATAC_143_standalone.r\n\n")

cat("This will:\n")
cat("  - Load atac_07_labeled.rds (from Step 8)\n")
cat("  - Run faithful LR test with latent.vars='nCount_ATAC'\n")
cat("  - Save results to same directories\n")
cat("  - Generate same plots\n\n")

cat("Proceeding to Step 10 (Coverage plots)...\n")

###############################################################################
# STEP 10: COVERAGE PLOTS AT KEY GENES
###############################################################################
cat("\n=== STEP 10: Coverage plots ===\n")

if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  levels(pbmc) <- unique(pbmc$predicted.id)
}

key_genes <- c("MS4A1", "CD19", "MGAT5", "ST6GAL1", "FUT8", "B4GALT1", 
               "IMPDH2", "IFI44L", "CD3D", "CD4")

for (gene in key_genes) {
  tryCatch({
    p <- CoveragePlot(pbmc, region = gene, extend.upstream = 40000, 
                      extend.downstream = 20000)
    save_plot(p, paste0("coverage_", gene), w=14, h=8)
  }, error = function(e) cat("Could not plot", gene, "\n"))
}

save_rds(pbmc, "atac_final")

###############################################################################
# PART 2: CITE-seq INTEGRATION
###############################################################################
cat("\n\n===== CITE-seq INTEGRATION =====\n")

cite_base <- "/lustre/project/crosslin/crosslin_team/irizarry/Bhargava/CITEseq/Anita_scCITE_GEX_hBCR_8_data/251002_CITE5hBCRGEXv3_Anita_8_cellranger"

cite_samples <- c("Anita_5'BCR_8samples", "Anita_CITE_GEX_8samples", 
                  "HighACE_BCR_4samples", "HighACE_CITE_GEX_4samples")

cat("\nAvailable CITE-seq samples:\n")
for (sample in cite_samples) {
  sample_path <- file.path(cite_base, sample)
  if (dir.exists(sample_path)) {
    cat("  ✓", sample, "\n")
  }
}

cat("\n=== CITE-seq Integration Instructions ===\n")
cat("1. Choose which CITE-seq sample to integrate\n")
cat("2. Load the filtered_feature_bc_matrix from that sample\n")
cat("3. Example for Anita_CITE_GEX_8samples:\n\n")

cat("cite_dir <- file.path(cite_base, 'Anita_CITE_GEX_8samples/outs/filtered_feature_bc_matrix')\n")
cat("cite_data <- Read10X(cite_dir)\n")
cat("pbmc_cite <- CreateSeuratObject(counts = cite_data$`Gene Expression`)\n")
cat("pbmc_cite[['ADT']] <- CreateAssayObject(counts = cite_data$`Antibody Capture`)\n\n")

template <- "
###############################################################################
# CITE-seq Integration Template
###############################################################################

cite_dir <- '/lustre/project/crosslin/crosslin_team/irizarry/Bhargava/CITEseq/Anita_scCITE_GEX_hBCR_8_data/251002_CITE5hBCRGEXv3_Anita_8_cellranger/Anita_CITE_GEX_8samples/outs/filtered_feature_bc_matrix'

cite_data <- Read10X(cite_dir)
pbmc_cite <- CreateSeuratObject(counts = cite_data\\$\\`Gene Expression\\`)
pbmc_cite[['ADT']] <- CreateAssayObject(counts = cite_data\\$\\`Antibody Capture\\`)

pbmc_cite <- NormalizeData(pbmc_cite)
pbmc_cite <- FindVariableFeatures(pbmc_cite)
DefaultAssay(pbmc_cite) <- 'ADT'
pbmc_cite <- NormalizeData(pbmc_cite, normalization.method = 'CLR', margin = 2)

pbmc_atac <- readRDS('atac_final.rds')

DefaultAssay(pbmc_cite) <- 'RNA'
DefaultAssay(pbmc_atac) <- 'RNA'

anchors <- FindTransferAnchors(
  reference = pbmc_cite, query = pbmc_atac,
  reference.assay = 'RNA', query.assay = 'RNA',
  reduction = 'pcaproject'
)

protein.pred <- TransferData(
  anchorset = anchors,
  refdata = GetAssayData(pbmc_cite, assay = 'ADT'),
  weight.reduction = pbmc_atac[['lsi']],
  dims = 2:30
)

pbmc_atac[['ADT']] <- CreateAssayObject(data = protein.pred)

DefaultAssay(pbmc_atac) <- 'ADT'
FeaturePlot(pbmc_atac, features = c('CD19', 'CD20', 'CD27', 'CD38'), 
            reduction = 'umap', ncol = 2)

b_cells <- subset(pbmc_atac, subset = CD19 > 1 & CD20 > 1)
"

writeLines(template, file.path(tabdir, "CITE_integration_template.R"))

###############################################################################
# SUMMARY
###############################################################################
cat("\n\n===== ANALYSIS COMPLETE (STEPS 1-8, 10) =====\n")
cat("ATAC-seq object:\n")
print(pbmc)
cat("\nOutputs in:", outs_dir, "\n")
cat("  Plots:", plotdir, "\n")
cat("  Tables:", tabdir, "\n")
cat("  RDS:", rdsdir, "\n")

cat("\nNext steps:\n")
cat("1. Run DA analysis separately: Rscript DA_ATAC_143_standalone.r\n")
cat("2. Review QC plots and marker gene expression\n")
cat("3. Choose CITE-seq sample to integrate\n")
cat("4. Add clinical metadata (nephritis status)\n")

cat("\nCompleted:", format(Sys.time()), "\n")
sink()
###############################################################################
