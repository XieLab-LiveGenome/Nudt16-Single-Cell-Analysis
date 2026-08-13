# =============================================================================
# MASTER SCRIPT — NUDT16 KO Female Spleen scRNA-seq  (F1 + F2 : 2 WT + 2 KO)
# =============================================================================
# ONE script to regenerate EVERY table and figure (main + supplementary) for
# the female NUDT16 spleen publication. Edit ONLY this file going forward.
#
# Run mode : ALWAYS RUN EVERYTHING (full pipeline recomputed each run; no cache).
# Primary annotation : celltype_marker  (canonical-marker cluster scoring).
#                      SingleR/ImmGen is kept as a SECONDARY validation label.
#
# SECTION MAP
#   0.  Config, libraries, helpers, palettes, marker sets, DDR gene sets
#   1.  Preprocessing  (load -> QC -> doublets -> MAD filter -> SCT ->
#                       Harmony -> cluster -> UMAP -> SingleR secondary ->
#                       composition)  --> saves female_obj_processed.rds
#   2.  Cluster QC diagnostics (flag doublet/junk clusters)
#   3.  Marker-score annotation (AddModuleScore z-score -> cluster labels)
#   4.  Finalize annotation  --> commits celltype_marker PRIMARY,
#                                saves female_obj_annotated.rds
#   5.  Marker dot plot + annotation-check dot plot
#   6.  Differential expression
#         6A. Pseudobulk DESeq2 (KO vs WT) : global + per celltype_marker
#         6B. Seurat FindMarkers (Wilcoxon): global + per celltype_marker
#   7.  Volcano plots  : BOTH methods x (global + each cell type)
#   8.  Figure panels  : UMAP highlights/panel, composition (dodged +
#                        diverging), NUDT16 validation, validation violins,
#                        DDR heatmap / WT-vs-KO / diagnostics
#   9.  Enrichment     : 9A pseudobulk ORA (GO/KEGG/Reactome/Hallmark)
#                        9B Hallmark GSEA + ORA on Seurat DE
#   10. Excel export (multi-sheet) + sessionInfo
#
# Author project path (Mac): /Users/budankm/Desktop/Sequencing/GONG
# Canonical output dir      : results/females/
# =============================================================================


# =============================================================================
# SECTION 0 — CONFIG, LIBRARIES, HELPERS, PALETTES, MARKER SETS
# =============================================================================

# ---- 0.1 Libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  # Core
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)

  # Normalization & integration
  library(glmGamPoi)
  library(harmony)

  # QC
  library(scDblFinder)
  library(scater)
  library(SingleCellExperiment)
  library(BiocParallel)

  # Annotation
  library(SingleR)
  library(celldex)

  # Differential expression
  library(DESeq2)

  # Enrichment
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(msigdbr)
  library(fgsea)

  # Visualization
  library(pheatmap)
  library(RColorBrewer)

  # Export
  library(openxlsx)
})

# Optional packages (loaded defensively; pipeline degrades gracefully) --------
HAVE_REACTOME  <- requireNamespace("ReactomePA",     quietly = TRUE)
HAVE_COMPLEXHM <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
                  requireNamespace("circlize",       quietly = TRUE)
HAVE_SPECKLE   <- requireNamespace("speckle",        quietly = TRUE)

# ---- 0.2 conflicted guards --------------------------------------------------
# ---- 0.2 conflicted guards --------------------------------------------------
# Register each base preference INDIVIDUALLY. Batching them into one
# conflicts_prefer() call is fragile: if any single symbol trips it, the whole
# call fails under try() and NONE register (incl. base::unname) -> Section 4
# aborts. One-by-one guarantees base::unname et al. are always set.
if (requireNamespace("conflicted", quietly = TRUE)) {
  try(conflicted::conflict_prefer_all("dplyr", quiet = TRUE), silent = TRUE)

  .prefer_base <- c(
    "rownames","colnames","intersect","setdiff","union","which","is.element",
    "table","which.max","which.min","as.data.frame","as.vector","as.factor",
    "unname","unlist","unique","duplicated","order","sort","rank","match",
    "append","rowSums","colSums","rowMeans","colMeans","scale","mean","sample",
    "apply","paste","paste0","Reduce","Map","Filter","Find","Position","get",
    "pmin","pmax","tapply","cbind","rbind","do.call","lengths","levels",
    "head","tail","nrow","ncol","t","expand.grid","mapply","vapply","sapply",
    "lapply","toupper","tolower","trimws","strsplit","gsub","sub","grepl",
    "grep","format","rev","factor","as.character","as.numeric","as.integer",
    "as.logical","is.na","na.omit"
  )
  for (.fn in unique(.prefer_base)) {
    suppressMessages(suppressWarnings(try(
      conflicted::conflict_prefer(.fn, "base", quiet = TRUE), silent = TRUE)))
  }
  for (.so in c("Assays","Layers")) {
    suppressMessages(suppressWarnings(try(
      conflicted::conflict_prefer(.so, "SeuratObject", quiet = TRUE), silent = TRUE)))
  }
  suppressWarnings(rm(list = intersect(c(".prefer_base",".fn",".so"), ls())))
}
set.seed(217)
options(Seurat.object.assay.version = "v5")

# ---- 0.3 Paths --------------------------------------------------------------
project_dir <- "/Users/budankm/Desktop/Sequencing/GONG"
RES_DIR     <- file.path(project_dir, "results", "females")

# All output subfolders used anywhere in this script (kept per binding decision
# "keep current subfolders"). Created up-front so every section can write.
OUT_SUBDIRS <- c(
  # preprocessing / core
  "qc", "umap", "composition", "volcano", "enrichment", "gene_plots", "heatmap",
  # annotation + qc
  "marker_annotation", "cluster_qc",
  # per-panel
  "panelE",
  # volcano method splits
  file.path("volcano", "deseq2"), file.path("volcano", "seurat"),
  # enrichment engines
  "enrichment_hallmark",
  # DDR
  "ddr_heatmap", "ddr_wt_vs_ko", "ddr_diagnostics",
  # tables
  "tables"
)
for (sub in OUT_SUBDIRS) {
  dir.create(file.path(RES_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}
if (dir.exists(project_dir)) setwd(project_dir)

# Convenience path helper
rf <- function(...) file.path(RES_DIR, ...)

# ---- 0.4 Global run knobs (single source of truth) --------------------------
# Primary cell-type column used by ALL downstream DE / figures / enrichment.
CELLTYPE_COL <- "celltype_marker"     # BINDING: marker-based primary
SECONDARY_COL <- "ImmGen.labels"      # SingleR secondary validation

# Manual annotation overrides (applied in Section 4). Leave empty by default.
#   DROP_CLUSTERS : character vector of cluster ids to remove after labelling
#   MANUAL_LABELS : named list  list("3" = "T cells", "7" = "DC")
DROP_CLUSTERS <- character(0)
MANUAL_LABELS <- list()

# QC / filtering
MAD_NMADS             <- 3
DOUBLET_SEED          <- 217
MIN_CELLS_PER_LABEL   <- 15      # ImmGen label retention threshold

# Embedding (harmony pipeline from preprocessfinal — NOT dims1:30/res0.5)
N_VAR_FEATURES        <- 3000
N_PCS                 <- 50
N_DIMS_USE            <- 20
CLUSTER_RES           <- 0.2

# Pseudobulk / Seurat DE
MIN_CELLS_PER_SAMPLE  <- 10
MIN_REPS_PER_GENOTYPE <- 2
MIN_COUNT             <- 10
MIN_SAMPLES           <- 2
LFC_CUT               <- 0.5
PADJ_CUT              <- 0.05

# Enrichment thresholds
ENR_PVAL              <- 0.05
ENR_QVAL              <- 0.25
SET_MIN               <- 10
SET_MAX               <- 500

# Volcano labelling
VOLC_N_LABEL          <- 20
VOLC_ANCHOR           <- "Nudt16"

# Figure output
FIG_DPI               <- 400

# ---- 0.5 Color anchors / palettes ------------------------------------------
COL_UP   <- "#B2182B"     # KO-up  (red)
COL_DOWN <- "#2166AC"     # WT-up  (blue)
COL_NS   <- "grey75"
COL_HM   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)

genotype_colors <- c(WT = "#3B6C9C", KO = "#D6604D")

# Okabe-Ito based cell-type palette (stable across all figures)
CELLTYPE_COLORS <- c(
  "B cells"           = "#0072B2",
  "T cells"           = "#D55E00",
  "NKT"               = "#00B0BE",
  "NK cells"          = "#009E73",
  "ILC"               = "#17BECF",
  "Tgd"               = "#BCBD22",
  "Monocytes"         = "#CC79A7",
  "Macrophages"       = "#E69F00",
  "DC"                = "#56B4E9",
  "Neutrophils"       = "#F0E442",
  "Basophils"         = "#8C564B",
  "Endothelial cells" = "#B279A2",
  "Stem cells"        = "#000000"
)

build_celltype_palette <- function(types) {
  types <- as.character(types)
  pal <- CELLTYPE_COLORS[types]
  miss <- is.na(pal)
  if (any(miss)) {
    extra <- scales::hue_pal()(sum(miss))
    pal[miss] <- extra
  }
  names(pal) <- types
  pal
}

lighten_hex <- function(hex, amt = 0.40) {
  rgb0 <- grDevices::col2rgb(hex) / 255
  rgb1 <- rgb0 + (1 - rgb0) * amt
  grDevices::rgb(rgb1[1, ], rgb1[2, ], rgb1[3, ])
}

# ---- 0.6 Canonical marker sets (13 immune/stromal types) --------------------
marker_list <- list(
  "B cells"           = c("Cd79a", "Ms4a1", "Cd79b", "Cd19"),
  "T cells"           = c("Cd3e", "Cd3d", "Cd3g", "Cd2"),
  "NKT"               = c("Zbtb16", "Klrb1c"),
  "NK cells"          = c("Ncr1", "Gzma", "Nkg7", "Klrk1"),
  "Tgd"               = c("Trdc", "Tcrg-C1"),
  "ILC"               = c("Il7r", "Gata3"),
  "DC"                = c("Itgax", "Flt3", "Xcr1"),
  "Monocytes"         = c("Ly6c2", "Ccr2"),
  "Macrophages"       = c("Adgre1", "C1qa", "C1qb"),
  "Neutrophils"       = c("S100a8", "Retnlg", "S100a9"),
  "Basophils"         = c("Mcpt8", "Cpa3", "Prss34"),
  "Endothelial cells" = c("Pecam1", "Cdh5"),
  "Stem cells"        = c("Cd34", "Kit")
)

# Compact 2-gene dot-plot marker set (block-diagonal display)
marker_list_dot <- list(
  "B cells"           = c("Cd79a", "Ms4a1"),
  "T cells"           = c("Cd3e", "Cd3d"),
  "NKT"               = c("Zbtb16", "Klrb1c"),
  "NK cells"          = c("Ncr1", "Nkg7"),
  "Tgd"               = c("Trdc", "Tcrg-C1"),
  "ILC"               = c("Il7r", "Gata3"),
  "DC"                = c("Itgax", "Flt3"),
  "Monocytes"         = c("Ly6c2", "Ccr2"),
  "Macrophages"       = c("Adgre1", "C1qa"),
  "Neutrophils"       = c("S100a8", "Retnlg"),
  "Basophils"         = c("Mcpt8", "Cpa3"),
  "Endothelial cells" = c("Pecam1", "Cdh5"),
  "Stem cells"        = c("Cd34", "Kit")
)

# ---- 0.7 DDR / DNA-repair gene sets (Knijnenburg et al., Cell Rep 2018) -----
# Each entry: pathway -> named list of gene -> synonym candidates.
DDR_CITATION <- paste(
  "DDR/DNA-repair gene sets adapted from Knijnenburg TA et al.,",
  "'Genomic and Molecular Landscape of DNA Damage Repair Deficiency across",
  "The Cancer Genome Atlas', Cell Reports 23(1):239-254 (2018)."
)

ddr_focused <- list(
  "BER"        = list(Ogg1 = "Ogg1", Mutyh = "Mutyh", Nthl1 = "Nthl1",
                      Apex1 = "Apex1", Parp1 = "Parp1", Xrcc1 = "Xrcc1"),
  "MMR"        = list(Msh2 = "Msh2", Msh6 = "Msh6", Mlh1 = "Mlh1", Pms2 = "Pms2"),
  "NER"        = list(Xpa = "Xpa", Xpc = "Xpc", Ercc1 = "Ercc1", Ercc4 = "Ercc4"),
  "HR"         = list(Brca1 = "Brca1", Brca2 = "Brca2", Rad51 = "Rad51",
                      Bard1 = "Bard1"),
  "NHEJ"       = list(Prkdc = "Prkdc", Xrcc6 = c("Xrcc6", "Ku70"),
                      Xrcc5 = c("Xrcc5", "Ku80"), Lig4 = "Lig4"),
  "FA"         = list(Fancd2 = "Fancd2", Fanca = "Fanca", Fancc = "Fancc"),
  "Direct-DR"  = list(Mgmt = "Mgmt"),
  "TLS"        = list(Rev1 = "Rev1", Polh = "Polh", Rev3l = "Rev3l"),
  "Checkpoint" = list(Atm = "Atm", Atr = "Atr", Chek1 = "Chek1",
                      Chek2 = "Chek2", Trp53 = c("Trp53", "Tp53")),
  "Nt.sanit."  = list(Nudt16 = "Nudt16", Nudt1 = c("Nudt1", "Mth1"),
                      Ung = "Ung", Dut = "Dut")
)

# ---- 0.8 Generic helpers ----------------------------------------------------
# Resolve a metadata column name from candidates present in an object.
pick_col <- function(obj, candidates, fallback = "seurat_clusters") {
  md <- colnames(obj@meta.data)
  hit <- candidates[candidates %in% md]
  if (length(hit)) hit[1] else fallback
}

# Resolve gene symbols/synonyms against rownames of an assay.
resolve_genes <- function(genes, universe) {
  out <- character(0)
  for (g in genes) {
    hit <- g[g %in% universe]
    if (length(hit)) out <- c(out, hit[1])
  }
  unique(out)
}

# Flatten a DDR pathway list to a resolved gene table given a feature universe.
resolve_ddr <- function(ddr, universe) {
  rows <- list()
  for (pw in names(ddr)) {
    for (canon in names(ddr[[pw]])) {
      hit <- resolve_genes(list(ddr[[pw]][[canon]]), universe)
      if (length(hit)) {
        rows[[length(rows) + 1]] <- data.frame(
          pathway = pw, gene = hit[1], canonical = canon,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(data.frame(pathway = character(0),
                                       gene = character(0),
                                       canonical = character(0)))
  do.call(rbind, rows)
}

# Save a ggplot to pdf + png (+ optional svg) at consistent DPI.
save_plot <- function(plot, path_noext, width = 8, height = 6, svg = FALSE) {
  ggsave(paste0(path_noext, ".pdf"), plot, width = width, height = height,
         device = grDevices::cairo_pdf)
  ggsave(paste0(path_noext, ".png"), plot, width = width, height = height,
         dpi = FIG_DPI)
  if (svg) {
    tryCatch(
      ggsave(paste0(path_noext, ".svg"), plot, width = width, height = height),
      error = function(e) message("  svg skipped: ", conditionMessage(e))
    )
  }
  invisible(path_noext)
}

# Resolve a normalized expression assay (prefer RNA 'data', else normalize,
# else fall back to SCT). Returns the modified object + the assay name.
ensure_norm_assay <- function(obj) {
  assay_use <- NULL
  if ("RNA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "RNA"
    obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)
    has_data <- tryCatch({
      d <- GetAssayData(obj, assay = "RNA", layer = "data")
      !is.null(d) && length(d@x) > 0
    }, error = function(e) FALSE)
    if (!has_data) obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
    assay_use <- "RNA"
  } else if ("SCT" %in% Assays(obj)) {
    assay_use <- "SCT"
    DefaultAssay(obj) <- "SCT"
  } else {
    assay_use <- DefaultAssay(obj)
  }
  list(obj = obj, assay = assay_use)
}

# Wrap a figure-panel block so one failing panel never aborts the whole run.
safe_panel <- function(name, expr) {
  message("\n>>> PANEL: ", name)
  tryCatch(
    force(expr),
    error = function(e)
      message("  [SKIP] panel '", name, "' failed: ", conditionMessage(e))
  )
}

# Excel-safe sheet name (<=31 chars, unique-ish)
safe_sheet <- function(x) substr(gsub("[^A-Za-z0-9_]", "_", x), 1, 31)

message("\n================  MASTER NUDT16 FEMALES — START  ================\n")
message("Output dir: ", RES_DIR)
message("Primary annotation column: ", CELLTYPE_COL,
        "  |  Secondary: ", SECONDARY_COL)


# =============================================================================
# SECTION 1 — PREPROCESSING
#   load -> QC -> doublets -> MAD filter -> SCT -> Harmony -> cluster ->
#   UMAP -> SingleR (secondary) -> composition -> save processed RDS
# =============================================================================

# Sample manifest (female only — 2 WT, 2 KO). PIPseeker filtered matrices.
sample_dirs <- list(
  Nudt16Spleen_WTF1 = "Z_pipseeker_output/Nudt16Spleen_WTF1/filtered_matrix/sensitivity_4",
  Nudt16Spleen_KOF1 = "Z_pipseeker_output/Nudt16Spleen_KOF1/filtered_matrix/sensitivity_3",
  Nudt16Spleen_WTF2 = "Z_pipseeker_output/Nudt16Spleen_WTF2/filtered_matrix/sensitivity_2",
  Nudt16Spleen_KOF2 = "Z_pipseeker_output/Nudt16Spleen_KOF2/filtered_matrix/sensitivity_3"
)

# ---- 1.1 Load and merge -----------------------------------------------------
message("\n[1.1] Loading PIPseeker matrices...")

seurat_list <- lapply(names(sample_dirs), function(s) {
  mtx <- Read10X(data.dir = sample_dirs[[s]], gene.column = 2)
  obj <- CreateSeuratObject(counts = mtx, project = s,
                            min.cells = 3, min.features = 200)
  obj$sample     <- s
  obj$orig.ident <- s
  obj
})
names(seurat_list) <- names(sample_dirs)

seu <- merge(seurat_list[[1]], y = seurat_list[-1],
             add.cell.ids = names(seurat_list))

seu$genotype  <- factor(ifelse(grepl("KO", seu$sample), "KO", "WT"),
                        levels = c("WT", "KO"))
seu$replicate <- ifelse(grepl("F1", seu$sample), "F1", "F2")

cat("Cells loaded:", ncol(seu), "\n")

# ---- 1.2 QC metrics ---------------------------------------------------------
message("\n[1.2] Computing QC metrics...")

seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[ls]")

VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "sample", ncol = 4, pt.size = 0.1)
ggsave(rf("qc", "QC_violins_prefilter.pdf"), width = 12, height = 4)

# ---- 1.3 Doublet removal (scDblFinder, per sample) --------------------------
message("\n[1.3] Running scDblFinder per sample...")

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

sce_dbl <- as.SingleCellExperiment(seu, assay = "RNA")
sce_dbl <- scDblFinder(sce_dbl, samples = sce_dbl$sample,
                       BPPARAM = SerialParam(RNGseed = DOUBLET_SEED))

seu$scDblFinder_class <- sce_dbl$scDblFinder.class
seu$scDblFinder_score <- sce_dbl$scDblFinder.score

cat("Doublet rates per sample:\n")
print(round(prop.table(table(seu$sample, seu$scDblFinder_class), 1), 3))

seu <- subset(seu, subset = scDblFinder_class == "singlet")
cat("Cells after doublet removal:", ncol(seu), "\n")

# ---- 1.4 MAD-based QC (nmads = 3, per sample) -------------------------------
message("\n[1.4] MAD-based QC filtering...")

discard <-
  isOutlier(seu$nFeature_RNA, nmads = MAD_NMADS, type = "both",   log = TRUE,  batch = seu$sample) |
  isOutlier(seu$nCount_RNA,   nmads = MAD_NMADS, type = "higher", log = TRUE,  batch = seu$sample) |
  isOutlier(seu$percent.mt,   nmads = MAD_NMADS, type = "higher", log = FALSE, batch = seu$sample) |
  isOutlier(seu$percent.ribo, nmads = MAD_NMADS, type = "higher", log = FALSE, batch = seu$sample)

cat("Discarded:", sum(discard), " | Retained:", sum(!discard), "\n")
seu <- seu[, !discard]

cat("\nCells per sample (post-QC):\n");   print(table(seu$sample))
cat("\nCells per genotype (post-QC):\n"); print(table(seu$genotype))

VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "sample", ncol = 4, pt.size = 0.1)
ggsave(rf("qc", "QC_violins_postfilter.pdf"), width = 12, height = 4)

# ---- 1.5 SCTransform (per sample) -------------------------------------------
message("\n[1.5] SCTransform per sample...")

obj_list <- SplitObject(seu, split.by = "sample")
obj_list <- lapply(obj_list, SCTransform,
                   vars.to.regress = c("percent.mt", "percent.ribo"),
                   method     = "glmGamPoi",
                   vst.flavor = "v2",
                   verbose    = FALSE)

hvgs       <- SelectIntegrationFeatures(obj_list, nfeatures = N_VAR_FEATURES)
female_obj <- merge(obj_list[[1]], y = obj_list[-1], merge.data = TRUE)
VariableFeatures(female_obj) <- hvgs
DefaultAssay(female_obj)     <- "SCT"

# ---- 1.6 PCA + Harmony + UMAP -----------------------------------------------
message("\n[1.6] PCA, Harmony, clustering, UMAP...")

female_obj <- RunPCA(female_obj, npcs = N_PCS, verbose = FALSE)
ggsave(rf("umap", "elbow.pdf"),
       ElbowPlot(female_obj, ndims = N_PCS) + ggtitle("Females — PCA elbow"),
       width = 7, height = 5)

female_obj <- RunHarmony(female_obj,
                         group.by.vars  = "sample",
                         assay.use      = "SCT",
                         reduction      = "pca",
                         reduction.save = "harmony",
                         verbose        = FALSE)

female_obj <- FindNeighbors(female_obj, reduction = "harmony",
                            dims = 1:N_DIMS_USE, verbose = FALSE)
female_obj <- FindClusters(female_obj, resolution = CLUSTER_RES, verbose = FALSE)
female_obj <- RunUMAP(female_obj, reduction = "harmony",
                      dims = 1:N_DIMS_USE, verbose = FALSE)

female_obj <- PrepSCTFindMarkers(female_obj)

p1 <- DimPlot(female_obj, group.by = "sample")          + ggtitle("By sample")
p2 <- DimPlot(female_obj, group.by = "genotype")        + ggtitle("By genotype")
p3 <- DimPlot(female_obj, group.by = "replicate")       + ggtitle("By replicate")
p4 <- DimPlot(female_obj, group.by = "seurat_clusters",
              label = TRUE)                              + ggtitle("Clusters")

ggsave(rf("umap", "UMAP_sample.pdf"),    p1, width = 8, height = 6)
ggsave(rf("umap", "UMAP_genotype.pdf"),  p2, width = 8, height = 6)
ggsave(rf("umap", "UMAP_replicate.pdf"), p3, width = 8, height = 6)
ggsave(rf("umap", "UMAP_clusters.pdf"),  p4, width = 8, height = 6)
ggsave(rf("umap", "UMAP_all.pdf"), (p1 | p2) / (p3 | p4),
       width = 14, height = 12)

# ---- 1.7 SingleR annotation (SECONDARY / validation) ------------------------
message("\n[1.7] SingleR annotation against ImmGen (secondary label)...")

DefaultAssay(female_obj) <- "RNA"
female_obj <- JoinLayers(female_obj, assay = "RNA")
# Ensure an RNA 'data' layer exists for all downstream figure panels.
female_obj <- NormalizeData(female_obj, assay = "RNA", verbose = FALSE)

cts <- GetAssayData(female_obj, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(assays  = list(counts = cts),
                            colData = female_obj@meta.data)
logcounts(sce) <- log1p(counts(sce))

immgen <- celldex::ImmGenData()
sr <- SingleR(test = sce, ref = immgen,
              labels = immgen$label.main,
              BPPARAM = MulticoreParam(workers = 4))

female_obj$ImmGen.labels <- sr$labels

lab_tab     <- table(female_obj$ImmGen.labels)
keep_labels <- names(lab_tab)[lab_tab >= MIN_CELLS_PER_LABEL]
female_obj  <- subset(female_obj, subset = ImmGen.labels %in% keep_labels)

cat("\nCells per ImmGen label (post-filter):\n")
print(sort(table(female_obj$ImmGen.labels), decreasing = TRUE))

p_ct <- DimPlot(female_obj, group.by = "ImmGen.labels",
                label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Females — SingleR/ImmGen (secondary)")
ggsave(rf("umap", "UMAP_ImmGen_secondary.pdf"), p_ct, width = 10, height = 8)

# ---- 1.8 Composition (SingleR secondary) ------------------------------------
message("\n[1.8] Composition analysis (secondary label)...")

celltype_genotype <- table(female_obj$ImmGen.labels, female_obj$genotype)
ct_counts <- as.data.frame.matrix(celltype_genotype) %>%
  tibble::rownames_to_column("CellType") %>%
  mutate(Total = KO + WT) %>% arrange(desc(Total))
write.csv(ct_counts, rf("composition", "ImmGen_by_Genotype_counts.csv"),
          row.names = FALSE)

ct_prop <- as.data.frame.matrix(prop.table(celltype_genotype, 2) * 100) %>%
  tibble::rownames_to_column("CellType")
write.csv(ct_prop, rf("composition", "ImmGen_by_Genotype_proportions.csv"),
          row.names = FALSE)

# ---- 1.9 Save processed object ----------------------------------------------
message("\n[1.9] Saving female_obj_processed.rds ...")
saveRDS(female_obj, rf("female_obj_processed.rds"))
cat("Saved:", rf("female_obj_processed.rds"), "\n")

# =============================================================================
# SECTION 2 — CLUSTER QC DIAGNOSTICS
#   Flags doublet / junk clusters WITHOUT modifying the object. Writes tables +
#   figures to cluster_qc/. Review cluster_qc_summary.csv (flag column) and set
#   DROP_CLUSTERS in Section 0 to remove any flagged clusters in Section 4.
# =============================================================================
message("\n[2] Cluster QC diagnostics...")

safe_panel("cluster_qc", {
  cqc_dir <- rf("cluster_qc")
  RUN_CLUSTER_MARKERS <- FALSE   # TRUE = also FindAllMarkers (slower)
  COEXPR_MIN    <- 0
  PURITY_DROP   <- 0.60
  PURITY_REVIEW <- 0.75
  MITO_HI       <- 5
  NFEAT_LO_FRAC <- 0.6
  DBL_HI        <- 0.30
  COEXPR_HI     <- 0.10

  lineage_markers <- list(
    B     = c("Cd79a","Cd79b","Ms4a1","Cd19"),
    T     = c("Cd3e","Cd3d","Cd3g","Cd8a","Cd4"),
    NK    = c("Ncr1","Klrb1c","Gzma","Klrk1"),
    Myelo = c("Lyz2","Csf1r","Itgam","Fcgr3"),
    DC    = c("Flt3","Itgax","Xcr1"),
    Neut  = c("S100a8","S100a9","Retnlg"),
    Baso  = c("Mcpt8","Cpa3","Prss34"),
    Ery   = c("Hba-a1","Hbb-bs","Alas2"),
    Stem  = c("Cd34","Kit")
  )

  obj <- female_obj
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)

  md <- obj[[]]
  clu_col  <- pick_col(obj, c("seurat_clusters"))
  lab_col  <- pick_col(obj, c("ImmGen.labels","CellType","celltype"), fallback = NA)
  geno_col <- pick_col(obj, c("genotype","condition","group"))
  clu  <- factor(as.character(md[[clu_col]]))
  lab  <- if (!is.na(lab_col))  as.character(md[[lab_col]])  else rep(NA_character_, nrow(md))
  geno <- as.character(md[[geno_col]])
  nFeat <- md$nFeature_RNA; nCnt <- md$nCount_RNA
  pmt <- md$percent.mt; prb <- md$percent.ribo
  dbl <- if (!is.null(md$scDblFinder_score)) md$scDblFinder_score else NA

  present_mod <- lapply(lineage_markers, function(g) intersect(unique(g), rownames(obj)))
  present_mod <- present_mod[sapply(present_mod, length) > 0]
  obj <- AddModuleScore(obj, features = present_mod, name = "lin_", seed = 1)
  mod_cols <- grep("^lin_[0-9]+$", colnames(obj[[]]), value = TRUE)
  names(mod_cols) <- names(present_mod)[seq_along(mod_cols)]
  mod_df <- obj[[]][, mod_cols, drop = FALSE]; colnames(mod_df) <- names(mod_cols)

  cnt <- GetAssayData(obj, "RNA", layer = "counts")
  b_g <- intersect(c("Cd79a","Cd79b","Ms4a1"), rownames(cnt))
  t_g <- intersect(c("Cd3e","Cd3d","Cd3g"),     rownames(cnt))
  b_pos <- if (length(b_g)) Matrix::colSums(cnt[b_g, , drop = FALSE] > COEXPR_MIN) > 0 else rep(FALSE, ncol(cnt))
  t_pos <- if (length(t_g)) Matrix::colSums(cnt[t_g, , drop = FALSE] > COEXPR_MIN) > 0 else rep(FALSE, ncol(cnt))
  bt_co <- b_pos & t_pos

  gmed_nfeat <- median(nFeat)
  entropy <- function(p) { p <- p[p > 0]; -sum(p * log2(p)) }
  clusters <- levels(clu)
  summ <- lapply(clusters, function(cl) {
    idx <- clu == cl
    lt  <- sort(table(lab[idx]), decreasing = TRUE); frac <- lt / sum(lt)
    top3 <- paste(sprintf("%s(%.0f%%)", names(frac)[1:min(3,length(frac))],
                          100*frac[1:min(3,length(frac))]), collapse = "; ")
    data.frame(cluster = cl, n_cells = sum(idx),
      n_KO = sum(idx & geno == "KO"), n_WT = sum(idx & geno == "WT"),
      majority = if (length(frac)) names(frac)[1] else NA,
      purity = if (length(frac)) as.numeric(frac[1]) else NA,
      label_entropy = if (length(frac)) entropy(as.numeric(frac)) else NA,
      med_nFeature = median(nFeat[idx]), med_nCount = median(nCnt[idx]),
      med_pct_mt = if (all(is.na(pmt))) NA else median(pmt[idx]),
      med_pct_ribo = if (all(is.na(prb))) NA else median(prb[idx]),
      mean_dblscore = if (all(is.na(dbl))) NA else mean(dbl[idx]),
      pct_BT_coexpr = mean(bt_co[idx]), top3_labels = top3, stringsAsFactors = FALSE)
  })
  summ <- do.call(rbind, summ)
  lin_means <- t(sapply(clusters, function(cl) colMeans(mod_df[clu == cl, , drop = FALSE])))
  lin_means <- as.data.frame(lin_means); colnames(lin_means) <- paste0("mod_", colnames(mod_df))
  on_thresh <- apply(lin_means, 2, function(x) quantile(x, 0.66, na.rm = TRUE))
  n_lineages_on <- rowSums(sweep(lin_means, 2, on_thresh, ">"), na.rm = TRUE)
  summ <- cbind(summ, lin_means, n_lineages_on = n_lineages_on)

  qc_redflag <- with(summ,
    (!is.na(med_pct_mt) & med_pct_mt > MITO_HI) |
    (med_nFeature < NFEAT_LO_FRAC * gmed_nfeat) |
    (!is.na(mean_dblscore) & mean_dblscore > DBL_HI) |
    (pct_BT_coexpr > COEXPR_HI) | (n_lineages_on >= 3))
  summ$flag <- with(summ, ifelse(purity < PURITY_DROP & qc_redflag, "drop",
                          ifelse(purity < PURITY_REVIEW | qc_redflag, "review", "keep")))
  summ$flag_reason <- with(summ, paste0(
    ifelse(purity < PURITY_REVIEW, sprintf("low-purity(%.2f) ", purity), ""),
    ifelse(!is.na(med_pct_mt) & med_pct_mt > MITO_HI, sprintf("high-mt(%.1f) ", med_pct_mt), ""),
    ifelse(med_nFeature < NFEAT_LO_FRAC * gmed_nfeat, sprintf("low-nFeat(%.0f) ", med_nFeature), ""),
    ifelse(!is.na(mean_dblscore) & mean_dblscore > DBL_HI, sprintf("dblscore(%.2f) ", mean_dblscore), ""),
    ifelse(pct_BT_coexpr > COEXPR_HI, sprintf("BT-coexpr(%.0f%%) ", 100*pct_BT_coexpr), ""),
    ifelse(n_lineages_on >= 3, sprintf("multi-lineage(%d) ", n_lineages_on), "")))
  summ <- summ[order(summ$flag != "drop", summ$flag != "review", -summ$n_cells), ]
  write.csv(summ, file.path(cqc_dir, "cluster_qc_summary.csv"), row.names = FALSE)
  message("  cluster flags:")
  print(summ[, c("cluster","n_cells","majority","purity","flag","flag_reason")], row.names = FALSE)

  Idents(obj) <- clu
  qc_long <- data.frame(cluster = clu, nFeature = nFeat, nCount = nCnt, pct_mt = pmt, dblscore = dbl)
  qc_m <- tidyr::pivot_longer(qc_long, -cluster, names_to = "metric", values_to = "value")
  qc_m <- qc_m[is.finite(qc_m$value), ]
  pv <- ggplot(qc_m, aes(cluster, value, fill = cluster)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 11) + theme(legend.position = "none") +
    labs(x = "cluster (res 0.2)", y = NULL, title = "Per-cluster QC")
  save_plot(pv, file.path(cqc_dir, "cluster_qc_violins"), 8, 9)

  lin_long <- data.frame(cluster = rep(rownames(lin_means), ncol(lin_means)),
                         lineage = rep(sub("^mod_","",colnames(lin_means)), each = nrow(lin_means)),
                         score   = as.vector(as.matrix(lin_means)))
  pl <- ggplot(lin_long, aes(lineage, cluster, fill = score)) +
    geom_tile(color = "grey85") +
    scale_fill_gradient2(low = "#4575B4", mid = "grey95", high = "#B2182B", midpoint = 0) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Lineage module score by cluster (2+ hot lineages = doublet suspect)",
         x = NULL, y = "cluster")
  save_plot(pl, file.path(cqc_dir, "cluster_lineage_heatmap"), 7, 6)

  summ$cluster <- factor(summ$cluster, levels = summ$cluster[order(summ$purity)])
  pp <- ggplot(summ, aes(cluster, purity, fill = flag)) +
    geom_col() + coord_flip() +
    geom_text(aes(label = majority), hjust = -0.05, size = 3) +
    scale_fill_manual(values = c(keep="#4DAF4A", review="#FF7F00", drop="#E41A1C")) +
    ylim(0, 1.15) + theme_bw(base_size = 11) +
    labs(title = "SingleR purity per cluster (label = majority type)", y = "purity", x = NULL)
  save_plot(pp, file.path(cqc_dir, "cluster_purity_bar"), 7, 5)

  co <- data.frame(cluster = summ$cluster, pct = 100*summ$pct_BT_coexpr, flag = summ$flag)
  pc <- ggplot(co, aes(cluster, pct, fill = flag)) + geom_col() + coord_flip() +
    scale_fill_manual(values = c(keep="#4DAF4A", review="#FF7F00", drop="#E41A1C")) +
    geom_hline(yintercept = 100*COEXPR_HI, linetype = 2, color = "grey40") +
    theme_bw(base_size = 11) +
    labs(title = "Heterotypic doublet signal: % cells co-expressing B & T markers",
         y = "% B+T co-expressing", x = NULL)
  save_plot(pc, file.path(cqc_dir, "cluster_BT_coexpr"), 7, 5)

  if (RUN_CLUSTER_MARKERS) {
    mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
    top <- mk %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 15)
    write.csv(top, file.path(cqc_dir, "cluster_top_markers.csv"), row.names = FALSE)
  }
})


# =============================================================================
# SECTION 3 — MARKER-SCORE ANNOTATION
#   AddModuleScore per canonical type -> mean per cluster -> column z-score ->
#   argmax = marker_label. Writes cluster_assignment.csv (consumed by Section 4)
#   + comparison tables/figures. Does NOT modify the object.
# =============================================================================
message("\n[3] Marker-score cluster annotation...")

MARKER_ASSIGN_MODE <- "zscore"   # "zscore" (recommended) or "raw"
MARKER_MIN_MARGIN  <- 0.15

safe_panel("marker_annotation", {
  ma_dir <- rf("marker_annotation")
  obj <- female_obj
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)

  md <- obj[[]]
  clu_col  <- pick_col(obj, c("seurat_clusters"))
  lab_col  <- pick_col(obj, c("ImmGen.labels","CellType","celltype"), fallback = NA)
  geno_col <- pick_col(obj, c("genotype","condition","group"))
  clu  <- factor(as.character(md[[clu_col]]))
  immg <- if (!is.na(lab_col)) as.character(md[[lab_col]]) else rep(NA_character_, nrow(md))
  geno <- as.character(md[[geno_col]])
  clusters <- levels(clu)

  feats  <- rownames(obj); lookup <- setNames(feats, tolower(feats))
  present <- lapply(marker_list, function(gs) unname(lookup[tolower(gs)]))
  present <- lapply(present, function(x) x[!is.na(x)])
  present <- present[sapply(present, length) > 0]
  types   <- names(present)

  obj <- AddModuleScore(obj, features = present, name = "mk_", seed = 1)
  mk_cols <- grep("^mk_[0-9]+$", colnames(obj[[]]), value = TRUE)
  score_cell <- obj[[]][, mk_cols, drop = FALSE]; colnames(score_cell) <- types
  clu_score  <- t(sapply(clusters, function(cl) colMeans(score_cell[clu == cl, , drop = FALSE])))
  clu_score  <- as.data.frame(clu_score); colnames(clu_score) <- types
  write.csv(cbind(cluster = clusters, round(clu_score, 4)),
            file.path(ma_dir, "cluster_x_celltype_score.csv"), row.names = FALSE)

  zmat <- scale(as.matrix(clu_score))
  zmat_plot <- zmat; zmat_plot[!is.finite(zmat_plot)] <- NA
  zmat[!is.finite(zmat)] <- -Inf
  assign_mat <- if (MARKER_ASSIGN_MODE == "zscore") zmat else as.matrix(clu_score)

  best   <- types[apply(assign_mat, 1, which.max)]
  ord    <- t(apply(assign_mat, 1, function(r) order(r, decreasing = TRUE)))
  second <- types[ord[, 2]]
  best_s <- apply(assign_mat, 1, function(r) sort(r, decreasing = TRUE)[1])
  scnd_s <- apply(assign_mat, 1, function(r) sort(r, decreasing = TRUE)[2])
  margin <- best_s - scnd_s
  marker_lab <- setNames(best, clusters)
  ambiguous  <- margin < MARKER_MIN_MARGIN

  override_csv <- file.path(ma_dir, "cluster_manual_labels.csv")
  if (file.exists(override_csv)) {
    ov <- read.csv(override_csv, stringsAsFactors = FALSE); ov$cluster <- as.character(ov$cluster)
    hit <- ov$cluster %in% names(marker_lab); marker_lab[ov$cluster[hit]] <- ov$label[hit]
    message("  applied ", sum(hit), " manual override(s)")
  }

  immg_major  <- sapply(clusters, function(cl) { t <- sort(table(immg[clu == cl]), decreasing = TRUE)
    if (length(t)) names(t)[1] else NA })
  immg_purity <- sapply(clusters, function(cl) { t <- table(immg[clu == cl]); if (sum(t)) max(t)/sum(t) else NA })
  cmp <- data.frame(cluster = clusters, n_cells = as.integer(table(clu)[clusters]),
    immgen_majority = immg_major, immgen_purity = round(immg_purity, 3),
    marker_label = marker_lab[clusters], marker_2nd = setNames(second, clusters)[clusters],
    marker_margin = round(margin, 3), ambiguous = ambiguous,
    agrees = tolower(gsub("[^a-z]", "", tolower(immg_major))) ==
             tolower(gsub("[^a-z]", "", tolower(marker_lab[clusters]))),
    stringsAsFactors = FALSE)
  cmp <- cmp[order(cmp$agrees, -cmp$n_cells), ]
  write.csv(cmp, file.path(ma_dir, "cluster_assignment.csv"), row.names = FALSE)

  marker_cell <- marker_lab[as.character(clu)]
  pct <- function(x) { t <- table(x); round(100 * as.numeric(t) / sum(t), 2) |> setNames(names(t)) }
  all_types <- sort(unique(c(immg, unname(marker_cell))))
  fill0 <- function(v) { out <- setNames(rep(0, length(all_types)), all_types); out[names(v)] <- v; out }
  comp_immg <- fill0(pct(immg)); comp_mark <- fill0(pct(marker_cell))
  comp_overall <- data.frame(cell_type = all_types,
    pct_ImmGen = comp_immg[all_types], pct_marker = comp_mark[all_types],
    delta_marker_minus_ImmGen = round(comp_mark[all_types] - comp_immg[all_types], 2), row.names = NULL)
  comp_overall <- comp_overall[order(-comp_overall$pct_ImmGen), ]
  write.csv(comp_overall, file.path(ma_dir, "composition_overall.csv"), row.names = FALSE)

  have_geno <- !all(is.na(geno)) && length(unique(na.omit(geno))) >= 2
  comp_g <- NULL
  if (have_geno) {
    gl <- unique(na.omit(geno))
    comp_g <- do.call(rbind, lapply(gl, function(g) {
      idx <- geno == g & !is.na(geno)
      ci <- fill0(pct(immg[idx])); cm <- fill0(pct(marker_cell[idx]))
      rbind(data.frame(genotype = g, method = "ImmGen", cell_type = all_types, pct = ci[all_types]),
            data.frame(genotype = g, method = "marker", cell_type = all_types, pct = cm[all_types]))
    }))
    write.csv(comp_g, file.path(ma_dir, "composition_by_genotype.csv"), row.names = FALSE)
    if (all(c("WT","KO") %in% gl)) {
      wide <- comp_g |> tidyr::pivot_wider(names_from = genotype, values_from = pct)
      wide$KO_minus_WT <- round(wide$KO - wide$WT, 2)
      wide <- wide[order(wide$method, -abs(wide$KO_minus_WT)), ]
      write.csv(wide, file.path(ma_dir, "composition_KO_minus_WT.csv"), row.names = FALSE)
    }
  }

  conf <- as.data.frame.matrix(table(ImmGen = immg, marker = marker_cell))
  write.csv(cbind(ImmGen = rownames(conf), conf),
            file.path(ma_dir, "confusion_immgen_vs_marker.csv"), row.names = FALSE)

  # figures
  zl <- data.frame(cluster = rep(clusters, length(types)),
                   type    = rep(types, each = length(clusters)),
                   z       = as.vector(as.matrix(zmat_plot)))
  zl$called <- mapply(function(cl, ty) identical(unname(marker_lab[cl]), ty), zl$cluster, zl$type)
  ph <- ggplot(zl, aes(type, factor(cluster, levels = rev(clusters)), fill = z)) +
    geom_tile(color = "grey88") +
    geom_point(data = subset(zl, called), shape = 8, size = 2, color = "black") +
    scale_fill_gradient2(low = "#4575B4", mid = "grey95", high = "#B2182B", midpoint = 0,
                         name = "z(score)", na.value = "grey85") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Marker module score per cluster (* = assigned type)", x = NULL, y = "cluster")
  save_plot(ph, file.path(ma_dir, "cluster_score_heatmap"), 8, 6)

  cl2 <- comp_overall |> tidyr::pivot_longer(c(pct_ImmGen, pct_marker),
                                             names_to = "method", values_to = "pct")
  cl2$method <- sub("pct_", "", cl2$method)
  cl2$cell_type <- factor(cl2$cell_type, levels = rev(comp_overall$cell_type))
  pcb <- ggplot(cl2, aes(cell_type, pct, fill = method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) + coord_flip() +
    scale_fill_manual(values = c(ImmGen = "#7570B3", marker = "#1B9E77")) +
    theme_bw(base_size = 12) +
    labs(title = "Cell composition: ImmGen/SingleR vs canonical-marker annotation",
         x = NULL, y = "% of cells")
  save_plot(pcb, file.path(ma_dir, "composition_method_bars"), 8, 6)

  if (have_geno && !is.null(comp_g)) {
    cg <- comp_g; cg$cell_type <- factor(cg$cell_type, levels = rev(comp_overall$cell_type))
    pg <- ggplot(cg, aes(cell_type, pct, fill = genotype)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.72) +
      coord_flip() + facet_wrap(~method) +
      scale_fill_manual(values = c(WT = "#3182BD", KO = "#DE2D26")) +
      theme_bw(base_size = 12) +
      labs(title = "Composition by genotype, per annotation method", x = NULL, y = "% of cells")
    save_plot(pg, file.path(ma_dir, "composition_genotype_bars"), 9, 6)
    if (all(c("WT","KO") %in% unique(na.omit(geno)))) {
      wide2 <- comp_g |> tidyr::pivot_wider(names_from = genotype, values_from = pct) |>
        dplyr::mutate(KO_minus_WT = KO - WT)
      wide2$cell_type <- factor(wide2$cell_type, levels = rev(comp_overall$cell_type))
      pd <- ggplot(wide2, aes(cell_type, KO_minus_WT, fill = method)) +
        geom_col(position = position_dodge(width = 0.75), width = 0.7) +
        coord_flip() + geom_hline(yintercept = 0, color = "grey40") +
        scale_fill_manual(values = c(ImmGen = "#7570B3", marker = "#1B9E77")) +
        theme_bw(base_size = 12) +
        labs(title = "KO - WT composition change (percentage points)", x = NULL, y = "KO - WT (pp)")
      save_plot(pd, file.path(ma_dir, "composition_delta_bars"), 8, 6)
    }
  }

  cf <- as.matrix(conf); rn <- sweep(cf, 1, pmax(rowSums(cf), 1), "/")
  cl_long <- data.frame(ImmGen = rep(rownames(rn), ncol(rn)),
                        marker = rep(colnames(rn), each = nrow(rn)), frac = as.vector(rn))
  pcf <- ggplot(cl_long, aes(marker, ImmGen, fill = frac)) +
    geom_tile(color = "grey90") +
    scale_fill_gradient(low = "white", high = "#08519C", name = "row frac", limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Where each ImmGen label lands in marker annotation (row-normalized)",
         x = "marker label", y = "ImmGen label")
  save_plot(pcf, file.path(ma_dir, "confusion_heatmap"), 8, 7)
})


# =============================================================================
# SECTION 4 — FINALIZE ANNOTATION (commit celltype_marker as PRIMARY)
#   Verifies marker labels vs data-driven FindAllMarkers, commits
#   celltype_marker, drops DROP_CLUSTERS, keeps ImmGen.labels as secondary,
#   saves female_obj_annotated.rds. Updates in-memory female_obj.
# =============================================================================
message("\n[4] Finalize annotation -> celltype_marker PRIMARY...")

FINALIZE_TOP_N       <- 15
FINALIZE_MARKER_PADJ <- 0.05

{
  ma_dir  <- rf("marker_annotation")
  obj <- female_obj
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)

  md <- obj[[]]
  clu_col  <- pick_col(obj, c("seurat_clusters"))
  lab_col  <- pick_col(obj, c("ImmGen.labels","CellType","celltype"), fallback = NA)
  geno_col <- pick_col(obj, c("genotype","condition","group"))
  clu  <- factor(as.character(md[[clu_col]]))
  geno <- as.character(md[[geno_col]])
  clusters <- levels(clu)

  assign_csv <- file.path(ma_dir, "cluster_assignment.csv")
  if (!file.exists(assign_csv))
    stop("Missing ", assign_csv, " -- Section 3 must run first.")
  asg <- read.csv(assign_csv, stringsAsFactors = FALSE); asg$cluster <- as.character(asg$cluster)
  lab_map <- setNames(asg$marker_label, asg$cluster)
  override_csv <- file.path(ma_dir, "cluster_manual_labels.csv")
  if (file.exists(override_csv)) {
    ov <- read.csv(override_csv, stringsAsFactors = FALSE); ov$cluster <- as.character(ov$cluster)
    lab_map[ov$cluster] <- ov$label
  }
  for (cl in names(MANUAL_LABELS)) lab_map[[cl]] <- MANUAL_LABELS[[cl]]
  missing_lab <- setdiff(clusters, names(lab_map))
  if (length(missing_lab)) lab_map[missing_lab] <- "Unassigned"

  lookup <- setNames(rownames(obj), tolower(rownames(obj)))
  canon  <- lapply(marker_list, function(g) unname(lookup[tolower(g)]))
  canon  <- lapply(canon, function(x) x[!is.na(x)])

  Idents(obj) <- clu
  mk <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
                       min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
  mk <- mk[mk$p_val_adj < FINALIZE_MARKER_PADJ, ]
  `%||%` <- function(a, b) if (is.null(a)) b else a
  ver <- lapply(clusters, function(cl) {
    sub <- mk[mk$cluster == cl, ]; sub <- sub[order(-sub$avg_log2FC), ]
    top <- head(sub$gene, FINALIZE_TOP_N); assigned <- unname(lab_map[[cl]])
    votes <- sapply(canon, function(cg) length(intersect(cg, top)))
    marker_vote <- if (max(votes) > 0) names(votes)[which.max(votes)] else NA
    assigned_hits <- intersect(canon[[assigned]] %||% character(0), top)
    data.frame(cluster = cl, assigned_label = assigned, marker_vote = marker_vote %||% NA,
      votes_assigned = if (!is.null(canon[[assigned]])) length(assigned_hits) else 0L,
      votes_top = as.integer(max(votes)),
      support = (length(assigned_hits) > 0) && (is.na(marker_vote) || marker_vote == assigned),
      assigned_canon_in_top = paste(assigned_hits, collapse = ", "),
      top_markers = paste(head(top, 10), collapse = ", "), stringsAsFactors = FALSE)
  })
  ver <- do.call(rbind, ver); ver <- ver[order(ver$support, ver$cluster), ]
  write.csv(ver, file.path(ma_dir, "annotation_verification.csv"), row.names = FALSE)
  bad <- ver[!ver$support, ]
  if (nrow(bad)) message("  ** ", nrow(bad),
    " cluster(s) not supported by top markers -- review MANUAL_LABELS/DROP_CLUSTERS: ",
    paste(bad$cluster, "->", bad$assigned_label, collapse = " | "))

  obj$celltype_marker <- factor(unname(lab_map[as.character(clu)]))
  if (length(DROP_CLUSTERS)) {
    keep_cells <- !(as.character(clu) %in% as.character(DROP_CLUSTERS))
    message("  dropping clusters ", paste(DROP_CLUSTERS, collapse = ", "),
            " (", sum(!keep_cells), " cells removed)")
    obj <- subset(obj, cells = colnames(obj)[keep_cells])
    obj$celltype_marker <- droplevels(obj$celltype_marker)
  }
  Idents(obj) <- "celltype_marker"

  final_clu <- factor(as.character(obj[[]][[clu_col]]))
  write.csv(data.frame(cluster = levels(final_clu),
                       final_label = unname(lab_map[levels(final_clu)])),
            file.path(ma_dir, "final_labels.csv"), row.names = FALSE)
  ct   <- as.character(obj$celltype_marker); gvec <- as.character(obj[[]][[geno_col]])
  pctf <- function(x) { t <- table(x); round(100 * as.numeric(t) / sum(t), 2) |> setNames(names(t)) }
  co   <- pctf(ct)
  write.csv(data.frame(cell_type = names(co), pct = as.numeric(co),
                       n = as.integer(table(ct)[names(co)])),
            file.path(ma_dir, "final_composition_overall.csv"), row.names = FALSE)
  if (!all(is.na(gvec)) && length(unique(na.omit(gvec))) >= 2) {
    gl <- unique(na.omit(gvec))
    cg <- do.call(rbind, lapply(gl, function(g) { p <- pctf(ct[gvec == g & !is.na(gvec)])
      data.frame(genotype = g, cell_type = names(p), pct = as.numeric(p)) }))
    write.csv(cg, file.path(ma_dir, "final_composition_by_genotype.csv"), row.names = FALSE)
  }

  feats_present <- unique(unlist(lapply(canon[names(canon) %in% ct], identity), use.names = FALSE))
  if (length(feats_present)) {
    dp <- tryCatch(DotPlot(obj, features = feats_present, cluster.idents = FALSE, dot.scale = 6) +
      scale_color_gradient2(low = "#4575B4", mid = "grey90", high = "#B2182B", midpoint = 0) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
      labs(title = "Canonical markers on final marker-based labels", x = NULL, y = NULL),
      error = function(e) NULL)
    if (!is.null(dp)) {
      wd <- max(8, 0.4*length(feats_present)+3); hd <- max(4.5, 0.4*length(unique(ct))+2)
      save_plot(dp, file.path(ma_dir, "annotation_check_dotplot"), wd, hd)
    }
  }

  saveRDS(obj, rf("female_obj_annotated.rds"))
  female_obj <- obj              # promote annotated object for all downstream work
  message("  saved female_obj_annotated.rds ; Idents = celltype_marker")
}

# Cell types actually present, ordered by abundance (used across all downstream)
CT_LEVELS <- names(sort(table(as.character(female_obj$celltype_marker)), decreasing = TRUE))
CT_PALETTE <- build_celltype_palette(CT_LEVELS)


# =============================================================================
# SECTION 5 — MARKER DOT PLOT + PRIMARY-LABEL UMAP
#   Block-diagonal canonical-marker dot plot on celltype_marker + a labelled
#   UMAP of the primary annotation. Written to panelE/ and umap/.
# =============================================================================
message("\n[5] Marker dot plot + primary-label UMAP...")

safe_panel("marker_dotplot", {
  pe_dir <- rf("panelE")
  na <- ensure_norm_assay(female_obj); obj <- na$obj; DefaultAssay(obj) <- na$assay

  feats  <- rownames(obj); lookup <- setNames(feats, tolower(feats))
  resolved <- lapply(marker_list_dot, function(gs) unname(lookup[tolower(gs)]))
  present  <- lapply(resolved, function(x) x[!is.na(x)])
  ct_in_data <- unique(as.character(obj[[CELLTYPE_COL]][, 1]))
  keep_ct    <- names(present)[sapply(present, length) > 0 & names(present) %in% ct_in_data]
  present    <- present[keep_ct]
  features   <- unique(unlist(present, use.names = FALSE))
  write.csv(data.frame(cell_type = rep(names(present), lengths(present)),
                       marker = unlist(present, use.names = FALSE)),
            file.path(pe_dir, "marker_dotplot_genes_used.csv"), row.names = FALSE)

  obj$.ct <- factor(as.character(obj[[CELLTYPE_COL]][, 1]), levels = rev(keep_ct))
  Idents(obj) <- ".ct"
  obj <- subset(obj, idents = keep_ct)

  p <- DotPlot(obj, features = features, cluster.idents = FALSE,
               dot.scale = 7, col.min = -1.5, col.max = 2.5) +
    scale_color_gradient2(low = "#4575B4", mid = "grey90", high = "#B2182B",
                          midpoint = 0, name = "Scaled\nexpression") +
    guides(size = guide_legend(title = "% expressing")) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          legend.title = element_text(size = 12))
  w <- max(7, 0.42 * length(features) + 3); h <- max(4.5, 0.42 * length(keep_ct) + 2)
  save_plot(p, file.path(pe_dir, "marker_dotplot"), w, h)
})

safe_panel("primary_label_umap", {
  pal <- CT_PALETTE
  p_ct <- DimPlot(female_obj, group.by = CELLTYPE_COL, cols = pal,
                  label = TRUE, repel = TRUE) +
    ggtitle("Females — cell types (marker-based, primary)")
  save_plot(p_ct, rf("umap", "UMAP_celltype_marker"), 10, 8)
  p_split <- DimPlot(female_obj, group.by = CELLTYPE_COL, split.by = "genotype",
                     cols = pal, label = TRUE, repel = TRUE) +
    ggtitle("Females — cell types by genotype")
  save_plot(p_split, rf("umap", "UMAP_celltype_marker_by_genotype"), 16, 6)
})

# =============================================================================
# SECTION 6 — DIFFERENTIAL EXPRESSION (KO vs WT)
#   6A. Pseudobulk DESeq2 (mice as replicates) : global + per celltype_marker
#   6B. Seurat FindMarkers (Wilcoxon, cells)    : global + per celltype_marker
#   Results held in memory for volcanoes (7), enrichment (9), Excel (10).
# =============================================================================
message("\n[6] Differential expression (KO vs WT)...")

DefaultAssay(female_obj) <- "RNA"
female_obj <- JoinLayers(female_obj, assay = "RNA")
# guarantee an RNA 'data' layer for Seurat FindMarkers / figures
if (!("data" %in% tryCatch(SeuratObject::Layers(female_obj[["RNA"]]),
                           error = function(e) character())))
  female_obj <- NormalizeData(female_obj, assay = "RNA", verbose = FALSE)

counts_all <- GetAssayData(female_obj, assay = "RNA", layer = "counts")
meta_all   <- female_obj@meta.data
sample_geno <- meta_all %>% dplyr::distinct(sample, genotype) %>% dplyr::arrange(sample)

# ---- 6A. Pseudobulk DESeq2 engine ------------------------------------------
detection_rates <- function(counts_sub, meta_sub) {
  ko <- rownames(meta_sub)[meta_sub$genotype == "KO"]
  wt <- rownames(meta_sub)[meta_sub$genotype == "WT"]
  pct1 <- if (length(ko)) Matrix::rowMeans(counts_sub[, ko, drop = FALSE] > 0) else NA_real_
  pct2 <- if (length(wt)) Matrix::rowMeans(counts_sub[, wt, drop = FALSE] > 0) else NA_real_
  data.frame(gene = rownames(counts_sub),
             pct.1 = round(as.numeric(pct1), 3),
             pct.2 = round(as.numeric(pct2), 3), stringsAsFactors = FALSE)
}

make_pseudobulk <- function(counts_sub, meta_sub) {
  samp <- factor(meta_sub$sample)
  ind <- Matrix::sparse.model.matrix(~ 0 + samp); colnames(ind) <- levels(samp)
  as.matrix(counts_sub %*% ind)
}

run_pseudobulk_deseq <- function(counts_sub, meta_sub, label) {
  message(sprintf("  [pseudobulk] %s", label))
  cps <- table(meta_sub$sample)
  usable <- names(cps)[cps >= MIN_CELLS_PER_SAMPLE]
  meta_sub <- meta_sub[meta_sub$sample %in% usable, , drop = FALSE]
  if (nrow(meta_sub) == 0) { message("    skip: no usable samples"); return(NULL) }
  geno_tab <- sample_geno %>% dplyr::filter(sample %in% usable) %>% dplyr::count(genotype)
  if (!all(c("KO","WT") %in% geno_tab$genotype) || any(geno_tab$n < MIN_REPS_PER_GENOTYPE)) {
    message("    skip: insufficient replicates per genotype"); return(NULL)
  }
  counts_sub2 <- counts_sub[, rownames(meta_sub), drop = FALSE]
  pb <- make_pseudobulk(counts_sub2, meta_sub)[, usable, drop = FALSE]
  coldata <- sample_geno[match(colnames(pb), sample_geno$sample), ]
  coldata$genotype <- factor(coldata$genotype, levels = c("WT","KO"))
  rownames(coldata) <- coldata$sample
  keep <- rowSums(pb >= MIN_COUNT) >= MIN_SAMPLES
  pb <- pb[keep, , drop = FALSE]
  if (nrow(pb) < 10) { message("    skip: <10 genes pass count filter"); return(NULL) }
  res_df <- tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = pb, colData = coldata, design = ~ genotype)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("genotype","KO","WT"))
    res <- tryCatch(lfcShrink(dds, contrast = c("genotype","KO","WT"), type = "ashr", res = res),
                    error = function(e) { message("    (no shrinkage: ", e$message, ")"); res })
    as.data.frame(res)
  }, error = function(e) { message("    DESeq2 failed: ", e$message); NULL })
  if (is.null(res_df)) return(NULL)
  res_df$gene <- rownames(res_df)
  det <- detection_rates(counts_sub2, meta_sub)
  res_df %>%
    dplyr::left_join(det, by = "gene") %>%
    dplyr::transmute(gene, log2FC = log2FoldChange, padj = padj, p_val = pvalue,
                     pct.1, pct.2, avg_log2FC = log2FoldChange, p_val_adj = padj, baseMean) %>%
    dplyr::filter(!is.na(log2FC)) %>%
    dplyr::arrange(dplyr::desc(log2FC))
}

# Shared ORA helper (used by Section 9A). Kept here so it is defined before use.
run_enrichment_analysis <- function(gene_list, gene_universe, analysis_name) {
  message(sprintf("  enrichment: %s (%d genes)", analysis_name, length(gene_list)))
  if (length(gene_list) == 0) return(list())
  gene_entrez <- tryCatch(suppressMessages(suppressWarnings(
    bitr(gene_list, "SYMBOL", "ENTREZID", org.Mm.eg.db, drop = TRUE))),
    error = function(e) data.frame(SYMBOL=character(0), ENTREZID=character(0)))
  universe_entrez <- tryCatch(suppressMessages(suppressWarnings(
    bitr(gene_universe, "SYMBOL", "ENTREZID", org.Mm.eg.db, drop = TRUE))),
    error = function(e) data.frame(SYMBOL=character(0), ENTREZID=character(0)))
  if (nrow(gene_entrez) == 0 || nrow(universe_entrez) == 0) return(list())
  results <- list()
  for (ont in c("BP","MF","CC")) {
    tryCatch({
      eg <- enrichGO(gene_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = ont,
                     universe = universe_entrez$ENTREZID, pAdjustMethod = "BH",
                     pvalueCutoff = ENR_PVAL, qvalueCutoff = ENR_PVAL, readable = TRUE)
      if (!is.null(eg) && nrow(eg@result) > 0) results[[paste0("GO_", ont)]] <- as.data.frame(eg)
    }, error = function(e) message("    GO ", ont, " failed: ", e$message))
  }
  tryCatch({
    kegg <- enrichKEGG(gene_entrez$ENTREZID, organism = "mmu",
                       universe = universe_entrez$ENTREZID, pAdjustMethod = "BH",
                       pvalueCutoff = ENR_PVAL, qvalueCutoff = ENR_PVAL)
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      kegg <- setReadable(kegg, org.Mm.eg.db, keyType = "ENTREZID")
      results$KEGG <- as.data.frame(kegg)
    }
  }, error = function(e) message("    KEGG failed: ", e$message))
  if (HAVE_REACTOME) tryCatch({
    rc <- ReactomePA::enrichPathway(gene_entrez$ENTREZID, organism = "mouse",
                        universe = universe_entrez$ENTREZID, pAdjustMethod = "BH",
                        pvalueCutoff = ENR_PVAL, qvalueCutoff = ENR_PVAL, readable = TRUE)
    if (!is.null(rc) && nrow(rc@result) > 0) results$Reactome <- as.data.frame(rc)
  }, error = function(e) message("    Reactome failed: ", e$message))
  tryCatch({
    h <- msigdbr(species = "Mus musculus", category = "H")
    t2g <- h %>% dplyr::select(gs_name, entrez_gene)
    hm <- enricher(gene_entrez$ENTREZID, universe = universe_entrez$ENTREZID,
                   TERM2GENE = t2g, pAdjustMethod = "BH",
                   pvalueCutoff = ENR_PVAL, qvalueCutoff = ENR_PVAL)
    if (!is.null(hm) && nrow(hm@result) > 0) {
      hr <- as.data.frame(hm)
      hr$geneID <- sapply(strsplit(hr$geneID, "/"), function(ids) {
        s <- gene_entrez$SYMBOL[match(ids, gene_entrez$ENTREZID)]
        paste(s[!is.na(s)], collapse = "/") })
      results$Hallmark <- hr
    }
  }, error = function(e) message("    Hallmark failed: ", e$message))
  results
}

message("\n[6A] Pseudobulk DESeq2: global...")
de_global <- run_pseudobulk_deseq(counts_all, meta_all, "Global")
if (is.null(de_global)) stop("Global pseudobulk DE returned no results.")
deg_global_filtered <- de_global %>%
  dplyr::filter(!is.na(padj), abs(log2FC) >= LFC_CUT, padj <= PADJ_CUT) %>%
  dplyr::arrange(dplyr::desc(log2FC))
gene_universe <- de_global$gene

message("\n[6A] Pseudobulk DESeq2: per cell type (celltype_marker)...")
celltypes <- CT_LEVELS
de_by_celltype <- list()
for (ct in celltypes) {
  ct_cells <- rownames(meta_all)[as.character(meta_all[[CELLTYPE_COL]]) == ct]
  if (!length(ct_cells)) next
  res_ct <- run_pseudobulk_deseq(counts_all[, ct_cells, drop = FALSE],
                                 meta_all[ct_cells, , drop = FALSE], ct)
  if (is.null(res_ct) || nrow(res_ct) == 0) next
  de_by_celltype[[ct]] <- res_ct %>% dplyr::mutate(celltype = ct) %>%
    dplyr::select(gene, celltype, log2FC, padj, p_val, pct.1, pct.2, avg_log2FC, p_val_adj, baseMean)
}
deg_by_celltype_filtered <- lapply(de_by_celltype, function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% dplyr::filter(!is.na(padj), abs(log2FC) >= LFC_CUT, padj <= PADJ_CUT) %>%
    dplyr::arrange(dplyr::desc(log2FC))
})
deg_by_celltype_filtered <- deg_by_celltype_filtered[!sapply(deg_by_celltype_filtered, is.null)]
deg_by_celltype_filtered <- deg_by_celltype_filtered[sapply(deg_by_celltype_filtered, nrow) > 0]
deg_ct_combined <- if (length(deg_by_celltype_filtered) > 0) {
  dplyr::bind_rows(deg_by_celltype_filtered) %>% dplyr::arrange(dplyr::desc(log2FC))
} else {
  data.frame(gene=character(0), celltype=character(0), log2FC=numeric(0), padj=numeric(0),
             p_val=numeric(0), pct.1=numeric(0), pct.2=numeric(0),
             avg_log2FC=numeric(0), p_val_adj=numeric(0), baseMean=numeric(0))
}
# annotate global tables with DE_in_celltypes
gene_to_DE_celltypes <- deg_ct_combined %>%
  dplyr::filter(!is.na(padj), abs(log2FC) >= LFC_CUT, padj <= PADJ_CUT) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(DE_in_celltypes = paste(sort(unique(celltype)), collapse = "; ")) %>%
  dplyr::ungroup()
de_global <- de_global %>% dplyr::left_join(gene_to_DE_celltypes, by = "gene") %>%
  dplyr::mutate(DE_in_celltypes = ifelse(is.na(DE_in_celltypes), "None", DE_in_celltypes)) %>%
  dplyr::select(gene, DE_in_celltypes, everything())
deg_global_filtered <- deg_global_filtered %>% dplyr::left_join(gene_to_DE_celltypes, by = "gene") %>%
  dplyr::mutate(DE_in_celltypes = ifelse(is.na(DE_in_celltypes), "None", DE_in_celltypes)) %>%
  dplyr::select(gene, DE_in_celltypes, everything())
# write DESeq2 tables
write.csv(de_global, rf("volcano", "deseq2", "DESeq2_Global_all.csv"), row.names = FALSE)
write.csv(deg_global_filtered, rf("volcano", "deseq2", "DESeq2_Global_filtered.csv"), row.names = FALSE)
for (ct in names(de_by_celltype))
  write.csv(de_by_celltype[[ct]],
            rf("volcano", "deseq2", paste0("DESeq2_", safe_sheet(ct), "_all.csv")), row.names = FALSE)

# ---- 6B. Seurat FindMarkers engine -----------------------------------------
SEURAT_MIN_CELLS <- 3
run_seurat_de_group <- function(cells, label) {
  message(sprintf("  [seurat] %s", label))
  o <- if (length(cells) == ncol(female_obj)) female_obj else subset(female_obj, cells = cells)
  DefaultAssay(o) <- "RNA"
  gvec <- as.character(o$genotype)
  Idents(o) <- factor(gvec, levels = c("WT", "KO"))
  tab <- table(as.character(Idents(o)))
  if (!all(c("WT","KO") %in% names(tab)) || any(tab[c("WT","KO")] < SEURAT_MIN_CELLS)) {
    message("    skip: too few cells per genotype"); return(NULL)
  }
  m <- tryCatch(FindMarkers(o, ident.1 = "KO", ident.2 = "WT",
                   test.use = "wilcox", logfc.threshold = 0, min.pct = 0, verbose = FALSE),
                error = function(e) { message("    FindMarkers failed: ", e$message); NULL })
  if (is.null(m) || !nrow(m)) return(NULL)
  m$gene <- rownames(m)
  p <- m$p_val; pos <- p[p > 0 & is.finite(p)]
  floor_p <- if (length(pos)) min(pos) * 0.1 else 1e-300
  p[!is.finite(p) | p == 0] <- floor_p
  stat <- sign(m$avg_log2FC) * pmin(-log10(p), 300) + m$avg_log2FC * 1e-6
  data.frame(gene = m$gene, log2FoldChange = m$avg_log2FC,
             pvalue = m$p_val, padj = m$p_val_adj, stat = stat,
             pct.1 = m$pct.1, pct.2 = m$pct.2,
             stringsAsFactors = FALSE)[order(m$p_val_adj), ]
}

message("\n[6B] Seurat FindMarkers: global + per cell type...")
seurat_de_global <- run_seurat_de_group(colnames(female_obj), "Global")
if (!is.null(seurat_de_global))
  write.csv(seurat_de_global, rf("volcano", "seurat", "Seurat_Global_all.csv"), row.names = FALSE)
seurat_de_by_celltype <- list()
for (ct in CT_LEVELS) {
  cells <- colnames(female_obj)[as.character(female_obj[[CELLTYPE_COL]][, 1]) == ct]
  res <- run_seurat_de_group(cells, ct)
  if (is.null(res) || !nrow(res)) next
  seurat_de_by_celltype[[ct]] <- res
  write.csv(res, rf("volcano", "seurat", paste0("Seurat_", safe_sheet(ct), "_all.csv")),
            row.names = FALSE)
}


# =============================================================================
# SECTION 7 — VOLCANO PLOTS (BOTH methods x global + each cell type)
#   DESeq2  -> volcano/deseq2/   ;  Seurat -> volcano/seurat/
# =============================================================================
message("\n[7] Volcano plots (DESeq2 + Seurat, global + per cell type)...")

# Generic volcano from a DE data frame. Column names supplied by caller.
make_volcano_plot <- function(df, gene_col, lfc_col, fdr_col, title, out_stem,
                              lfc_cut = LFC_CUT, fdr_cut = PADJ_CUT,
                              n_label = VOLC_N_LABEL, anchor = VOLC_ANCHOR) {
  if (is.null(df) || !nrow(df)) { message("  [skip volcano] ", title); return(invisible()) }
  d <- data.frame(gene = as.character(df[[gene_col]]),
                  lfc  = suppressWarnings(as.numeric(df[[lfc_col]])),
                  fdr  = suppressWarnings(as.numeric(df[[fdr_col]])))
  d <- d[is.finite(d$lfc) & is.finite(d$fdr) & d$fdr > 0, ]
  if (!nrow(d)) { message("  [skip volcano, no finite rows] ", title); return(invisible()) }
  d$neglog <- -log10(d$fdr)
  YCAP <- stats::quantile(d$neglog, 0.999, na.rm = TRUE)
  d$neglog_c <- pmin(d$neglog, YCAP)
  d$dir <- "NS"
  d$dir[d$fdr < fdr_cut & d$lfc >=  lfc_cut] <- "Up in KO"
  d$dir[d$fdr < fdr_cut & d$lfc <= -lfc_cut] <- "Down in KO"
  d$dir <- factor(d$dir, levels = c("Down in KO", "NS", "Up in KO"))
  dir_colors <- c("Down in KO" = COL_DOWN, "NS" = COL_NS, "Up in KO" = COL_UP)
  named <- d[!grepl("^ENSMUSG", d$gene) & d$dir != "NS", ]
  lab <- named %>% arrange(fdr) %>% slice_head(n = n_label)
  if (anchor %in% d$gene)
    lab <- bind_rows(lab, d[d$gene == anchor, ]) %>% distinct(gene, .keep_all = TRUE)
  n_up <- sum(d$dir == "Up in KO"); n_dn <- sum(d$dir == "Down in KO")
  p <- ggplot(d, aes(lfc, neglog_c, colour = dir)) +
    geom_point(size = 1.1, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey55") +
    scale_colour_manual(values = dir_colors, name = NULL) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 4, fontface = "italic",
      max.overlaps = Inf, box.padding = 0.4, segment.size = 0.3,
      min.segment.length = 0, colour = "black", seed = 1) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
             label = paste0("Up in KO: ", n_up), colour = COL_UP, size = 4.5) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
             label = paste0("Down in KO: ", n_dn), colour = COL_DOWN, size = 4.5) +
    labs(title = title, x = expression(log[2]~"fold change (KO / WT)"),
         y = expression(-log[10]~"(adjusted "*italic(p)*")")) +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          axis.text = element_text(colour = "black"), legend.position = "top")
  save_plot(p, out_stem, 7, 6.5)
}

# --- DESeq2 volcanoes ---
safe_panel("volcano_deseq2_global", {
  make_volcano_plot(de_global, "gene", "log2FC", "padj",
                    "Pseudobulk DESeq2 — Global (KO vs WT)",
                    rf("volcano", "deseq2", "volcano_DESeq2_Global"))
})
for (ct in names(de_by_celltype)) {
  local({
    ctt <- ct
    safe_panel(paste0("volcano_deseq2_", ctt), {
      make_volcano_plot(de_by_celltype[[ctt]], "gene", "log2FC", "padj",
                        paste0("Pseudobulk DESeq2 — ", ctt, " (KO vs WT)"),
                        rf("volcano", "deseq2", paste0("volcano_DESeq2_", safe_sheet(ctt))))
    })
  })
}

# --- Seurat volcanoes ---
safe_panel("volcano_seurat_global", {
  make_volcano_plot(seurat_de_global, "gene", "log2FoldChange", "padj",
                    "Seurat FindMarkers — Global (KO vs WT)",
                    rf("volcano", "seurat", "volcano_Seurat_Global"))
})
for (ct in names(seurat_de_by_celltype)) {
  local({
    ctt <- ct
    safe_panel(paste0("volcano_seurat_", ctt), {
      make_volcano_plot(seurat_de_by_celltype[[ctt]], "gene", "log2FoldChange", "padj",
                        paste0("Seurat FindMarkers — ", ctt, " (KO vs WT)"),
                        rf("volcano", "seurat", paste0("volcano_Seurat_", safe_sheet(ctt))))
    })
  })
}

# =============================================================================
# SECTION 8 — PUBLICATION FIGURE PANELS
# -----------------------------------------------------------------------------
# Every figure-panel script from the project is folded in here as a self-
# contained safe_panel() block. All panels operate on the in-memory,
# marker-annotated `female_obj` (celltype_marker = PRIMARY labels), route
# outputs through rf() to the existing results/females/ subfolders, and use the
# shared palettes/params from SECTION 0. One failing panel never aborts the run.
# Edit THIS section to regenerate any main/supplementary figure.
# =============================================================================

# ---- 8.0 Shared setup for figure panels -------------------------------------
{
  .norm       <- ensure_norm_assay(female_obj)
  female_obj  <- .norm$obj
  ASSAY_FIG   <- .norm$assay                 # normalized assay for expression plots
  CT_COL      <- CELLTYPE_COL                # "celltype_marker" (BINDING #1)
  GENO_COL    <- pick_col(female_obj, c("genotype", "Genotype", "condition", "group"), "genotype")
  SAMPLE_COL  <- pick_col(female_obj, c("sample", "Sample", "orig.ident", "mouse", "library"), "sample")
  message("SECTION 8 setup: assay=", ASSAY_FIG, " celltype=", CT_COL,
          " genotype=", GENO_COL, " sample=", SAMPLE_COL)
  # Ensure celltype factor is ordered by abundance (stable legend/axis order)
  female_obj@meta.data[[CT_COL]] <- factor(as.character(female_obj@meta.data[[CT_COL]]),
                                           levels = CT_LEVELS)
}

# ---- 8.1 Split-violin geom (used by the Nudt16 validation panel) ------------
# Standard split-violin implementation (jan-glx / vdb style), defined once.
GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin", ggplot2::GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- plyr_arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
                            if (grp %% 2 == 1) data$y else -data$y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)),
                         setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
    }
  })
# Small helper to avoid a hard plyr dependency for the arrange used above.
plyr_arrange <- function(df, ord) df[order(ord), , drop = FALSE]
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., draw_quantiles = NULL,
                              trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles,
                  na.rm = na.rm, ...))
}

# ---- 8.2 Nudt16 KO-validation panel (panelE) --------------------------------
safe_panel("nudt16_validation", {
  gene    <- VOLC_ANCHOR                     # "Nudt16"
  out     <- function(...) rf("panelE", ...)
  univ    <- rownames(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data"))
  if (!(gene %in% univ)) stop(gene, " not present in assay ", ASSAY_FIG)
  DefaultAssay(female_obj) <- ASSAY_FIG

  point_cols <- c(WT = "#1B3F5C", KO = "#7A2A20")

  # (a) FeaturePlot split by genotype + combined
  p_split <- FeaturePlot(female_obj, features = gene, split.by = GENO_COL,
                         order = TRUE, pt.size = 0.35) &
    theme(plot.title = element_text(face = "italic"))
  save_plot(p_split, out(paste0(gene, "_UMAP_split")), width = 10, height = 5)
  p_comb <- FeaturePlot(female_obj, features = gene, order = TRUE, pt.size = 0.35) +
    ggtitle(gene) + theme(plot.title = element_text(face = "italic"))
  save_plot(p_comb, out(paste0(gene, "_UMAP_combined")), width = 6, height = 5.5)

  # Expression vector + genotype
  expr <- as.numeric(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")[gene, ])
  geno <- factor(as.character(female_obj@meta.data[[GENO_COL]]), levels = c("WT", "KO"))
  ct   <- factor(as.character(female_obj@meta.data[[CT_COL]]), levels = CT_LEVELS)
  df   <- data.frame(expr = expr, genotype = geno, celltype = ct)

  # (b) Overall violin + jitter with % expressing labels
  pct_lab <- df %>% group_by(genotype) %>%
    summarise(pct = 100 * mean(expr > 0), y = max(df$expr) * 1.05, .groups = "drop")
  p_ov <- ggplot(df, aes(genotype, expr, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, colour = "grey25", linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.25, alpha = 0.25, colour = "grey20") +
    geom_text(data = pct_lab, aes(genotype, y, label = sprintf("%.1f%%", pct)),
              inherit.aes = FALSE, size = 3.4, fontface = "bold") +
    scale_fill_manual(values = genotype_colors, guide = "none") +
    labs(x = NULL, y = paste0(gene, " (log-norm)"),
         title = paste0(gene, " expression by genotype")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "italic"))
  save_plot(p_ov, out(paste0(gene, "_violin_overall")), width = 4.5, height = 5)

  # (c) Expressing-cells-only violin
  df_e <- df[df$expr > 0, ]
  if (nrow(df_e) > 5) {
    p_e <- ggplot(df_e, aes(genotype, expr, fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, colour = "grey25", linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.3, alpha = 0.3, colour = "grey20") +
      scale_fill_manual(values = genotype_colors, guide = "none") +
      labs(x = NULL, y = paste0(gene, " (log-norm, expressing only)"),
           title = paste0(gene, ", expressing cells")) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "italic"))
    save_plot(p_e, out(paste0(gene, "_violin_overall_expressing")), width = 4.5, height = 5)
  }

  # (d) Per-celltype split violin (WT|KO)
  p_ct <- ggplot(df, aes(celltype, expr, fill = genotype)) +
    geom_split_violin(scale = "width", trim = TRUE, colour = "grey25", linewidth = 0.25) +
    scale_fill_manual(values = genotype_colors, name = NULL) +
    labs(x = NULL, y = paste0(gene, " (log-norm)"),
         title = paste0(gene, " by cell type and genotype")) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          plot.title = element_text(face = "italic"), legend.position = "top")
  save_plot(p_ct, out(paste0(gene, "_violin_bycelltype")), width = 9, height = 5)

  # (e) Detection statistics CSV
  det <- df %>% group_by(celltype, genotype) %>%
    summarise(n_cells = dplyr::n(), pct_expressing = 100 * mean(expr > 0),
              mean_lognorm = mean(expr), .groups = "drop")
  write.csv(det, out(paste0(gene, "_detection_stats.csv")), row.names = FALSE)
  message("  panelE: wrote UMAPs, violins and detection stats for ", gene)
})

# ---- 8.3 Validation violins for KO-response genes (gene_plots) --------------
safe_panel("validation_violins", {
  validation_genes <- c("Ccnb2", "Cdk1", "E2f8", "Ermap", "Marchf8", "Pdzk1ip1", "Ung")
  DefaultAssay(female_obj) <- ASSAY_FIG
  univ <- rownames(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data"))
  genes_use <- resolve_genes(as.list(validation_genes), univ)
  if (!length(genes_use)) stop("None of the validation genes found in assay.")
  geno <- factor(as.character(female_obj@meta.data[[GENO_COL]]), levels = c("WT", "KO"))
  ct   <- factor(as.character(female_obj@meta.data[[CT_COL]]), levels = CT_LEVELS)

  violin_layers <- function(p) {
    p +
      geom_violin(scale = "width", trim = TRUE, colour = "grey25", linewidth = 0.25,
                  position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = genotype_colors, name = NULL) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            plot.title = element_text(face = "italic"), legend.position = "top")
  }

  # Per-gene dodged violins (cell type x genotype)
  for (g in genes_use) {
    local({
      gg <- g
      ex <- as.numeric(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")[gg, ])
      d  <- data.frame(expr = ex, genotype = geno, celltype = ct)
      p  <- violin_layers(ggplot(d, aes(celltype, expr, fill = genotype))) +
        labs(x = NULL, y = paste0(gg, " (log-norm)"), title = gg)
      save_plot(p, rf("gene_plots", paste0("violin_celltype_by_genotype_", gg)),
                width = 9, height = 4.5)
    })
  }

  # Faceted multi-gene panel
  long <- do.call(rbind, lapply(genes_use, function(g) {
    ex <- as.numeric(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")[g, ])
    data.frame(gene = g, expr = ex, genotype = geno, celltype = ct)
  }))
  long$gene <- factor(long$gene, levels = genes_use)
  p_panel <- violin_layers(ggplot(long, aes(celltype, expr, fill = genotype))) +
    facet_wrap(~ gene, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "log-normalized expression",
         title = "KO-response validation genes by cell type and genotype")
  save_plot(p_panel, rf("gene_plots", "validation_violins_PANEL"),
            width = 12, height = 2.6 * ceiling(length(genes_use) / 2))
  message("  gene_plots: ", length(genes_use), " validation-gene violins + panel")
})

# ---- 8.4 Per-cell-type highlight UMAPs (umap) -------------------------------
safe_panel("celltype_highlight_umaps", {
  PT_SIZE <- 0.30
  ct_vec  <- factor(as.character(female_obj@meta.data[[CT_COL]]), levels = CT_LEVELS)
  n_total <- length(ct_vec)
  for (ctype in CT_LEVELS) {
    local({
      ctype   <- ctype
      cells   <- colnames(female_obj)[which(ct_vec == ctype)]
      if (!length(cells)) return(invisible(NULL))
      hex     <- unname(CT_PALETTE[ctype]); if (is.na(hex)) hex <- "#D55E00"
      pct     <- 100 * length(cells) / n_total
      ttl     <- sprintf("%s (n = %s, %.1f%%)", ctype, format(length(cells), big.mark = ","), pct)
      p <- DimPlot(female_obj, cells.highlight = cells, sizes.highlight = PT_SIZE,
                   cols = "grey85", cols.highlight = hex, pt.size = PT_SIZE) +
        ggtitle(ttl) + NoLegend() +
        theme(plot.title = element_text(size = 11, face = "bold"))
      stem <- gsub("[^A-Za-z0-9]+", "_", ctype)
      save_plot(p, rf("umap", paste0("umap_highlight_", stem)), width = 5.5, height = 5)
    })
  }
  message("  umap: highlight UMAPs for ", length(CT_LEVELS), " cell types")
})

# ---- 8.5 Combined cell-type UMAP panel grid (umap) --------------------------
safe_panel("celltype_umap_panel", {
  ct_vec <- factor(as.character(female_obj@meta.data[[CT_COL]]), levels = CT_LEVELS)
  emb    <- Embeddings(female_obj, "umap")
  base_df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], celltype = ct_vec)
  make_panel <- function(ctype) {
    hex <- unname(CT_PALETTE[ctype]); if (is.na(hex)) hex <- "#D55E00"
    d   <- base_df; d$fg <- d$celltype == ctype
    ggplot() +
      geom_point(data = d[!d$fg, ], aes(UMAP_1, UMAP_2), colour = "grey88", size = 0.15) +
      geom_point(data = d[d$fg, ],  aes(UMAP_1, UMAP_2), colour = hex,      size = 0.30) +
      ggtitle(ctype) +
      theme_void(base_size = 10) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  }
  panels <- lapply(CT_LEVELS, make_panel)
  grid   <- patchwork::wrap_plots(panels, ncol = 4)
  n_row  <- ceiling(length(CT_LEVELS) / 4)
  save_plot(grid, rf("umap", "celltype_umap_panel_grid"),
            width = 14, height = 3.2 * n_row)
  message("  umap: combined ", length(CT_LEVELS), "-type panel grid")
})

# ---- 8.6 Figure 1 composition: dodged bars + differential abundance --------
safe_panel("composition_figure1", {
  KO_LIGHTEN <- 0.40
  md <- female_obj@meta.data
  md$celltype <- factor(as.character(md[[CT_COL]]), levels = CT_LEVELS)
  md$genotype <- factor(as.character(md[[GENO_COL]]), levels = c("WT", "KO"))
  md$sample   <- as.character(md[[SAMPLE_COL]])

  # Per-genotype cell-type fractions
  frac <- md %>% group_by(genotype, celltype) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    group_by(genotype) %>% mutate(frac = 100 * n / sum(n)) %>% ungroup()
  write.csv(frac, rf("composition", "CellType_fraction_by_genotype.csv"), row.names = FALSE)

  # PanelA: dodged bars, WT solid + KO lightened (fill identity so hue = cell type)
  frac$bar_fill <- ifelse(frac$genotype == "KO",
                          lighten_hex(unname(CT_PALETTE[as.character(frac$celltype)]), KO_LIGHTEN),
                          unname(CT_PALETTE[as.character(frac$celltype)]))
  pA <- ggplot(frac, aes(celltype, frac, group = genotype, fill = bar_fill)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             colour = "grey25", linewidth = 0.2) +
    scale_fill_identity() +
    labs(x = NULL, y = "% of cells", title = "Cell-type composition (WT solid, KO lighter)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          plot.title = element_text(face = "bold", size = 12))
  save_plot(pA, rf("composition", "Fig1A_composition_bars"), width = 8.5, height = 5)

  # PanelB: differential abundance (propeller if available + >=3 samples/group, else prop.test)
  n_samp <- length(unique(md$sample))
  da_tbl <- NULL
  if (HAVE_SPECKLE && n_samp >= 3) {
    da <- tryCatch(
      speckle::propeller(clusters = md$celltype, sample = md$sample, group = md$genotype),
      error = function(e) NULL)
    if (!is.null(da)) {
      da_tbl <- data.frame(celltype = rownames(da), da, row.names = NULL)
      names(da_tbl) <- tolower(names(da_tbl))
      pcol <- grep("fdr|adj|p.value|pval", names(da_tbl), value = TRUE)[1]
      lcol <- grep("logfc|log2", names(da_tbl), value = TRUE)[1]
      da_tbl$stat <- if (!is.na(lcol)) da_tbl[[lcol]] else NA
      da_tbl$padj <- if (!is.na(pcol)) da_tbl[[pcol]] else NA
    }
  }
  if (is.null(da_tbl)) {
    # prop.test per cell type: KO vs WT counts
    tab <- md %>% group_by(genotype, celltype) %>%
      summarise(n = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = genotype, values_from = n, values_fill = 0)
    tot <- md %>% group_by(genotype) %>% summarise(N = dplyr::n(), .groups = "drop")
    NWT <- tot$N[tot$genotype == "WT"]; NKO <- tot$N[tot$genotype == "KO"]
    res <- lapply(seq_len(nrow(tab)), function(i) {
      wt <- ifelse("WT" %in% names(tab), tab$WT[i], 0)
      ko <- ifelse("KO" %in% names(tab), tab$KO[i], 0)
      pt <- tryCatch(prop.test(c(ko, wt), c(NKO, NWT)), error = function(e) NULL)
      data.frame(celltype = tab$celltype[i],
                 log2FC = log2(((ko / NKO) + 1e-6) / ((wt / NWT) + 1e-6)),
                 p = if (is.null(pt)) NA else pt$p.value)
    })
    da_tbl <- do.call(rbind, res)
    da_tbl$padj <- p.adjust(da_tbl$p, "BH")
    da_tbl$stat <- da_tbl$log2FC
  }
  write.csv(da_tbl, rf("composition", "CellType_differential_abundance.csv"), row.names = FALSE)
  da_tbl$celltype <- factor(as.character(da_tbl$celltype), levels = CT_LEVELS)
  da_tbl$dir <- ifelse(da_tbl$stat >= 0, "KO up", "WT up")
  da_tbl$star <- ifelse(!is.na(da_tbl$padj) & da_tbl$padj < 0.05, "*", "")
  pB <- ggplot(da_tbl, aes(celltype, stat, fill = dir)) +
    geom_col(colour = "grey25", linewidth = 0.2, width = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = star, vjust = ifelse(stat >= 0, -0.1, 1.1)), size = 5) +
    scale_fill_manual(values = c("KO up" = COL_UP, "WT up" = COL_DOWN), name = NULL) +
    labs(x = NULL, y = "log2 abundance ratio (KO/WT)",
         title = "Differential cell-type abundance") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          plot.title = element_text(face = "bold", size = 12), legend.position = "top")
  save_plot(pB, rf("composition", "Fig1B_differential_abundance"), width = 8.5, height = 5)

  combined <- pA / pB
  save_plot(combined, rf("composition", "Fig1_composition_combined"), width = 9, height = 10)
  message("  composition: Fig1A/B + combined + CSVs (method: ",
          if (HAVE_SPECKLE && n_samp >= 3) "propeller" else "prop.test", ")")
})

# ---- 8.7 Diverging cell-fraction pyramid + per-sample Wilcoxon (composition)
safe_panel("cellfraction_diverging", {
  BASE_SIZE <- 16
  md <- female_obj@meta.data
  md$celltype <- factor(as.character(md[[CT_COL]]), levels = CT_LEVELS)
  md$genotype <- factor(as.character(md[[GENO_COL]]), levels = c("WT", "KO"))
  md$sample   <- as.character(md[[SAMPLE_COL]])

  # Per-sample fractions, then per-genotype mean
  comp_sample <- md %>% group_by(sample, genotype, celltype) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    group_by(sample) %>% mutate(freq = 100 * n / sum(n)) %>% ungroup()
  comp_geno <- comp_sample %>% group_by(genotype, celltype) %>%
    summarise(mean_freq = mean(freq), .groups = "drop")

  # PanelA: population pyramid (WT to the left/negative, KO to the right/positive)
  comp_geno$signed <- ifelse(comp_geno$genotype == "WT",
                             -comp_geno$mean_freq, comp_geno$mean_freq)
  pyr <- ggplot(comp_geno, aes(celltype, signed, fill = celltype)) +
    geom_col(colour = "grey25", linewidth = 0.2, width = 0.78) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = CT_PALETTE, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(abs(x), "%")) +
    coord_flip() +
    annotate("text", x = length(CT_LEVELS) + 0.5, y = -max(comp_geno$mean_freq) * 0.6,
             label = "WT", fontface = "bold", size = 5) +
    annotate("text", x = length(CT_LEVELS) + 0.5, y =  max(comp_geno$mean_freq) * 0.6,
             label = "KO", fontface = "bold", size = 5) +
    labs(x = NULL, y = "Mean % of cells", title = "Cell-type composition (WT | KO)") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(face = "bold"))
  save_plot(pyr, rf("composition", "Fig1_cellfraction_diverging"), width = 9, height = 8)

  # PanelB: per-sample Wilcoxon (KO vs WT) on per-sample fractions, BH-adjusted
  wilx <- comp_sample %>% group_by(celltype) %>%
    summarise(
      p = tryCatch(suppressWarnings(
            wilcox.test(freq[genotype == "KO"], freq[genotype == "WT"])$p.value),
            error = function(e) NA_real_),
      wt_mean = mean(freq[genotype == "WT"]),
      ko_mean = mean(freq[genotype == "KO"]), .groups = "drop") %>%
    mutate(padj = p.adjust(p, "BH"),
           log2FC = log2((ko_mean + 1e-6) / (wt_mean + 1e-6)),
           dir = ifelse(log2FC >= 0, "KO up", "WT up"),
           star = ifelse(!is.na(padj) & padj < 0.05, "*", ""))
  write.csv(wilx, rf("composition", "Fig1_diffabundance_wilcox_stats.csv"), row.names = FALSE)
  wilx$celltype <- factor(as.character(wilx$celltype), levels = CT_LEVELS)
  pW <- ggplot(wilx, aes(celltype, log2FC, fill = dir)) +
    geom_col(colour = "grey25", linewidth = 0.2, width = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = star, vjust = ifelse(log2FC >= 0, -0.1, 1.1)), size = 6) +
    scale_fill_manual(values = c("KO up" = COL_UP, "WT up" = COL_DOWN), name = NULL) +
    coord_flip() +
    labs(x = NULL, y = "log2 fraction ratio (KO/WT)",
         title = "Per-sample differential abundance (Wilcoxon, BH)") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  save_plot(pW, rf("composition", "Fig1_diffabundance_wilcox_padj"), width = 9, height = 8)

  combined <- pyr / pW
  save_plot(combined, rf("composition", "Fig1_cellfraction_combined"), width = 9, height = 14)
  message("  composition: diverging pyramid + per-sample Wilcoxon + CSVs")
})

# ---- 8.8 DDR / repair heatmap (ddr_heatmap) ---------------------------------
safe_panel("ddr_heatmap", {
  PADJ_STAR <- PADJ_CUT      # 0.05
  LFC_SAT   <- 2
  DefaultAssay(female_obj) <- ASSAY_FIG
  univ <- rownames(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data"))
  gtab <- resolve_ddr(ddr_focused, univ)     # pathway / gene / canonical
  if (!nrow(gtab)) stop("No DDR genes resolved against the data.")
  gtab <- gtab[!duplicated(gtab$gene), ]

  geno <- as.character(female_obj@meta.data[[GENO_COL]])
  names(geno) <- colnames(female_obj)
  dat  <- GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")[gtab$gene, , drop = FALSE]
  wt_c <- names(geno)[geno == "WT"]; ko_c <- names(geno)[geno == "KO"]
  wt_mean <- Matrix::rowMeans(dat[, wt_c, drop = FALSE])
  ko_mean <- Matrix::rowMeans(dat[, ko_c, drop = FALSE])
  expr_mat <- cbind(WT = wt_mean, KO = ko_mean)
  z_mat <- t(scale(t(expr_mat))); z_mat[!is.finite(z_mat)] <- 0

  # log2FC + padj from in-memory global pseudobulk DESeq2 (SECTION 6)
  gtab$log2FC <- NA_real_; gtab$padj <- NA_real_
  if (exists("de_global") && all(c("gene", "log2FC", "padj") %in% names(de_global))) {
    m <- match(tolower(gtab$gene), tolower(de_global$gene))
    gtab$log2FC <- de_global$log2FC[m]
    gtab$padj   <- de_global$padj[m]
  }
  need_lfc <- is.na(gtab$log2FC)
  if (any(need_lfc)) {
    eps <- 1e-9
    gtab$log2FC[need_lfc] <- log2((expm1(ko_mean[need_lfc]) + eps) /
                                  (expm1(wt_mean[need_lfc]) + eps))
  }
  gtab$WT_mean <- wt_mean; gtab$KO_mean <- ko_mean
  gtab$WT_high <- !is.na(gtab$log2FC) & gtab$log2FC < 0
  write.csv(gtab, rf("ddr_heatmap", "DDR_gene_table.csv"), row.names = FALSE)

  # Keep WT-high genes for the focused heatmap
  keep <- gtab[gtab$WT_high, ]
  if (!nrow(keep)) { message("  ddr_heatmap: no WT-high DDR genes; skipping heatmap."); return(invisible()) }
  keep <- keep[order(factor(keep$pathway, levels = names(ddr_focused))), ]
  z_use    <- z_mat[keep$gene, , drop = FALSE]
  lfc_vec  <- keep$log2FC
  star_vec <- ifelse(!is.na(keep$padj) & keep$padj < PADJ_STAR, "*", "")
  row_split <- droplevels(factor(keep$pathway, levels = names(ddr_focused)))
  Hh <- max(5, 0.17 * nrow(z_use) + 1.6); Wh <- 5.6

  if (HAVE_COMPLEXHM) {
    suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })
    z_col <- colorRamp2(c(-2, 0, 2), c(COL_DOWN, "white", COL_UP))
    l_col <- colorRamp2(c(-LFC_SAT, 0, LFC_SAT), c(COL_DOWN, "white", COL_UP))
    lfc_m <- matrix(lfc_vec, ncol = 1, dimnames = list(rownames(z_use), "log2FC"))
    ht1 <- Heatmap(z_use, name = "row z", col = z_col,
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   row_split = row_split, row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   row_names_gp = gpar(fontsize = 9, fontface = "italic"),
                   column_names_rot = 0, border = TRUE,
                   rect_gp = gpar(col = "grey85", lwd = 0.4),
                   column_title = "Expression", width = unit(20, "mm"))
    ht2 <- Heatmap(lfc_m, name = "log2FC\nKO/WT", col = l_col,
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = FALSE, column_names_rot = 0, border = TRUE,
                   rect_gp = gpar(col = "grey85", lwd = 0.4),
                   column_title = "KO vs WT", width = unit(12, "mm"),
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if (star_vec[i] != "")
                       grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "bold"))
                   })
    render <- function(dev, path) {
      dev(path)
      draw(ht1 + ht2, merge_legend = TRUE,
           column_title = "DNA damage response / repair genes (KO vs WT)",
           column_title_gp = gpar(fontsize = 12, fontface = "bold"))
      grid::grid.text(DDR_CITATION, x = unit(1, "npc") - unit(2, "mm"),
                      y = unit(1.5, "mm"), just = c("right", "bottom"),
                      gp = gpar(fontsize = 6, col = "grey40"))
      dev.off()
    }
    render(function(p) grDevices::cairo_pdf(p, width = Wh, height = Hh),
           rf("ddr_heatmap", "DDR_heatmap.pdf"))
    render(function(p) grDevices::png(p, width = Wh, height = Hh, units = "in",
                                      res = FIG_DPI, bg = "white"),
           rf("ddr_heatmap", "DDR_heatmap.png"))
    message("  ddr_heatmap: rendered with ComplexHeatmap (", nrow(z_use), " genes)")
  } else {
    suppressPackageStartupMessages(library(pheatmap))
    ann_row <- data.frame(pathway = row_split, log2FC = lfc_vec)
    rownames(ann_row) <- rownames(z_use)
    gaps <- cumsum(as.integer(table(row_split)))
    render <- function(dev, path) {
      dev(path)
      print(pheatmap(z_use, cluster_rows = FALSE, cluster_cols = FALSE,
                     gaps_row = gaps, annotation_row = ann_row,
                     color = COL_HM, fontsize_row = 9,
                     main = "DDR/repair genes (WT vs KO, row z-score)"))
      dev.off()
    }
    render(function(p) grDevices::cairo_pdf(p, width = Wh, height = Hh),
           rf("ddr_heatmap", "DDR_heatmap.pdf"))
    render(function(p) grDevices::png(p, width = Wh, height = Hh, units = "in",
                                      res = FIG_DPI, bg = "white"),
           rf("ddr_heatmap", "DDR_heatmap.png"))
    message("  ddr_heatmap: rendered with pheatmap fallback (", nrow(z_use), " genes)")
  }
})

# ---- 8.9 DDR WT-vs-KO panels: lollipop / log2FC bars / %-expressing ---------
safe_panel("ddr_wt_vs_ko", {
  PADJ_STAR <- PADJ_CUT; LFC_SAT <- 2
  DefaultAssay(female_obj) <- ASSAY_FIG
  univ <- rownames(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data"))
  gtab <- resolve_ddr(ddr_focused, univ)
  gtab <- gtab[!duplicated(gtab$gene), ]
  if (!nrow(gtab)) stop("No DDR genes resolved.")

  geno <- as.character(female_obj@meta.data[[GENO_COL]]); names(geno) <- colnames(female_obj)
  wt_c <- names(geno)[geno == "WT"]; ko_c <- names(geno)[geno == "KO"]
  dat  <- GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")[gtab$gene, , drop = FALSE]
  wt_m <- dat[, wt_c, drop = FALSE]; ko_m <- dat[, ko_c, drop = FALSE]
  gtab$WT_mean <- Matrix::rowMeans(wt_m); gtab$KO_mean <- Matrix::rowMeans(ko_m)
  gtab$pct_WT  <- Matrix::rowMeans(wt_m > 0); gtab$pct_KO <- Matrix::rowMeans(ko_m > 0)
  eps <- 1e-9
  gtab$log2FC <- log2((expm1(gtab$KO_mean) + eps) / (expm1(gtab$WT_mean) + eps))
  # Per-gene Wilcoxon (cells as replicates -> descriptive; documented caveat)
  wt_d <- as.matrix(wt_m); ko_d <- as.matrix(ko_m)
  gtab$wilcox_p <- vapply(seq_len(nrow(dat)), function(i) {
    a <- ko_d[i, ]; b <- wt_d[i, ]
    if (all(a == a[1]) && all(b == b[1]) && a[1] == b[1]) return(NA_real_)
    tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_)
  }, numeric(1))
  gtab$wilcox_padj <- p.adjust(gtab$wilcox_p, "BH")
  gtab$higher_in <- ifelse(gtab$log2FC < 0, "WT", "KO")
  write.csv(gtab, rf("ddr_wt_vs_ko", "DDR_WT_vs_KO_stats.csv"), row.names = FALSE)

  gtab$pathway <- factor(gtab$pathway, levels = names(ddr_focused))
  gtab <- gtab[order(gtab$pathway, gtab$log2FC), ]
  gtab$y <- rev(seq_len(nrow(gtab))); gtab$label <- gtab$gene
  gtab$star <- ifelse(!is.na(gtab$wilcox_padj) & gtab$wilcox_padj < PADJ_STAR, "*", "")
  runs <- rle(as.character(gtab$pathway))
  pw <- data.frame(pathway = runs$values, n = runs$lengths)
  pw$ymax <- rev(cumsum(rev(pw$n))); pw$ymin <- pw$ymax - pw$n + 1
  pw$ymid <- (pw$ymin + pw$ymax) / 2
  WT_COL <- COL_DOWN; KO_COL <- COL_UP
  band_fill <- rep(c("grey96", "white"), length.out = nrow(pw))
  base_bands <- function(p) p +
    geom_rect(data = pw, aes(xmin = -Inf, xmax = Inf, ymin = ymin - 0.5, ymax = ymax + 0.5),
              fill = band_fill, inherit.aes = FALSE) +
    geom_text(data = pw, aes(x = -Inf, y = ymid, label = pathway),
              hjust = -0.05, size = 3.1, fontface = "bold", colour = "grey30",
              inherit.aes = FALSE)
  Wd <- 6.4; Hd <- max(5, 0.24 * nrow(gtab) + 1.8)

  # (1) lollipop
  long <- rbind(data.frame(y = gtab$y, grp = "WT", val = gtab$WT_mean),
                data.frame(y = gtab$y, grp = "KO", val = gtab$KO_mean))
  seg <- data.frame(y = gtab$y, x0 = gtab$WT_mean, x1 = gtab$KO_mean,
                    down = gtab$KO_mean < gtab$WT_mean)
  p1 <- base_bands(ggplot()) +
    geom_segment(data = seg, aes(x0, y, xend = x1, yend = y, colour = down),
                 linewidth = 0.7, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = KO_COL, `FALSE` = "grey60"), guide = "none") +
    geom_point(data = long, aes(val, y, fill = grp), shape = 21, size = 2.6,
               stroke = 0.3, colour = "grey20") +
    scale_fill_manual(name = NULL, values = c(WT = WT_COL, KO = KO_COL), breaks = c("WT", "KO")) +
    scale_y_continuous(breaks = gtab$y, labels = gtab$label, expand = expansion(add = 0.8)) +
    labs(x = "Mean log-normalized expression", y = NULL,
         title = "DDR / repair genes: WT vs NUDT16 KO", caption = DDR_CITATION) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(face = "italic", size = 8), legend.position = "top",
          plot.title = element_text(face = "bold"), plot.caption = element_text(size = 6, colour = "grey45"))
  save_plot(p1, rf("ddr_wt_vs_ko", "DDR_WT_vs_KO_lollipop"), width = Wd, height = Hd)

  # (2) log2FC bars
  gtab$lfc_clip <- pmax(pmin(gtab$log2FC, LFC_SAT), -LFC_SAT)
  p2 <- base_bands(ggplot()) +
    geom_col(data = gtab, aes(lfc_clip, y, fill = higher_in), width = 0.62,
             colour = "grey20", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(data = gtab, aes(lfc_clip, y, label = star,
              hjust = ifelse(lfc_clip >= 0, -0.2, 1.2)), size = 4.4, vjust = 0.75) +
    scale_fill_manual(name = "Higher in", values = c(WT = WT_COL, KO = KO_COL), breaks = c("KO", "WT")) +
    scale_y_continuous(breaks = gtab$y, labels = gtab$label, expand = expansion(add = 0.8)) +
    labs(x = "log2 fold change (KO / WT)", y = NULL,
         title = "DDR / repair genes: KO vs WT fold change",
         subtitle = "* Wilcoxon BH-padj < 0.05 (cells as replicates; descriptive)",
         caption = DDR_CITATION) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(face = "italic", size = 8), legend.position = "top",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8.5, colour = "grey35"),
          plot.caption = element_text(size = 6, colour = "grey45"))
  save_plot(p2, rf("ddr_wt_vs_ko", "DDR_WT_vs_KO_log2FC"), width = Wd, height = Hd)

  # (3) % expressing
  longp <- rbind(data.frame(y = gtab$y, grp = "WT", val = 100 * gtab$pct_WT),
                 data.frame(y = gtab$y, grp = "KO", val = 100 * gtab$pct_KO))
  segp <- data.frame(y = gtab$y, x0 = 100 * gtab$pct_WT, x1 = 100 * gtab$pct_KO,
                     down = gtab$pct_KO < gtab$pct_WT)
  p3 <- base_bands(ggplot()) +
    geom_segment(data = segp, aes(x0, y, xend = x1, yend = y, colour = down),
                 linewidth = 0.7, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = KO_COL, `FALSE` = "grey60"), guide = "none") +
    geom_point(data = longp, aes(val, y, fill = grp), shape = 21, size = 2.6,
               stroke = 0.3, colour = "grey20") +
    scale_fill_manual(name = NULL, values = c(WT = WT_COL, KO = KO_COL), breaks = c("WT", "KO")) +
    scale_y_continuous(breaks = gtab$y, labels = gtab$label, expand = expansion(add = 0.8)) +
    labs(x = "% of cells expressing", y = NULL,
         title = "DDR / repair genes: fraction expressing, WT vs KO", caption = DDR_CITATION) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(face = "italic", size = 8), legend.position = "top",
          plot.title = element_text(face = "bold"), plot.caption = element_text(size = 6, colour = "grey45"))
  save_plot(p3, rf("ddr_wt_vs_ko", "DDR_WT_vs_KO_pctcells"), width = Wd, height = Hd)
  message("  ddr_wt_vs_ko: lollipop + log2FC + pctcells (", nrow(gtab), " genes)")
})

# ---- 8.10 DDR diagnostics: is the "DDR-down-in-KO" pattern real? ------------
safe_panel("ddr_diagnostics", {
  MIN_PCT <- 0.10; N_PERM <- 2000; N_EXPR_BINS <- 20
  set.seed(217)
  DefaultAssay(female_obj) <- ASSAY_FIG
  univ <- rownames(GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data"))
  gtab <- resolve_ddr(ddr_focused, univ); gtab <- gtab[!duplicated(gtab$gene), ]
  ddr_genes <- gtab$gene
  if (!length(ddr_genes)) stop("No DDR genes resolved.")

  md   <- female_obj@meta.data
  geno <- factor(as.character(md[[GENO_COL]]), levels = c("WT", "KO")); names(geno) <- rownames(md)
  samp <- as.character(md[[SAMPLE_COL]]); names(samp) <- rownames(md)
  wt_c <- names(geno)[geno == "WT"]; ko_c <- names(geno)[geno == "KO"]
  samp_geno <- tapply(as.character(geno), samp, function(x) names(sort(table(x), decreasing = TRUE))[1])

  cnt <- GetAssayData(female_obj, assay = ASSAY_FIG, layer = "counts")
  nCount   <- if ("nCount_RNA" %in% colnames(md))   md$nCount_RNA   else Matrix::colSums(cnt)
  nFeature <- if ("nFeature_RNA" %in% colnames(md)) md$nFeature_RNA else Matrix::colSums(cnt > 0)

  # (1) DEPTH
  depth_cell <- data.frame(cell = rownames(md), sample = samp, genotype = geno,
                           nCount = as.numeric(nCount), nFeature = as.numeric(nFeature))
  depth_mouse <- depth_cell %>% group_by(sample) %>%
    summarise(genotype = dplyr::first(genotype), mean_nCount = mean(nCount),
              mean_nFeature = mean(nFeature), n_cells = dplyr::n(), .groups = "drop")
  write.csv(depth_cell,  rf("ddr_diagnostics", "depth_per_cell.csv"),  row.names = FALSE)
  write.csv(depth_mouse, rf("ddr_diagnostics", "depth_per_mouse.csv"), row.names = FALSE)
  dt_count <- tryCatch(t.test(mean_nCount ~ genotype, depth_mouse)$p.value, error = function(e) NA)
  dt_feat  <- tryCatch(t.test(mean_nFeature ~ genotype, depth_mouse)$p.value, error = function(e) NA)
  wt_md <- median(depth_cell$nCount[depth_cell$genotype == "WT"])
  ko_md <- median(depth_cell$nCount[depth_cell$genotype == "KO"])
  depth_ratio <- ko_md / wt_md
  grDevices::png(rf("ddr_diagnostics", "depth_WT_vs_KO.png"), width = 6.5, height = 5,
                 units = "in", res = FIG_DPI, bg = "white")
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (v in c("nCount", "nFeature")) {
    vals <- depth_cell[[v]]
    boxplot(log10(vals + 1) ~ depth_cell$genotype, col = c(COL_DOWN, COL_UP),
            outline = FALSE, ylab = paste0("log10(", v, ")"), xlab = "", main = v)
  }
  par(op); dev.off()

  # (2) PSEUDOBULK per mouse
  agg <- tryCatch(
    AggregateExpression(female_obj, assays = ASSAY_FIG, group.by = SAMPLE_COL,
                        slot = "counts")[[ASSAY_FIG]],
    error = function(e)
      sapply(unique(samp), function(s) Matrix::rowSums(cnt[, names(samp)[samp == s], drop = FALSE])))
  agg <- as.matrix(agg)
  clean <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
  sg_clean <- setNames(as.character(samp_geno), clean(names(samp_geno)))
  pb_geno <- unname(sg_clean[clean(colnames(agg))])
  logcpm <- log2(t(t(agg) / colSums(agg)) * 1e6 + 1)
  wt_s <- colnames(agg)[pb_geno == "WT"]; ko_s <- colnames(agg)[pb_geno == "KO"]
  lc <- logcpm[ddr_genes, , drop = FALSE]
  pb <- data.frame(gene = ddr_genes, pathway = gtab$pathway,
                   WT_logcpm = rowMeans(lc[, wt_s, drop = FALSE]),
                   KO_logcpm = rowMeans(lc[, ko_s, drop = FALSE]))
  pb$pb_log2FC <- pb$KO_logcpm - pb$WT_logcpm
  dat <- GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")[ddr_genes, , drop = FALSE]
  pb$cell_log2FC <- log2((expm1(Matrix::rowMeans(dat[, ko_c, drop = FALSE])) + 1e-9) /
                         (expm1(Matrix::rowMeans(dat[, wt_c, drop = FALSE])) + 1e-9))
  concord    <- suppressWarnings(cor(pb$cell_log2FC, pb$pb_log2FC, use = "complete.obs"))
  sign_agree <- mean(sign(pb$cell_log2FC) == sign(pb$pb_log2FC), na.rm = TRUE)
  write.csv(pb, rf("ddr_diagnostics", "pseudobulk_DDR.csv"), row.names = FALSE)
  grDevices::png(rf("ddr_diagnostics", "pseudobulk_concordance.png"), width = 6.5, height = 6,
                 units = "in", res = FIG_DPI, bg = "white")
  plot(pb$cell_log2FC, pb$pb_log2FC, pch = 21, bg = "#4477AA",
       xlab = "cell-level log2FC (KO/WT)", ylab = "pseudobulk log2FC (per-mouse)",
       main = sprintf("cell vs pseudobulk (r=%.2f, sign agree=%.0f%%)", concord, 100 * sign_agree))
  abline(0, 1, lty = 2, col = "grey50"); abline(h = 0, v = 0, col = "grey80")
  text(pb$cell_log2FC, pb$pb_log2FC, pb$gene, pos = 3, cex = 0.45, col = "grey30"); dev.off()

  # (3) EXPRESSION FILTER
  pct_wt <- Matrix::rowMeans(cnt[ddr_genes, wt_c, drop = FALSE] > 0)
  pct_ko <- Matrix::rowMeans(cnt[ddr_genes, ko_c, drop = FALSE] > 0)
  expr_tbl <- data.frame(gene = ddr_genes, pathway = gtab$pathway,
                         pct_WT = pct_wt, pct_KO = pct_ko,
                         max_pct = pmax(pct_wt, pct_ko),
                         expressed = pmax(pct_wt, pct_ko) >= MIN_PCT)
  write.csv(expr_tbl, rf("ddr_diagnostics", "expression_filter.csv"), row.names = FALSE)
  n_expr <- sum(expr_tbl$expressed)

  # (4) SPECIFICITY vs expression-matched random gene sets
  data_all <- GetAssayData(female_obj, assay = ASSAY_FIG, layer = "data")
  gene_mean_all <- Matrix::rowMeans(data_all)
  pool <- names(gene_mean_all)[gene_mean_all > 0]
  bins <- cut(rank(gene_mean_all[pool]), breaks = N_EXPR_BINS, labels = FALSE); names(bins) <- pool
  km <- Matrix::rowMeans(data_all[pool, ko_c, drop = FALSE])
  wm <- Matrix::rowMeans(data_all[pool, wt_c, drop = FALSE])
  lfc_all <- log2((expm1(km) + 1e-9) / (expm1(wm) + 1e-9)); names(lfc_all) <- pool
  ddr_in_pool <- ddr_genes[ddr_genes %in% pool]
  obs_mean_lfc <- mean(lfc_all[ddr_in_pool], na.rm = TRUE)
  ddr_bins <- bins[ddr_in_pool]
  null_mean <- replicate(N_PERM, {
    picks <- vapply(ddr_bins, function(b) {
      cand <- setdiff(names(bins)[bins == b], ddr_genes)
      if (length(cand)) sample(cand, 1) else NA_character_
    }, character(1))
    mean(lfc_all[picks[!is.na(picks)]], na.rm = TRUE)
  })
  emp_p <- (1 + sum(null_mean <= obs_mean_lfc)) / (N_PERM + 1)
  grDevices::png(rf("ddr_diagnostics", "specificity_vs_random.png"), width = 6.5, height = 5,
                 units = "in", res = FIG_DPI, bg = "white")
  hist(null_mean, breaks = 40, col = "grey85", border = "white",
       xlab = "mean log2FC(KO/WT) of matched random set",
       main = sprintf("DDR vs %d random sets (obs=%.3f, emp p=%.3f)", N_PERM, obs_mean_lfc, emp_p),
       xlim = range(c(null_mean, obs_mean_lfc)))
  abline(v = obs_mean_lfc, col = COL_UP, lwd = 2); dev.off()

  # (5) MODULE SCORE tested per mouse
  female_obj <- tryCatch(
    AddModuleScore(female_obj, features = list(DDR = ddr_genes), name = "DDRscore",
                   ctrl = min(100, floor(length(pool) / 10)), seed = 1),
    error = function(e) female_obj)
  score_col <- grep("^DDRscore", colnames(female_obj@meta.data), value = TRUE)[1]
  mod_p <- NA
  if (!is.na(score_col)) {
    sc <- female_obj@meta.data[[score_col]]
    mod_cell <- data.frame(sample = samp, genotype = geno, score = sc)
    mod_tbl <- mod_cell %>% group_by(sample) %>%
      summarise(genotype = dplyr::first(genotype), mean_score = mean(score), .groups = "drop")
    write.csv(mod_tbl, rf("ddr_diagnostics", "module_score_per_mouse.csv"), row.names = FALSE)
    mod_p <- tryCatch(t.test(mean_score ~ genotype, mod_tbl)$p.value, error = function(e) NA)
    grDevices::png(rf("ddr_diagnostics", "module_score_WT_vs_KO.png"), width = 5, height = 5,
                   units = "in", res = FIG_DPI, bg = "white")
    boxplot(score ~ genotype, mod_cell, col = c(COL_DOWN, COL_UP), outline = FALSE,
            ylab = "DDR module score (cell)", xlab = "",
            main = sprintf("DDR module (per-mouse t p = %s)",
                           ifelse(is.na(mod_p), "NA", sprintf("%.3f", mod_p)))); dev.off()
  }

  # VERDICT
  verdict <- c(
    "================ DDR WT-vs-KO DIAGNOSTIC SUMMARY ================",
    sprintf("DEPTH        : median nCount KO/WT = %.3f | mouse t p(nCount)=%s p(nFeature)=%s",
            depth_ratio, ifelse(is.na(dt_count), "NA", sprintf("%.3f", dt_count)),
            ifelse(is.na(dt_feat), "NA", sprintf("%.3f", dt_feat))),
    if (depth_ratio < 0.9) "   -> KO cells notably SHALLOWER: global downshift EXPECTED; regress/downsample depth."
    else if (depth_ratio > 1.1) "   -> KO cells DEEPER; a global downshift is NOT explained by depth."
    else "   -> Depth comparable (within 10%): a uniform downshift is NOT simply depth.",
    sprintf("PSEUDOBULK   : cell-vs-mouse log2FC r=%.2f, sign agreement=%.0f%%", concord, 100 * sign_agree),
    if (sign_agree < 0.6) "   -> Poor concordance: cell-level FCs NOT reproduced at mouse level."
    else "   -> Direction largely survives pseudobulk.",
    sprintf("EXPRESSION   : %d of %d DDR genes detected in >%.0f%% of a genotype. Trust only these.",
            n_expr, nrow(expr_tbl), 100 * MIN_PCT),
    sprintf("SPECIFICITY  : DDR mean log2FC=%.3f vs random %.3f; emp p(one-sided down)=%.3f",
            obs_mean_lfc, mean(null_mean), emp_p),
    if (emp_p > 0.05) "   -> DDR set does NOT shift more than matched-random: effect is GLOBAL, not DDR-specific."
    else "   -> DDR set shifts more than matched-random: effect has DDR specificity.",
    sprintf("MODULE SCORE : per-mouse t-test p = %s", ifelse(is.na(mod_p), "NA", sprintf("%.3f", mod_p))),
    "----------------------------------------------------------------",
    "PAPER GUIDANCE: n = 2 WT + 2 KO mice -> hypothesis-generating.",
    "  * depth KO/WT < ~0.9 OR specificity emp p > 0.05 -> likely technical; not a main claim.",
    "  * keep only EXPRESSED genes; replace cell-Wilcoxon stars with per-mouse pseudobulk t-p.",
    "================================================================")
  writeLines(verdict, rf("ddr_diagnostics", "VERDICT.txt"))
  message("  ddr_diagnostics: 5 checks + VERDICT.txt (depth KO/WT=",
          sprintf("%.2f", depth_ratio), ", spec emp p=", sprintf("%.3f", emp_p), ")")
})

# =============================================================================
# SECTION 9 — FUNCTIONAL ENRICHMENT
#   9A. Pseudobulk DESeq2 DEGs -> ORA (GO/KEGG/Reactome/Hallmark)  -> enrichment/
#   9B. Seurat FindMarkers -> MSigDB Hallmark GSEA + directional ORA
#                                                            -> enrichment_hallmark/
# =============================================================================
message("\n[9] Functional enrichment (pseudobulk ORA + Hallmark GSEA/ORA)...")

# ---- 9A. Over-representation on pseudobulk DEGs (up / down separately) -------
safe_panel("enrichment_pseudobulk", {
  write_enr <- function(res_list, prefix) {
    if (!length(res_list)) { message("    no enrichment terms for ", prefix); return(invisible()) }
    for (nm in names(res_list)) {
      df <- res_list[[nm]]
      if (is.null(df) || !nrow(df)) next
      write.csv(df, rf("enrichment", paste0(prefix, "_", nm, ".csv")), row.names = FALSE)
    }
  }
  bar_top <- function(res_list, prefix, title, top = 15) {
    combo <- do.call(rbind, lapply(names(res_list), function(nm) {
      df <- res_list[[nm]]
      if (is.null(df) || !nrow(df)) return(NULL)
      data.frame(source = nm, Description = df$Description,
                 padj = df$p.adjust, Count = df$Count, stringsAsFactors = FALSE)
    }))
    if (is.null(combo) || !nrow(combo)) return(invisible())
    combo <- combo[order(combo$padj), ]
    combo <- head(combo, top)
    combo$Description <- factor(combo$Description, levels = rev(combo$Description))
    p <- ggplot(combo, aes(-log10(padj), Description, fill = source)) +
      geom_col(colour = "grey25", linewidth = 0.2) +
      labs(x = expression(-log[10]~"(adjusted "*italic(p)*")"), y = NULL, title = title) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11))
    save_plot(p, rf("enrichment", paste0(prefix, "_topterms")),
              width = 9, height = 0.35 * nrow(combo) + 1.5)
  }

  run_dir <- function(deg_df, label) {
    up <- deg_df$gene[deg_df$log2FC > 0]
    dn <- deg_df$gene[deg_df$log2FC < 0]
    if (length(up)) {
      r <- run_enrichment_analysis(up, gene_universe, paste0(label, " UP"))
      write_enr(r, paste0(safe_sheet(label), "_UP")); bar_top(r, paste0(safe_sheet(label), "_UP"),
                paste0(label, " — up in KO"))
    }
    if (length(dn)) {
      r <- run_enrichment_analysis(dn, gene_universe, paste0(label, " DOWN"))
      write_enr(r, paste0(safe_sheet(label), "_DOWN")); bar_top(r, paste0(safe_sheet(label), "_DOWN"),
                paste0(label, " — down in KO"))
    }
  }

  # Global
  if (nrow(deg_global_filtered)) run_dir(deg_global_filtered, "Global")
  # Per cell type
  for (ct in names(deg_by_celltype_filtered)) {
    local({ ctt <- ct; run_dir(deg_by_celltype_filtered[[ctt]], ctt) })
  }
  message("  enrichment: pseudobulk ORA written to enrichment/")
})

# ---- 9B. MSigDB Hallmark GSEA + directional ORA on Seurat FindMarkers --------
safe_panel("enrichment_hallmark", {
  # Build Hallmark gene sets in mouse symbols (once)
  hm <- tryCatch(msigdbr(species = "Mus musculus", category = "H"),
                 error = function(e) NULL)
  if (is.null(hm) || !nrow(hm)) stop("msigdbr Hallmark set unavailable.")
  hm_list <- split(hm$gene_symbol, hm$gs_name)              # for fgsea
  t2g_sym <- unique(hm[, c("gs_name", "gene_symbol")])      # for enricher (symbol)

  up_fill <- COL_UP; dn_fill <- COL_DOWN

  run_gsea_group <- function(de_df, label) {
    if (is.null(de_df) || !nrow(de_df)) return(invisible())
    rk <- de_df$stat; names(rk) <- de_df$gene
    rk <- rk[is.finite(rk)]; rk <- rk[!duplicated(names(rk))]
    rk <- sort(rk, decreasing = TRUE)
    if (length(rk) < 10) return(invisible())
    set.seed(1)
    fg <- tryCatch(fgsea::fgsea(pathways = hm_list, stats = rk,
                                minSize = SET_MIN, maxSize = SET_MAX, eps = 0),
                   error = function(e) { message("    fgsea failed: ", e$message); NULL })
    if (is.null(fg) || !nrow(fg)) return(invisible())
    fg <- fg[order(fg$padj), ]
    out <- as.data.frame(fg)
    out$leadingEdge <- vapply(out$leadingEdge, function(x) paste(x, collapse = "/"), character(1))
    write.csv(out, rf("enrichment_hallmark", paste0("GSEA_Hallmark_", safe_sheet(label), ".csv")),
              row.names = FALSE)
    sig <- out[!is.na(out$padj) & out$padj < 0.25, ]
    if (!nrow(sig)) return(invisible())
    sig <- sig[order(sig$NES), ]
    sig$pathway <- factor(gsub("^HALLMARK_", "", sig$pathway),
                          levels = gsub("^HALLMARK_", "", sig$pathway))
    p <- ggplot(sig, aes(NES, pathway, fill = NES > 0)) +
      geom_col(colour = "grey25", linewidth = 0.2) +
      geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
      scale_fill_manual(values = c(`TRUE` = up_fill, `FALSE` = dn_fill),
                        labels = c(`TRUE` = "Up in KO", `FALSE` = "Down in KO"), name = NULL) +
      labs(x = "Normalized enrichment score (KO vs WT)", y = NULL,
           title = paste0("Hallmark GSEA — ", label), subtitle = "FDR < 0.25") +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11))
    save_plot(p, rf("enrichment_hallmark", paste0("GSEA_Hallmark_", safe_sheet(label))),
              width = 8.5, height = 0.32 * nrow(sig) + 1.6)
  }

  run_ora_group <- function(de_df, label) {
    if (is.null(de_df) || !nrow(de_df)) return(invisible())
    universe <- de_df$gene
    up <- de_df$gene[de_df$padj < PADJ_CUT & de_df$log2FoldChange >=  LFC_CUT]
    dn <- de_df$gene[de_df$padj < PADJ_CUT & de_df$log2FoldChange <= -LFC_CUT]
    do_dir <- function(genes, dir_lab, fill) {
      genes <- genes[!is.na(genes)]
      if (length(genes) < 3) return(NULL)
      e <- tryCatch(enricher(genes, universe = universe, TERM2GENE = t2g_sym,
                             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
                    error = function(e) NULL)
      if (is.null(e) || !nrow(e@result)) return(NULL)
      df <- as.data.frame(e)
      write.csv(df, rf("enrichment_hallmark",
                       paste0("ORA_Hallmark_", safe_sheet(label), "_", dir_lab, ".csv")),
                row.names = FALSE)
      sig <- df[df$p.adjust < 0.25, ]
      if (!nrow(sig)) return(NULL)
      sig <- head(sig[order(sig$p.adjust), ], 15)
      sig$Description <- factor(gsub("^HALLMARK_", "", sig$Description),
                                levels = rev(gsub("^HALLMARK_", "", sig$Description)))
      p <- ggplot(sig, aes(-log10(p.adjust), Description)) +
        geom_col(fill = fill, colour = "grey25", linewidth = 0.2) +
        labs(x = expression(-log[10]~"(adjusted "*italic(p)*")"), y = NULL,
             title = paste0("Hallmark ORA — ", label, " (", dir_lab, ")")) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 11))
      save_plot(p, rf("enrichment_hallmark",
                      paste0("ORA_Hallmark_", safe_sheet(label), "_", dir_lab)),
                width = 8.5, height = 0.32 * nrow(sig) + 1.6)
    }
    do_dir(up, "UP", up_fill); do_dir(dn, "DOWN", dn_fill)
  }

  # Global + every celltype_marker group
  groups <- c(list(Global = seurat_de_global), seurat_de_by_celltype)
  for (label in names(groups)) {
    local({
      lab <- label; de <- groups[[lab]]
      safe_panel(paste0("hallmark_", lab), {
        run_gsea_group(de, lab); run_ora_group(de, lab)
      })
    })
  }
  message("  enrichment_hallmark: GSEA + directional ORA for ", length(groups), " groups")
})

# =============================================================================
# SECTION 10 — CONSOLIDATED EXCEL WORKBOOK + SESSION INFO
#   One multi-sheet .xlsx bundling every publication table, written to tables/.
# =============================================================================
message("\n[10] Consolidated Excel workbook + sessionInfo...")

safe_panel("excel_export", {
  wb <- openxlsx::createWorkbook()
  hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1",
                              halign = "left", border = "bottom")
  add_sheet <- function(name, df) {
    if (is.null(df) || !nrow(df)) return(invisible())
    sn <- safe_sheet(name)
    # openxlsx requires unique sheet names <=31 chars
    existing <- openxlsx::sheets(wb)
    if (sn %in% existing) sn <- substr(paste0(sn, "_", length(existing)), 1, 31)
    openxlsx::addWorksheet(wb, sn)
    openxlsx::writeData(wb, sn, df, headerStyle = hs)
    openxlsx::freezePane(wb, sn, firstRow = TRUE)
    openxlsx::setColWidths(wb, sn, cols = seq_len(ncol(df)), widths = "auto")
  }

  # Differential expression (pseudobulk DESeq2)
  add_sheet("DESeq2_Global_all", de_global)
  add_sheet("DESeq2_Global_DEG", deg_global_filtered)
  add_sheet("DESeq2_byCelltype_DEG", deg_ct_combined)
  # Differential expression (Seurat FindMarkers)
  if (exists("seurat_de_global") && !is.null(seurat_de_global))
    add_sheet("Seurat_Global_all", seurat_de_global)
  if (exists("seurat_de_by_celltype") && length(seurat_de_by_celltype)) {
    seurat_ct_combined <- dplyr::bind_rows(lapply(names(seurat_de_by_celltype), function(ct) {
      d <- seurat_de_by_celltype[[ct]]; d$celltype <- ct; d }))
    add_sheet("Seurat_byCelltype_all", seurat_ct_combined)
  }

  # Composition / abundance tables (read back the CSVs the panels wrote)
  read_if <- function(path) if (file.exists(path))
    tryCatch(read.csv(path, check.names = FALSE), error = function(e) NULL) else NULL
  add_sheet("CellType_fractions", read_if(rf("composition", "CellType_fraction_by_genotype.csv")))
  add_sheet("CellType_diffAbund",  read_if(rf("composition", "CellType_differential_abundance.csv")))
  add_sheet("CellType_Wilcox",     read_if(rf("composition", "Fig1_diffabundance_wilcox_stats.csv")))
  add_sheet("Cluster_QC",          read_if(rf("cluster_qc", "cluster_qc_summary.csv")))
  add_sheet("Cluster_assignment",  read_if(rf("marker_annotation", "cluster_assignment.csv")))
  add_sheet("Nudt16_detection",    read_if(rf("panelE", paste0(VOLC_ANCHOR, "_detection_stats.csv"))))
  add_sheet("DDR_gene_table",      read_if(rf("ddr_heatmap", "DDR_gene_table.csv")))
  add_sheet("DDR_WT_vs_KO_stats",  read_if(rf("ddr_wt_vs_ko", "DDR_WT_vs_KO_stats.csv")))
  add_sheet("DDR_pseudobulk",      read_if(rf("ddr_diagnostics", "pseudobulk_DDR.csv")))

  xlsx_path <- rf("tables", "NUDT16_females_publication_tables.xlsx")
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  message("  excel: ", length(openxlsx::sheets(wb)), " sheets -> ", xlsx_path)
})

# ---- Session info (provenance) ----------------------------------------------
safe_panel("session_info", {
  si_path <- rf("tables", "sessionInfo.txt")
  writeLines(capture.output(sessionInfo()), si_path)
  message("  sessionInfo -> ", si_path)
})

message("\n========================================================")
message("MASTER PIPELINE COMPLETE.")
message("All outputs under: ", RES_DIR)
message("Primary cell-type labels: ", CELLTYPE_COL,
        " (secondary: ", SECONDARY_COL, ")")
message("========================================================")

# >>> MASTER_END <<<
