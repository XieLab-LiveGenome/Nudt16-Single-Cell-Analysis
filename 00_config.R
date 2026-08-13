# =============================================================================
# 00_config.R  —  Shared configuration for the NUDT16 KO female spleen
#                 scRNA-seq publication pipeline.

# options(nudt16.config.loaded = TRUE).
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

# Optional packages  --------
HAVE_REACTOME  <- requireNamespace("ReactomePA",     quietly = TRUE)
HAVE_COMPLEXHM <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
                  requireNamespace("circlize",       quietly = TRUE)
HAVE_SPECKLE   <- requireNamespace("speckle",        quietly = TRUE)

# ---- 0.2 conflicted guards --------------------------------------------------

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
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB for parallel steps

# ---- 0.3 Paths --------------------------------------------------------------
project_dir <- "/Users/budankm/Desktop/Sequencing/GONG"
RES_DIR     <- file.path(project_dir, "results", "females")
PIPE_DIR    <- file.path(RES_DIR, "pipeline")   # location of these scripts
SCENIC_DIR  <- file.path(RES_DIR, "scenic")     # pySCENIC inputs/outputs (stage 05)

# All output subfolders used anywhere in the pipeline (kept per binding decision
# "keep current subfolders"). Created up-front so every stage can write.
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
  "tables",
  # pySCENIC 
  "scenic", file.path("scenic", "tcell"), file.path("scenic", "bcell")
)
for (sub in OUT_SUBDIRS) {
  dir.create(file.path(RES_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}
if (dir.exists(project_dir)) setwd(project_dir)


rf <- function(...) file.path(RES_DIR, ...)


RDS_PROCESSED <- rf("female_obj_processed.rds")   # after stage 01
RDS_ANNOTATED <- rf("female_obj_annotated.rds")   # after stage 02
RDS_DE        <- rf("de_results.rds")             # after stage 03

# ---- 0.4 Global run knobs--------------------------
# Primary cell-type column used by ALL downstream DE / figures / enrichment.
CELLTYPE_COL  <- "celltype_marker"     # BINDING: marker-based primary
SECONDARY_COL <- "ImmGen.labels"       # SingleR secondary validation

# Manual annotation overrides 
#   DROP_CLUSTERS : character vector of cluster ids to remove after labelling
#   MANUAL_LABELS : named list, e.g. list("3" = "T cells", "7" = "DC")
DROP_CLUSTERS <- character(0)
MANUAL_LABELS <- list()

# QC / filtering
MAD_NMADS             <- 3
DOUBLET_SEED          <- 217
MIN_CELLS_PER_LABEL   <- 15      # ImmGen label retention threshold

# Embedding (harmony pipeline from preprocessfinal — NOT dims 1:30 / res 0.5)
N_VAR_FEATURES        <- 3000
N_PCS                 <- 50
N_DIMS_USE            <- 20      # paper STAR Methods uses 30
CLUSTER_RES           <- 0.2     # paper STAR Methods uses 0.5

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


RUN_MARKER_VERIFICATION <- FALSE

# ---- 0.4b The 7 validated genes (drive the pySCENIC regulator analysis) -----
VALIDATED_GENES <- c("Ccnb2", "Cdk1", "E2f8", "Ermap", "Marchf8", "Pdzk1ip1", "Ung")

# ---- 0.5 Color anchors / palettes ------------------------------------------
COL_UP   <- "#B2182B"     # KO-up  (red)
COL_DOWN <- "#2166AC"     # WT-up  (blue)
COL_NS   <- "grey75"
COL_HM   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)

genotype_colors <- c(WT = "#3B6C9C", KO = "#D6604D")


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

pick_col <- function(obj, candidates, fallback = "seurat_clusters") {
  md <- colnames(obj@meta.data)
  hit <- candidates[candidates %in% md]
  if (length(hit)) hit[1] else fallback
}


resolve_genes <- function(genes, universe) {
  out <- character(0)
  for (g in genes) {
    hit <- g[g %in% universe]
    if (length(hit)) out <- c(out, hit[1])
  }
  unique(out)
}


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

safe_panel <- function(name, expr) {
  message("\n>>> PANEL: ", name)
  tryCatch(
    force(expr),
    error = function(e)
      message("  [SKIP] panel '", name, "' failed: ", conditionMessage(e))
  )
}

# Excel-safe sheet name (<=31 chars)
safe_sheet <- function(x) substr(gsub("[^A-Za-z0-9_]", "_", x), 1, 31)

# load_stage_rds(): read a prior stage's .rds with a clear error if missing.
load_stage_rds <- function(path, produced_by) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path,
         "\n  Run ", produced_by, " first.", call. = FALSE)
  }
  message("[00_config] loading ", basename(path))
  readRDS(path)
}

# ---- 0.9 Done ---------------------------------------------------------------
options(nudt16.config.loaded = TRUE)
message("\n================  NUDT16 FEMALES config loaded  ================")
message("Output dir: ", RES_DIR)
message("Primary annotation column: ", CELLTYPE_COL, "  |  Secondary: ", SECONDARY_COL)
message("CLUSTER_RES = ", CLUSTER_RES, " | N_DIMS_USE = ", N_DIMS_USE,
        " | RUN_MARKER_VERIFICATION = ", RUN_MARKER_VERIFICATION)
