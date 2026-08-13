# =============================================================================
# 05a_export_for_scenic.R  —  NUDT16 KO female spleen scRNA-seq
#   Export the annotated Seurat object into the inputs pySCENIC needs.
#
# INPUT  : female_obj_annotated.rds   (RDS_ANNOTATED, from stage 02)
# OUTPUT : scenic/expr_counts.mtx        (genes x cells, raw counts, MatrixMarket)
#          scenic/genes.txt              (gene symbols, matrix row order)
#          scenic/barcodes.txt           (cell barcodes, matrix col order)
#          scenic/cell_metadata.csv      (barcode, celltype_marker, genotype, sample, cluster)
#          scenic/validated_genes.txt    (the 7 validated genes present in the data)
#          scenic/female_counts.loom     (ONLY if SCopeLoomR is installed; optional)
#
# The .mtx trio + metadata is the portable, dependency-light export. 05b builds
# the loom pySCENIC consumes from these files (via loompy) if no loom exists.
# =============================================================================

# --- source shared config ----------------------------------------------------
if (!isTRUE(getOption("nudt16.config.loaded"))) {
  cand <- c(
    "00_config.R",
    file.path("pipeline", "00_config.R"),
    file.path("results", "females", "pipeline", "00_config.R"),
    "/Users/budankm/Desktop/Sequencing/GONG/results/females/pipeline/00_config.R"
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Cannot locate 00_config.R — run from the pipeline folder.")
  source(hit[1])
}

message("\n=========  STAGE 05a — EXPORT FOR pySCENIC  =========")

female_obj <- load_stage_rds(RDS_ANNOTATED, "02_cell_annotation.R")

DefaultAssay(female_obj) <- "RNA"
female_obj <- JoinLayers(female_obj, assay = "RNA")

# Raw counts (genes x cells) — SCENIC/GRNBoost2 expects raw UMI counts.
counts <- GetAssayData(female_obj, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
message("  counts matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

# Write MatrixMarket trio.
Matrix::writeMM(counts, rf("scenic", "expr_counts.mtx"))
writeLines(rownames(counts), rf("scenic", "genes.txt"))
writeLines(colnames(counts), rf("scenic", "barcodes.txt"))

# Cell metadata (barcode order matches barcodes.txt / matrix columns).
md <- female_obj@meta.data
clu_col  <- pick_col(female_obj, c("seurat_clusters"))
cell_meta <- data.frame(
  barcode        = colnames(counts),
  celltype_marker = as.character(md[[CELLTYPE_COL]]),
  genotype       = as.character(md$genotype),
  sample         = as.character(md$sample),
  cluster        = as.character(md[[clu_col]]),
  stringsAsFactors = FALSE
)
write.csv(cell_meta, rf("scenic", "cell_metadata.csv"), row.names = FALSE)

# The 7 validated genes present in the data (targets of interest for 05c).
vg_present <- VALIDATED_GENES[VALIDATED_GENES %in% rownames(counts)]
vg_missing <- setdiff(VALIDATED_GENES, vg_present)
if (length(vg_missing))
  message("  NOTE: validated genes absent from data: ", paste(vg_missing, collapse = ", "))
writeLines(vg_present, rf("scenic", "validated_genes.txt"))

# Optional: write a loom directly if SCopeLoomR is available (skips 05b's loom build).
if (requireNamespace("SCopeLoomR", quietly = TRUE)) {
  message("  SCopeLoomR found — writing female_counts.loom directly")
  loom_path <- rf("scenic", "female_counts.loom")
  if (file.exists(loom_path)) file.remove(loom_path)
  SCopeLoomR::build_loom(
    file.name       = loom_path,
    dgem            = counts,                # genes x cells
    title           = "NUDT16 females spleen",
    default.embedding = NULL
  )
} else {
  message("  SCopeLoomR not installed — 05b will build the loom from the .mtx trio via loompy.")
}

message("\n  Exported to: ", rf("scenic"))
message("  Next: run 05b_run_pyscenic.sh (edit the RESOURCES paths at the top first).")
message("\n=========  STAGE 05a complete  =========")
