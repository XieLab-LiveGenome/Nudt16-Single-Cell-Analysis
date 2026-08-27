# =============================================================================
# run_all.R  —  Driver for the NUDT16 KO female spleen scRNA-seq pipeline
#
#   This sources the numbered stage scripts in order.
#     00_config.R              shared config / libs / helpers / palettes / markers
#     01_preprocessing.R       QC -> SCT -> Harmony -> clusters -> UMAP  (-> processed.rds)
#     02_cell_annotation.R     cluster QC + marker-score annotation      (-> annotated.rds)
#     03_deg_seurat_deseq2.R   pseudobulk DESeq2 + Seurat FindMarkers DE  (-> de_results.rds)
#     04_gsea_hallmark.R       GO/KEGG/Reactome ORA + Hallmark GSEA/ORA
#     05a_export_for_scenic.R  export counts trio + metadata for pySCENIC
#     05b_run_pyscenic.sh      <SHELL>  GRNBoost2 -> cisTarget -> AUCell
#     05c_regulons_tcell_bcell.R  regulators of the 7 validated genes (T/B cells)

t0 <- Sys.time()

if (!isTRUE(getOption("nudt16.config.loaded"))) {
  cand <- c(
    "00_config.R",
    file.path("pipeline", "00_config.R"),
    file.path("results", "females", "pipeline", "00_config.R"),
    "/Users/budankm/Desktop/Sequencing/GONG/results/females/pipeline/00_config.R"
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("run_all.R: cannot locate 00_config.R — run from the pipeline folder.")
  CONFIG_PATH <- normalizePath(hit[1])
  source(CONFIG_PATH)
} else {
  CONFIG_PATH <- getOption("nudt16.config.path", "00_config.R")
}

PIPE_SCRIPTS <- if (exists("PIPE_DIR") && dir.exists(PIPE_DIR)) PIPE_DIR else dirname(CONFIG_PATH)

run_stage <- function(fname) {
  path <- file.path(PIPE_SCRIPTS, fname)
  if (!file.exists(path)) stop("run_all.R: stage script not found: ", path)
  message("\n#############################################################")
  message("### RUN_ALL  ->  ", fname)
  message("#############################################################")
  ts <- Sys.time()
  source(path, local = FALSE)   # stages define objects in .GlobalEnv
  message("### DONE ", fname, "  (", round(difftime(Sys.time(), ts, units = "mins"), 2),
          " min)")
  invisible(TRUE)
}
R_STAGES <- c(
  "01_preprocessing.R",
  "02_cell_annotation.R",
  "03_deg_seurat_deseq2.R",
  "04_gsea_hallmark.R",
  "05a_export_for_scenic.R"
)

for (s in R_STAGES) run_stage(s)

message("\n=============================================================")
message(" run_all.R: R stages 01–05a complete in ",
        round(difftime(Sys.time(), t0, units = "mins"), 2), " min.")
message("=============================================================")
