# =============================================================================
# 05c_regulons_tcell_bcell.R  —  NUDT16 KO female spleen scRNA-seq
#   pySCENIC downstream: within the T-cell and B-cell compartments, identify and
#   plot the regulons (TFs) that regulate the 7 validated genes, and test
#   KO-vs-WT differential regulon activity.
#
# INPUT  : female_obj_annotated.rds        (RDS_ANNOTATED, from stage 02)
#          scenic/auc_mtx.csv              (cells x regulons; from 05b)
#          scenic/regulon_targets.csv      (regulon,target long table; from 05b)
#          scenic/validated_genes.txt      (7 validated genes present; from 05a)
# OUTPUT : scenic/regulons_targeting_validated_genes.csv
#          scenic/tcell/ and scenic/bcell/:
#             regulon_activity_by_genotype.csv
#             differential_regulon_activity.csv       (Wilcoxon KO vs WT, all regulons)
#             focus_regulons_differential.csv         (regulators of the 7 genes)
#             heatmap_focus_regulon_activity.pdf/.png
#             barplot_focus_regulon_deltaAUC.pdf/.png
#             regulator_target_map.pdf/.png           (regulon x 7-gene edge map)
#
# Differential testing is a plain Wilcoxon on the AUCell matrix (no presto), so
# it cannot trigger the Seurat v5 marker crash.
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

message("\n=========  STAGE 05c — SCENIC REGULONS (T & B cells, 7 genes)  =========")

# --- Compartment definitions (edit to widen the T-lineage if desired) --------
TCELL_TYPES <- c("T cells")            # add "NKT","Tgd" to include those lineages
BCELL_TYPES <- c("B cells")

# --- Load metadata + SCENIC outputs ------------------------------------------
female_obj <- load_stage_rds(RDS_ANNOTATED, "02_cell_annotation.R")
md <- female_obj@meta.data
meta <- data.frame(
  barcode        = colnames(female_obj),
  celltype_marker = as.character(md[[CELLTYPE_COL]]),
  genotype       = as.character(md$genotype),
  sample         = as.character(md$sample),
  stringsAsFactors = FALSE
)

auc_path <- rf("scenic", "auc_mtx.csv")
tgt_path <- rf("scenic", "regulon_targets.csv")
vg_path  <- rf("scenic", "validated_genes.txt")
for (p in c(auc_path, tgt_path))
  if (!file.exists(p)) stop("Missing SCENIC output: ", p, "\n  Run 05b_run_pyscenic.sh first.")

auc <- read.csv(auc_path, row.names = 1, check.names = FALSE)
message("  AUCell matrix: ", nrow(auc), " cells x ", ncol(auc), " regulons")
reg_targets <- read.csv(tgt_path, stringsAsFactors = FALSE)
vg <- if (file.exists(vg_path)) readLines(vg_path) else VALIDATED_GENES
vg <- vg[nzchar(vg)]
message("  validated genes: ", paste(vg, collapse = ", "))

# Align AUC rows to metadata barcodes.
common <- intersect(rownames(auc), meta$barcode)
if (!length(common)) stop("No overlap between AUCell cell IDs and Seurat barcodes.")
auc  <- auc[common, , drop = FALSE]
meta <- meta[match(common, meta$barcode), ]
stopifnot(identical(rownames(auc), meta$barcode))

# =============================================================================
# Identify the FOCUS regulons: those that regulate any of the 7 validated genes
#   (target-based) OR whose TF is itself one of the 7 (e.g. E2f8 is a TF).
# =============================================================================
targets_hit <- reg_targets[reg_targets$target %in% vg, , drop = FALSE]
focus_from_targets <- unique(targets_hit$regulon)
# regulon names look like "E2f8(+)"; TF-in-7 check by stripping the suffix.
reg_tf <- sub("\\(.*$", "", colnames(auc))
focus_from_tf <- colnames(auc)[reg_tf %in% vg]
focus_regulons <- intersect(unique(c(focus_from_targets, focus_from_tf)), colnames(auc))

# Table: which validated gene(s) each focus regulon targets.
reg_vg_map <- targets_hit %>%
  dplyr::group_by(regulon) %>%
  dplyr::summarize(targeted_validated_genes = paste(sort(unique(target)), collapse = "; "),
                   n_validated_targets = dplyr::n_distinct(target),
                   .groups = "drop")
focus_tbl <- data.frame(regulon = focus_regulons, stringsAsFactors = FALSE) %>%
  dplyr::left_join(reg_vg_map, by = "regulon") %>%
  dplyr::mutate(tf = sub("\\(.*$", "", regulon),
                tf_is_validated = tf %in% vg,
                targeted_validated_genes = ifelse(is.na(targeted_validated_genes),
                                                  "", targeted_validated_genes))
write.csv(focus_tbl, rf("scenic", "regulons_targeting_validated_genes.csv"), row.names = FALSE)
message("  focus regulons (regulators of the 7 genes): ", length(focus_regulons))
if (!length(focus_regulons))
  message("  NOTE: no regulon targets the 7 genes; the compartment figures will still ",
          "produce genome-wide differential tables.")

# =============================================================================
# Per-compartment analysis
# =============================================================================
# Wilcoxon KO-vs-WT differential regulon activity across ALL regulons.
diff_regulon_activity <- function(auc_sub, geno_sub) {
  ko_i <- geno_sub == "KO"; wt_i <- geno_sub == "WT"
  res <- lapply(colnames(auc_sub), function(r) {
    x <- auc_sub[ko_i, r]; y <- auc_sub[wt_i, r]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    pv <- tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    data.frame(regulon = r, mean_KO = mean(x), mean_WT = mean(y),
               delta_AUC = mean(x) - mean(y),
               log2FC = log2((mean(x) + 1e-6) / (mean(y) + 1e-6)),
               p = pv, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res) || !nrow(res)) return(res)
  res$padj <- stats::p.adjust(res$p, method = "BH")
  res[order(res$padj, -abs(res$delta_AUC)), ]
}

analyze_compartment <- function(types, out_sub, comp_label) {
  cells <- meta$barcode[meta$celltype_marker %in% types]
  if (length(cells) < 10) {
    message("  [skip] ", comp_label, ": <10 cells"); return(invisible())
  }
  out_dir <- rf("scenic", out_sub)
  auc_sub  <- auc[cells, , drop = FALSE]
  geno_sub <- meta$genotype[match(cells, meta$barcode)]
  message("  ", comp_label, ": ", length(cells), " cells (",
          sum(geno_sub == "KO"), " KO / ", sum(geno_sub == "WT"), " WT)")

  # ---- mean regulon activity by genotype (all regulons) --------------------
  act_by_geno <- data.frame(
    regulon = colnames(auc_sub),
    mean_KO = colMeans(auc_sub[geno_sub == "KO", , drop = FALSE]),
    mean_WT = colMeans(auc_sub[geno_sub == "WT", , drop = FALSE]),
    row.names = NULL, stringsAsFactors = FALSE
  )
  write.csv(act_by_geno, file.path(out_dir, "regulon_activity_by_genotype.csv"),
            row.names = FALSE)

  # ---- differential regulon activity (Wilcoxon KO vs WT) -------------------
  diff_all <- diff_regulon_activity(auc_sub, geno_sub)
  if (is.null(diff_all) || !nrow(diff_all)) {
    message("  [skip figures] ", comp_label, ": no testable regulons"); return(invisible())
  }
  write.csv(diff_all, file.path(out_dir, "differential_regulon_activity.csv"),
            row.names = FALSE)

  # ---- focus regulons (regulators of the 7 genes) --------------------------
  foc <- intersect(focus_regulons, colnames(auc_sub))
  if (length(foc)) {
    diff_focus <- diff_all[diff_all$regulon %in% foc, , drop = FALSE]
    diff_focus <- dplyr::left_join(diff_focus, focus_tbl[, c("regulon",
                     "targeted_validated_genes", "tf_is_validated")], by = "regulon")
    write.csv(diff_focus, file.path(out_dir, "focus_regulons_differential.csv"),
              row.names = FALSE)

    # (a) heatmap of focus-regulon mean activity by genotype (z per regulon) --
    safe_panel(paste0("scenic_heatmap_", out_sub), {
      hm <- as.matrix(act_by_geno[match(foc, act_by_geno$regulon), c("mean_WT", "mean_KO")])
      rownames(hm) <- foc
      hm_z <- t(scale(t(hm)))              # z-score each regulon across WT/KO
      hm_z[!is.finite(hm_z)] <- 0
      pheatmap::pheatmap(
        hm_z, cluster_cols = FALSE, cluster_rows = nrow(hm_z) > 1,
        color = COL_HM, border_color = "grey85",
        main = paste0(comp_label, " — regulators of the 7 genes (row z of mean AUCell)"),
        filename = file.path(out_dir, "heatmap_focus_regulon_activity.png"),
        width = 5.5, height = max(3, 0.28 * nrow(hm_z) + 1.5))
      pheatmap::pheatmap(
        hm_z, cluster_cols = FALSE, cluster_rows = nrow(hm_z) > 1,
        color = COL_HM, border_color = "grey85",
        main = paste0(comp_label, " — regulators of the 7 genes (row z of mean AUCell)"),
        filename = file.path(out_dir, "heatmap_focus_regulon_activity.pdf"),
        width = 5.5, height = max(3, 0.28 * nrow(hm_z) + 1.5))
    })

    # (b) delta-AUC bar plot for focus regulons, colored by significance ------
    safe_panel(paste0("scenic_deltabar_", out_sub), {
      dfb <- diff_focus
      dfb$sig <- ifelse(!is.na(dfb$padj) & dfb$padj < 0.05, "FDR<0.05", "ns")
      dfb$dir <- ifelse(dfb$delta_AUC >= 0, "Up in KO", "Down in KO")
      dfb <- dfb[order(dfb$delta_AUC), ]
      dfb$regulon <- factor(dfb$regulon, levels = dfb$regulon)
      p <- ggplot(dfb, aes(delta_AUC, regulon, fill = dir, alpha = sig)) +
        geom_col(colour = "grey25", linewidth = 0.2) +
        geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
        scale_fill_manual(values = c("Up in KO" = COL_UP, "Down in KO" = COL_DOWN), name = NULL) +
        scale_alpha_manual(values = c("FDR<0.05" = 1, "ns" = 0.45), name = NULL) +
        labs(x = "Δ mean AUCell (KO − WT)", y = NULL,
             title = paste0(comp_label, " — differential activity of 7-gene regulators")) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 11))
      save_plot(p, file.path(out_dir, "barplot_focus_regulon_deltaAUC"),
                width = 8, height = 0.32 * nrow(dfb) + 1.8)
    })

    # (c) regulon x validated-gene edge map (weights), regulon dir annotated --
    safe_panel(paste0("scenic_regmap_", out_sub), {
      edges <- reg_targets[reg_targets$regulon %in% foc & reg_targets$target %in% vg, ,
                           drop = FALSE]
      if (nrow(edges)) {
        dir_map <- setNames(ifelse(diff_focus$delta_AUC >= 0, "Up in KO", "Down in KO"),
                            diff_focus$regulon)
        edges$dir <- dir_map[edges$regulon]
        edges$target  <- factor(edges$target, levels = vg)
        edges$regulon <- factor(edges$regulon, levels = rev(sort(unique(edges$regulon))))
        p <- ggplot(edges, aes(target, regulon, fill = weight)) +
          geom_tile(colour = "grey80") +
          scale_fill_gradient(low = "grey92", high = "#08519C", name = "GRN\nweight") +
          labs(x = "validated target gene", y = "regulon (TF)",
               title = paste0(comp_label, " — 7-gene regulator → target map")) +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
                panel.grid = element_blank())
        save_plot(p, file.path(out_dir, "regulator_target_map"),
                  width = max(5, 0.8 * length(vg) + 2),
                  height = 0.3 * length(unique(edges$regulon)) + 2)
      } else {
        message("  [skip regmap] ", comp_label, ": no focus regulon→7-gene edges")
      }
    })
  } else {
    message("  ", comp_label, ": no focus regulons present in this compartment; ",
            "wrote genome-wide differential table only.")
  }
}

analyze_compartment(TCELL_TYPES, "tcell", "T cells")
analyze_compartment(BCELL_TYPES, "bcell", "B cells")

message("\n=========  STAGE 05c complete  =========")
