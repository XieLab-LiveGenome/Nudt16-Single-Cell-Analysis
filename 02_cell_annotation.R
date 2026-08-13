# =============================================================================
# 02_cell_annotation.R  —  NUDT16 KO female spleen scRNA-seq
#   cluster QC diagnostics -> marker-score annotation -> finalize (commit
#   celltype_marker PRIMARY) -> marker dot plot + primary-label UMAP
#
# INPUT  : female_obj_processed.rds   (RDS_PROCESSED, )
# OUTPUT : female_obj_annotated.rds   (RDS_ANNOTATED)
#          cluster_qc/, marker_annotation/, panelE/, umap/ figures + tables

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

message("\n=========  STAGE 02 — CELL ANNOTATION  =========")

female_obj <- load_stage_rds(RDS_PROCESSED, "01_preprocessing.R")

# =============================================================================
# SECTION 2 — CLUSTER QC DIAGNOSTICS
# =============================================================================
message("\n[2] Cluster QC diagnostics...")

safe_panel("cluster_qc", {
  cqc_dir <- rf("cluster_qc")
  RUN_CLUSTER_MARKERS <- FALSE  
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
  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (isTRUE(RUN_MARKER_VERIFICATION)) {
    message("  [verify] running FindAllMarkers cross-check (RUN_MARKER_VERIFICATION = TRUE)...")
    Idents(obj) <- clu
    mk <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
                         min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
    mk <- mk[mk$p_val_adj < FINALIZE_MARKER_PADJ, ]
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
  } else {
    message("  [skip] FindAllMarkers verification (RUN_MARKER_VERIFICATION = FALSE). ",
            "Committing marker-score labels directly.")
  }


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

  saveRDS(obj, RDS_ANNOTATED)
  female_obj <- obj              
  message("  saved female_obj_annotated.rds ; Idents = celltype_marker")
}

# Cell types actually present, ordered by abundance (used across all downstream)
CT_LEVELS <- names(sort(table(as.character(female_obj$celltype_marker)), decreasing = TRUE))
CT_PALETTE <- build_celltype_palette(CT_LEVELS)


# =============================================================================
# SECTION 5 — MARKER DOT PLOT + PRIMARY-LABEL UMAP
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

message("\n=========  STAGE 02 complete  =========")
