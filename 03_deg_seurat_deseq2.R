# =============================================================================
# 03_deg_seurat_deseq2.R  —  NUDT16 KO female spleen scRNA-seq
# =============================================================================
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

message("\n=========  STAGE 03 — DIFFERENTIAL EXPRESSION  =========")

female_obj <- load_stage_rds(RDS_ANNOTATED, "02_cell_annotation.R")
CT_LEVELS  <- names(sort(table(as.character(female_obj$celltype_marker)), decreasing = TRUE))
CT_PALETTE <- build_celltype_palette(CT_LEVELS)

# =============================================================================
# DIFFERENTIAL EXPRESSION (KO vs WT)
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
    dds <- DESeqDataSetFromMatrix(countData = round(pb), colData = coldata, design = ~ genotype)
    
    dds <- DESeq(dds, sfType = "poscounts", fitType = "local", quiet = TRUE)
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
# VOLCANO PLOTS 
#   DESeq2  -> volcano/deseq2/   ;  Seurat -> volcano/seurat/
# =============================================================================
message("\n[7] Volcano plots (DESeq2 + Seurat, global + per cell type)...")

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
de_results <- list(
  de_global               = de_global,
  deg_global_filtered     = deg_global_filtered,
  de_by_celltype          = de_by_celltype,
  deg_by_celltype_filtered = deg_by_celltype_filtered,
  gene_universe           = gene_universe,
  seurat_de_global        = seurat_de_global,
  seurat_de_by_celltype   = seurat_de_by_celltype,
  CT_LEVELS               = CT_LEVELS
)
saveRDS(de_results, RDS_DE)
message("\n  saved de_results.rds  (", RDS_DE, ")")

message("\n=========  STAGE 03 complete  =========")
