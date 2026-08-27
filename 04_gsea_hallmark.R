# =============================================================================
# 04_gsea_hallmark.R  —  NUDT16 KO female spleen scRNA-seq
# =============================================================================

if (!isTRUE(getOption("nudt16.config.loaded"))) {
  cand <- c(
    "00_config.R",
    file.path("pipeline", "00_config.R"),
    file.path("results", "females", "pipeline", "00_config.R"),
    "/00_config.R"
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Cannot locate 00_config.R — run from the pipeline folder.")
  source(hit[1])
}

message("\n=========  STAGE 04 — FUNCTIONAL ENRICHMENT (GSEA + Hallmark)  =========")

de <- load_stage_rds(RDS_DE, "03_deg_seurat_deseq2.R")
de_global                <- de$de_global
deg_global_filtered      <- de$deg_global_filtered
de_by_celltype           <- de$de_by_celltype
deg_by_celltype_filtered <- de$deg_by_celltype_filtered
gene_universe            <- de$gene_universe
seurat_de_global         <- de$seurat_de_global
seurat_de_by_celltype    <- de$seurat_de_by_celltype

# ---- ORA helper 
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

# =============================================================================
# FUNCTIONAL ENRICHMENT
# =============================================================================
message("\n[9] Functional enrichment (pseudobulk ORA + Hallmark GSEA/ORA)...")

# ---- 9A. Over-representation on pseudobulk DEGs -------
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
  if (!is.null(deg_global_filtered) && nrow(deg_global_filtered)) run_dir(deg_global_filtered, "Global")
  # Per cell type
  for (ct in names(deg_by_celltype_filtered)) {
    local({ ctt <- ct; run_dir(deg_by_celltype_filtered[[ctt]], ctt) })
  }
  message("  enrichment: pseudobulk ORA written to enrichment/")
})

# ---- 9B. MSigDB Hallmark GSEA  --------
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
                    
  groups <- c(list(Global = seurat_de_global), seurat_de_by_celltype)
  for (label in names(groups)) {
    local({
      lab <- label; de_grp <- groups[[lab]]
      safe_panel(paste0("hallmark_", lab), {
        run_gsea_group(de_grp, lab); run_ora_group(de_grp, lab)
      })
    })
  }
  message("  enrichment_hallmark: GSEA + directional ORA for ", length(groups), " groups")
})

message("\n=========  STAGE 04 complete  =========")
