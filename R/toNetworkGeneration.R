build_network_by_group_generalized <- function(
    seurat_obj,
    group_col,          # Metadata column name for grouping, e.g., "WT_DSPKD"
    group_value,        # Group value, e.g., "WT" or "KO"
    network_obj,        # Network object (contains Symbol.1, Symbol.2, MI_WT / MI_KO etc.)
    random_network_obj, # Random network object (contains source,target,MI_WT / MI_KO etc.)
    housekeeping_gmt,   # Path to housekeeping genes file
    cluster_col = NULL, # If clustering column exists, can specify
    run_clustering = TRUE,
    resolution = 0.5,
    dims = 1:20,
    top_n_marker = 300,
    save_prefix = "./data_save"
) {

  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(purrr)
  library(COSG)
  library(GSEABase)

  #----------------------------------------
  message(">>> Step 1.1: Subset Seurat object and preprocess")
  seurat_sub <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_col]] %in% group_value, ]))

  if (run_clustering) {
    seurat_sub <- NormalizeData(seurat_sub)
    seurat_sub <- FindVariableFeatures(seurat_sub)
    seurat_sub <- ScaleData(seurat_sub)
    seurat_sub <- RunPCA(seurat_sub, features = VariableFeatures(seurat_sub))
    seurat_sub <- FindNeighbors(seurat_sub, dims = dims)
    seurat_sub <- FindClusters(seurat_sub, resolution = resolution)
    seurat_sub <- RunUMAP(seurat_sub, dims = dims)
    seurat_sub$clusters <- seurat_sub$seurat_clusters
  } else if (!is.null(cluster_col)) {
    seurat_sub$clusters <- seurat_sub@meta.data[[cluster_col]]
  } else {
    stop("Please provide cluster_col or set run_clustering=TRUE")
  }

  #----------------------------------------
  message(">>> Step 1.2: Prepare network column mapping")
  suffix <- toupper(group_value) # Automatically identify WT / KO
  mi_col <- paste0("MI_", suffix)
  d1_col <- paste0("DREMI_1_", suffix)
  d2_col <- paste0("DREMI_2_", suffix)

  if (!all(c(mi_col, d1_col, d2_col) %in% colnames(network_obj))) {
    stop(paste0("Error: Cannot find columns ", mi_col, ", ", d1_col, ", ", d2_col, " in network_obj"))
  }

  #----------------------------------------
  message(">>> Step 1.3: Build main network")
  Network_x <- network_obj %>%
    dplyr::select(Symbol.1, Symbol.2,
                  !!sym(mi_col), !!sym(d1_col), !!sym(d2_col)) %>%
    dplyr::mutate(DREMI = pmax(!!sym(d1_col), !!sym(d2_col))) %>%
    dplyr::rename(MI = !!sym(mi_col)) %>%
    dplyr::select(Symbol.1, Symbol.2, MI, DREMI) %>%
    dplyr::mutate(edge = paste0(Symbol.1, "_", Symbol.2),
                  random = "true")

  dir.create(save_prefix, showWarnings = FALSE, recursive = TRUE)


  #----------------------------------------
  message(">>> Step 1.4: Build random network")
  random_Network <- random_network_obj %>%
    dplyr::select(source, target,
                  !!sym(mi_col), !!sym(d1_col), !!sym(d2_col)) %>%
    dplyr::rename(Symbol.1 = source, Symbol.2 = target) %>%
    dplyr::mutate(DREMI = pmax(!!sym(d1_col), !!sym(d2_col)),
                  MI = !!sym(mi_col)) %>%
    dplyr::select(Symbol.1, Symbol.2, MI, DREMI) %>%
    dplyr::mutate(edge = paste0(Symbol.1, "_", Symbol.2))

  gene_map <- data.frame(gene = unique(c(Network_x$Symbol.1, Network_x$Symbol.2)))
  gene_map$index <- seq_len(nrow(gene_map)) - 1

  random_clean <- random_Network %>%
    filter(Symbol.1 %in% gene_map$gene & Symbol.2 %in% gene_map$gene)
  dup <- intersect(Network_x$edge, random_clean$edge)
  random_clean <- random_clean[!random_clean$edge %in% dup, ]
  random_clean$random <- "random"

  Network_x <- rbind(Network_x, random_clean)

  #----------------------------------------
  # message(">>> Step 1.5: Sigmoid transform")
  # sigmold_x <- function(x, mu, k) 1 / (1 + exp(-k * (x - mu)))
  # a1 <- median(log(Network_x$MI), na.rm = TRUE)
  # Network_x$sigmold_MI <- sigmold_x(log(Network_x$MI), mu = a1, k = 0.8)
  # a2 <- median(log(Network_x$DREMI), na.rm = TRUE)
  # Network_x$sigmold_DREMI <- sigmold_x(log(Network_x$DREMI), mu = a2, k = 0.8)
  message(">>> Step 5: Sigmoid transform")
  sigmold_x <- function(x, mu, k) 1 / (1 + exp(-k * (x - mu)))
  # MI
  mi_log <- log(Network_x$MI)
  a1 <- median(mi_log, na.rm = TRUE)
  k1 <- 2 / IQR(mi_log, na.rm = TRUE)
  Network_x$sigmold_MI <- sigmold_x(mi_log, mu = a1, k = k1)
  # DREMI
  dre_log <- log(Network_x$DREMI)
  a2 <- median(dre_log, na.rm = TRUE)
  k2 <- 2 / IQR(dre_log, na.rm = TRUE)
  Network_x$sigmold_DREMI <- sigmold_x(dre_log, mu = a2, k = k2)

  save(Network_x, file = file.path(save_prefix, paste0("Network_", group_value, ".Rdata")))
  #----------------------------------------
  message(">>> Step 1.6: Identify COSG markers")
  Idents(seurat_sub) <- seurat_sub$clusters
  COSG_markers <- cosg(seurat_sub, groups = 'all', assay = 'RNA', slot = 'data', mu = 1, n_genes_user = top_n_marker)
  celltype_marker <- COSG_markers[["names"]]

  #----------------------------------------
  message(">>> Step 1.7: Integrate housekeeping genes")
  gene_sets <- getGmt(housekeeping_gmt)
  House_keep_gene <- geneIds(gene_sets[[1]])
  cell_types <- colnames(celltype_marker)
  positive_marker_list <- vector("list", length(cell_types))
  names(positive_marker_list) <- gsub("[ +]", "_", cell_types)

  for (ct in cell_types) {
    marker_gene <- na.omit(celltype_marker[, ct])
    obj_ct <- subset(seurat_sub, clusters %in% ct)
    expr_matrix <- GetAssayData(obj_ct, slot = "data")
    gene_pct <- Matrix::rowMeans(expr_matrix > 0)
    gene_mean <- Matrix::rowMeans(expr_matrix)
    gene_stats <- data.frame(gene = rownames(expr_matrix),
                             pct_cells_expressed = gene_pct,
                             mean_expression = gene_mean)
    hk_stats <- gene_stats[House_keep_gene, ] %>% na.omit()
    hk_good <- hk_stats %>%
      filter(pct_cells_expressed > 0.75, mean_expression > 1) %>%
      pull(gene)
    positive_marker_list[[gsub("[ +]", "_", ct)]] <- unique(c(marker_gene, hk_good))
  }

  #----------------------------------------
  message(">>> Step 1.8: Build subnetwork per cell type")
  network_sub_list <- vector("list", length(positive_marker_list))
  names(network_sub_list) <- names(positive_marker_list)
  for (i in seq_along(positive_marker_list)) {
    ct_name <- names(positive_marker_list)[i]
    markers <- positive_marker_list[[i]]
    sub_net <- Network_x[Network_x$Symbol.1 %in% markers & Network_x$Symbol.2 %in% markers, ]
    sub_net$cell_type <- ct_name
    network_sub_list[[ct_name]] <- sub_net
  }
  network_all <- imap_dfr(network_sub_list, ~mutate(.x, cell_type = .y))
  Positive_Network_all <- subset(network_all, random == "true")

  save(Positive_Network_all, file = file.path(save_prefix, paste0("Positive_Network_all_", group_value, ".Rdata")))
  #----------------------------------------
  message(">>> Step 1.9: Compute mean expression per cell type")
  network_list <- vector("list", length(cell_types))
  names(network_list) <- gsub("[ +]", "_", cell_types)
  for (ct in cell_types) {
    obj_ct <- subset(seurat_sub, clusters == ct)
    expr_matrix <- GetAssayData(obj_ct, slot = "data")
    gene_stats <- data.frame(gene = rownames(expr_matrix),
                             mean_expression = Matrix::rowMeans(expr_matrix))
    gene2expr <- setNames(gene_stats$mean_expression, gene_stats$gene)
    net_tmp <- Network_x
    expr1 <- gene2expr[net_tmp$Symbol.1]
    expr2 <- gene2expr[net_tmp$Symbol.2]
    net_tmp$geom_mean_expr <- sqrt(expr1 * expr2)
    net_tmp$celltype <- ct
    network_list[[gsub("[ +]", "_", ct)]] <- net_tmp
  }
  network_all_celltype <- imap_dfr(network_list, ~mutate(.x, celltype = .y))
  network_all_celltype$cluster <- paste0("cluster_", network_all_celltype$celltype)
  network_all_celltype$cluster <- factor(network_all_celltype$cluster)
  one_hot <- model.matrix(~ cluster - 1, data = network_all_celltype)
  network_all_celltype <- cbind(network_all_celltype, one_hot)

  #----------------------------------------
  message(">>> Step 1.10: Save and return results")
  save(network_all_celltype, file = file.path(save_prefix, paste0("network_all_celltype_", group_value, ".Rdata")))

  return(list(
    seurat_obj = seurat_sub,
    Network = Network_x,
    Positive_Network_all = Positive_Network_all,
    network_all_celltype = network_all_celltype
  ))
}

#' run_network_pipeline
#'
#' @param seurat_obj seurat_obj
#' @param group_col
#' @param group_values
#' @param network_obj
#' @param random_network_obj
#' @param housekeeping_gmt
#' @param output_dir
#' @param run_clustering
#' @param cluster_col
#' @param dims
#' @param resolution
#' @param top_n_marker
#' @param negative_sample_size
#' @param save_intermediate
#'
#' @return
#' @export
#'
#' @examples
run_network_pipeline <- function(
    seurat_obj,
    group_col,
    group_values = c("WT", "KO"),
    network_obj,
    random_network_obj,
    housekeeping_gmt,
    output_dir = "./output",
    run_clustering = TRUE,
    cluster_col = NULL,
    dims = 1:20,
    resolution = 0.5,
    top_n_marker = 300,
    negative_sample_size = 50000,
    save_intermediate = TRUE
) {
  library(dplyr)
  library(purrr)
  library(Seurat)
  library(Matrix)
  library(GSEABase)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("========== NETWORK PIPELINE START ==========")
  message("Groups: ", paste(group_values, collapse = ", "))

  results_list <- list()

  # ==========================================================
  # 🧩 Step 1: Build individual group networks
  # ==========================================================
  for (grp in group_values) {
    message(">>> Running group: ", grp)

    save_prefix <- file.path(output_dir, grp)
    dir.create(save_prefix, recursive = TRUE, showWarnings = FALSE)

    #⚙️ Suppress warnings for this step
    result <- suppressWarnings({
      build_network_by_group_generalized(
        seurat_obj = seurat_obj,
        group_col = group_col,
        group_value = grp,
        network_obj = network_obj,
        random_network_obj = random_network_obj,
        housekeeping_gmt = housekeeping_gmt,
        run_clustering = run_clustering,
        cluster_col = cluster_col,
        dims = dims,
        resolution = resolution,
        top_n_marker = top_n_marker,
        save_prefix = save_prefix
      )
    })

    # Save individual group results
    saveRDS(result, file = file.path(save_prefix, paste0(grp, "_network_result.rds")))
    results_list[[grp]] <- result
  }

  # Save all group results
  saveRDS(results_list, file = file.path(output_dir, "all_groups_network_results.rds"))

  # ==========================================================
  # 🧩 Step 2–5：Merge + Generate training data
  # ==========================================================
  suppressWarnings({
    # ------------------------------
    # Step 2: Merge WT/KO networks
    # ------------------------------
    message(">>> Step 2: Merge WT/KO networks")

    network_all_celltype <- bind_rows(
      results_list[[group_values[1]]]$network_all_celltype[, 1:11] %>%
        mutate(stim = group_values[1]),
      results_list[[group_values[2]]]$network_all_celltype[, 1:11] %>%
        mutate(stim = group_values[2])
    )

    Positive_Network_all <- bind_rows(
      results_list[[group_values[1]]]$Positive_Network_all %>%
        mutate(stim = group_values[1]),
      results_list[[group_values[2]]]$Positive_Network_all %>%
        mutate(stim = group_values[2])
    )

    # ------------------------------
    # Step 3: Gene index mapping
    # ------------------------------
    message(">>> Step 3: Assign gene indices")

    gene_map <- data.frame(gene = unique(c(
      network_all_celltype$Symbol.1,
      network_all_celltype$Symbol.2
    )))
    gene_map$index <- seq_len(nrow(gene_map)) - 1
    colnames(gene_map) <- c("gene", "index")

    # Create export directory
    export_dir <- file.path(output_dir, "data_export")
    dir.create(export_dir, showWarnings = FALSE)

    orignal_network <- network_obj[, c("Symbol.1", "Symbol.2")]
    orignal_network$index1 <- gene_map$index[match(orignal_network$Symbol.1, gene_map$gene)]
    orignal_network$index2 <- gene_map$index[match(orignal_network$Symbol.2, gene_map$gene)]

    # 💾 Save gene_map
    write.table(
      orignal_network[, c("index1", "index2")],
      file = file.path(export_dir, "gene_map.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Add index information
    network_all_celltype$index1 <- gene_map$index[match(network_all_celltype$Symbol.1, gene_map$gene)]
    network_all_celltype$index2 <- gene_map$index[match(network_all_celltype$Symbol.2, gene_map$gene)]
    network_all_celltype$index  <- seq_len(nrow(network_all_celltype)) - 1

    network_all_celltype$cell_index <- paste0(
      network_all_celltype$stim, "_",
      network_all_celltype$Symbol.1, "_",
      network_all_celltype$Symbol.2, "_",
      network_all_celltype$cluster
    )

    network_all_celltype$celltype_stim <- paste0(network_all_celltype$stim, "_", network_all_celltype$cluster)

    # Ensure cluster is factor (not character)
    network_all_celltype$celltype_stim <- factor(network_all_celltype$celltype_stim)

    # Generate one-hot matrix (remove intercept)
    one_hot <- model.matrix(~ celltype_stim - 1, data = network_all_celltype)
    # Bind one-hot encoding back to original dataframe
    network_all_celltype <- cbind(network_all_celltype, one_hot)

    # save(network_all_celltype, file = "./network_all_celltype.Rdata")
    # ------------------------------
    # Step 4:  Export training data
    # ------------------------------
    message(">>> Step 4: Export training data")

    write.table(
      network_all_celltype[, c("index1", "index2")],
      file = file.path(export_dir, "index_pairs_all.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    cols <- c(
      grep("^celltype_stim", names(network_all_celltype), value = TRUE),
      "sigmold_MI", "sigmold_DREMI", "geom_mean_expr"
    )
    cols <- setdiff(cols, "celltype_stim")
    write.table(
      network_all_celltype[, cols],
      file = file.path(export_dir, "edge_feature_all.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )



    # ------------------------------
    # Step 5: Positive / Negative edges
    # ------------------------------
    message(">>> Step 5: Generate Positive / Negative edges")

    Positive_Network_all$cluster <- paste0("cluster_", Positive_Network_all$cell_type)
    Positive_Network_all$cell_index <- paste0(
      Positive_Network_all$stim, "_",
      Positive_Network_all$Symbol.1, "_",
      Positive_Network_all$Symbol.2, "_",
      Positive_Network_all$cluster
    )

    Positive_Network_all1 <- network_all_celltype %>%
      filter(cell_index %in% Positive_Network_all$cell_index)

    write.table(
      Positive_Network_all1[, c("index1", "index2", "index")],
      file = file.path(export_dir, "Positive_edge.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Negative samples: random + no expression
    network_all_celltype_random <- network_all_celltype[network_all_celltype$random == "random", ]
    network_all_celltype_true   <- network_all_celltype[network_all_celltype$random == "true", ]

    Negative_network_all_celltype1 <- network_all_celltype_true[
      network_all_celltype_true$geom_mean_expr == 0, ]
    Negative_network_all_celltype2 <- network_all_celltype_random

    Negative_network_all_celltype <- rbind(
      Negative_network_all_celltype1,
      Negative_network_all_celltype2
    )

    # Random sampling of negative samples
    set.seed(123)
    Negative_network_all_celltype <- Negative_network_all_celltype %>%
      sample_n(min(negative_sample_size, nrow(network_all_celltype_random)))

    # Save negative samples
    write.table(
      Negative_network_all_celltype[, c("index1", "index2", "index")],
      file = file.path(export_dir, "Negative_edge.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # Save true network index pairs
    write.table(
      network_all_celltype_true[, c("index1", "index2")],
      file = file.path(export_dir, "index_pairs_true.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    write.table(
      network_all_celltype_true[, cols],
      file = file.path(export_dir, "edge_feature_true.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    write.table(
      network_all_celltype_true[, c("edge", "celltype_stim")],
      file = file.path(export_dir, "symbol_pair.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
  })

  message("========== PIPELINE FINISHED ==========")
}
