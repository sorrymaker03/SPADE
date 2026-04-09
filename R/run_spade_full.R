#' Spatial Proximity Analysis of Differential Expression (SPADE)
#'
#' This function analyzes gene expression differences between cells near vs far
#' from a reference cell type in spatial transcriptomics data.
#'
#' @param spe A Seurat object containing spatial transcriptomics data.
#' @param cell_a_list List of cell type names for the query cells.
#' @param cell_b_list List of cell type names for the reference cells.
#' @param celltype_col Column name in metadata containing cell type annotations.
#' @param output_base Base directory for output files.
#' @param near_dist Distance threshold (pixels/um) defining "near" cells.
#' @param far_dist Distance threshold (pixels/um) defining "far" cells.
#' @param spatial_width Width of output spatial plots (inches).
#' @param spatial_height Height of output spatial plots (inches).
#' @param pt.size.factor Point size factor for spatial plots.
#'
#' @return Invisibly returns a list of results for each comparison.
#' @export
#'
#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import dplyr
#' @import FNN
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import scales
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' run_spade_full(
#'     spe = seurat_object,
#'     cell_a_list = list(c("Endo")),
#'     cell_b_list = list(c("Mono")),
#'     celltype_col = "predicted.id",
#'     output_base = "results/",
#'     near_dist = 60,
#'     far_dist = 200
#' )
#' }
run_spade_full <- function(
    spe,
    cell_a_list,
    cell_b_list,
    celltype_col,
    output_base,
    near_dist,
    far_dist,
    spatial_width = 10,
    spatial_height = 10,
    pt.size.factor = 1
) {


  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }

  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("FNN package is required")
  }


  coords <- tryCatch({
    Seurat::GetTissueCoordinates(spe)
  }, error = function(e) {
    spe@images[[1]]@coordinates
  })
  coords <- as.data.frame(coords)


  if(!"x" %in% colnames(coords)) {
    colnames(coords)[1:2] <- c("x", "y")
  }

  meta <- spe@meta.data
  if(!celltype_col %in% colnames(meta)) {
    stop("celltype_col not found in metadata")
  }

  celltype <- meta[[celltype_col]]
  names(celltype) <- colnames(spe)

  results_list <- list()

  for(cell_a in cell_a_list){
    for(cell_b in cell_b_list){

      tag <- paste0(paste(cell_a, collapse = "_"), "_vs_", paste(cell_b, collapse = "_"))
      outdir <- file.path(output_base, tag)
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

      cells_a <- names(celltype)[celltype %in% cell_a]
      cells_b <- names(celltype)[celltype %in% cell_b]

      if(length(cells_a) < 3 | length(cells_b) < 3) next

      coords_a <- coords[cells_a, c("x","y"), drop = FALSE]
      coords_b <- coords[cells_b, c("x","y"), drop = FALSE]


      dist_b <- FNN::get.knnx(data = coords_b, query = coords_a, k = 1)$nn.dist

      spe$dist_b <- NA
      spe$dist_b[cells_a] <- dist_b

      spe$group <- NA
      spe$group[cells_a] <- ifelse(dist_b < near_dist, "near",
                                   ifelse(dist_b > far_dist, "far", "mid"))

      obj <- subset(spe, cells = cells_a)
      Seurat::Idents(obj) <- obj$group

      group_counts <- table(obj$group)
      n_near <- ifelse("near" %in% names(group_counts), group_counts["near"], 0)
      n_far  <- ifelse("far"  %in% names(group_counts), group_counts["far"], 0)

      if (n_near < 3 | n_far < 3) next


      markers <- Seurat::FindMarkers(obj, ident.1 = "near", ident.2 = "far",
                                     logfc.threshold = 0.25, min.pct = 0.1)


      ctrl_cells <- names(celltype)[celltype %in% unlist(cell_b_list)]

      near_cells <- colnames(obj)[obj$group == "near"]
      far_cells  <- colnames(obj)[obj$group == "far"]

      coords_spe <- coords[colnames(spe), c("x","y"), drop = FALSE]
      coords_near <- coords_spe[near_cells, , drop = FALSE]
      coords_ctrl <- coords_spe[ctrl_cells, , drop = FALSE]

      if(nrow(coords_near) > 0 && nrow(coords_ctrl) > 0) {
        dist_ctrl_to_far <- FNN::get.knnx(
          data = coords_near,
          query = coords_ctrl,
          k = 1
        )$nn.dist

        ctrl_cells <- ctrl_cells[dist_ctrl_to_far < near_dist]
      }

      if(length(far_cells) >= 3 & length(ctrl_cells) >= 3){
        obj_ctrl <- subset(spe, cells = c(far_cells, ctrl_cells))
        obj_ctrl$group_for_filter <- ifelse(colnames(obj_ctrl) %in% far_cells, "far", "ctrl")
        Seurat::Idents(obj_ctrl) <- obj_ctrl$group_for_filter

        ctrl_markers <- Seurat::FindMarkers(obj_ctrl, ident.1 = "ctrl", ident.2 = "far",
                                            logfc.threshold = 0.25, min.pct = 0.1)

        markers_filtered <- markers[!rownames(markers) %in%
                                      rownames(ctrl_markers[ctrl_markers$avg_log2FC > 0,]), ]
      } else {
        markers_filtered <- markers
      }


      write.csv(markers_filtered, file = file.path(outdir, paste0(tag, "_markers_filtered.csv")))


      p1 <- Seurat::SpatialFeaturePlot(obj, features = "dist_b",
                                       pt.size.factor = pt.size.factor,
                                       image.alpha = 0.1) +
        ggplot2::scale_fill_gradientn(
          colors = c("#440154", "#3B528B", "#21908C", "#5DC963", "#FDE725"),
          limits = c(0, far_dist),
          oob = scales::squish
        ) +
        ggplot2::coord_fixed()

      ggplot2::ggsave(file.path(outdir, paste0(tag, "_spatial_distance.pdf")), p1,
                      width = spatial_width, height = spatial_height, units = "in", dpi = 300)

      cols <- c("near"= "#D62728","mid"="grey80","far"="#1F77B4")
      p2 <- Seurat::SpatialDimPlot(obj, group.by = "group",
                                   pt.size.factor = pt.size.factor,
                                   cols = cols,
                                   image.alpha = 0.1) +
        ggplot2::coord_fixed()

      ggplot2::ggsave(file.path(outdir, paste0(tag, "_spatial_group.pdf")), p2,
                      width = spatial_width, height = spatial_height, units = "in", dpi = 300)


      if(requireNamespace("clusterProfiler", quietly = TRUE) &&
         requireNamespace("org.Hs.eg.db", quietly = TRUE)) {

        res <- markers_filtered
        res$gene <- rownames(res)

        res <- dplyr::mutate(res,
                             log10_p = -log10(p_val_adj + 1e-300),
                             change = dplyr::case_when(
                               p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up",
                               p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
                               TRUE ~ "NS"
                             ))

        up_genes <- res %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% dplyr::pull(gene)
        down_genes <- res %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>% dplyr::pull(gene)

        if(length(up_genes) > 0) {
          up_df <- clusterProfiler::bitr(up_genes, fromType="SYMBOL", toType="ENTREZID",
                                         OrgDb=org.Hs.eg.db)
          ego_up <- clusterProfiler::enrichGO(up_df$ENTREZID, OrgDb=org.Hs.eg.db,
                                              ont="BP", readable=TRUE)
        } else {
          ego_up <- data.frame()
        }

        if(length(down_genes) > 0) {
          down_df <- clusterProfiler::bitr(down_genes, fromType="SYMBOL", toType="ENTREZID",
                                           OrgDb=org.Hs.eg.db)
          ego_down <- clusterProfiler::enrichGO(down_df$ENTREZID, OrgDb=org.Hs.eg.db,
                                                ont="BP", readable=TRUE)
        } else {
          ego_down <- data.frame()
        }

        ego_all <- rbind(
          as.data.frame(ego_up) %>% dplyr::mutate(direction="Up"),
          as.data.frame(ego_down) %>% dplyr::mutate(direction="Down")
        )

        if (!is.null(ego_all) && nrow(ego_all) > 0) {
          write.csv(ego_all, file = file.path(outdir, paste0(tag, "_GO_results.csv")))
        }

        # GSEA
        gene_list <- res$avg_log2FC
        names(gene_list) <- res$gene
        gene_list <- gene_list[!is.na(names(gene_list))]

        if(length(gene_list) > 0) {
          gene_df <- tryCatch({
            clusterProfiler::bitr(names(gene_list),
                                  fromType="SYMBOL", toType="ENTREZID",
                                  OrgDb=org.Hs.eg.db)
          }, error = function(e) data.frame())

          if(nrow(gene_df) > 0) {
            gene_list_entrez <- gene_list[gene_df$SYMBOL]
            names(gene_list_entrez) <- gene_df$ENTREZID
            gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

            gsea_go <- clusterProfiler::gseGO(
              geneList = gene_list_entrez,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              verbose = FALSE
            )

            if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
              write.csv(as.data.frame(gsea_go),
                        file = file.path(outdir, paste0(tag, "_GSEA_results.csv")))
            }
          }
        }
      }


      results_list[[tag]] <- list(
        markers = markers_filtered,
        near_cells = near_cells,
        far_cells = far_cells,
        output_dir = outdir
      )
    }
  }

  invisible(results_list)
}
