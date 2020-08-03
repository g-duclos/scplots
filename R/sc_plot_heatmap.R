#' Function for overlaying gene expression data onto dimensionality reduction results with the option of outlining cells based on known phenotype assignments
#'

#' This function plots dimensionality reduction results and highlights cells based on known phenotype assignments
#' @param dataset ExpressionSet object
#' @param dimred 2-D output of t-SNE, UMAP, or PCA: Matrix with 2 columns corrersponding to dimensions 1 (x-axis) & 2 (y-axis), row numbers equivalent to the number of cells, and rownames corresponding to cell IDs
#' @param dimred.type Options include "tsne", "umap", or "pca"
#' @param dimred.type Options include "tsne", "umap", or "pca"
#' @param gene.id.type Options include "ID" or "ENS" (ensembl ID)
#' @param pheno.var.sub Character; Cells are highlighted according to specified phenotype variable
#' @param pheno.var.sub.target Logical or Character; FALSE if default order; Vector of pheno.var sub-groups to be plotted
#' @param pheno.var.bg Logical or Character; FALSE if no additional subsetting, otherwise specify secondary phenotype variable to subset from pheno.var pheno.order
#' @param pheno.var.bg.target Character; specify pheno.var.2 sub-group(s) to be visualized
#' @param pheno.var.bg.colors Character; Colors to highlight pheno.var sub-groups ordered according to pheno.order
#' @param ax.labels Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @param log Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @param max.val Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @param min.val Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @param title Logical or Character; FALSE if no title, otherwise specify title in quotes
#' @param size Numeric; point size, default is 6
#' @param bg.size.mult Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @param seed Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @return Scatter plot of dimensionality reduction results
#' @export
#' @examples
#' plot_dimred_pheno(dataset=dataset, dimred=dimred, type="umap, pheno.var=pheno.var, pheno.order=FALSE, pheno.colors=pheno.colors, pheno.var.2=FALSE, pheno.var.2.target=pheno.var.2.target, ax.labels=TRUE, title="Title", size=6)
#

sc_plot_heatmap <- function(
  dataset=dataset,
  heatmap.matrix=heatmap.matrix,
  cluster.var=cluster.var,
  cluster.colors = cluster.colors,
  cluster.order=cluster.order,
  pheno=FALSE,
  pheno.var=pheno.var,
  pheno.colors=pheno.colors,
  pheno.order=pheno.order,
  show.rownames=FALSE) {
  #
  if (cluster.order == FALSE) {
    cluster.sort = as.character(na.omit(unique(dataset[[cluster.var]])))
      } else {
        cluster.sort = as.character(cluster.order)
      }
  #
  if (pheno == TRUE) {
    if (pheno.order == FALSE) {
      pheno.sort = as.character(unique(dataset[[pheno.var]]))
    } else {
      pheno.sort = as.character(pheno.order)
    }
  }
  # Clusters in pData slot
  # Data.frame containing phenotype information for cells to be plotted
  col.annot = data.frame(row.names=unlist(lapply(cluster.sort, function(cluster) {
    return(dataset$ID[which(dataset[[cluster.var]] == cluster & dataset$Cell_Filter == "Cell")])
  })))
  #
  if (pheno == TRUE) {
    col.annot[[pheno.var]] <- unlist(lapply(cluster.sort, function(cluster) {
      return(dataset[[pheno.var]][which(dataset[[cluster.var]] == cluster & dataset$Cell_Filter == "Cell")])
    }))
  }
  #
  col.annot[[cluster.var]] <- unlist(lapply(cluster.sort, function(cluster) {
       return(dataset[[cluster.var]][which(dataset[[cluster.var]] == cluster & dataset$Cell_Filter == "Cell")])
  }))
  #
  if (pheno == TRUE) {
    col.colors[1] <- list(colorRampPalette(pheno.colors)(length(pheno.colors)))
    col.colors[2] <- list(colorRampPalette(cluster.colors)(length(cluster.colors)))
    names(col.colors) <- c(pheno.var, cluster.var)
    names(col.colors[[pheno.var]]) <- unique(pheno.sort)
    names(col.colors[[cluster.var]]) <- unique(cluster.sort)
    } else {
      col.colors[1] <- list(colorRampPalette(cluster.colors)(length(cluster.colors)))
      names(col.colors) <- cluster.var
      names(col.colors[[cluster.var]]) <- unique(cluster.sort)
      }
  #
  hmap.colors = colorRampPalette(c("white", "grey80", "grey80", "red"))(256)
  #
  pheatmap(mat=heatmap.matrix,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = col.annot, # specify gene categories in data.frame
           annotation_colors = col.colors, # specify colors corresponding to those categories
           color = hmap.colors,
           annotation_legend = FALSE, legend=FALSE,
           annotation_names_row = FALSE, annotation_names_col = FALSE,
           show_rownames = show.rownames, show_colnames = FALSE,
           fontsize = 13)
}
