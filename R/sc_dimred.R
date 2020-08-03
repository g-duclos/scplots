#' Function for plotting dimensionality reduction results and highlighting cells based on known phenotype assignments
#'

#' This function plots dimensionality reduction results and highlights cells based on known phenotype assignments
#' @param dataset ExpressionSet object
#' @param dimred 2-D output of t-SNE, UMAP, or PCA: Matrix with 2 columns corrersponding to dimensions 1 (x-axis) & 2 (y-axis), row numbers equivalent to the number of cells, and rownames corresponding to cell IDs
#' @param dimred.type Options include "tsne", "umap", or "pca"
#' @param pheno.var Character; Cells are highlighted according to specified phenotype variable
#' @param pheno.order Logical or Character; FALSE if default order; Vector of pheno.var sub-groups to be plotted
#' @param pheno.colors Character; Colors to highlight pheno.var sub-groups ordered according to pheno.order
#' @param pheno.var.2 Logical or Character; FALSE if no additional subsetting, otherwise specify secondary phenotype variable to subset from pheno.var pheno.order
#' @param pheno.var.2.target Character; specify pheno.var.2 sub-group(s) to be visualized
#' @param ax.labels Logical; RUE if x- and y-axis labeled according to specified type, FALSE if no axis labels
#' @param title Logical or Character; FALSE if no title, otherwise specify title in quotes
#' @param size Numeric; point size, default is 6
#' @return Scatter plot of dimensionality reduction results
#' @export
#' @examples
#' sc_dimred(dataset=dataset, dimred=dimred, dimred.type="umap", pheno.var=pheno.var, pheno.order=FALSE, pheno.colors=pheno.colors, pheno.var.2=FALSE, pheno.var.2.target=pheno.var.2.target, ax.labels=TRUE, title="Title", size=6)
#

sc_dimred <- function(
  dataset=dataset,
  dimred=dimred,
  dimred.type=dimred.type,
  pheno.var=pheno.var,
  pheno.order=FALSE,
  pheno.colors=pheno.colors,
  pheno.var.2=FALSE,
  pheno.var.2.target=pheno.var.2.target,
  ax.labels=TRUE,
  title=FALSE,
  size=6) {
  #
  if (ncol(dimred) != 2) {
    stop("Dimensionality reduction 'dimred' input must have 2 dimensions (2-column matrix)")
  }
  #
  point.size=(0.35*(size/6))
  #
  if (pheno.var.2 == FALSE) {
    cells = dataset$ID[which(dataset$Cell_Filter == "Cell")]
    } else {
      cells = unlist(lapply(pheno.var.2.target, function(pheno) {
      return(dataset$ID[which(dataset[[pheno.var.2]] == pheno & dataset$Cell_Filter == "Cell")])
    }))
  }
  #
  if (pheno.order == FALSE) {
    pheno.order = as.character(na.omit(unique(dataset[[pheno.var]])))
  }
  #
  cell.colors = unlist(lapply(1:length(pheno.order), function(i) {
    id = intersect(cells, dataset$ID[which(dataset[[pheno.var]] == pheno.order[i])])
    id.colors = rep(pheno.colors[i], length(id))
    names(id.colors) <- id
    #
    return(id.colors)
    }))

  # randomize the cells
  cell.colors.sample <- sample(cell.colors)
  #
  if (ax.labels == FALSE) {
    x.lab=""
    y.lab=""
    }
  #
  if (ax.labels == TRUE) {
    if (dimred.type == "tsne") {
      x.lab="t-SNE Dim 1"
      y.lab="t-SNE Dim 2"
    } else if (dimred.type == "umap") {
      x.lab="UMAP Dim 1"
      y.lab="UMAP Dim 2"
      } else if (dimred.type == "pca") {
        # Edit this to generalize
        x.lab=paste(unlist(strsplit(colnames(dimred)[1],"_")), collapse=" ")
        y.lab=paste(unlist(strsplit(colnames(dimred)[2],"_")), collapse=" ")
      }
  }
  #
  if (ax.labels == TRUE) {
    par(mar=c(3,3,2,1), mgp=c(1.8, 0.6, 0), las=1)
    } else {
      par(mar=c(2,2,2,1), mgp=c(1.8, 0.6, 0), las=1)
      }
  #
  plot(NA, NA,
    xlim = c(min(dimred[cells,1]), max(dimred[cells,1])),
    ylim = c(min(dimred[cells,2]), max(dimred[cells,2])),
    cex.lab=1.3, xlab=x.lab, ylab=y.lab)
  #
  abline(v=axis(1), h=axis(2), lwd=0.25, col="grey80")
  #
  points(dimred[names(cell.colors.sample),1], dimred[names(cell.colors.sample),2],
    pch=21, cex=point.size,
    bg=cell.colors.sample, col=cell.colors.sample)
  #
  if (title == FALSE) {
    title("", line=0.5, cex.main=1.5)
    } else {
      title(title, line=0.5, cex.main=1.5)
      }
  #
  box(lwd=3)
}
