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

sc_dimred_express <- function(
  dataset=dataset,
  dimred=dimred,
  dimred.type=dimred.type,
  gene.target=gene.target,
  gene.id.type="ID",
  pheno.var.sub=FALSE,
  pheno.var.sub.target=pheno.var.sub.target,
  pheno.var.bg=FALSE,
  pheno.var.bg.target=pheno.var.bg.target,
  pheno.var.bg.colors=pheno.var.bg.colors,
  log=FALSE,
  max.val=2,
  min.val=-2,
  ax.labels=TRUE,
  title=FALSE,
  size=6,
  bg.size.mult=8,
  seed=FALSE) {
  #
  if (ncol(dimred) != 2) {
    stop("Dimensionality reduction 'dimred' input must have 2 dimensions (2-column matrix)")
  }
  #
  point.size=(0.35*(size/6))
  #
  bg.point.size=bg.size.mult*point.size
  #
  if (pheno.var.sub == FALSE) {
    cells = dataset$ID[which(dataset$Cell_Filter == "Cell")]
    } else {
      cells = unlist(lapply(pheno.var.2.target, function(pheno) {
        return(dataset$ID[which(dataset[[pheno.var.2]] == pheno & dataset$Cell_Filter == "Cell")])
        }))
      }
  #
  if (gene.id.type == "ID") {
    ens.target = unlist(lapply(1:length(gene.target), function(x) { rownames(fData(dataset)[grep(paste("^", gene.target[x], "$",sep=""), fData(dataset)[,"ID"]),]) }))
    } else if (gene.id.type == "ENS") {
      ens.target = gene.target
      }
  #
  umis = dataset$UMIs
  names(umis) <- dataset$ID
  #
  if (length(gene.target) > 1) {
    #
    gene.tpm = sweep(exprs(dataset)[ens.target, cells], MARGIN=2, STATS=umis[cells]/1E6, FUN="/")
    colnames(gene.tpm) <- cells
    #
    } else if (length(gene.target == 1)) {
      #
      gene.tpm = as.numeric(sweep(t(exprs(dataset)[ens.target, cells]), MARGIN=2, STATS=umis[cells]/1E6, FUN="/"))
      names(gene.tpm) <- cells
      }
  #
  range = round(seq(min.val, max.val, 0.01), 2)
  #
  if (log == TRUE) {
    meta.gene = log2(gene.tpm+1)
    colors = colorRampPalette(c("grey80", "grey90", "lightpink", "pink1", "salmon", "indianred2", "red"))(401)
    } else {
      meta.gene = gene.tpm
      colors = colorRampPalette(c("grey80","grey80","red"))(401)
      }
  # Calculate MetaGene (mean Z-score) for each gene set
  if (length(gene.target) > 1) {
    #
    meta.gene = colMeans(t(scale(t(meta.gene))))
    } else if (length(gene.target) == 1 & log == FALSE) {
      #
      meta.gene = as.numeric(scale(gene.tpm))
      names(meta.gene) <- names(gene.tpm)
      }
  #
  meta.gene[meta.gene > max.val] <- max.val
  meta.gene[meta.gene < min.val] <- min.val

  # round z-normalized expression values to 2 decimal places
  norm.mat.r = unlist(list(round(meta.gene, 2)))
  names(norm.mat.r) <- names(meta.gene)

  # assign values from -2 to 2 to grey (low) to red (high) color spectrum
  color.values = as.numeric(as.character(factor(round(seq(min.val, max.val, 0.01),2))))
  names(color.values) <- colors

  # assign colors to gene expression values
  color.assign = names(color.values)[match(norm.mat.r, color.values)]
  names(color.assign) <- names(norm.mat.r)

  #
  if (seed != FALSE & is.numeric(seed) == TRUE) {
    set.seed(seed)
    }
  # randomize points for plotting
  color.assign = sample(color.assign)

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
          #
          x.lab=paste("PC", colnames(dimred)[1])
          y.lab=paste("PC", colnames(dimred)[2])
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
         xlim=c(min(dimred[,1]), max(dimred[,1])),
         ylim=c(min(dimred[,2]), max(dimred[,2])),
         cex.lab=1.3, ylab=y.lab, xlab=x.lab)
    #
    abline(v=axis(1), h=axis(2), lwd=0.25, col="grey80")
    #
    if (pheno.var.bg != FALSE) {
      #
      bg.cells = unlist(lapply(pheno.var.bg.target, function(pheno) {
        return(dataset$ID[which(dataset[[pheno.var.bg]] == pheno & dataset$Cell_Filter == "Cell")])
      }))
      #
      bg.colors = unlist(lapply(1:length(pheno.var.bg.target), function(pheno) {
        id = dataset$ID[which(dataset[[pheno.var.bg]] == pheno.var.bg.target[pheno])]
        id.colors = rep(adjustcolor(pheno.var.bg.colors[pheno], alpha.f = 0.03), length(id))
        names(id.colors) <- id
        #
        return(id.colors)
        #
      }))
      #
      bg.cells = intersect(bg.cells, cells)
      #
      bg.colors.sample = sample(bg.colors[bg.cells])
      #
      points(dimred[names(bg.colors.sample), 1], dimred[names(bg.colors.sample), 2],
             pch=21, cex=bg.point.size,
             bg=bg.colors.sample, col=NA)
      }
    #
    points(dimred[names(color.assign), 1], dimred[names(color.assign), 2],
      pch=21, cex=point.size,
      bg=color.assign, col=color.assign)
    #
    if (label == FALSE) {
      title("", line=0.5, cex.main=1.5)
      } else {
        title(label, line=0.5, cex.main=1.5)
        }
    #
    box(lwd=3)
    #
}
