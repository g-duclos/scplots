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

sc_build_heatmap <- function(
  dataset=dataset,
  cluster.var=cluster.var,
  cluster.order=FALSE,
  gene.sets=gene.sets,
  col.space=200,
  row.space=1,
  rank=25) {
  #
  if (cluster.order == FALSE) {
    cluster.sort = as.character(na.omit(unique(dataset[[cluster.var]])))
      } else {
        cluster.sort = as.character(cluster.order)
        }
  #
  clusters = unlist(lapply(cluster.sort, function(cluster) {
    cluster.tmp == dataset[[cluster.var]][which(dataset[[cluster.var]] == cluster & dataset$Cell_Filter == "Cell")]
    names(cluster.tmp) <-  dataset$ID[which(dataset[[cluster.var]] == cluster & dataset$Cell_Filter == "Cell")]
    #
    return(cluster.tmp)
    #
  }))
  #
  umis = dataset$UMIs[which(dataset$Cell_Filter == "Cell")]
  names(umis) = dataset$ID[which(dataset$Cell_Filter == "Cell")]

  # Gene Sets in fData slot
  gene.sets <- fData(dataset)[[gene.var]]
  names(gene.sets) <- rownames(fData(dataset))
  # Remove unnassigned genes
  gene.sets <- na.omit(gene.sets)

  # Mean Cluster Z score mat
  #cluster.res = assayData(dataset)$Secondary[["Mean_Cluster_Z"]][[cluster.var]]
  cluster.res = sapply(cluster.sort, function(cluster) {
    if (length(names(clusters)[which(clusters == cluster)]) > 1) {
      return(rowMeans(tpm.z[names(gene.sets), names(clusters)[which(clusters == cluster)]]))
      } else {
        return(tpm.z[names(gene.sets), names(clusters)[which(clusters == cluster)]])
      }
  })

  # select top X expressed genes
  gene.sets.rank = unlist(lapply(cluster.sort, function(cluster) {
    if (rank > 0 & length(names(gene.sets)[which(gene.sets == cluster)]) > 1) {
      return(na.omit(gene.sets[names(sort(cluster.res[names(gene.sets)[which(gene.sets == cluster)], cluster], decreasing=TRUE)[1:rank])]))
      #
      } else if (rank == 0 & length(names(gene.sets)[which(gene.sets == cluster)]) > 1) {
        return(na.omit(gene.sets[names(sort(cluster.res[names(gene.sets)[which(gene.sets == cluster)], cluster], decreasing=TRUE))]))
        #
        } else if (length(names(gene.sets)[which(gene.sets == cluster)]) > 1) {
          return(gene.sets[names(gene.sets)[which(gene.sets == cluster)]])
          #
          } else if (length(names(gene.sets)[which(gene.sets == cluster)]) == 0) {
            cat(paste("No genes assigned to:", cluster), "\n")
            }
  }))
  #
  gene.sets.sort = unique(gene.sets.rank)

  # Subset counts matrix based on cell and gene input
  counts = exprs(dataset)[names(gene.sets.rank), names(clusters)]

  # Convert to TPM - select only 'expressed' genes & targeted cells
  tpm.z = t(scale(t(sweep(counts[names(gene.sets.rank), names(clusters)], MARGIN=2, STATS=umis[names(clusters)]/1E6, FUN="/"))))

  #
  meta.gene = t(sapply(gene.sets.sort, function(gene.set) {
    if (length(names(gene.sets.rank)[which(gene.sets.rank == gene.set)]) == 1) {
      return(tpm.z[names(gene.sets.rank)[which(gene.sets.rank == gene.set)],])
      } else {
        return(colMeans(tpm.z[names(gene.sets.rank)[which(gene.sets.rank == gene.set)],]))
        }
  }))

  #####
  # what if no assigned genes to a cluster?
  rank=hclust(dist(t(mat)), method="single")$order
  ####

  # Assign each Gene to its affiliated mean-target-Z-score
  target.gene.z = unlist(lapply(gene.sets.sort, function(gene.set) {
      target.gene = cluster.res[names(gene.sets.rank)[which(gene.sets.rank == gene.set)], as.character(gene.set)]
      names(target.gene) <- names(gene.sets.rank)[which(gene.sets.rank == gene.set)]
      return(target.gene)
  }))


  #### Added correction for situations where no genes assigned to a cluster
  # Assign each Cell to its affiliated target-specific MetaGene score
  target.meta = unlist(lapply(cluster.sort, function(cluster) {
    if (length(names(gene.sets.rank)[which(gene.sets.rank == gene.set)]) > 0) {
      return(meta.gene[as.character(cluster), names(clusters)[which(clusters == cluster)]])
      } else {
        cells = names(clusters)[which(clusters == cluster)]
        cell.order = cells[hclust(dist(t(tpm.z[names(gene.sets.rank), cells])), method="single")$order]
        meta.gene.hclust = length(cell.order):1
        names(meta.gene.hclust) <- cell.order
        return(meta.gene.hclust)
        }
    }))

  # initiate heatmap.matrix and cap highs and lows at +2 and -2
  tpm.z[tpm.z < -2] <- -2
  tpm.z[tpm.z > 2] <- 2
  # Specify numbers of spaces in between targets
  # Re-organize matrix so there are spaces in between cell targets
  # Order cells in each target based on MetaGene expression of target-specific gene set
  for (cluster in 1:length(cluster.sort)) {
    cell.order = names(sort(target.meta[names(clusters)[which(clusters == cluster.sort[cluster])]], decreasing=TRUE))
    if (cluster == 1) {
      heatmap.cell.space = tpm.z[,cell.order]
    } else {
      heatmap.cell.space = cbind(heatmap.cell.space, matrix(NA, nrow=nrow(heatmap.cell.space), ncol=col.space), tpm.z[,cell.order])
    }
  }

  # Specify numbers of spaces in between gene sets
  # Re-organize matrix so there are spaces in between target-specific gene sets
  # Order genes in each gene set based on target-specific Z-score
  for (gene.set in 1:length(gene.sets.sort)) {
      gene.order = names(sort(target.gene.z[names(gene.sets.rank)[which(gene.sets.rank == gene.set)]], decreasing=TRUE))
      if (gene.set == 1) {
        heatmap.gene.space = heatmap.cell.space[gene.order,]
      } else {
        heatmap.gene.space = rbind(heatmap.gene.space, matrix(NA, ncol=ncol(heatmap.gene.space), nrow=row.space), heatmap.cell.space[gene.order,])
      }
    }
  # Convert all NAs to -3 (display as white spaces)
  heatmap.gene.space[is.na(heatmap.gene.space)] <- -3
  #
  return(heatmap.gene.space)
  #
}
