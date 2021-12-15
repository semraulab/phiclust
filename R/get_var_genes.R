#' Extract variance driving genes
#'
#' This function lists the variance driving genes for each significant singular vector
#'
#' @param measure_obj the measure-object as the output of phiclust
#' @param cluster the name of the cluster of interest, e.g. "Pod"
#'
#' @return A matrix, with the columns as the number of significant gene singular vecotrs and the rows contain the 1% genes with highest and lowest values in the gene singular vector
#'
#'@examples
#' data("out")
#' get_var_genes(out, "Group2")
#'
#' @export

get_var_genes <- function(measure_obj, cluster){

  gene_info <- measure_obj$genes[[cluster]]

  return(gene_info)
}
