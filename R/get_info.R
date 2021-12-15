#' Get informaton about each singular vectos
#'
#' This function extracts information about sigma in a cluster.
#'
#' @param measure_obj the measure-object as the output of phiclust
#' @param cluster the name of the cluster of interest, e.g. "Pod"
#'
#' @return A data.frame with the following information:
#' \itemize{
#' \item phiclust: the value of sigma for each singular value
#' \item g_phiclust: the value of g-sigma for each sigular value
#' \item theta: the singular values of the signal matrix
#' \item r2vals: the adjusted r squared obtained during the regression
#' \item singular_value: the index of the singular values (from highest to lowest after the first calculation)
#' \item celltype: the selected cluster
#' }
#'
#' @examples
#' data("out")
#' get_info(out, "Group2")
#'
#' @export

get_info <- function(measure_obj, cluster){

  angle_info <- measure_obj$all_info[grep(cluster, measure_obj$all_info$celltype),]

  return(angle_info)
}
