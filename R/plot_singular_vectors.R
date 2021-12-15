
#' Plot singular vectors for a cluster
#'
#' scatter plot of two singular vectors for a selected cluster
#'
#' @param measure_obj the measure-object as the output of phiclust
#' @param cluster the name of the cluster of interest, e.g. "Pod"
#' @param v1 index of vector to plot on the x-axis, default is 1
#' @param v2 index of vector to plot on the y-axis, default is 1
#' @param colour either a gene name or a vector the same length as the number of cells, default is NULL (no colour)
#' @param scaled TRUE/FALSE, if the scaled expression of a gene should be shown, default is FALSE
#'
#' @import ggplot2
#' @importFrom grDevices rgb
#'
#' @return ggplot plot
#'
#' @examples
#' data("out")
#' plot_singular_vectors(out, "Group2")
#'
#' @export

plot_singular_vectors <- function(measure_obj, cluster, v1 = 1, v2 = 2, colour = NULL, scaled = FALSE){

  L <- measure_obj$rmt_out[[cluster]]

  v.plot <- data.frame(X1 = L$eigen$vectors[,v1], X2 = L$eigen$vectors[,v2])

  if(length(colour) == 1){

    if(scaled){
      v.plot$gene <- L$input_parameters$expr[colour,]
    }else{
      if(length(measure_obj$cell.index[[cluster]]) > 0){
        v.plot$gene <- unlist(measure_obj$input_parameters$expr[colour, measure_obj$input_parameters$clusters == cluster][-measure_obj$cell.index[[cluster]]])
      }else{
        v.plot$gene <- measure_obj$input_parameters$expr[colour, measure_obj$input_parameters$clusters == cluster]
      }
    }

    v.plot$gene[v.plot$gene == 0] <- NA

    p <- ggplot(v.plot, aes(x = .data$X1, y = .data$X2, colour = .data$gene)) + geom_point() + theme_light(base_size = 14) +
      xlab(paste0("V", v1)) + ylab(paste0("V", v2)) +
      scale_color_gradientn(colours = c("navy", "yellow"),
                            na.value = rgb(0.75,0.75,0.75,alpha = 0.2), name = colour)

  }else if(length(colour) == nrow(v.plot)){
    v.plot$colour <- colour
    p <- ggplot(v.plot, aes(x = .data$X1, y = .data$X2, colour = .data$colour)) + geom_point() + theme_light(base_size = 14) +
      xlab(paste0("V", v1)) + ylab(paste0("V", v2))

  }else if(length(colour) == 0){
    p <- ggplot(v.plot, aes(x = .data$X1, y = .data$X2)) + geom_point() + theme_light(base_size = 14) +
      xlab(paste0("V", v1)) + ylab(paste0("V", v2))

  }else{
    cat("Not sure what to do with this colour information")
    p <- NA
  }

  print(p)
}
