
#' Plot all g-phiclusts
#'
#' plot of all g-phiclusts obtained per cluster
#'
#' @param measure_obj the measure-object as the output of phiclust
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#'
#' @return ggplot plot
#'
#' @examples
#' data("out")
#' plot_all_g_phiclusts(out)
#'
#' @export



plot_all_g_phiclusts <- function(measure_obj){

  lvls <- paste0(c(sort(unique(as.numeric(measure_obj$all_info$singular_value))), 0))

  if(lvls[1] == "0"){
    lvls <- lvls[-1]
  }

  measure_obj$all_info$singular_value <- factor(measure_obj$all_info$singular_value, levels = lvls)
  measure_obj$all_info$celltype <- factor(measure_obj$all_info$celltype, levels = names(measure_obj$phiclust)[order(measure_obj$phiclust)])

  cls <- colorRampPalette(c("orchid2", "steelblue1", "olivedrab2", "tomato"))

  if(!is.na(match("0", measure_obj$all_info$singular_value)[1]) > 0){
    cls <- c(cls(length(unique(measure_obj$all_info$singular_value)) - 1), "lightgrey")
  }else{
    cls <- c(cls(length(unique(measure_obj$all_info$singular_value))))
  }

  p <- ggplot(measure_obj$all_info, aes(x = .data$celltype, y = .data$g_phiclust, colour = .data$singular_value)) +
    geom_point(size = 2) + theme_light(base_size = 12) + xlab("") +
    ylab(expression(phi["clust"]^"g")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_manual(values = cls) + guides(colour=guide_legend(title = "singular value"))


  print(p)
}
