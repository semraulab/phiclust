#' Plot phiclust
#'
#' plot the value of phiclust for every cluster
#'
#' @param measure_obj the measure-object as the output of phiclust
#'
#' @import ggplot2
#'
#' @return ggplot plot
#'
#' @examples
#' data("out")
#' plot_phiclust(out)
#'
#' @export

plot_phiclust <- function(measure_obj){

  if(is.data.frame(measure_obj$phiclust)){
    m.plot <- measure_obj$phiclust
    m.plot$celltype = rownames(m.plot)

    m.plot$celltype <- factor(m.plot$celltype, levels = m.plot$celltype[order(m.plot$EV)])

    p <- ggplot(m.plot, aes(x = .data$celltype, y = .data$EV, colour = .data$celltype)) + geom_point(size = 1) + theme_light(base_size = 12) + xlab("") +
      ylab(expression(phi["clust"])) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.3)

    print(p)
  }else{
    m.plot <- data.frame(celltype = names(measure_obj$phiclust), values = measure_obj$phiclust)
    m.plot$celltype <- factor(m.plot$celltype, levels = m.plot$celltype[order(m.plot$values)])

    p <- ggplot(m.plot, aes(x = .data$celltype, y = .data$values, colour = .data$celltype)) + geom_point(size = 2) + theme_light(base_size = 12) + xlab("") +
      ylab(expression(phi["clust"])) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

    print(p)
  }

}
