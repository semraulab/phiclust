#' Plot fit of MP
#'
#' singular values density plot and MP density as a solid line
#'
#' @param measure_obj the measure-object as the output of phiclust
#' @param cluster the name of the cluster of interest, e.g. "Pod"
#' @param plot.title optional title of the plot
#'
#' @import ggplot2
#' @importFrom S4Vectors isEmpty
#' @importFrom stats density
#'
#' @return ggplot plot
#'
#' @examples
#' data("out")
#' plot_MP(out, "Group2")
#'
#' @export


plot_MP <- function(measure_obj, cluster, plot.title = ""){

  L <- measure_obj$rmt_out[[cluster]]

  Q <- L$N/L$M

  if(Q > 1){
    singvals <- c(L$svd$d/(sqrt(L$M-1)), rep(0, L$N - L$M))
  }else{
    singvals <- L$svd$d/(sqrt(L$M-1))
  }


  a <- sqrt(abs(L$minEigen))
  b <- sqrt(L$maxEigen)

  sigma <- (1 - max(singvals)/sum(singvals))
  n <- seq(a, b, length.out = 10000)
  distr <- matrix(0, ncol = length(n))
  #Plug in results sigma, ...
  for (i in 1:length(n)){
    distr[1,i] <- sv_density(x = n[i], sigma = sigma, M = L$M, N = L$N, RMTminEig = L$minEigen, RMTmaxEig = L$maxEigen)
  }

  plot.dat <- data.frame("val" = singvals)

  #MP data frame
  plot.mp <- data.frame(val.mp=distr[1,])
  plot.mp$x <- n

  p <- ggplot() + geom_histogram(data=plot.dat, aes(x = .data$val, y=..density..), binwidth = max(singvals[-1])/150, fill = "skyblue") + theme_bw() +
    geom_line(data=plot.mp, aes(x = .data$x, y = .data$val.mp), col = "tomato", size = 1, alpha = 0.8) + xlab("Singular Values") + ylab("Density")+
    theme(axis.text=element_text(size = 18), axis.title = element_text(size=16, face="bold"), legend.position = "right",
          legend.background = element_rect(fill = "transparent", size=0.7), legend.key = element_rect(fill = "transparent"))+
    scale_color_manual(guide = guide_legend(
      label.theme = element_text(size=16, family = "Helvetica"),
      title.theme = element_text(size=16, face = "bold",family = "Helvetica"),
      title = "Methods",
      title.position = "top",
      label.position = "top",
      keywidth = 1,
      label.hjust = 0,
      title.hjust = 0,
      override.aes = list(alpha = 1, size=1.7),
      label.vjust = 0.5)) + ggtitle(plot.title)

  if(!S4Vectors::isEmpty(L$sig_vectors)){
    plot.sig.vec <- data.frame(sig = singvals[L$sig_vectors])
    plot.sig.vec$y <- 0.15
    plot.sig.vec$label <- "*"
    p <- p + geom_text(data = plot.sig.vec, aes(x = .data$sig, y = .data$y, label = .data$label), size = 8)
  }

  print(p)
}

#Auxillary function
sv_density <- function(x, sigma, RMTmaxEig, RMTminEig, M, N, singvals = NULL){
  Q <- N/M
  # a <- (1 - sqrt(Q))^2
  # b <- (1 + sqrt(Q))^2

  if(is.null(singvals)){
    if(x^2 > RMTminEig & x^2 < RMTmaxEig){
      out <- sqrt((x^2 - RMTminEig)*(RMTmaxEig - x^2))/(Q*pi*x)
    }else{
      out <- 0
    }
  }else{
    if(x^2 > RMTminEig & x^2 < RMTmaxEig){
      out <- sqrt((x^2 - RMTminEig)*(RMTmaxEig - x^2))/(Q*pi*x*mean(singvals^2))
    }else{
      out <- 0
    }
  }


  return(out)
}

