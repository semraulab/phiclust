#' Fit the Marchenko-Pastur distribution
#'
#' This function takes a data matrix as input and out puts the parameters asscociated to the MP distribution.
#'
#' @param expr a data matrix with cells in the columns and genes in the rows, preferably standardized gene-wise
#' @param sample TRUE/FALSE, if the parameters should be estimated by random sampling or not, default is FALSE
#' @param cor TRUE/FALSE, if the svd should be calculated on the correlation matrix or not (covariance matrix), default is TRUE
#' @param nu the number of gene singular vectors to calculate in the process (the more, the more time expansive), default is 50
#' @param p.val the p-value to be used in the test of normality for the singular vectors, default is 0.01
#'
#' @return A MP object:
#' \itemize{
#' \item eigen: eigenvalues and vectors of the cell-cell correlation matrix
#' \item maxEigen: maximum eigenvalue of the MP distribution
#' \item minEigen: minimum eigenvalue of the MP distribution
#' \item sig_vectors: singular vectors that lie significantly above the MP distribution
#' \item M: the number of genes
#' \item N: the number of cells
#' \item svd: the singular value decomposition
#' \item genes.used: genes that have been used for the calculation of the MP distribution
#' \item p.value_mp_fit: p-value for the similarity of a MP distribution to the eigenvalues found as noise
#' \item transcriptome_mode: the index of the transcriptome mode
#' \item input_parameters: the inputs to the function
#' }
#'
#' @importFrom stats quantile shapiro.test cor ks.test
#' @importFrom utils tail
#' @importFrom RMTstat rmp
#'
#' @examples
#' library(splatter)
#' data("splatO")
#' expr <- counts(splatO)
#' expr <- expr[rowSums(expr)>0,]
#'
#' #Normalize and log-transform the data
#' expr.norm <- t(t(expr)/colSums(expr))*10000
#' expr.norm.log <- log(expr.norm + 1)
#'
#' expr.scale <- t(scale(t(expr.norm.log)))
#' L <- fit_mp(expr.scale)
#'
#' @export


fit_mp <- function(expr, sample = FALSE, cor = TRUE, nu = 50, p.val = 0.01){

  N <- ncol(expr)
  p.val.thres <- p.val

  cat("Calculating svd ...", "\n")
  if(cor){
    svd.expr <- svd(scale(expr, scale = TRUE), nv = N, nu = nu)
  }else{
    svd.expr <- svd(scale(expr, scale = FALSE), nv = N, nu = nu)
  }
  M <- nrow(expr)


  s <- list()
  s[["values"]] <- c((svd.expr$d^2)/(M-1), rep(0, N-length(svd.expr$d)))
  s[["vectors"]] <- svd.expr$v
  s[["scores"]] <- svd.expr$u

  ind <- order(s$values)
  eigvals <- s$values[ind]
  V <- s$vectors[,ind]

  sigma <- 1 - (max(s$values)/sum(s$values))
  #cat("Scaled by: ", sigma, "\n")

  Q <- M/N

  if(!sample){

    #Maximum eigenval
    RMTmaxEig <- (mean(s$values))*sigma*(1 + (1/Q) + 2*sqrt(1/Q))
    RMTminEig <- (mean(s$values))*sigma*(1 + (1/Q) - 2*sqrt(1/Q))

  }else{
    out <- unlist(lapply(1:20, function(x) random_sampling(expr)))
    RMTmaxEig <- quantile(out[grep("max", names(out))], prob = 0.99)
    RMTminEig <- quantile(out[grep("min", names(out))], prob = 0.01)
  }

  tw <- tracy_widom(N, M)

  #Minimum eigenval
  RMTminIndex <- which(eigvals < RMTminEig)
  RMTminIndex <- tail(RMTminIndex, 1)

  #Calculating transcriptome mode
  transcriptome_mode <- which(colSums(V>0) == 0 | colSums(V<0) == 0)
  cat("Market Mode: ", N + 1 - transcriptome_mode, "\n")

  #Calculating p-value
  if(N > 50){
    r2 <- s$values[ !(colSums(V>0) == 0 | colSums(V<0) == 0) & s$values <= RMTmaxEig]
    r <- sigma*rmp(1000, ndf = M, pdim = N)
    p.val.mp <- ks.test(r, r2)$p.value
  }else{
    p.val.mp <- NA
  }


  #Critical eigenvalue
  RMTmaxIndex <- which(eigvals > tw)

  if (isEmpty(RMTmaxIndex)){

    out <- list(eigen = s, maxEigen = RMTmaxEig, minEigen = RMTminEig, sig_vectors = c(),
                M = M, N = N, svd = svd.expr, genes.used = rownames(expr), p.value_mp_fit = p.val.mp, transcriptome_mode = transcriptome_mode,
                input_parameters = list(sample = sample, cor = cor, expr = expr))
    return(out)

  }else{
    RMTmaxIndex <- RMTmaxIndex[1]}


  #Check loacalisation of eigenvectors above Marchenko-Pastur
  to.test <- N:RMTmaxIndex

  if(isEmpty(to.test)){
    out <- list(eigen = s, maxEigen = RMTmaxEig, minEigen = RMTminEig, sig_vectors = c(),
                M = M, N = N, svd = svd.expr, genes.used = rownames(expr), p.value_mp_fit = p.val.mp, transcriptome_mode = transcriptome_mode,
                input_parameters = list(sample = sample, cor = cor, expr = expr))
    return(out)
  }

  if(length(to.test) <= 5){
    passed_test <- c(0)
  }else if(length(to.test) <= 10 & length(to.test) > 5){
    passed_test <- c(1 : (length(to.test) - 5))
  }else{
    passed_test <- c(1 : (round((length(to.test))*0.5) - 1))
  }

  vec.to.test <- V[,to.test, drop = FALSE]
  rmt_indices <- to.test

  for(i in (max(passed_test)+1):dim(vec.to.test)[2]){

    if(length(vec.to.test[,i]) > 5000){
      p.val <- shapiro.test(vec.to.test[sample(1:nrow(vec.to.test), size = 5000),i])$p.value
      passed_test <- c(passed_test, if(p.val < p.val.thres) i )
    }else{
      p.val <- shapiro.test(vec.to.test[,i])$p.value
      passed_test <- c(passed_test, if(p.val < p.val.thres) i )}

  }


  out <- list(eigen = s, maxEigen = RMTmaxEig, minEigen = RMTminEig, sig_vectors = (N - c(rmt_indices[passed_test])+1), M = M,
              N = N, svd = svd.expr, genes.used = rownames(expr), p.value_mp_fit = p.val.mp,
              transcriptome_mode = N - transcriptome_mode + 1, input_parameters = list(sample = sample, cor = cor, expr = expr))
  return(out)
}




#Auxillary functions
#Tracy WIdom
tracy_widom <- function(N, M, k = 1, singvals = NULL){

  Q <- N/M

  #L <- eigen_values[(eigen_values > minEig) & (eigen_values < maxEig)]

  if(is.null(singvals)){
    lambda_max <- (1+sqrt(Q))^2
  }else{
    lambda_max <- mean(singvals^2)*(1+sqrt(Q))^2
  }

  gamma <- sqrt(Q) * lambda_max^(2/3)

  lambda_c <- lambda_max + (gamma*N^(-2/3))*k
  return(lambda_c)
}

#Sampling the data
random_sampling <- function(expr.norm.log){
  expr.norm.log.sample <- apply(expr.norm.log, 2, function(x) x[sample(1:length(x), length(x))])
  C <- cor(expr.norm.log.sample)
  s <- eigen(C)

  return(list("min" = tail(s$values, 1),"max" = s$values[1]))
}
