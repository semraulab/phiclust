#' Main function to calculate phiclust
#'
#' This function will calculate the value of phiclust for each cluster and output a measure object which can be used with all the plotting functions
#'
#' @param expr a data matrix with cells in the columns and genes in the rows, preferably normalized and log-transformed
#' @param clusters a vector of the same length as the number of cells, indicating which cell type they belong to
#' @param exclude a data.frame of variables to reduce in the measure, e.g. total number of transcripts, average expression of MT, Rb or stress genes.
#'  The data frame should have the same number of rows as cells (also in the same order), and each column corresponds to a different variable.
#' @param confidence TRUE/FALSE if the confidence interval for phiclust should be calculated. Caution: Increases computational time significantly.
#' @param exp_genes percentage of variance driving genes to extract per sigificant singular vector
#' @param exclude_outlier_cells TRUE/FALSE if outlier cells should be excluded, default is FALSE (this functions is not fully tested)
#' @param outlier_value cutoff for outlier cells
#' @param p.val the p-value to be used in the test of normality for the singular vectors, default is 0.01
#' @param nu number of left singular vectors to calculate, default is 50. High values increase computational time.
#'
#' @return A measure object:
#' \itemize{
#' \item phiclust: phiclust for each cluster
#' \item g_phiclust: g-phiclust for each cluster
#' \item all_info: detailed information about each singular value (see function get_info)
#' \item genes: list of variance driving genes (see function get_var_genes)
#' \item rmt_out: MP object for each cluster (see function fit_mp)
#' \item cell.index: if exclude_outlier_cells = TRUE, then the index of the excluded cells can be found here
#' \item input_parameters: the inputs to the function
#' }
#'
#' @importFrom stats median mad lm p.adjust rnorm
#' @importFrom parallel mclapply
#'
#' @examples
#' #Load sample data simulated with splatter
#' library(splatter)
#' data("splatO")
#' expr <- counts(splatO)
#' expr <- expr[rowSums(expr)>0,]
#'
#' #Normalize and log-transform the data
#' expr.norm <- t(t(expr)/colSums(expr))*10000
#' expr.norm.log <- log(expr.norm + 1)
#'
#' #Create toy example of a data set
#' test.cluster <- as.character(splatO$Group)
#' test.cluster[test.cluster == "Group3"] <- "Group2"
#' test.cluster[test.cluster == "Group4"] <- "Group2"
#'
#' #Main funcion that calculates the clusterability
#' out <- phiclust(expr = expr.norm.log, clusters = test.cluster,
#'               exclude = data.frame(clsm = log(colSums(expr) + 1)))
#'
#' @export


phiclust <- function(expr, clusters, exclude = NULL, confidence = F, exp_genes = 0.01, exclude_outlier_cells = F, outlier_value = 10, p.val = 0.01, nu = 50){

  celltype <- as.character(clusters)
  uclt <- unique(celltype)

  thetas.add <- list()
  lambdas.add <- list()
  lambdas.corr <- list()
  angles.add <- list()
  angles.add2 <- list()
  angles.add_u <- list()
  all_r2vals <- list()
  gene.list <- list()
  sd.up <- list()
  sd.down <- list()
  sd.up.u <- list()
  sd.down.u <- list()

  rmt.list <- list()
  cell.index <- list()

  for(i in 1:length(uclt)){
    if(sum(celltype == uclt[i]) >= 10){

      cat("Calculating values for cluster ", uclt[i], "\n")

      data.small2 <- expr[,celltype == uclt[i]]

      filter_genes <- apply(data.small2, 1,
                            function(x) length(x[x > 1]) >= 2)

      out.gene <- t(apply(data.small2, 1, robust_scale))
      ind <- rowSums(out.gene>50) <= 10 & rowSums(out.gene>50) > 0

      data.small.scale <- t(scale(t(data.small2[filter_genes & !ind,])))

      if(exclude_outlier_cells){

        s <- svd(scale(data.small.scale), nu = 0, nv = min(50, ncol(data.small.scale)))
        out <- apply(s$v, 2, robust_scale)
        ind.cells <- unique(unlist(apply(out, 2, function(y) which(abs(y) > outlier_value))))

        if(length(ind.cells) > 0 & length(ind.cells) <= 5){

          cell.index[[uclt[i]]] <- ind.cells

          data.small2 <- data.small2[,-ind.cells]

          filter_genes <- apply(data.small2, 1,
                                function(x) length(x[x > 1]) >= 2)

          out.gene <- t(apply(data.small2, 1, robust_scale))
          ind <- rowSums(out.gene>20) <= 10 & rowSums(out.gene>20) > 0

          data.small.scale <- t(scale(t(data.small2[filter_genes & !ind,])))
        }
      }

      data.small.scale[data.small.scale < - sqrt(ncol(data.small.scale))] <- - sqrt(ncol(data.small.scale))
      data.small.scale[data.small.scale > sqrt(ncol(data.small.scale))] <- sqrt(ncol(data.small.scale))


      cat("Dim: ",dim(data.small.scale), "\n")

      if(ncol(data.small.scale) >= 50){
        L <- fit_mp(expr = data.small.scale, sample = FALSE, cor = T, p.val = p.val, nu = nu)
      }else{
        L <- fit_mp(expr = data.small.scale, sample = TRUE, cor = T, p.val = p.val, nu = nu)
      }

      rmt.list[[uclt[i]]] <- L

      nn <- max(round((exp_genes)*nrow(data.small.scale)), 10)

      if(length(L$sig_vectors) > 0){

        sel.genes <- data.frame(apply(L$svd$u[,L$sig_vectors, drop = F], 2, function(x){
          as.character(c(rownames(data.small.scale)[order(x, decreasing = T)[1:nn]], rownames(data.small.scale)[order(x, decreasing = F)[1:nn]]))
        }), stringsAsFactors = F)

        rownames(sel.genes) <- c(paste0("Highest-", 1:nn), paste0("Lowest-", 1:nn))
        colnames(sel.genes) <- paste0("Singular vector-", L$sig_vectors)
        sel.genes <- data.frame(sel.genes)

        gene.list[[uclt[i]]] <- sel.genes

        if(!is.null(exclude)){

          #regress out unwanted variations
          df.responds <- data.frame(x = L$eigen$vectors[,L$sig_vectors])
          df.fit <- data.frame(exclude[celltype == uclt[i],])

          if(length(cell.index[[uclt[i]]]) > 0){
            df.fit <- df.fit[-ind.cells,, drop = F]
          }

          df.fit <- data.frame(t(t(df.fit)/sqrt(colSums(df.fit^2))))

          r2vals <- apply(df.responds,2, function(k){
            lr <- lm(k ~ ., data = df.fit[,, drop = FALSE])
            lrs <- summary(lr)
            return(lrs$adj.r.squared)
          })
        }else{
          r2vals <- rep(0, length(L$sig_vectors))
        }

        r2vals[r2vals < 0] <- 0

        #Adjustements
        perc.lambda <- (1 - r2vals)*(L$eigen$values[L$sig_vectors])
        perc.sigma <- sqrt(perc.lambda*(L$M - 1))

        #save resulting angle
        all_r2vals[[uclt[i]]] <- r2vals
        names(all_r2vals[[uclt[i]]]) <- L$sig_vectors
        thetas.add[[uclt[i]]] <- obtain_theta_sv(lambda = perc.sigma/sqrt(L$M-1), L = L)
        names(thetas.add[[uclt[i]]]) <- L$sig_vectors

        lambdas.add[[uclt[i]]] <- L$svd$d[L$sig_vectors]/sqrt(L$M-1)
        names(lambdas.add[[uclt[i]]]) <- L$sig_vectors
        lambdas.corr[[uclt[i]]] <- perc.sigma/sqrt(L$M-1)
        names(lambdas.corr[[uclt[i]]]) <- L$sig_vectors

        angles.add[[uclt[i]]] <- vec_norm_sv(theta = thetas.add[[uclt[i]]], L = L)
        names(angles.add[[uclt[i]]]) <- L$sig_vectors
        angles.add_u[[uclt[i]]] <- vec_norm_sv_u(theta = thetas.add[[uclt[i]]], L = L)
        names(angles.add_u[[uclt[i]]]) <- L$sig_vectors

        #Calculating confidence intervals
        if(confidence){
          sig.ind <- L$sig_vectors
          s.signal <- L$svd$u[,sig.ind]%*%diag(sqrt(L$M - 1)*thetas.add[[uclt[i]]], nrow = length(L$sig_vectors), ncol = length(L$sig_vectors))%*%t(L$svd$v[,sig.ind])

          conf.out <- unlist(mclapply(1:50, function(x){
            set.seed(x)
            mat.nois <- matrix(rnorm(L$M*L$N, mean = 0, sd = 1), ncol = L$N)
            new_mat <- mat.nois + s.signal
            L2 <- svd(new_mat, nu = 0, nv = 0)

            return(max(L2$d/sqrt(L$M - 1)))
          }))

          thetas.conf <- obtain_theta_sv(lambda = conf.out, L = L)
          vec.conf <- vec_norm_sv(theta = thetas.conf, L = L)
          up.ind <- vec.conf >= max(angles.add[[uclt[i]]])

          sd.up[[uclt[i]]] <- sqrt((1/(sum(up.ind) - 1))*sum((vec.conf[up.ind] - max(angles.add[[uclt[i]]]))^2))
          sd.down[[uclt[i]]] <- sqrt((1/(sum(!up.ind) - 1))*sum((vec.conf[!up.ind] - max(angles.add[[uclt[i]]]))^2))

          #For g-phiclust
          vec.conf.u <- vec_norm_sv_u(theta = thetas.conf, L = L)
          up.ind <- vec.conf.u >= max(angles.add_u[[uclt[i]]])

          sd.up.u[[uclt[i]]] <- sqrt((1/(sum(up.ind) - 1))*sum((vec.conf.u[up.ind] - max(angles.add_u[[uclt[i]]]))^2))
          sd.down.u[[uclt[i]]] <- sqrt((1/(sum(!up.ind) - 1))*sum((vec.conf.u[!up.ind] - max(angles.add_u[[uclt[i]]]))^2))
        }

      }else{
        gene.list[[uclt[i]]] <- "No further clusters"

        #save resulting angle
        all_r2vals[[uclt[i]]] <- 0
        names(all_r2vals[[uclt[i]]]) <- 0
        thetas.add[[uclt[i]]] <- 0
        names(thetas.add[[uclt[i]]]) <- 0

        lambdas.add[[uclt[i]]] <- 0
        names(lambdas.add[[uclt[i]]]) <- 0
        lambdas.corr[[uclt[i]]] <- 0
        names(lambdas.corr[[uclt[i]]]) <- 0

        angles.add[[uclt[i]]] <- 0
        names(angles.add[[uclt[i]]]) <-0
        angles.add_u[[uclt[i]]] <- 0
        names(angles.add_u[[uclt[i]]]) <- 0

        if(confidence){
          sd.up[[uclt[i]]] <- 0
          sd.down[[uclt[i]]] <- 0
          sd.up.u[[uclt[i]]] <- 0
          sd.down.u[[uclt[i]]] <- 0
        }

      }


    }

  }

  pvals.mp <- c()
  u.cl <- unique(clusters)
  for(i in u.cl){
    pvals.mp <- c(pvals.mp, rmt.list[[i]]$p.value_mp_fit)
  }
  names(pvals.mp) <- u.cl
  pvals.mp <- c(p.adjust(pvals.mp[!is.na(pvals.mp)], "BH"), pvals.mp[is.na(pvals.mp)])

  if(confidence){
    maximum_measure <- data.frame( EV = unlist(lapply(angles.add, max)), upper = unlist(lapply(angles.add, max)) + unlist(sd.up),
                                   lower = unlist(lapply(angles.add, max)) - unlist(sd.down))

    maximum_measure.u <- data.frame( EV = unlist(lapply(angles.add_u, max)), upper = unlist(lapply(angles.add_u, max)) + unlist(sd.up.u),
                                   lower = unlist(lapply(angles.add_u, max)) - unlist(sd.down.u))
  }else{
    maximum_measure <- unlist(lapply(angles.add, max))
    maximum_measure.u <- unlist(lapply(angles.add_u, max))
  }

  all_info <- data.frame(phiclust = unlist(angles.add), g_phiclust = unlist(angles.add_u), lambda = unlist(lambdas.add),
                         r2vals = unlist(all_r2vals), lambda_corrected = unlist(lambdas.corr), theta = unlist(thetas.add))
  all_info$singular_value <- unlist(lapply(angles.add, names))
  all_info$celltype <- rep(names(angles.add), lengths(angles.add))
  rownames(all_info) <- NULL

  return(list(phiclust = maximum_measure, g_phiclust = maximum_measure.u, all_info = all_info, genes = gene.list, rmt_out = rmt.list, cell.index = cell.index, p.val.mp.fit = pvals.mp,
              input_parameters = list(expr = expr, clusters = clusters, exclude = exclude)))
}


#Auxillary functions

robust_scale <- function(x) {
  return((x - median(x)) / (mad(x) + .Machine$double.eps))
}

max_eig_sv <- function(theta, L){ # function that calculates lambda_max

  Q <- L$N/L$M
  v <- c() #where to save the output

  gbplus <- Q^(1/4) # 1/G(b+) threshold

  for(i in 1:length(theta)){
    if(abs(theta[i]) > gbplus ){
      v <- c(v, sqrt( ((1 + theta[i]^2)*(Q + theta[i]^2))/theta[i]^2)  )  # G inverse of 1/theta for (v): see above
    }else{
      v <- c(v, 1 + sqrt(Q) ) # maximum MP otherwise
    }
  }
  return(v)
}


vec_norm_sv <- function(theta, L){ # function that calculates lambda_max

  Q <- L$N/L$M
  v <- c() #where to save the output

  gbplus <- Q^(1/4) # 1/G(b+) threshold

  for(i in 1:length(theta)){
    if(abs(theta[i]) > gbplus ){
      v <- c(v,  1 - (Q*(1 + theta[i]^2)/( theta[i]^2*(theta[i]^2 + Q))))  # G inverse of 1/theta for (v): see above   1 - ( (Q+theta[i]^2)/(theta[i]^2*(theta[i]^2 + 1))  )
    }else{
      v <- c(v, 0 ) # maximum MP otherwise
    }
  }
  return(v)
}

vec_norm_sv_u <- function(theta, L){ # function that calculates lambda_max

  Q <- L$N/L$M
  v <- c() #where to save the output

  gbplus <- Q^(1/4) # 1/G(b+) threshold

  for(i in 1:length(theta)){
    if(abs(theta[i]) > gbplus ){
      v <- c(v,  1 - ((Q + theta[i]^2)/( theta[i]^2*(theta[i]^2 + 1))))  # G inverse of 1/theta for (v): see above   1 - ( (Q+theta[i]^2)/(theta[i]^2*(theta[i]^2 + 1))  )
    }else{
      v <- c(v, 0 ) # maximum MP otherwise
    }
  }
  return(v)
}


obtain_theta_sv <- function(lambda, L){

  Q <- L$N/L$M
  v <- c()

  aa <- 1 - sqrt(Q)
  bb <- 1 + sqrt(Q)

  for(i in 1:length(lambda)){
    if(abs(lambda[i]) > bb){
      v <- c( v, sqrt(2*Q/( lambda[i]^2 - (Q+1) - sqrt( (lambda[i]^2 - (Q+1))^2 - 4*Q ) )) )
    }else{
      v <- c(v, 0)
    }
  }


  return(v)

}
