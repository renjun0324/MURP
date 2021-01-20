#' LogLikly
#'
#' @description
#' multi-threads likelyhood function for multivariate distribution.
#' (Prior probability)
#'
#' @param GeneNum: gene number
#' @param CellNum: cell number
#' @param K: the number of clusters
#' @param Res: the residual of cells and center in each cluster
#' @param sigma2: the coariance matrix of the multivariate gaussian distribution
#' @param cores: the number of threads
#'
#' @importFrom parallel mclapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom mvnfast dmvn
#' @import grDevices
#'
#' @noRd
#'
LogLikly <- function(Data = NULL,
                     centers = NULL,
                     cluster = NULL,
                     K = NULL,
                     seed = NULL,
                     cores = 1) {
  # K <- K_series[i]
  # cl <- cls[[i]]

  CellNum = nrow(Data)
  GeneNum = ncol(Data)

  Yt_filt <- centers[cluster, ] #Replace all points with all current center stores
  Yt_res <- Data - Yt_filt
  #Yt_res <- MaxIfMatrix(Yt_res)

  SStot <- (Data - t(colMeans(Data) %*% matrix(1, 1, CellNum))) ^ 2
  SSpam <- (Yt_filt - t(colMeans(Data) %*% matrix(1, 1, CellNum))) ^ 2
  SSres <- sapply(1:GeneNum, function(i) (sum(SStot[, i]) - sum(SSpam[, i]))) #the res of every gene

  SSres <- MaxIfVector(SSres)
  sigma2res <- SSres / (CellNum - K)
  #sigma2res[which(sigma2res<0)] <- 1e-32

  LL <- pbmclapply(1:GeneNum,
                 function(g) {
                   set.seed(seed)
                   mvnfast::dmvn(Yt_res[, g], mu = rep(0, CellNum), sigma = diag(sigma2res[g], CellNum), log = TRUE)},
                 mc.cores = cores )

  loglike <- sum(unlist(LL))

  #write.table(sigma2res, file = 'var.txt')
  return(loglike)
}


#' MaxIfMatrix
#'
#' @description
#' determine whether the expression level of genes is all 0
#'
#' @param dat: orig data
#' @noRd
#'
MaxIfMatrix <- function(dat) {
  a <- apply(dat,2,max)
  index <- which(a==0)
  if( length(index)==0 ) {
    return(dat)
  } else {
    dat[nrow(dat),index] <- 1e-32
    return(dat)
  }
}

#' MaxIfVector
#'
#' @description
#' determine whether the expression level of genes is all 0
#'
#' @param dat: orig data
#' @noRd
#'
MaxIfVector <- function(dat) {
  index <- which(dat==0)
  if( length(index)==0 ) {
    return(dat)
  } else {
    dat[index] <- 1e-32
    return(dat)
  }
}
