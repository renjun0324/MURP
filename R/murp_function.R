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
#' @importFrom pbapply pblapply
#' @importFrom NPflow mvnpdfC
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @import grDevices
#'
#' @noRd
#'
LogLikly <- function(Data = NULL,
                     centers = NULL,
                     cluster = NULL,
                     K = NULL,
                     cores = 1) {

  CellNum = nrow(Data)
  GeneNum = ncol(Data)

  Yt_filt <- centers[cluster, ] #Replace all points with all current center stores
  Yt_res <- Data - Yt_filt

  SStot <- (Data - t(colMeans(Data) %*% matrix(1, 1, CellNum))) ^ 2
  SSpam <- (Yt_filt - t(colMeans(Data) %*% matrix(1, 1, CellNum))) ^ 2
  SSres <- sapply(1:GeneNum, function(i) (sum(SStot[, i]) - sum(SSpam[, i]))) #the res of every gene

  SSres <- MaxIfVector(SSres)
  sigma2res <- SSres / (CellNum - K)

  cl <- makeCluster(cores)
  clusterExport(cl, c("Data", "centers", "cluster", "mvnpdfC",
                      "K", "CellNum", "GeneNum", "Yt_res"), envir = environment())
  LL <- parLapply(cl, 1:ncol(Data), function(g){
    yr = Yt_res[, g]
    if(!is.matrix(yr)){
      yr = as.matrix(yr)
    }
    mvnpdfC(yr, mean = rep(0, CellNum), varcovM = diag(sigma2res[g], CellNum), Log = TRUE)
  })
  stopCluster(cl)


  loglike <- sum(unlist(LL))

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
