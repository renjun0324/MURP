#' LogLikly
#'
#' @description
#' multi-threads likelyhood function for multivariate distribution.
#' (Prior probability)
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

  # SSres <- MaxIfVector(SSres)
  sigma2res <- SSres / (CellNum - K)

  cl <- makeCluster(cores)
  clusterExport(cl, c("Data", "centers", "cluster", "mvnpdfC",
                      "K", "CellNum", "GeneNum", "Yt_res"), envir = environment())
  LL <- parLapply(cl, 1:ncol(Data), function(g){
    yr = Yt_res[, g]
    if(!is.matrix(yr)){
      yr = as.matrix(yr)
    }
    sig = sigma2res[g]
    if(sig==0){sig=1e-16}
    mvnpdfC(yr, mean = rep(0, CellNum), varcovM = diag(sig, CellNum), Log = TRUE)
  })
  stopCluster(cl)


  loglike <- sum(unlist(LL))

  return(loglike)
}


#' FastLogLikly
#'
#' @description
#' multi-threads likelyhood function for multivariate distribution.
#' (Prior probability)
#'
#' @importFrom pbapply pblapply
#' @importFrom NPflow mvnpdfC
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @import grDevices
#'
#' @noRd
#'
FastLogLikly <- function(Data = NULL,
                         centers = NULL,
                         cluster = NULL,
                         K = NULL) {

  CellNum = nrow(Data)
  GeneNum = ncol(Data)

  Yt_filt <- centers[cluster, ] #Replace all points with all current center stores
  Yt_res <- Data - Yt_filt

  SStot <- (Data - t(colMeans(Data) %*% matrix(1, 1, CellNum))) ^ 2
  SSpam <- (Yt_filt - t(colMeans(Data) %*% matrix(1, 1, CellNum))) ^ 2
  SSres <- sapply(1:GeneNum, function(i) (sum(SStot[, i]) - sum(SSpam[, i]))) #the res of every gene

  # SSres <- MaxIfVector(SSres)
  sigma2res <- SSres / (CellNum - K)

  # cl <- makeCluster(cores)
  # clusterExport(cl, c("Data", "centers", "cluster", "mvnpdfC",
  #                     "K", "CellNum", "GeneNum", "Yt_res"), envir = environment())
  # LL <- parLapply(cl, 1:ncol(Data), function(g){
  #   yr = Yt_res[, g]
  #   if(!is.matrix(yr)){
  #     yr = as.matrix(yr)
  #   }
  #   mvnpdfC(yr, mean = rep(0, CellNum), varcovM = diag(sigma2res[g], CellNum), Log = TRUE)
  # })
  # stopCluster(cl)
  # loglike <- sum(unlist(LL))

  yr = apply(Yt_res, 1, sum)
  if(!is.matrix(yr)){
    yr = as.matrix(yr)
  }
  loglike = mvnpdfC(yr, mean = rep(0, CellNum), varcovM = diag(sum(sigma2res), CellNum), Log = TRUE)

  return(loglike)
}

