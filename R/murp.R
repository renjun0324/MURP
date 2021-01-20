#' MURP
#'
#' @description
#' This new method is optimized from the MURP algorithm, it's more efficient.
#' Loglike: Prior probability estimated by likelyhood function;
#' MURP: Bayesian Information Criterions value;
#' Recommended_K: Recommended K value.
#'
#' @param Data expression matrix with Cells * Genes
#' @param omega omega value in pseudo-BIC calculating
#' @param iter the number of iterations
#' @param cores the number of cores to use in multi-threads calculating, default is 1.
#' Multi-threaded calculation is not recommended except for users with server support
#' @param seed random seed, default is 723
#'
#' @return a list include recommend_K, object of clustering result and so on
#'
#' @export
#'

MURP <- function(Data,
                 cores = 1,
                 iter = 20,
                 omega = 1/6,
                 seed = 723){

  f <- MURP.Kmeans(Data = Data, cores = cores, iter = iter, omega = omega, seed = seed)

  return(f)
}

#' BicSpan
#'
#' @description
#' Calculate the pseu-BIC value under different K.
#' Loglike: Prior probability estimated by likelyhood function;
#' MURP: Bayesian Information Criterions value;
#' Recommended_K: Recommended K value.
#'
#' @param Data expression matrix with Cells * Genes
#' @param omega omega value in pseudo-BIC calculating
#' @param h hclust object of orig data, which is not needed without hierarchical clustering
#' @param d dist matrix of orig data
#' @param K the number of specified clusters
#' @param cores the number of cores to use in multi-threads calculating,
#' default is 1. I do not recommend using multi-threaded calculations on
#' a local computer, but if you have server support, you can try.
#' @param seed random seed, default is 723
#'
#' @importFrom graphics abline plot
#' @importFrom stats dbeta kmeans median quantile sd var
#'
#' @return a list include recommend_K, object of clustering result and so on
#'
#' @export
#'
murp_specific <- function(Data = NULL,
                          h = NULL,
                          d = NULL,
                          cores = NULL,
                          omega = NULL,
                          K = NULL,
                          seed = 723){

  f <- murp_specific_Kmeans(Data = Data, cores = cores, K = K, omega = omega, seed = seed)

  return(f)
}






