#' MURP.Kmeans
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom NPflow mvnpdfC
#'
#' @export
#' @noRd
#'
MURP.Kmeans <- function(Data,
                        cores = 1,
                        iter = 20,
                        omega = 1/50,
                        seed = 723,
                        fast = FALSE,
                        cluster_iter_max = 1000){

  cat('Initiating ...        ',date(),'\n')

  CellNum <- dim(Data)[1]
  GeneNum <- dim(Data)[2]

  cat('K-Means calculating   ',date(),'\n')
  clsSave <- list()
  k <- c()
  Loglike <- c()
  BIC <- c()
  Step <- c()
  Penalty <- c()

  maxK <- round(CellNum * 1)
  K_series <- ceiling(quantile(1:maxK, 0.5))
  min_k <- round(K_series * 1)
  step <- 0

  for (i in 1:iter) {
    cat('Iter',i,'              ',date(),'\n')

    if (i != 1 & step <= 1) {
      cat('Step less than 1,stop iter \n')
      break()
    }

    if (i != 1) {
      min_k <- k[which(BIC == min(BIC))]
    }

    step <- (maxK*1) * 1/2^(i + 1)
    K_series <- ceiling(c(max(min_k - step, 2),
                          min(min_k + step,maxK - 1)))
    cat('Step:',step,'  k:',min_k,'\n')
    cat(K_series,'\n')

    cl <- makeCluster(cores)
    clusterExport(cl, c("seed", "Data", "K_series"), envir = environment())
    cls <- parLapply(cl, K_series, function(Centers){
      set.seed(seed)
      kmeans(Data, centers = Centers, iter.max = cluster_iter_max)
    })
    stopCluster(cl)

    gc()

    for(r in 1:length(cls) ){
      clsSave[[dim(cls[[r]]$centers)[1]]] <- cls[[r]]
    }

    K_series_len <- length(K_series)
    if(fast){
      Loglike_iter <-  unlist(pblapply(1:K_series_len,
                                       function(l){FastLogLikly(Data = Data,
                                                                centers = cls[[l]]$centers,
                                                                cluster = cls[[l]]$cluster,
                                                                K = K_series[l]) }))
    }else{
      Loglike_iter <-  unlist(pblapply(1:K_series_len,
                                       function(l){LogLikly(Data = Data,
                                                            centers = cls[[l]]$centers,
                                                            cluster = cls[[l]]$cluster,
                                                            K = K_series[l],
                                                            cores = cores ) }))
    }

    a <- 1 + 1e-9
    gamma <- sqrt(a * K_series * K_series * log(GeneNum))
    if(fast){
      penalty <- (1/100 * omega) * (1 + gamma) * log(GeneNum) * K_series
    }else{
      penalty <- (6 * omega) * (1 + gamma) * log(GeneNum) * K_series
    }

    BIC_iter <- penalty - 2 * Loglike_iter
    k <- append(k, K_series)
    Penalty <- append(Penalty, penalty)
    Loglike <- append(Loglike, Loglike_iter)
    BIC <- append(BIC, BIC_iter)
    Step <- append(Step, step)

    if ((min_k - step < 1 || min_k + step > CellNum - 1)) {
      break
    }
  }

  Recommended_K <- k[which(BIC == min(BIC))[1]]

  command <- paste('Data =',dim(Data)[1],'*',dim(Data)[2],
                   'omega =', omega,
                   'cores =',cores,
                   'iter =',iter)

  final <- list(k = k,
                Step = Step,
                Loglike = Loglike,
                BIC = BIC,
                Penalty = Penalty,
                Recommended_K = Recommended_K,
                Recommended_K_cl = clsSave[[Recommended_K]],
                Command = command,
                cl = clsSave,
                rawdata = Data)

  cat('BIC completed         ',date(),'\n')
  return(final)
}

#' BicSpan.Kmeans
#'
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom NPflow mvnpdfC
#'
#' @export
#' @noRd
#'
murp_specific_Kmeans <- function(Data,
                                 cores = 1,
                                 K = NULL,
                                 omega = 1/50,
                                 seed = 723,
                                 fast = FALSE,
                                 cluster_iter_max = 1000){

  cat('Initiating ...        ',date(),'\n')

  a = 1 + 1e-9

  CellNum <- dim(Data)[1]
  GeneNum <- dim(Data)[2]

  cat('K-Means calculating   ',date(),'\n')
  cat('K:',K,'\n')

  set.seed(seed)
  cls <- kmeans(Data, centers = K, iter.max = cluster_iter_max)

  if(fast){
    Loglike_iter = FastLogLikly(Data = Data,
                                centers = cls$centers,
                                cluster = cls$cluster,
                                K = K)
  }else{
    Loglike_iter = LogLikly(Data = Data,
                            centers = cls$centers,
                            cluster = cls$cluster,
                            K = K,
                            cores = cores )
  }


  gamma <- sqrt(a * K * K * log(GeneNum))
  penalty <- (6 * omega) * (1 + gamma) * log(GeneNum) * K
  BIC_iter <- penalty - 2 * Loglike_iter

  command <- paste('Data =', dim(Data)[1],'*',dim(Data)[2],
                   'omega = ', omega,
                   'cores =', cores,
                   'K =', K)

  final <- list(k = K,
                Loglike = Loglike_iter,
                BIC = BIC_iter,
                Penalty = penalty,
                Recommended_K = K,
                Command = command,
                cl = cls)

  cat('BIC completed         ',date(),'\n')

  gc()
  return(final)
}
