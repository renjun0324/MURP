

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
                        cluster_iter_max = 1000,
                        max_murp = round(nrow(Data)*0.2) ){

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

  max_k <- max_murp
  K_series <- ceiling(quantile(1:max_k, 0.5))
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

    step <- (max_k*1) * 1/2^(i + 1)
    K_series <- ceiling(c(max(min_k - step, 2),
                          min(min_k + step,max_k - 1)))
    cat('Step:',step,'  k:',min_k,'\n')
    cat(K_series,'\n')

    cl <- makeCluster(cores)
    clusterExport(cl, c("seed", "Data", "K_series","cluster_iter_max"), envir = environment())
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

  command <- paste("Function: MURP.Kmeans",
                   'Data =',dim(Data)[1],'*',dim(Data)[2],
                   '| omega =', omega,
                   '| cores =',cores,
                   '| iter =',iter,
                   '| fast = ', fast,
                   '| seed = ', seed,
                   '| cluster_iter_max = ', cluster_iter_max,
                   '| max_murp = ', max_murp)

  final <- list(k = k,
                Step = Step,
                Loglike = Loglike,
                BIC = BIC,
                Penalty = Penalty,
                Recommended_K = Recommended_K,
                Recommended_K_cl = clsSave[[Recommended_K]],
                Command = command,
                cl = clsSave,
                rawdata = Data,
                date = date())

  cat('BIC completed         ',date(),'\n')
  return(final)
}


#' MURP.Kmeans.pca
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom NPflow mvnpdfC irlba
#'
#' @export
#' @noRd
#'
MURP.Kmeans.pca <- function(Data,
                            cores = 1,
                            iter = 20,
                            omega = 1/50,
                            seed = 723,
                            fast = FALSE,
                            cluster_iter_max = 1000,
                            center = F,
                            scale = F,
                            pc_number = 30,
                            max_murp = round(nrow(Data)*0.2)){

  cat('Initiating ...        ',date(),'\n')
  CellNum <- dim(Data)[1]
  GeneNum <- dim(Data)[2]

  cat('PCA ...        ',date(),'\n')
  set.seed(seed)
  pca = irlba::prcomp_irlba(Data, center = center, scale. = center, n = pc_number)
  pca_dat = pca$x[,1:pc_number]

  cat('K-Means calculating   ',date(),'\n')
  clsSave <- list()
  k <- c()
  Loglike <- c()
  BIC <- c()
  Step <- c()
  Penalty <- c()

  max_k <- max_murp
  K_series <- ceiling(quantile(1:max_k, 0.5))
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

    step <- (max_k*1) * 1/2^(i + 1)
    K_series <- ceiling(c(max(min_k - step, 2),
                          min(min_k + step,max_k - 1)))
    cat('Step:',step,'  k:',min_k,'\n')
    cat(K_series,'\n')

    cl <- makeCluster(cores)
    clusterExport(cl, c("seed", "Data", "pca_dat", "K_series", "cluster_iter_max"), envir = environment())
    cls <- parLapply(cl, K_series, function(Centers){
      set.seed(seed)
      tmp = kmeans(pca_dat, centers = Centers, iter.max = cluster_iter_max)
      cluster = tmp$cluster
      cluster_factor = factor(cluster, levels = 1:length(unique(cluster)))
      tapply(1:CellNum, cluster_factor,
             function(x, y) {
               # Data[x,,drop=F]
               apply(Data[x,,drop=F], 2, mean)
             }) -> tmplist
      centers = do.call(rbind, tmplist)
      return(list(cluster = cluster, centers = centers))
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

  command <- paste("Function: MURP.Kmeans.pca",
                   '| Data =',dim(Data)[1],'*',dim(Data)[2],
                   '| omega =', omega,
                   '| cores =',cores,
                   '| iter =',iter,
                   '| fast = ', fast,
                   '| seed = ', seed,
                   '| cluster_iter_max = ', cluster_iter_max,
                   '| max_murp = ', max_murp,
                   '| center =', center,
                   '| scale =', scale,
                   '| pc_number =', pc_number)

  final <- list(k = k,
                Step = Step,
                Loglike = Loglike,
                BIC = BIC,
                Penalty = Penalty,
                Recommended_K = Recommended_K,
                Recommended_K_cl = clsSave[[Recommended_K]],
                Command = command,
                cl = clsSave,
                rawdata = Data,
                date = date())

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

  command <- paste("Function: murp_specific_Kmeans",
                   '| Data =', dim(Data)[1],'*',dim(Data)[2],
                   '| omega = ', omega,
                   '| cores =', cores,
                   '| K =', K,
                   '| fast =', fast,
                   '| cluster_iter_max = ',cluster_iter_max)

  final <- list(k = K,
                Loglike = Loglike_iter,
                BIC = BIC_iter,
                Penalty = penalty,
                Recommended_K = K,
                Command = command,
                cl = cls,
                date = date())

  cat('BIC completed         ',date(),'\n')

  gc()
  return(final)
}
