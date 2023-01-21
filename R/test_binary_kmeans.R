#
# library(Seurat)
# library(MURP)
# library(dplyr)
# data(sdata)
#
# data = sdata
# CellNum = nrow(data)
# GeneNum = ncol(data)
#
# kmeans_binary <- function(dat,
#                           iter.max = 1000,
#                           nstart = 1){
#   if(nrow(dat)<3){
#     k_dat = tapply(1:nrow(dat), as.factor(rep(1, nrow(dat))),
#                    function(x, y) dat[x, , drop = FALSE])
#     return(k_dat)
#   }
#   k_r = kmeans(dat, centers = 2, iter.max = iter.max, nstart = nstart)
#   k_dat = tapply(1:nrow(dat), as.factor(k_r$cluster),
#                  function(x, y) dat[x, , drop = FALSE])
#   return(k_dat)
# }
#
# append_list <- function(dat_list){
#   t = list()
#   for(i in 1:length(dat_list)){
#     t = append(t, dat_list[[i]])
#   }
#   return(t)
# }
#
#
# # 1. binary kmeans clustering
# binary_list = list()
# dat_list = list(sdata)
# for(i in 1:10){
#   cat(i, "\n")
#
#   # binary kmeans clustering
#   dat_list = lapply(dat_list, kmeans_binary)
#   dat_list = append_list(dat_list)
#
#   # calculate cluster
#   names(dat_list) <- paste0("c",1:length(dat_list))
#   nr = sapply(dat_list, nrow)
#   cluster = rep(names(nr), times = nr)
#   cluster = as.numeric(gsub("c","",cluster))
#
#   # calculate centers
#   dat = do.call(rbind, dat_list) %>% data.frame
#   centers = do.call(rbind, lapply(dat_list, colMeans) )
#   binary_list = append(binary_list,
#                        list(list(cluster = cluster,
#                                  centers = centers,
#                                  k = nrow(centers))) )
#   names(binary_list)[i] = paste0("c", nrow(centers))
# }
#
# # 2. loglikly
# c = names(binary_list)
# llk = sapply(1:length(binary_list), function(i){
#   cat(i, "\n")
#   r = binary_list[[i]]
#   centers = r$centers
#   cluster = r$cluster
#   k = r$k
#   if(fast){
#     llk <- FastLogLikly(Data = data, centers = centers, cluster = cluster, K = k)
#   }else{
#     llk <- LogLikly(Data = data, centers = centers, cluster = cluster, K = k)
#   }
#
#   a <- 1 + 1e-9
#   gamma <- sqrt(a * K_series * K_series * log(GeneNum))
#   if(fast){
#     penalty <- (1/100 * omega) * (1 + gamma) * log(GeneNum) * K_series
#   }else{
#     penalty <- (6 * omega) * (1 + gamma) * log(GeneNum) * K_series
#   }
#   BIC_iter <- penalty - 2 * llk
#   return(BIC_iter)
# })
#
