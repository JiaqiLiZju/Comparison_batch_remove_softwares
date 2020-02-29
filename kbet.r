#library(devtools)
#install_github('theislab/kBET')
library("kBET")
require('FNN')
library("ggplot2")

run_kbet <- function(data, batch, subset_size=0.2, k_fractions=c(10, 15, 20, 25, 30, 40, 50), nPCs=20,
                     show_fig=FALSE, save_fig=FALSE) {
  #data: a matrix (rows: samples, columns: features (genes))
  #batch: vector or factor with batch label of each cell 
  
  # data sample
  # subsample to 10% default
  if (subset_size){
    subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
    data <- data[subset_id,]
    batch <- batch[subset_id]
  }
  
  # repeated runs of kBET
  kBET_result_list <- list()
  mean_kBET_list <- list()
  for (fr in k_fractions){
    # k <- floor(fr*length(batch)) #neighbourhood size: mean batch size 
    k <- floor(fr)
    fr <- as.character(fr)
    knn <- get.knn(data, k=k, algorithm = 'cover_tree')
    #now run kBET with pre-defined nearest neighbours.
    batch.estimate <- kBET(data, batch, k=k, knn = knn, testSize=NULL, do.pca=TRUE, dim.pca=nPCs)
    kBET_result_list[[fr]] <- batch.estimate
    mean_kBET <- batch.estimate$summary$kBET.observed[1]
    mean_kBET_list[[fr]] <- mean_kBET
  }
  
  # plot 
  if (show_fig||save_fig){
    plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                      each=length(batch.estimate$stats$kBET.observed)), 
                            data =  c(batch.estimate$stats$kBET.observed,
                                      batch.estimate$stats$kBET.expected))
    g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
      labs(x='Test', y='Rejection rate',title='kBET test results') +
      theme_bw() +  
      scale_y_continuous(limits=c(0,1))
    if (save_fig) {pdf(save_fig); g; dev.off()}
  }
  
  return(list(mean_kBET_list, kBET_result_list))
}

run_silhouette <- function(data, batch, repeat_nb=10, subset_size=0.8, nPCs=20) {
  batch_silh_list = list()
  data_orig <- data
  batch_orig <- batch
  for (i in 1:repeat_nb){
    # data sample
    # subsample to 10% default
    if (subset_size){
      subset_id <- sample.int(n = length(batch_orig), size = floor(subset_size * length(batch_orig)), replace=FALSE)
      data <- data_orig[subset_id,]
      batch <- batch_orig[subset_id]
    }
    
    # Compute a silhouette width and PCA-based measure:
    pca.data <- prcomp(data, center=TRUE) #compute PCA representation of the data
    batch.silhouette <- batch_sil(pca.data, batch, nPCs=nPCs)
    
    batch_silh_list[[as.character(i)]] <- batch.silhouette
  }
  
  return(batch_silh_list)
}

# 前面准备文件就用下面的代码
# python
#import scanpy as sc
#import pandas as pd

#adata=sc.read_h5ad("adata_combat.h5ad")
#df=pd.DataFrame(adata.X)
#df.index = adata.obs.index
#df.columns = adata.var.index
#df.to_csv("adata_combat.dge.txt")
#adata.obs.to_csv("adata_combat.dge_obs.txt")

# batch这里都一样随便选一个
# 这个dge_obs就是python scanpy adata.obs写出的文件
# batch <- "cluster_Baseline.txt"
# batch <- read.table(batch, sep=',', header=T)$batch #, row.name='index')

# # 这里加上其他的软件，
# # 重复代码片段就可以了
# # 这个dge就是python scanpy adata.X写出的文件
# data <- "adata_base.dge.csv"
# data_base <- read.table(data, sep=',', header=T, row.name='index')
# data <- "adata_combat.dge.csv"
# data_combat <- read.table(data, sep=',', header=T, row.name='index')
# data <- "adata_reg.dge.csv"
# data_reg <- read.table(data, sep=',', header=T, row.name='index')
# data <- "adata_scanorama.dge.csv"
# data_scano <- read.table(data, sep=',', header=T, row.name='index')
# data <- "adata_mnns.dge.csv"
# data_mnns <- read.table(data, sep=',', header=T, row.name='index')
# save.image("workspace.rdata")

# k_fractions = c(10, 15, 20, 25, 30, 40, 50)
# mean_kBET_list <- run_kbet(data_base, batch, subset_size=0.1, k_fractions=k_fractions)[[1]]
# combat_mean_kBET_list <- run_kbet(data_combat, batch, subset_size=0.1, k_fractions=k_fractions)[[1]]
# reg_mean_kBET_list <- run_kbet(data_reg, batch, subset_size=0.1, k_fractions=k_fractions)[[1]]
# scano_mean_kBET_list <- run_kbet(data_scano, batch, subset_size=0.1, k_fractions=k_fractions)[[1]]
# mnn_mean_kBET_list <- run_kbet(data_mnns, batch, subset_size=0.1, k_fractions=k_fractions)[[1]]

# # Wilcoxon statistical test
# # 如果检验结果差异显著p<0.05，就手动标注出来s
# wilcox.test(x=unlist(mean_kBET_list), y=unlist(combat_mean_kBET_list), alternative='greater', correct=TRUE)
# #p=0.9714
# wilcox.test(x=unlist(mean_kBET_list), y=unlist(reg_mean_kBET_list), alternative='greater', correct=TRUE)
# #p=0.1
# wilcox.test(x=unlist(mean_kBET_list), y=unlist(scano_mean_kBET_list), alternative='greater', correct=TRUE)
# #p=0.9
# wilcox.test(x=unlist(mean_kBET_list), y=unlist(mnn_mean_kBET_list), alternative='greater', correct=TRUE)
# #p=1

# # 这里在做图，先做出一个dataframe，后做折线图
# # rbind of kBET
# kBET_df = data.frame(class=rep(c('base', 'combat', 'regress_out', 'scanorama', 'mnns'), 
#                                each=length(mean_kBET_list)),
#                      fr=k_fractions,
#                      data=unlist(cbind(mean_kBET_list, combat_mean_kBET_list, reg_mean_kBET_list,
#                                        scano_mean_kBET_list, mnn_mean_kBET_list)))
# # do plot
# p <- ggplot(kBET_df, aes(x=factor(fr), y=data, group=class)) + 
#   geom_line(aes(colour=class)) + 
#   geom_point(size=6, aes(shape=class, colour=class)) + 
#   labs(x='K', y='Rejection rate', title='kBET results')
# p
# pdf("kBET_batch.pdf");p;dev.off()

load("workspace.rdata")

silh_list <- run_silhouette(data_base, batch, repeat_nb=10, subset_size=0.02)
combat_silh_list <- run_silhouette(data_combat, batch, repeat_nb=10, subset_size=0.02)
reg_silh_list <- run_silhouette(data_reg, batch, repeat_nb=10, subset_size=0.02)
scano_silh_list <- run_silhouette(data_scano, batch, repeat_nb=10, subset_size=0.02)
mnn_silh_list <- run_silhouette(data_mnns, batch, repeat_nb=10, subset_size=0.02)

wilcox.test(x=unlist(silh_list), y=unlist(combat_silh_list), alternative='greater', correct=TRUE)
#p=0.6579
wilcox.test(x=unlist(silh_list), y=unlist(reg_silh_list), alternative='greater', correct=TRUE)
#p=0.2406
wilcox.test(x=unlist(silh_list), y=unlist(scano_silh_list), alternative='greater', correct=TRUE)
#p=0.03151
wilcox.test(x=unlist(silh_list), y=unlist(mnn_silh_list), alternative='greater', correct=TRUE)
#p=0.3153

# 这里在做图，先做出一个dataframe，后做箱型图
# rbind of silh
silh_df = data.frame(class=rep(c('base', 'combat', 'regress_out', 'scanorama', 'mnns'), 
                               each=length(silh_list)),
                     data=unlist(cbind(silh_list, combat_silh_list, reg_silh_list,
                                       scano_silh_list, mnn_silh_list)))
# do plot
# 这里可能这个limits=c(0, 0.1)上限要调整
p <- ggplot(silh_df, aes(class, data)) + geom_boxplot() + 
  theme_bw() +  
  scale_y_continuous(limits=c(-0.5, 0.5)) +
  labs(x='Methods', y='Average Silhouette Width', title='silh results')
p
pdf("ASW_batch.pdf");p;dev.off()

save.image("workspace.rdata")

cluster <- read.table("cluster_Baseline.txt", sep=',', header=T)$louvain
cluster_silh_list <- run_silhouette(data_base, cluster, repeat_nb=10, subset_size=0.02)
cluster_combat <- read.table("cluster_Combat.txt", sep=',', header=T)$louvain
cluster_combat_silh_list <- run_silhouette(data_combat, cluster_combat, repeat_nb=10, subset_size=0.02)
cluster_reg <- read.table("cluster_reg.txt", sep=',', header=T)$louvain
cluster_reg_silh_list <- run_silhouette(data_reg, cluster_reg, repeat_nb=10, subset_size=0.02)
cluster_scano <- read.table("cluster_Scanorama.txt", sep=',', header=T)$louvain
cluster_scano_silh_list <- run_silhouette(data_scano, cluster_scano, repeat_nb=10, subset_size=0.02)
cluster_mnn <- read.table("cluster_mnn.txt", sep=',', header=T)$louvain
cluster_mnn_silh_list <- run_silhouette(data_mnns, cluster_mnn, repeat_nb=10, subset_size=0.02)

wilcox.test(x=unlist(cluster_silh_list), y=unlist(cluster_combat_silh_list), alternative='greater', correct=TRUE)
#p-value = 5.413e-06
wilcox.test(x=unlist(cluster_silh_list), y=unlist(cluster_reg_silh_list), alternative='greater', correct=TRUE)
#p=5.413e-06
wilcox.test(x=unlist(cluster_silh_list), y=unlist(cluster_scano_silh_list), alternative='greater', correct=TRUE)
#p=1
wilcox.test(x=unlist(cluster_silh_list), y=unlist(cluster_mnn_silh_list), alternative='greater', correct=TRUE)
#p-value = 0.009272

cluster_silh_df = data.frame(class=rep(c('base', 'combat', 'regress_out', 'scanorama', 'mnns'), 
                               each=length(cluster_silh_list)),
                     data=unlist(cbind(cluster_silh_list, cluster_combat_silh_list, 
                                       cluster_reg_silh_list,
                                       cluster_scano_silh_list, cluster_mnn_silh_list)))
# do plot
# 这里可能这个limits=c(0, 0.1)上限要调整
p <- ggplot(cluster_silh_df, aes(class, data)) + geom_boxplot() + 
  theme_bw() +  
  scale_y_continuous(limits=c(-0.5, 0.5)) +
  labs(x='Methods', y='Average Silhouette Width', title='silh results')
p
pdf("ASW_cluster.pdf");p;dev.off()

save.image("workspace.rdata")