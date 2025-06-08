# import R package
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(SingleCellExperiment)
library(ggthemes)
library(stringr)
library(philentropy)
options(warn=-1)
setwd("D:/academic_relate_code_two/Nessie-main/DeepTXcopy5/DeepTX-main/inferObservedData")
# setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData")

source("utils.R")

result_dir = "result/IdUDMSO"
data_dir = "data/IdUDMSO"

# load data dmso
scRNA_dmso = read.table(sprintf("%s/DMSO_norm_filter.csv",data_dir),header=T, sep = ',')

ssa_dmso_dir =sprintf("%s/dmsoSSADistribution/",result_dir)
dmso_kl_df = calculate_KL(scRNA_dmso,ssa_dmso_dir)
dmso_kl_df["geneName"] = names(scRNA_dmso)

mean_var_res_dmso = calculate_mean_var(scRNA_dmso,ssa_dmso_dir)
dmso_kl_df["mean_val"] = mean_var_res_dmso$mean_val
dmso_kl_df["var_val"] = mean_var_res_dmso$var_val
write.csv(dmso_kl_df, sprintf("%s/dmso_kl.csv",result_dir),quote = F,row.names=FALSE) 

# load data IdU
scRNA_IdU = read.table(sprintf("%s/IdU_norm_filter.csv",data_dir),header=T, sep = ',')
ssa_IdU_dir =sprintf("%s/iduSSADistribution/",result_dir)
IdU_kl_df = calculate_KL(scRNA_IdU,ssa_IdU_dir)
IdU_kl_df["geneName"] = names(scRNA_IdU)

mean_var_res_IdU = calculate_mean_var(scRNA_IdU,ssa_IdU_dir)
IdU_kl_df["mean_val"] = mean_var_res_IdU$mean_val
IdU_kl_df["var_val"] = mean_var_res_IdU$var_val
write.csv(IdU_kl_df, sprintf("%s/IdU_kl.csv",result_dir),quote = F,row.names=FALSE) 

gene_estimated_matrix_idu = read.csv(file = sprintf("%s/idu_estimated_model_stats_prob.csv",result_dir), row.names =
                                       6)
gene_estimated_matrix_dmso = read.csv(file = sprintf("%s/dmso_estimated_model_stats_prob.csv",result_dir), row.names =
                                        6)
keep_gene_IdU = filter_gene_by_kl(scRNA_IdU, IdU_kl_df)
keep_gene_dmso = filter_gene_by_kl(scRNA_dmso, dmso_kl_df)
keep_intersect_gene = intersect(x = keep_gene_IdU, y = keep_gene_dmso)
gene_estimated_matrix_dmso = calculateBurstIndicator(gene_estimated_matrix_dmso)
gene_estimated_matrix_idu = calculateBurstIndicator(gene_estimated_matrix_idu)
gene_estimated_matrix_dmso = combineBurstIndicator(gene_estimated_matrix_dmso, gene_estimated_matrix_idu)
gene_estimated_matrix_dmso = calculateBurstDistance(gene_estimated_matrix_dmso)
gene_estimated_matrix_dmso = na.omit(gene_estimated_matrix_dmso[keep_intersect_gene,])
write.csv(gene_estimated_matrix_dmso, sprintf("%s/gene_estimated_stats_matrix.csv",result_dir),quote = F,row.names=TRUE) 
