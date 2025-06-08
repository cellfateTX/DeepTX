# import R package
options(warn=-1)
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(ggthemes)
library(stringr)
# library(philentropy)
library(SingleCellExperiment)
setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData/data/IdUDMSO")
theme_pubr_new = function(){
  theme_pubr(base_family = "Helvetica", base_size = 14, legend = "right") %+replace%
    theme(
      axis.title = element_text(face = "bold", size = 18), 
      legend.text = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 16)
    )
}
color.feedback <- c("#3AA438", "#547CBE", "#E71419")
dmso_idu_cell_gene_matrix = read.csv(file = "DMSO_IdU_scRNA_counts.csv", row.names = 1)
seurat_scal_dmso_idu = matrix(data = NA, nrow = 27998, ncol = 1556)
# dmso_idu_cell_gene_matrix=log10(dmso_idu_cell_gene_matrix+1)
for (i in 1:1556){
  seurat_scal_dmso_idu[,i] = round((((dmso_idu_cell_gene_matrix[,i]*5000)/sum(dmso_idu_cell_gene_matrix[,i]))*10))
}

dmso_idu_filter = rowMeans(seurat_scal_dmso_idu) > 1
rownames(seurat_scal_dmso_idu) = rownames(dmso_idu_cell_gene_matrix)
seurat_dmso_mean_var = data.frame(dmso_mean = rowMeans(seurat_scal_dmso_idu[dmso_idu_filter,1:812]), dmso_var = (rowSds(as.matrix(seurat_scal_dmso_idu[dmso_idu_filter,1:812])))^2, dmso_cv = ((rowSds(seurat_scal_dmso_idu[dmso_idu_filter,1:812]))^2)/((rowMeans(seurat_scal_dmso_idu[dmso_idu_filter,1:812]))^2))
seurat_idu_mean_var = data.frame(idu_mean = rowMeans(seurat_scal_dmso_idu[dmso_idu_filter,813:1556]), idu_var = (rowSds(as.matrix(seurat_scal_dmso_idu[dmso_idu_filter,813:1556])))^2, idu_cv = ((rowSds(seurat_scal_dmso_idu[dmso_idu_filter,813:1556]))^2)/((rowMeans(seurat_scal_dmso_idu[dmso_idu_filter,813:1556]))^2))
seurat_idu_dmso_mean_var = cbind(seurat_idu_mean_var,seurat_dmso_mean_var)

dmso_scal_data = seurat_scal_dmso_idu[dmso_idu_filter,1:812]
idu_scal_data = seurat_scal_dmso_idu[dmso_idu_filter,813:1556]
rownames(dmso_scal_data) =rownames(seurat_scal_dmso_idu[dmso_idu_filter,])
rownames(idu_scal_data) =rownames(seurat_scal_dmso_idu[dmso_idu_filter,])

dmso_data = as.data.frame(dmso_scal_data)
write.table(t(dmso_data),"DMSO_norm_filter.csv",row.names=FALSE,col.names=TRUE,sep=",")

# save the normlize and filter data
idu_data = as.data.frame(idu_scal_data)
write.table(t(idu_data),"IdU_norm_filter.csv",row.names=FALSE,col.names=TRUE,sep=",")
