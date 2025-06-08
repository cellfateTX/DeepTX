# import R package
options(warn=-1)
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")
library(Matrix)
library(ggthemes)
library(stringr)
library(philentropy)
library(SingleCellExperiment)
setwd("D:/academic_relate_code_two/Nessie-main/DeepTX/inferObservedData/data")

filter_data = function(seq_data){
  seq_data_sub = seq_data[,colSums(seq_data>=1)>=1200]
  seq_data_sub = seq_data_sub[rowSums(seq_data_sub>1)>=40,]
  return(seq_data_sub)
}
rko.5FU_only.dose.zero = read.csv(file = "FUcommon/rkoDoseZero.csv", row.names = 1)
rko.5FU_only.dose.ten = read.csv(file = "FUcommon/rkoDoseTen.csv", row.names = 1)
rko.5FU_only.dose.fifty = read.csv(file = "FUcommon/rkoDoseFifty.csv", row.names = 1)
rko.5FU_only.dose.hundred = read.csv(file = "FUcommon/rkoDoseHundred.csv", row.names = 1)

rko.5FU_only.dose.zero = filter_data(rko.5FU_only.dose.zero)
rko.5FU_only.dose.ten = filter_data(rko.5FU_only.dose.ten)
rko.5FU_only.dose.fifty = filter_data(rko.5FU_only.dose.fifty)

zero_names = row.names(rko.5FU_only.dose.zero)
ten_names = row.names(rko.5FU_only.dose.ten)
zero_ten_intersect_name = intersect(zero_names,ten_names)

fifty_names = row.names(rko.5FU_only.dose.fifty)
data_type_li = c("doseTen","doseFifty","doseHundred")
data_type = data_type_li[1]
if(data_type == data_type_li[1]) {
  zero_type = "doseZero"
  zero_names = row.names(rko.5FU_only.dose.zero)
  ten_names = row.names(rko.5FU_only.dose.ten)
  zero_ten_intersect_name = intersect(zero_names,ten_names)
  rko.5FU_only.dose.zero = rko.5FU_only.dose.zero[zero_ten_intersect_name,]
  rko.5FU_only.dose.ten = rko.5FU_only.dose.ten[zero_ten_intersect_name,]
  rko.dose.nonZero = rko.5FU_only.dose.ten
}else if (data_type == data_type_li[2]) {
  zero_type = "doseOne"
  zero_names = row.names(rko.5FU_only.dose.zero)
  fifty_names = row.names(rko.5FU_only.dose.fifty)
  zero_fifty_intersect_name = intersect(zero_names,fifty_names)
  rko.5FU_only.dose.zero = rko.5FU_only.dose.zero[zero_fifty_intersect_name,]
  rko.5FU_only.dose.fifty = rko.5FU_only.dose.fifty[zero_fifty_intersect_name,]
  rko.dose.nonZero = rko.5FU_only.dose.fifty
}else if (data_type == data_type_li[3]) {
  rko.dose.nonZero = rko.5FU_only.dose.hundred
}
rko_dose_zero_nonZero = cbind(rko.5FU_only.dose.zero,rko.dose.nonZero)

end_i = dim(rko.5FU_only.dose.zero)[2]
end_j = dim(rko.5FU_only.dose.zero)[2] + dim(rko.dose.nonZero)[2]
seurat_rko_dose_zero_nonZero = matrix(data = NA, nrow = dim(rko_dose_zero_nonZero)[1], ncol = end_j)

for (i in 1:end_j){
  seurat_rko_dose_zero_nonZero[,i] = round((((rko_dose_zero_nonZero[,i]*2500)/sum(rko_dose_zero_nonZero[,i]))))
  # seurat_rko_dose_zero_nonZero[,i] =rko_dose_zero_nonZero[,i]
}

dose_zero_nonZero_filter = (rowMeans(seurat_rko_dose_zero_nonZero) > 1)
rownames(seurat_rko_dose_zero_nonZero) =rownames(rko_dose_zero_nonZero)
file_name_zero = sprintf("Ten5FU/%s_norm_filter.csv",zero_type)
file_name_nonZero = sprintf("Ten5FU/%s_norm_filter.csv",data_type)

write.table(t(seurat_rko_dose_zero_nonZero[dose_zero_nonZero_filter,1:end_i]),file_name_zero,row.names=FALSE,col.names=TRUE,sep=",")
write.table(t(seurat_rko_dose_zero_nonZero[dose_zero_nonZero_filter,(end_i+1):end_j]),file_name_nonZero,row.names=FALSE,col.names=TRUE,sep=",")
print("done")
