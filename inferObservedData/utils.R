
filter_gene_by_kl <- function(gene_df,kl_list_data,
threshold=0.14){
  gene_names= names(gene_df)
  kl_keep= kl_list_data<threshold
  keep_gene= gene_names[kl_keep]
  return(keep_gene)
}

calculateBurstIndicator <- function(paramter_data){
  

   paramter_data$tau = paramter_data[,9]/paramter_data[,5] 
  paramter_data$tau_on = paramter_data[,1]/paramter_data[,2] 
  paramter_data$tau_off =paramter_data[,3]/paramter_data[,4] 
  return(paramter_data)
  
}

combineBurstIndicator<- function(gene_estimated_matrix_dmso,gene_estimated_matrix_idu){
  
  gene_estimated_matrix_dmso$tau_on_idu = gene_estimated_matrix_idu$tau_on
  gene_estimated_matrix_dmso$tau_off_idu = gene_estimated_matrix_idu$tau_off
  gene_estimated_matrix_dmso$bs_idu = gene_estimated_matrix_idu$bs
  gene_estimated_matrix_dmso$bf_idu = gene_estimated_matrix_idu$bf
  gene_estimated_matrix_dmso$x5_idu = gene_estimated_matrix_idu$x5
  gene_estimated_matrix_dmso$mean_true_idu = gene_estimated_matrix_idu$mean_true
  gene_estimated_matrix_dmso$tau_idu = gene_estimated_matrix_idu$tau
  gene_estimated_matrix_dmso$var_true_idu = gene_estimated_matrix_idu$var_true
  return(gene_estimated_matrix_dmso)
  
}

calculateBurstDistance <- function(gene_estimated_matrix_dmso){
  gene_estimated_matrix_dmso$tau_distance = (gene_estimated_matrix_dmso[,"tau"]-gene_estimated_matrix_dmso[,"tau_idu"])/gene_estimated_matrix_dmso[,"tau"]
  gene_estimated_matrix_dmso$bs_distance = (gene_estimated_matrix_dmso[,"bs"]-gene_estimated_matrix_dmso[,"bs_idu"])/gene_estimated_matrix_dmso[,"bs"]
  gene_estimated_matrix_dmso$bf_distance = (gene_estimated_matrix_dmso[,"bf"]-gene_estimated_matrix_dmso[,"bf_idu"])/gene_estimated_matrix_dmso[,"bf"]
  gene_estimated_matrix_dmso$tau_off_distance = (gene_estimated_matrix_dmso[,"tau_off"]-gene_estimated_matrix_dmso[,"tau_off_idu"])/gene_estimated_matrix_dmso[,"tau_off"]
  gene_estimated_matrix_dmso$tau_on_distance = (gene_estimated_matrix_dmso[,"tau_on"]-gene_estimated_matrix_dmso[,"tau_on_idu"])/gene_estimated_matrix_dmso[,"tau_on"]
  gene_estimated_matrix_dmso$var_distance = (gene_estimated_matrix_dmso[,"var_true"]-gene_estimated_matrix_dmso[,"var_true_idu"])/gene_estimated_matrix_dmso[,"var_true"]
  gene_estimated_matrix_dmso$mean_distance = (gene_estimated_matrix_dmso[,"mean_true"]-gene_estimated_matrix_dmso[,"mean_true_idu"])/gene_estimated_matrix_dmso[,"mean_true"]
  return(gene_estimated_matrix_dmso)
}

filterBurstDistance <-function(burst_item,gene_estimated_matrix_dmso,burst_type){
  if(burst_type=="negative"){
    median_burst_item = median(gene_estimated_matrix_dmso[gene_estimated_matrix_dmso[,burst_item]<=0,][,burst_item],na.rm=TRUE)
    gene_estimated_matrix_dmso_sub = gene_estimated_matrix_dmso[gene_estimated_matrix_dmso[,burst_item]<=median_burst_item,]   
  }else{
    median_burst_item = median(gene_estimated_matrix_dmso[gene_estimated_matrix_dmso[,burst_item]>=0,][,burst_item],na.rm=TRUE)
    gene_estimated_matrix_dmso_sub = gene_estimated_matrix_dmso[gene_estimated_matrix_dmso[,burst_item]>=median_burst_item,]
  }
  return(gene_estimated_matrix_dmso_sub)
}

calculate_KL <- function(RNA_histogram_DATA,ssa_idu_dir) {
    kl_list = c()
    gene_names= names(RNA_histogram_DATA)
#     for (i in 1:10)   
    for (i in 1:length(gene_names))   
      {
          col_name= gene_names[i]
          NN_predicted_file_name = str_glue(ssa_idu_dir,"distribution_{i}.csv")
          if(file.exists(NN_predicted_file_name)) 
             {
                 RNA_prob_DATA=read.table(NN_predicted_file_name,header=T, sep = ',')    
                a = table(cut(RNA_histogram_DATA[,col_name],c(seq(-1,max(RNA_histogram_DATA[,col_name]),1))))
                a_prob = prop.table(a)
                RNA_prob_DATA[,"x1"]=RNA_prob_DATA[,"x1"]/sum(RNA_prob_DATA[,"x1"])
                dataTemp <- list(a_prob, RNA_prob_DATA[,"x1"])
                dataTemp2 <- do.call(rbind,lapply(lapply(dataTemp, unlist), `length<-`, max(lengths(dataTemp))))
                dataTemp2[is.na(dataTemp2)] = 0  
              kl_list[i]=KL(dataTemp2,unit='log')
             }
          else{
              kl_list[i]=-1
          }
      }
    kl_list_df <- data.frame(kl_list) 
    return(kl_list_df)
}

calculate_mean_var <- function(RNA_histogram_DATA,ssa_data_dir) {
        mean_val = c()
        var_val = c()
        gene_names= names(RNA_histogram_DATA)
        for (i in 1:length(gene_names)) {
             NN_predicted_file_name = str_glue(ssa_data_dir,"distribution_{i}.csv")
              if(file.exists(NN_predicted_file_name)) 
                 {
                     RNA_prob_DATA=read.table(NN_predicted_file_name,header=T, sep = ',')    
                  x = c(0:(dim(RNA_prob_DATA)-1))
                  Px  = RNA_prob_DATA["x1"]
                  mean_val[i] = sum(x*Px)
                  var_val[i] = sum((x-mean_val[i])^2*Px)  
                   }
              else {
                      mean_val[i]=0.001
                      var_val[i] = 0.001
                  }
        }
    return( list(mean_val=mean_val, var_val=var_val))
} 
