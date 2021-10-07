
#######
#Cdc25a
#######
rm(list = ls())

library(sqldf)
library(vegan)
library(xtable)
library(parallel)

setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")


sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")

expression_data = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Expression_dataset_corrected_for_celllines_difference_and_dox_nodox_HL_robustSD_07022019.csv", row.names = "X")
summary = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Step_1_summary_euclidean_distance_on_dox_corrected_data_21122018.csv")

na_count_expression_data = apply(expression_data,1,function(x){sum(is.na(x))})


expression_data = expression_data[which(na_count_expression_data==0),]

t_expression_data = t(expression_data)


for(cellline in as.character(unique(sample_annotation_file$Cell_Line)))
{
  setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")
  
  
  sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")
  sample_annotation_file = sample_annotation_file[which(sample_annotation_file$Cell_Line==cellline),]
  sample_annotation_file = sample_annotation_file[which(sample_annotation_file$Type1%in%c("Empty", "Cdc25a")),]
  match_vec = match(sample_annotation_file$Sample_ID,rownames(t_expression_data))
  
  sample_annotation_file_v1 = sample_annotation_file[!(is.na(match_vec)),]
  match_vec = match_vec[!(is.na(match_vec))]
  t_expression_data_v1 = t_expression_data[match_vec,]
  
  t_expression_data_v2 = t_expression_data_v1[which(!is.na(sample_annotation_file_v1$Cdc25a_analysis_v1)),]
  sample_annotation_file_v2 = sample_annotation_file_v1[which(!is.na(sample_annotation_file_v1$Cdc25a_analysis_v1)),]
  
  
  permanova_per_gene = function(i,data)
  {
    
    permanova_result_1_euclidean = adonis(t_expression_data_v2[,i] ~ Cdc25a_analysis_v1
                                          , data = sample_annotation_file_v2
                                          , method = "euclidean"
                                          , permutations = 5000
    )
    coef_mat = as.matrix(permanova_result_1_euclidean$aov.tab)
    
    HR_status_coefs_1 =c(coef_mat[which(rownames(coef_mat)%in%c("Cell_Line","Cdc25a_analysis_v1")),],permanova_result_1_euclidean$coefficients[which(rownames(permanova_result_1_euclidean$coefficients)=="Cdc25a_analysis_v1")])
    names(HR_status_coefs_1)[7] = "coefficient"
    return(HR_status_coefs_1)
  }
  
  time1 = proc.time()[3]
  no_cores <- 10
  cl <- makeCluster(no_cores, type = "FORK")
  
  time2 = proc.time()[3]
  set.seed(1235)
  # x = parLapply(cl, c(1:160),permanova_per_gene,sample_annotation_file_v2)
  x = parLapply(cl, c(1:dim(t_expression_data_v2)[2]),permanova_per_gene,sample_annotation_file_v2)
  
  stopCluster(cl)
  # df <- data.frame(matrix(unlist(x), nrow=dim(t_expression_data_v2)[2], byrow=T),stringsAsFactors=FALSE)
  df <- as.data.frame(do.call(rbind,lapply(x,matrix,ncol=7,byrow=FALSE)))
  colnames(df) = names(x[[1]])
  rownames(df) = rownames(expression_data)
  
  print((proc.time()[3]-time1)/60)
  
  write.csv(df, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Cdc25a_step1_euclidean_non_na_genes_",gsub("[/]","_",cellline),"_09022019.csv",sep = ""))
  
}

print("Cdc25a_euclidean")

#######
#CyclinE
#######
rm(list = ls())

library(sqldf)
library(vegan)
library(xtable)
library(parallel)

setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")


sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")

expression_data = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Expression_dataset_corrected_for_celllines_difference_and_dox_nodox_HL_robustSD_07022019.csv", row.names = "X")
summary = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Step_1_summary_euclidean_distance_on_dox_corrected_data_21122018.csv")

na_count_expression_data = apply(expression_data,1,function(x){sum(is.na(x))})


expression_data = expression_data[which(na_count_expression_data==0),]

t_expression_data = t(expression_data)


for(cellline in as.character(unique(sample_annotation_file$Cell_Line)))
{
  setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")
  
  
  sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")
  sample_annotation_file = sample_annotation_file[which(sample_annotation_file$Cell_Line==cellline),]
  sample_annotation_file = sample_annotation_file[which(sample_annotation_file$Type1%in%c("Empty", "CyclinE")),]
  match_vec = match(sample_annotation_file$Sample_ID,rownames(t_expression_data))
  
  sample_annotation_file_v1 = sample_annotation_file[!(is.na(match_vec)),]
  match_vec = match_vec[!(is.na(match_vec))]
  t_expression_data_v1 = t_expression_data[match_vec,]
  
  t_expression_data_v2 = t_expression_data_v1[which(!is.na(sample_annotation_file_v1$CyclinE_analysis_v1)),]
  sample_annotation_file_v2 = sample_annotation_file_v1[which(!is.na(sample_annotation_file_v1$CyclinE_analysis_v1)),]
  
  
  permanova_per_gene = function(i,data)
  {
    
    permanova_result_1_euclidean = adonis(t_expression_data_v2[,i] ~ CyclinE_analysis_v1
                                          , data = sample_annotation_file_v2
                                          , method = "euclidean"
                                          , permutations = 5000
    )
    coef_mat = as.matrix(permanova_result_1_euclidean$aov.tab)
    
    HR_status_coefs_1 =c(coef_mat[which(rownames(coef_mat)%in%c("Cell_Line","CyclinE_analysis_v1")),],permanova_result_1_euclidean$coefficients[which(rownames(permanova_result_1_euclidean$coefficients)=="CyclinE_analysis_v1")])
    names(HR_status_coefs_1)[7] = "coefficient"
    return(HR_status_coefs_1)
  }
  
  time1 = proc.time()[3]
  no_cores <- 10
  cl <- makeCluster(no_cores, type = "FORK")
  
  time2 = proc.time()[3]
  set.seed(1235)
  # x = parLapply(cl, c(1:160),permanova_per_gene,sample_annotation_file_v2)
  x = parLapply(cl, c(1:dim(t_expression_data_v2)[2]),permanova_per_gene,sample_annotation_file_v2)
  
  stopCluster(cl)
  # df <- data.frame(matrix(unlist(x), nrow=dim(t_expression_data_v2)[2], byrow=T),stringsAsFactors=FALSE)
  df <- as.data.frame(do.call(rbind,lapply(x,matrix,ncol=7,byrow=FALSE)))
  colnames(df) = names(x[[1]])
  rownames(df) = rownames(expression_data)
  
  print((proc.time()[3]-time1)/60)
  
  write.csv(df, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/CyclinE_step1_euclidean_non_na_genes_",gsub("[/]","_",cellline),"_09022019.csv",sep = ""))
  
}

print("CyclinE_euclidean")
#######
#Myc
#######
rm(list = ls())

library(sqldf)
library(vegan)
library(xtable)
library(parallel)

setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")


sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")

expression_data = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Expression_dataset_corrected_for_celllines_difference_and_dox_nodox_HL_robustSD_07022019.csv", row.names = "X")
summary = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Step_1_summary_euclidean_distance_on_dox_corrected_data_21122018.csv")

na_count_expression_data = apply(expression_data,1,function(x){sum(is.na(x))})


expression_data = expression_data[which(na_count_expression_data==0),]

t_expression_data = t(expression_data)


for(cellline in as.character(unique(sample_annotation_file$Cell_Line)))
{
  setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")
  
  
  sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")
  sample_annotation_file = sample_annotation_file[which(sample_annotation_file$Cell_Line==cellline),]
  sample_annotation_file = sample_annotation_file[which(sample_annotation_file$Type1%in%c("Empty", "Myc")),]
  match_vec = match(sample_annotation_file$Sample_ID,rownames(t_expression_data))
  
  sample_annotation_file_v1 = sample_annotation_file[!(is.na(match_vec)),]
  match_vec = match_vec[!(is.na(match_vec))]
  t_expression_data_v1 = t_expression_data[match_vec,]
  
  t_expression_data_v2 = t_expression_data_v1[which(!is.na(sample_annotation_file_v1$Myc_analysis_v1)),]
  sample_annotation_file_v2 = sample_annotation_file_v1[which(!is.na(sample_annotation_file_v1$Myc_analysis_v1)),]
  
  
  permanova_per_gene = function(i,data)
  {
    
    permanova_result_1_euclidean = adonis(t_expression_data_v2[,i] ~ Myc_analysis_v1
                                          , data = sample_annotation_file_v2
                                          , method = "euclidean"
                                          , permutations = 5000
    )
    coef_mat = as.matrix(permanova_result_1_euclidean$aov.tab)
    
    HR_status_coefs_1 =c(coef_mat[which(rownames(coef_mat)%in%c("Cell_Line","Myc_analysis_v1")),],permanova_result_1_euclidean$coefficients[which(rownames(permanova_result_1_euclidean$coefficients)=="Myc_analysis_v1")])
    names(HR_status_coefs_1)[7] = "coefficient"
    return(HR_status_coefs_1)
  }
  
  time1 = proc.time()[3]
  no_cores <- 10
  cl <- makeCluster(no_cores, type = "FORK")
  
  time2 = proc.time()[3]
  set.seed(1235)
  # x = parLapply(cl, c(1:160),permanova_per_gene,sample_annotation_file_v2)
  x = parLapply(cl, c(1:dim(t_expression_data_v2)[2]),permanova_per_gene,sample_annotation_file_v2)
  
  stopCluster(cl)
  # df <- data.frame(matrix(unlist(x), nrow=dim(t_expression_data_v2)[2], byrow=T),stringsAsFactors=FALSE)
  df <- as.data.frame(do.call(rbind,lapply(x,matrix,ncol=7,byrow=FALSE)))
  colnames(df) = names(x[[1]])
  rownames(df) = rownames(expression_data)
  
  print((proc.time()[3]-time1)/60)
  
  write.csv(df, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Myc_step1_euclidean_non_na_genes_",gsub("[/]","_",cellline),"_09022019.csv",sep = ""))
  
}

print("Myc_euclidean")

