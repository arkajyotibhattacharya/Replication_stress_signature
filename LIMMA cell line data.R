
library(gdata)

setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")



sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")
rownames(sample_annotation_file) = sample_annotation_file$Sample_ID
sample_annotation_file$Treatment = trim(sample_annotation_file$Treatment)

expression_data = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Expression_dataset_corrected_for_celllines_difference_and_dox_nodox_HL_robustSD_07022019.csv", row.names = "X")

geneswithid = read.csv("Genomic_mapping_Ensembl_to_geneifo.csv", row.names = "X")

# BiocManager::install("limma")
# BiocManager::install("statmod")
library(limma)



#MYC

sample_annotation_file_MYC = sample_annotation_file[which(!is.na(sample_annotation_file$Myc_analysis_v1)),]
match_vec = match(sample_annotation_file_MYC$Sample_ID,colnames(expression_data))
sample_annotation_file_MYC = sample_annotation_file_MYC[!(is.na(match_vec)),]

match_vec = match_vec[!(is.na(match_vec))]
expression_data_MYC = expression_data[,match_vec]
# sample_annotation_file_MYC$Treatment = trim(sample_annotation_file_MYC$Treatment)

sample_annotation_file_MYC = sample_annotation_file_MYC[,c("Treatment", "Cell_Line", "Myc_analysis_v1")]
sample_annotation_file_MYC$Treatment = as.factor(sample_annotation_file_MYC$Treatment)
sample_annotation_file_MYC$Cell_Line = as.factor(sample_annotation_file_MYC$Cell_Line)
# sample_annotation_file_MYC$Myc_analysis_v1 = as.factor(sample_annotation_file_MYC$Myc_analysis_v1)


all(rownames(sample_annotation_file_MYC) %in% colnames(expression_data_MYC))
all(rownames(sample_annotation_file_MYC) == colnames(expression_data_MYC))

fit <- lmFit(expression_data_MYC, design=sample_annotation_file_MYC$Myc_analysis_v1)
fit <- eBayes(fit)
a = topTable(fit,sort="none",n=Inf)

write.csv(as.data.frame(a), 
          file="/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Deseq2/Myc_limma_cell_line.csv")

#CCNE1

sample_annotation_file_CCNE1 = sample_annotation_file[which(!is.na(sample_annotation_file$CyclinE_analysis_v1)),]
match_vec = match(sample_annotation_file_CCNE1$Sample_ID,colnames(expression_data))
sample_annotation_file_CCNE1 = sample_annotation_file_CCNE1[!(is.na(match_vec)),]

match_vec = match_vec[!(is.na(match_vec))]
expression_data_CCNE1 = expression_data[,match_vec]
# sample_annotation_file_CCNE1$Treatment = trim(sample_annotation_file_CCNE1$Treatment)

sample_annotation_file_CCNE1 = sample_annotation_file_CCNE1[,c("Treatment", "Cell_Line", "CyclinE_analysis_v1")]
sample_annotation_file_CCNE1$Treatment = as.factor(sample_annotation_file_CCNE1$Treatment)
sample_annotation_file_CCNE1$Cell_Line = as.factor(sample_annotation_file_CCNE1$Cell_Line)
# sample_annotation_file_CCNE1$CyclinE_analysis_v1 = as.factor(sample_annotation_file_CCNE1$CyclinE_analysis_v1)


all(rownames(sample_annotation_file_CCNE1) %in% colnames(expression_data_CCNE1))
all(rownames(sample_annotation_file_CCNE1) == colnames(expression_data_CCNE1))

fit <- lmFit(expression_data_CCNE1, design=sample_annotation_file_CCNE1$CyclinE_analysis_v1)
fit <- eBayes(fit)
a = topTable(fit,sort="none",n=Inf)

write.csv(as.data.frame(a), 
          file="/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Deseq2/CCNE1_limma_cell_line.csv")

#CDC25A

sample_annotation_file_CDC25A = sample_annotation_file[which(!is.na(sample_annotation_file$Cdc25a_analysis_v1)),]
match_vec = match(sample_annotation_file_CDC25A$Sample_ID,colnames(expression_data))
sample_annotation_file_CDC25A = sample_annotation_file_CDC25A[!(is.na(match_vec)),]

match_vec = match_vec[!(is.na(match_vec))]
expression_data_CDC25A = expression_data[,match_vec]
# sample_annotation_file_CDC25A$Treatment = trim(sample_annotation_file_CDC25A$Treatment)

sample_annotation_file_CDC25A = sample_annotation_file_CDC25A[,c("Treatment", "Cell_Line", "Cdc25a_analysis_v1")]
sample_annotation_file_CDC25A$Treatment = as.factor(sample_annotation_file_CDC25A$Treatment)
sample_annotation_file_CDC25A$Cell_Line = as.factor(sample_annotation_file_CDC25A$Cell_Line)
# sample_annotation_file_CDC25A$Cdc25a_analysis_v1 = as.factor(sample_annotation_file_CDC25A$Cdc25a_analysis_v1)


all(rownames(sample_annotation_file_CDC25A) %in% colnames(expression_data_CDC25A))
all(rownames(sample_annotation_file_CDC25A) == colnames(expression_data_CDC25A))

fit <- lmFit(expression_data_CDC25A, design=sample_annotation_file_CDC25A$Cdc25a_analysis_v1)
fit <- eBayes(fit)
a = topTable(fit,sort="none",n=Inf)

write.csv(as.data.frame(a), 
          file="/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Deseq2/CDC25A_limma_cell_line.csv")


