
library(sqldf)
library(vegan)
library(xtable)
library(parallel)
library(DescTools)
library(jointseg)
setwd("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/")



sample_annotation_file = read.csv("Sergi_RNAseqSampleinformation.csv")


expression_data = read.csv("Expression_dataset.csv", row.names = "probe")

geneswithid = read.csv("Genomic_mapping_Ensembl_to_geneifo.csv", row.names = "X")

match_vec = match(sample_annotation_file$Sample_ID,colnames(expression_data))

sample_annotation_file = sample_annotation_file[!(is.na(match_vec)),]
match_vec = match_vec[!(is.na(match_vec))]
expression_data = expression_data[,match_vec]


count_low_expression_genes = apply(expression_data,1,function(x){length(which(x<10))})

length(which(count_low_expression_genes<= 60))


################
#cell-line correction whole cell line_HL_robustSD
################

BT549_samples = which(sample_annotation_file$Cell_Line=="BT549")
HCC1806_samples = which(sample_annotation_file$Cell_Line=="HCC1806")
MDAMB231_samples = which(sample_annotation_file$Cell_Line=="MDAMB231")
RPE1_samples = which(sample_annotation_file$Cell_Line=="RPE1")
RPE1P53_samples = which(sample_annotation_file$Cell_Line=="RPE1(P53/)")


expression_data_BT549_samples_robustsd = apply(expression_data[,BT549_samples],1, function(x){estimateSd(x[which(!is.na(x))])})
expression_data_HCC1806_samples_robustsd = apply(expression_data[,HCC1806_samples],1, function(x){estimateSd(x[which(!is.na(x))])})
expression_data_MDAMB231_samples_robustsd = apply(expression_data[,MDAMB231_samples],1, function(x){estimateSd(x[which(!is.na(x))])})
expression_data_RPE1_samples_robustsd = apply(expression_data[,RPE1_samples],1, function(x){estimateSd(x[which(!is.na(x))])})
expression_data_RPE1P53_samples_robustsd = apply(expression_data[,RPE1P53_samples],1, function(x){estimateSd(x[which(!is.na(x))])})

expression_data_BT549_samples_hl = apply(expression_data[,BT549_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})
expression_data_HCC1806_samples_hl = apply(expression_data[,HCC1806_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})
expression_data_MDAMB231_samples_hl = apply(expression_data[,MDAMB231_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})
expression_data_RPE1_samples_hl = apply(expression_data[,RPE1_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})
expression_data_RPE1P53_samples_hl = apply(expression_data[,RPE1P53_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})

expression_data_adjusted = expression_data
expression_data_adjusted[,BT549_samples] = (expression_data[,BT549_samples] - expression_data_BT549_samples_hl)/expression_data_BT549_samples_robustsd
expression_data_adjusted[,HCC1806_samples] = (expression_data[,HCC1806_samples] - expression_data_HCC1806_samples_hl)/expression_data_HCC1806_samples_robustsd
expression_data_adjusted[,MDAMB231_samples] = (expression_data[,MDAMB231_samples] - expression_data_MDAMB231_samples_hl)/expression_data_MDAMB231_samples_robustsd
expression_data_adjusted[,RPE1_samples] = (expression_data[,RPE1_samples] - expression_data_RPE1_samples_hl)/expression_data_RPE1_samples_robustsd
expression_data_adjusted[,RPE1P53_samples] = (expression_data[,RPE1P53_samples] - expression_data_RPE1P53_samples_hl)/expression_data_RPE1P53_samples_robustsd

#checkstart
summary = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Step_1_summary_euclidean_distance_on_dox_corrected_data_21122018.csv")
number_of_na_samples_in_each_gene_expression_data_adjusted = apply(expression_data_adjusted,1,function(x){sum(is.na(x))})
rownames(summary) = summary$Ensemble_ID
summary_na_genes = summary[names(number_of_na_samples_in_each_gene_expression_data_adjusted)[which(number_of_na_samples_in_each_gene_expression_data_adjusted>0)],]
summary_na_genes_qualified = summary_na_genes[which(summary_na_genes$quality_controlled_genes==1),]
# ENSG00000113430 -- expressed only in one cell line rest are almost zero
plot(as.numeric(expression_data[summary_na_genes_qualified$Ensemble_ID,]))
#checkend

write.csv(expression_data_adjusted, file = "/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Expression_dataset_corrected_for_celllines_difference_HL_robustSD_07022019.csv")

#dox no-dox correction

empty_dox_samples = intersect(which(sample_annotation_file$Empty==1), which(sample_annotation_file$Treatment=="Doxycycline "))

empty_no_dox_samples = intersect(which(sample_annotation_file$Empty==1), which(sample_annotation_file$Treatment=="NO Dox"))
sum(colnames(expression_data_adjusted)==sample_annotation_file$Sample_ID)
expression_data_empty_dox_samples_robustsd = apply(expression_data_adjusted[,empty_dox_samples],1, function(x){estimateSd(x[which(!is.na(x))])})
expression_data_empty_no_dox_samples_robustsd = apply(expression_data_adjusted[,empty_no_dox_samples],1, function(x){estimateSd(x[which(!is.na(x))])})

expression_data_empty_dox_samples_hl = apply(expression_data_adjusted[,empty_dox_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})
expression_data_empty_no_dox_samples_hl = apply(expression_data_adjusted[,empty_no_dox_samples],1, function(x){HodgesLehmann(x[which(!is.na(x))])})

  

all_dox_samples = which(sample_annotation_file$Treatment%in%c("Doxycycline","Doxycycline "))

all_no_dox_samples = which(sample_annotation_file$Treatment=="NO Dox")
expression_data_dox_nodox_adjusted = expression_data_adjusted
expression_data_dox_nodox_adjusted[,all_dox_samples] = (expression_data_adjusted[,all_dox_samples] - expression_data_empty_dox_samples_hl)/expression_data_empty_dox_samples_robustsd
expression_data_dox_nodox_adjusted[,all_no_dox_samples] = (expression_data_adjusted[,all_no_dox_samples] - expression_data_empty_no_dox_samples_hl)/expression_data_empty_no_dox_samples_robustsd


#checkstart
number_of_na_samples_in_each_gene_expression_data_dox_nodox_adjusted = apply(expression_data_dox_nodox_adjusted,1,function(x){sum(is.na(x))})
summary_dox_nodox_na_genes = summary[names(number_of_na_samples_in_each_gene_expression_data_dox_nodox_adjusted)[which(number_of_na_samples_in_each_gene_expression_data_dox_nodox_adjusted>0)],]
summary_dox_nodox_na_genes_qualified = summary_na_genes[which(summary_na_genes$quality_controlled_genes==1),]
# ENSG00000113430 -- expressed only in one cell line rest are almost zero
plot(as.numeric(expression_data[summary_na_genes_qualified$Ensemble_ID,]))
#checkend



write.csv(expression_data_dox_nodox_adjusted, file = "/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Expression_dataset_corrected_for_celllines_difference_and_dox_nodox_HL_robustSD_07022019.csv")


