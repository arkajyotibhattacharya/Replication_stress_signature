
library(sqldf)
library(DescTools)
library(ggplot2)
library(ggrepel)
library(data.table)
#load GEO data

tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GEO/Data/GPL570__Affy_hgu133plus2_NoDuplicatesamples_CleanedIdentifiers_RMA-sketch_NormalAndCancer_QCed.txt", sep = "\t" , header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
GEO_mrna_expression <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GEO/Data/GPL570__Affy_hgu133plus2_NoDuplicatesamples_CleanedIdentifiers_RMA-sketch_NormalAndCancer_QCed.txt", sep = "\t" , header = TRUE, colClasses = classes, row.names = "X")

genomic_mapping_file = fread("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.txt", sep = "\t" , header = TRUE)
genomic_mapping_file = genomic_mapping_file[which(genomic_mapping_file$top_probe_indicator==1),]

rownames(genomic_mapping_file) = genomic_mapping_file$PROBESET
GEO_mrna_expression = GEO_mrna_expression[rownames(genomic_mapping_file),]
rownames(GEO_mrna_expression) = genomic_mapping_file$ENTREZID

#load color codes
color_code = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/color_code_template.csv")
for(i in 1:45)
{
  sub_data = colMeans(color_code[sample(1:dim(color_code)[1],14),])
  color_code = rbind(color_code, sub_data)
}

#load sample annotation
sample_annotation = read.csv("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GEO/Data/GPL570__Sample_To_TumorType_with_common_cancer_type_mapping_GEO_TCGA.csv")
sample_annotation = sample_annotation[which(sample_annotation$GSM_IDENTIFIER%in%colnames(GEO_mrna_expression)),]
unique_tissue_type = unique(sample_annotation$TYPE2)
sample_annotation_normals = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/GSE7307_GEO_Sample_Info.csv")
sample_annotation_normals = sample_annotation_normals[which(as.character(sample_annotation_normals$Disease.Normal.or.Treatment..C.)=="normal"),]
length(which(sample_annotation$GSM_IDENTIFIER%in%as.character(sample_annotation_normals$Sample)))

sample_annotation = sample_annotation[-c(intersect(which(sample_annotation$TYPE3=="Normal")
                                                   ,which(!sample_annotation$GSM_IDENTIFIER%in%as.character(sample_annotation_normals$Sample)))),]

GEO_mrna_expression = GEO_mrna_expression[,as.character(sample_annotation$GSM_IDENTIFIER)]

replication_stress_signature = c("65265",
                                 "55661",
                                 "27304",
                                 "51678",
                                 "55226",
                                 "197407")

GEO_mrna_expression_replication_stress_signature = GEO_mrna_expression[replication_stress_signature,]
GEO_mrna_expression_replication_stress_signature = t(scale(t(GEO_mrna_expression_replication_stress_signature)))

GEO_mrna_expression_replication_stress_signature_score = colMeans(GEO_mrna_expression_replication_stress_signature)
Sample_identifier = names(GEO_mrna_expression_replication_stress_signature_score)
score_data = as.data.frame(GEO_mrna_expression_replication_stress_signature_score)
score_data$Sample_identifier = Sample_identifier
rownames(sample_annotation) = sample_annotation$GSM_IDENTIFIER
sample_annotation = sample_annotation[rownames(score_data),]
  score_data$Tumor_type = sample_annotation$TYPE3
  score_data$normalized_GEO_mrna_expression_replication_stress_signature_score = (score_data$GEO_mrna_expression_replication_stress_signature_score - min(score_data$GEO_mrna_expression_replication_stress_signature_score))/(max(score_data$GEO_mrna_expression_replication_stress_signature_score) - min(score_data$GEO_mrna_expression_replication_stress_signature_score))
  

write.table(score_data, file = "/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Replication_stress_signature_score_per_sample_GEO_18072019.txt", sep = "\t", row.names = FALSE)


unique_all_combinations = as.character(unique(score_data$Tumor_type))
color_sample<-rgb(red = color_code[c(1:length(unique_all_combinations)),1],green = color_code[c(1:length(unique_all_combinations)),2], blue = color_code[c(1:length(unique_all_combinations)),3])
color_code = as.data.frame(cbind(unique_all_combinations, color_sample))
sample_annotation_file = sqldf("select a.*, b.color_sample from score_data a
                               left join color_code b 
                               on a.Tumor_type = b.unique_all_combinations
                               order by Tumor_type")

pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Replication stress signatrue score distribution/GEO replication stress signature score_18072019.pdf", sep = "")
    , height = 8, width = 16)
  ggplot(sample_annotation_file, aes(x=Tumor_type, y=GEO_mrna_expression_replication_stress_signature_score)) + 
  geom_dotplot(binaxis="y",stackdir="center", dotsize=0.1, color=sample_annotation_file$color_sample, binwidth = 0.0015) + 
  geom_boxplot(fill=color_sample, color=color_sample, alpha=0.4) + scale_y_continuous(name = "VALUE") + 
  scale_x_discrete(name = "Cancer types") +  geom_hline(yintercept=0) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  ggtitle(paste("GEO replication stress signature score"))
dev.off()

pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Replication stress signatrue score distribution/GEO replication stress signature score normalized_18072019.pdf", sep = "")
    , height = 8, width = 16)
ggplot(sample_annotation_file, aes(x=Tumor_type, y=normalized_GEO_mrna_expression_replication_stress_signature_score)) + 
  geom_dotplot(binaxis="y",stackdir="center", dotsize=0.1, color=sample_annotation_file$color_sample, binwidth = 0.0015) + 
  geom_boxplot(fill=color_sample, color=color_sample, alpha=0.4) + scale_y_continuous(name = "VALUE") + 
  scale_x_discrete(name = "Cancer types") +  geom_hline(yintercept=0) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  ggtitle(paste("GEO replication stress signature score"))
dev.off()


score_summary = sqldf("select Tumor_type, avg(GEO_mrna_expression_replication_stress_signature_score) as GEO_mrna_expression_replication_stress_signature_mean_score,
median(GEO_mrna_expression_replication_stress_signature_score) as GEO_mrna_expression_replication_stress_signature_median_score,
max(GEO_mrna_expression_replication_stress_signature_score) as GEO_mrna_expression_replication_stress_signature_max_score,
min(GEO_mrna_expression_replication_stress_signature_score) as GEO_mrna_expression_replication_stress_signature_min_score,
avg(normalized_GEO_mrna_expression_replication_stress_signature_score) as normalized_GEO_mrna_expression_replication_stress_signature_mean_score,
median(normalized_GEO_mrna_expression_replication_stress_signature_score) as normalized_GEO_mrna_expression_replication_stress_signature_median_score,
max(normalized_GEO_mrna_expression_replication_stress_signature_score) as normalized_GEO_mrna_expression_replication_stress_signature_max_score,
                      min(normalized_GEO_mrna_expression_replication_stress_signature_score) as normalized_GEO_mrna_expression_replication_stress_signature_min_score
                      from score_data
                      group by 1")

#Should this be adjusted with respect to normals? like normal should have zero?

for(i in 1:dim(score_summary)[1])
{
  score_summary$GEO_mrna_expression_replication_stress_signature_HL_score[i] = HodgesLehmann(score_data$GEO_mrna_expression_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
  score_summary$GEO_mrna_expression_replication_stress_signature_sd_score[i] = sd(score_data$GEO_mrna_expression_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
  score_summary$normalized_GEO_mrna_expression_replication_stress_signature_HL_score[i] = HodgesLehmann(score_data$normalized_GEO_mrna_expression_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
  score_summary$normalized_GEO_mrna_expression_replication_stress_signature_sd_score[i] = sd(score_data$normalized_GEO_mrna_expression_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
}

write.table(score_summary, file = "/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Replication_stress_signature_score_summary_GEO_18072019.txt", sep = "\t", row.names = FALSE)

