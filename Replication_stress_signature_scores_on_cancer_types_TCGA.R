
library(sqldf)
library(DescTools)
library(ggplot2)
library(ggrepel)

#load TCGA RNA seq data

tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed.txt", sep = "\t" , header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
TCGA_rnaseq <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed.txt", sep = "\t" , header = TRUE, colClasses = classes)
colnames_tcga = gsub("\\.", "-", colnames(TCGA_rnaseq))
colnames(TCGA_rnaseq) = colnames_tcga
#load color codes
color_code = read.csv("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/color_code_template.csv")
for(i in 1:15)
{
  sub_data = colMeans(color_code[sample(1:dim(color_code)[1],14),])
  color_code = rbind(color_code, sub_data)
}

#load sample annotation
sample_annotation = read.csv("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__Sample_To_TumorType_with_common_cancer_type_mapping_GEO_TCGA.csv")
sample_annotation = sample_annotation[which(sample_annotation$Name%in%colnames(TCGA_rnaseq)),]
unique_tissue_type = unique(sample_annotation$TYPE3)
TCGA_rnaseq = TCGA_rnaseq[,as.character(sample_annotation$Name)]

sum(colnames(TCGA_rnaseq)==as.character(sample_annotation$Name))

replication_stress_signature = c("65265",
                                 "55661",
                                 "27304",
                                 "51678",
                                 "55226",
                                 "197407")

TCGA_rnaseq_replication_stress_signature = TCGA_rnaseq[replication_stress_signature,]
TCGA_rnaseq_replication_stress_signature = t(scale(t(TCGA_rnaseq_replication_stress_signature)))

TCGA_rnaseq_replication_stress_signature_score = colMeans(TCGA_rnaseq_replication_stress_signature)
Sample_identifier = names(TCGA_rnaseq_replication_stress_signature_score)
score_data = as.data.frame(TCGA_rnaseq_replication_stress_signature_score)
score_data$Sample_identifier = Sample_identifier
rownames(sample_annotation) = sample_annotation$Name
sample_annotation = sample_annotation[rownames(score_data),]
score_data$Tumor_type = sample_annotation$TYPE3
score_data$normalized_TCGA_rnaseq_replication_stress_signature_score = (score_data$TCGA_rnaseq_replication_stress_signature_score - min(score_data$TCGA_rnaseq_replication_stress_signature_score))/(max(score_data$TCGA_rnaseq_replication_stress_signature_score) - min(score_data$TCGA_rnaseq_replication_stress_signature_score))

write.table(score_data, file = "/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Replication_stress_signature_score_per_sample_TCGA.txt", sep = "\t", row.names = FALSE)


unique_all_combinations = as.character(unique(score_data$Tumor_type))
color_sample<-rgb(red = color_code[c(1:length(unique_all_combinations)),1],green = color_code[c(1:length(unique_all_combinations)),2], blue = color_code[c(1:length(unique_all_combinations)),3])
color_code = as.data.frame(cbind(unique_all_combinations, color_sample))
sample_annotation_file = sqldf("select a.*, b.color_sample from score_data a
                               left join color_code b 
                               on a.Tumor_type = b.unique_all_combinations
                               order by Tumor_type")

pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Replication stress signatrue score distribution/TCGA replication stress signature score.pdf", sep = "")
    , height = 8, width = 16)
  ggplot(sample_annotation_file, aes(x=Tumor_type, y=TCGA_rnaseq_replication_stress_signature_score)) + 
  geom_dotplot(binaxis="y",stackdir="center", dotsize=0.1, color=sample_annotation_file$color_sample, binwidth = 0.0015) + 
  geom_boxplot(fill=color_sample, color=color_sample, alpha=0.4) + scale_y_continuous(name = "VALUE") + 
  scale_x_discrete(name = "Cancer types") +  geom_hline(yintercept=0) + theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  ggtitle(paste("TCGA replication stress signature score"))
dev.off()

pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Replication stress signatrue score distribution/TCGA replication stress signature score normalized.pdf", sep = "")
    , height = 8, width = 16)
ggplot(sample_annotation_file, aes(x=Tumor_type, y=normalized_TCGA_rnaseq_replication_stress_signature_score)) + 
  geom_dotplot(binaxis="y",stackdir="center", dotsize=0.1, color=sample_annotation_file$color_sample, binwidth = 0.0015) + 
  geom_boxplot(fill=color_sample, color=color_sample, alpha=0.4) + scale_y_continuous(name = "VALUE") + 
  scale_x_discrete(name = "Cancer types") +  geom_hline(yintercept=0) + theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  ggtitle(paste("TCGA replication stress signature score"))
dev.off()

pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Replication stress signatrue score distribution/TCGA replication stress signature score normalized horizonntal.pdf", sep = "")
    , height = 12, width = 8)
ggplot(sample_annotation_file, aes(x=Tumor_type, y=normalized_TCGA_rnaseq_replication_stress_signature_score)) + 
  geom_dotplot(binaxis="y",stackdir="center", dotsize=0.1, color=sample_annotation_file$color_sample, binwidth = 0.0015) + 
  geom_boxplot(fill=color_sample, color=color_sample, alpha=0.4) + scale_y_continuous(name = "Replication stress signature score") + 
  scale_x_discrete(name = "") + 
  # geom_hline(yintercept=0) + 
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  ggtitle(paste("TCGA replication stress signature score"))+ coord_flip()
dev.off()

score_summary = sqldf("select Tumor_type, avg(TCGA_rnaseq_replication_stress_signature_score) as TCGA_rnaseq_replication_stress_signature_mean_score,
median(TCGA_rnaseq_replication_stress_signature_score) as TCGA_rnaseq_replication_stress_signature_median_score,
max(TCGA_rnaseq_replication_stress_signature_score) as TCGA_rnaseq_replication_stress_signature_max_score,
min(TCGA_rnaseq_replication_stress_signature_score) as TCGA_rnaseq_replication_stress_signature_min_score,
avg(normalized_TCGA_rnaseq_replication_stress_signature_score) as normalized_TCGA_rnaseq_replication_stress_signature_mean_score,
median(normalized_TCGA_rnaseq_replication_stress_signature_score) as normalized_TCGA_rnaseq_replication_stress_signature_median_score,
max(normalized_TCGA_rnaseq_replication_stress_signature_score) as normalized_TCGA_rnaseq_replication_stress_signature_max_score,
                      min(normalized_TCGA_rnaseq_replication_stress_signature_score) as normalized_TCGA_rnaseq_replication_stress_signature_min_score
                      from score_data
                      group by 1")

#Should this be adjusted with respect to normals? like normal should have zero?

for(i in 1:dim(score_summary)[1])
{
  score_summary$TCGA_rnaseq_replication_stress_signature_HL_score[i] = HodgesLehmann(score_data$TCGA_rnaseq_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
  score_summary$TCGA_rnaseq_replication_stress_signature_sd_score[i] = sd(score_data$TCGA_rnaseq_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
  score_summary$normalized_TCGA_rnaseq_replication_stress_signature_HL_score[i] = HodgesLehmann(score_data$normalized_TCGA_rnaseq_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
  score_summary$normalized_TCGA_rnaseq_replication_stress_signature_sd_score[i] = sd(score_data$normalized_TCGA_rnaseq_replication_stress_signature_score[which(score_data$Tumor_type==score_summary$Tumor_type[i])])
}

write.table(score_summary, file = "/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Replication_stress_signature_score_summary_TCGA.txt", sep = "\t", row.names = FALSE)

