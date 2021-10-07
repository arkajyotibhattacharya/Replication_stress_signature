library(sqldf)
library(vegan)
library(xtable)
library(parallel)
library(DescTools)
library(jointseg)
library(data.table)
genes_of_interest = c("1869","4342","3845","2064","4602","595")
genes_of_interest = "4609"

  time_full = proc.time()[3]
  time1 = proc.time()[3]
  #load TCGA SNP data
  
  tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA_non_adjusted_cn_data_ordered_combined.txt", sep = "\t" , header = TRUE, nrows = 50, row.names = "gene")
  classes <- sapply(tab5rows, class)
  TCGA_snp <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA_non_adjusted_cn_data_ordered_combined.txt", sep = "\t" , header = TRUE, colClasses = classes, row.names = "gene")
  colnames_tcga = gsub("\\.", "-", colnames(TCGA_snp))
  colnames(TCGA_snp) = colnames_tcga
  
  print(paste((proc.time()[3]-time1)/60,"mins taken to load TCGA SNP data "))
  
  time1 = proc.time()[3]
  #choose only genes_of_interest
  
  TCGA_snp_of_interest = TCGA_snp[genes_of_interest,]
  colnames_tcga = gsub("\\.", "-", colnames(TCGA_snp_of_interest))
  colnames(TCGA_snp_of_interest) = colnames_tcga
  
  rm(TCGA_snp)
  print(paste((proc.time()[3]-time1)/60,"mins taken to choose only genes_of_interest"))
  
  time1 = proc.time()[3]
  
  #load TCGA RNA seq data
  
  tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed.txt", sep = "\t" , header = TRUE, nrows = 5)
  classes <- sapply(tab5rows, class)
  TCGA_rnaseq <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed.txt", sep = "\t" , header = TRUE, colClasses = classes)
  colnames_tcga = gsub("\\.", "-", colnames(TCGA_rnaseq))
  colnames(TCGA_rnaseq) = colnames_tcga
  
  TCGA_rnaseq = TCGA_rnaseq[,which(!(colnames(TCGA_rnaseq)%in%colnames(TCGA_rnaseq)[which(substr(colnames(TCGA_rnaseq),1,15)%in%names(which(table(substr(colnames(TCGA_rnaseq),1,15))>1)))][20:38]))]
  
  which(table(substr(colnames(TCGA_rnaseq),1,15))>1)
  print(paste((proc.time()[3]-time1)/60,"mins taken to load TCGA RNA seq data"))
  
  time1 = proc.time()[3]
  #work on those samples which have corresponding RNASEQ data
  
  TCGA_snp_of_interest = TCGA_snp_of_interest[,which(colnames(TCGA_snp_of_interest)%in%substr(colnames(TCGA_rnaseq),1,15))]
  TCGA_rnaseq = TCGA_rnaseq[,which(substr(colnames(TCGA_rnaseq),1,15)%in%colnames(TCGA_snp_of_interest))]
  print(paste((proc.time()[3]-time1)/60,"mins to work on those samples which have corresponding RNASEQ data"))
  
  time1 = proc.time()[3]
  #load sample annotation
  sample_annotation = read.csv("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA_Sample_To_TumorType_20190920.csv")
  sample_annotation = sample_annotation[which(sample_annotation$Name%in%colnames(TCGA_rnaseq)),]
  unique_tissue_type = unique(sample_annotation$Type_updated)
  print(paste((proc.time()[3]-time1)/60,"mins to load sample annotation"))
  
  time1 = proc.time()[3]
  #load gene info
  gene_info = fread("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Data/Entrezid_mapping_using_org_Hs_eg_db_09052018.txt", sep = "\t", header = TRUE)
  print(paste((proc.time()[3]-time1)/60,"mins to load gene info"))
  
  # time1 = proc.time()[3]
  #Distribution of (Segment mean CN - 2) for oncogenes
  #################change location for pdf print##########################
  # pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Distribution_of_segment_mean_CN_for_chosen_oncogenes",Sys.Date(),".pdf"))
  # 
  # for( gene in genes_of_interest)
  # {
  #   print(plot(density(as.numeric(TCGA_snp_of_interest[gene,]), na.rm = TRUE), main = paste("Distribution of (Segment mean CN - 2) for ", gene_info$SYMBOL[which(gene_info$Entrezid==gene)], sep = "")))
  # }
  # 
  # dev.off()
  # print(paste((proc.time()[3]-time1)/60,"mins to plot distribution of (Segment mean CN - 2) for oncogenes"))
  
  time1 = proc.time()[3]
  #choose samples with normal and amplification separately for three oncogenes
  #http://mcr.aacrjournals.org/content/12/4/485.long
  
  samples_distribution_amplified = list()
  samples_distribution_normal = list()
  
  for( gene in genes_of_interest)
  {
    samples_distribution_amplified[[gene]] = colnames(TCGA_snp_of_interest)[which(TCGA_snp_of_interest[gene,]> 0.2)]
    samples_distribution_normal[[gene]] = colnames(TCGA_snp_of_interest)[intersect(which(TCGA_snp_of_interest[gene,]<= 0.2)
                                                                                   ,which(TCGA_snp_of_interest[gene,]> -0.2))]
  }
  samples_to_work_with = unique(c(unlist(samples_distribution_amplified)
                                  ,unlist(samples_distribution_normal)))
  
  print(paste((proc.time()[3]-time1)/60,"mins to choose samples with normal and amplification separately for three oncogenes"))
  
  time1 = proc.time()[3]
  # For each tumor type and all tumor types together
  # Remove samples for corresponding genes where cn < -0.2
  # Correlate CN of the gene with rnaseq of all genes and keep p-values
  # load and format TCGA rnaseq corrected for cancer types dataset
  
  tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_normalized_withrespectto_cancer_types.txt", sep = "\t" , header = TRUE, nrows = 5)
  classes <- sapply(tab5rows, class)
  TCGA_rnaseq_adjusted <- read.table("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_normalized_withrespectto_cancer_types.txt", sep = "\t" , header = TRUE, colClasses = classes)
  colnames_tcga = gsub("\\.", "-", colnames(TCGA_rnaseq_adjusted))
  colnames(TCGA_rnaseq_adjusted) = colnames_tcga
  
  TCGA_rnaseq$ENTREZID = rownames(TCGA_rnaseq)
  TCGA_rnaseq = TCGA_rnaseq[,c(dim(TCGA_rnaseq)[2],1:(dim(TCGA_rnaseq)[2]-1))]
  TCGA_rnaseq_adjusted_v1 = TCGA_rnaseq_adjusted
  TCGA_rnaseq_adjusted_v1$ENTREZID = rownames(TCGA_rnaseq_adjusted_v1)
  TCGA_rnaseq_adjusted_v1 = TCGA_rnaseq_adjusted_v1[,c(dim(TCGA_rnaseq_adjusted_v1)[2],1:(dim(TCGA_rnaseq_adjusted_v1)[2]-1))]
  print(paste((proc.time()[3]-time1)/60,"mins to load and format TCGA rnaseq corrected for cancer types dataset"))
  
  time1 = proc.time()[3]
  #load cancer type with sample info
  load("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_hl_robustsd_rnaseq_expression_for_all_tissue_types.RData")
  print(paste((proc.time()[3]-time1)/60,"mins to load load cancer type with sample info"))
  
  time1 = proc.time()[3]
  #Build datasets original and corrected for plot using Analyzer tool
  #sample annotation for each gene
  sample_annotation_for_plot = sample_annotation[,c("Name","Type_updated")]
  sample_annotation_for_plot$Type_updated = as.character(sample_annotation_for_plot$Type_updated)
  sample_annotation_for_plot_all_genes = list()
  for( gene in genes_of_interest)
  {
    sample_annotation_for_plot_all_genes[[gene]] = sample_annotation_for_plot[which(substr(sample_annotation_for_plot$Name,1,15)%in%c(samples_distribution_normal[[gene]],samples_distribution_amplified[[gene]])),]
    temp = sample_annotation_for_plot_all_genes[[gene]]
    temp$Type_updated = "1 All cancer types"
    sample_annotation_for_plot_all_genes[[gene]] = rbind(sample_annotation_for_plot_all_genes[[gene]], temp)
    sample_annotation_for_plot_all_genes[[gene]]$Type_updated[which(substr(sample_annotation_for_plot_all_genes[[gene]]$Name,1,15)%in%samples_distribution_amplified[[gene]])] = paste(gene_info$SYMBOL[which(gene_info$Entrezid==gene)],"///", as.character(sample_annotation_for_plot_all_genes[[gene]]$Type_updated[which(substr(sample_annotation_for_plot_all_genes[[gene]]$Name,1,15)%in%samples_distribution_amplified[[gene]])]), "/// amplified")
    sample_annotation_for_plot_all_genes[[gene]]$Type_updated[which(substr(sample_annotation_for_plot_all_genes[[gene]]$Name,1,15)%in%samples_distribution_normal[[gene]])] = paste(gene_info$SYMBOL[which(gene_info$Entrezid==gene)],"///", as.character(sample_annotation_for_plot_all_genes[[gene]]$Type_updated[which(substr(sample_annotation_for_plot_all_genes[[gene]]$Name,1,15)%in%samples_distribution_normal[[gene]])]), "/// neutral")

    temp_for_sort = sample_annotation_for_plot_all_genes[[gene]]
    sample_annotation_for_plot_all_genes[[gene]] = sqldf("select * from temp_for_sort order by 2 desc")
    sample_annotation_for_plot_all_genes[[gene]]$Type_updated = gsub(" /// 1", " ///",sample_annotation_for_plot_all_genes[[gene]]$Type_updated)
    #################change location for writing##########################
    # write.table(sample_annotation_for_plot_all_genes[[gene]], file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/sample_annotation_for_plot_",gene_info$SYMBOL[which(gene_info$Entrezid==gene)],".txt",sep = ""), sep = "\t", row.names = FALSE)

  }

  TCGA_rnaseq_for_plot_all_genes = list()
  TCGA_rnaseq_adjusted_for_plot_all_genes = list()

  for( gene in genes_of_interest)
  {
    TCGA_rnaseq_for_plot_all_genes[[gene]] = TCGA_rnaseq[,c("ENTREZID", as.character(sample_annotation_for_plot_all_genes[[gene]]$Name))]
    TCGA_rnaseq_adjusted_for_plot_all_genes[[gene]] = TCGA_rnaseq_adjusted_v1[,c("ENTREZID", as.character(sample_annotation_for_plot_all_genes[[gene]]$Name))]
    #################change location for writing##########################
    # write.table(TCGA_rnaseq_for_plot_all_genes[[gene]], file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_original_",gene_info$SYMBOL[which(gene_info$Entrezid==gene)],".txt",sep = ""), sep = "\t", row.names = FALSE)
    # write.table(TCGA_rnaseq_adjusted_for_plot_all_genes[[gene]], file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_corrected_",gene_info$SYMBOL[which(gene_info$Entrezid==gene)],".txt",sep = ""), sep = "\t", row.names = FALSE)

  }

  rm(TCGA_rnaseq_for_plot_all_genes)
  rm(TCGA_rnaseq_adjusted_for_plot_all_genes)
  rm(sample_annotation_for_plot_all_genes)
  rm(sample_annotation_for_plot)
  print(paste((proc.time()[3]-time1)/60,"mins to build datasets original and corrected for plot using Analyzer tool"))
  
  
  time1 = proc.time()[3]
  #build dataset to push to welch t test
  t_TCGA_rnaseq = as.data.frame(t(TCGA_rnaseq_adjusted))
  t_TCGA_rnaseq = t_TCGA_rnaseq[which(substr(rownames(t_TCGA_rnaseq),1,15)%in%samples_to_work_with),]
  
  #########assigning 1 and zero to samples for each gene's amplified or deleted state
  time3 = proc.time()[3]
  permanova_summary = list()
  # for( gene in genes_of_interest)
  # {
  gene = genes_of_interest
    t_TCGA_rnaseq$gene_amplified_indicator = NA
    t_TCGA_rnaseq$gene_amplified_indicator[which(substr(rownames(t_TCGA_rnaseq),1,15)%in%samples_distribution_amplified[[gene]])] = 1
    t_TCGA_rnaseq$gene_amplified_indicator[which(substr(rownames(t_TCGA_rnaseq),1,15)%in%samples_distribution_normal[[gene]])] = 0
    colnames(t_TCGA_rnaseq)[dim(t_TCGA_rnaseq)[2]] = gene_info$SYMBOL[which(gene_info$Entrezid==gene)]
  
  
  print(paste((proc.time()[3]-time1)/60,"mins to build dataset to push to permanova"))
  t_TCGA_rnaseq = t_TCGA_rnaseq[which(!is.na(t_TCGA_rnaseq$MYC)),]
  list_of_files = ls()
  list_of_files = list_of_files[-which(list_of_files == "t_TCGA_rnaseq")]
  rm(list = list_of_files)
  # TCGA_rnaseq_adjusted_for_permanova = t(t_TCGA_rnaseq)
  # 
  # design_vec= TCGA_rnaseq_adjusted_for_permanova[dim(TCGA_rnaseq_adjusted_for_permanova)[1],]
  # TCGA_rnaseq_adjusted_for_permanova = TCGA_rnaseq_adjusted_for_permanova[-dim(TCGA_rnaseq_adjusted_for_permanova)[1],]
  # TCGA_rnaseq_adjusted_for_permanova = TCGA_rnaseq_adjusted_for_permanova[,which(!is.na(design_vec))]
  # design_vec = design_vec[which(!is.na(design_vec))]
    time2 = proc.time()[3]
    
    
    
    permanova_per_gene = function(i,data)
    {
      time2 = proc.time()[3]
      
      permanova_result_1_euclidean = adonis(t_TCGA_rnaseq[,i] ~ MYC
                                            , data = t_TCGA_rnaseq
                                            , method = "euclidean"
                                            , permutations = 100
                                            , na.rm = TRUE
      )
      print((proc.time()[3]-time2)/60)
      
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
    
    
    
    fit <- lmFit(TCGA_rnaseq_adjusted_for_limma, design=design_vec)
    fit <- eBayes(fit)
    limma_summary[[gene]] = topTable(fit,sort="none",n=Inf)
    
    print((proc.time()[3]-time2)/60)
    print(gene)
    print("limma done")
    
  }
  
  print((proc.time()[3]-time3)/60)
  print("Whole welch t test analysis for all the genes is done")
  
  time1 = proc.time()[3]
  ######################change the location for saving the RData
  save(limma_summary, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/limma_summary_TCGA_for_all_tissue_types_and_combined_",paste(gene_info$SYMBOL[which(gene_info$Entrezid%in%genes_of_interest)], collapse = "_"),".RData", sep = ""))
  
  print(paste((proc.time()[3]-time1)/60,"mins to obtain summary for GSEA"))
  print(paste((proc.time()[3]-time_full)/60,"mins to run the whole analysis"))
  
  return(1)
}


name_the_function(genes_of_interest = c("4609","993", "898","1869","4342","3845","2064","4602","595"), sub_summary_request = c("65265", "55661", "79929", "27304", "51678", "55226", "79691", "197407", "22847"))
# genes_of_interest = c("4609","993", "898")

