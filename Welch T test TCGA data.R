library(sqldf)
library(vegan)
library(xtable)
library(parallel)
library(DescTools)
library(jointseg)
library(data.table)
genes_of_interest = c("1869","4342","3845","2064","4602","595")

name_the_function = function(genes_of_interest, sub_summary_request = NULL)
{
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
  
  time1 = proc.time()[3]
  #Distribution of (Segment mean CN - 2) for oncogenes
  #################change location for pdf print##########################
  pdf(paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Plots/Plots\ on\ results/Distribution_of_segment_mean_CN_for_chosen_oncogenes",Sys.Date(),".pdf"))
  
  for( gene in genes_of_interest)
  {
    print(plot(density(as.numeric(TCGA_snp_of_interest[gene,]), na.rm = TRUE), main = paste("Distribution of (Segment mean CN - 2) for ", gene_info$SYMBOL[which(gene_info$Entrezid==gene)], sep = "")))
  }
  
  dev.off()
  print(paste((proc.time()[3]-time1)/60,"mins to plot distribution of (Segment mean CN - 2) for oncogenes"))
  
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
    write.table(sample_annotation_for_plot_all_genes[[gene]], file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/sample_annotation_for_plot_",gene_info$SYMBOL[which(gene_info$Entrezid==gene)],".txt",sep = ""), sep = "\t", row.names = FALSE)

  }

  TCGA_rnaseq_for_plot_all_genes = list()
  TCGA_rnaseq_adjusted_for_plot_all_genes = list()

  for( gene in genes_of_interest)
  {
    TCGA_rnaseq_for_plot_all_genes[[gene]] = TCGA_rnaseq[,c("ENTREZID", as.character(sample_annotation_for_plot_all_genes[[gene]]$Name))]
    TCGA_rnaseq_adjusted_for_plot_all_genes[[gene]] = TCGA_rnaseq_adjusted_v1[,c("ENTREZID", as.character(sample_annotation_for_plot_all_genes[[gene]]$Name))]
    #################change location for writing##########################
    write.table(TCGA_rnaseq_for_plot_all_genes[[gene]], file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_original_",gene_info$SYMBOL[which(gene_info$Entrezid==gene)],".txt",sep = ""), sep = "\t", row.names = FALSE)
    write.table(TCGA_rnaseq_adjusted_for_plot_all_genes[[gene]], file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_corrected_",gene_info$SYMBOL[which(gene_info$Entrezid==gene)],".txt",sep = ""), sep = "\t", row.names = FALSE)

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
  for( gene in genes_of_interest)
  {
    t_TCGA_rnaseq$gene_amplified_indicator = NA
    t_TCGA_rnaseq$gene_amplified_indicator[which(substr(rownames(t_TCGA_rnaseq),1,15)%in%samples_distribution_amplified[[gene]])] = 1
    t_TCGA_rnaseq$gene_amplified_indicator[which(substr(rownames(t_TCGA_rnaseq),1,15)%in%samples_distribution_normal[[gene]])] = 0
    colnames(t_TCGA_rnaseq)[dim(t_TCGA_rnaseq)[2]] = gene_info$SYMBOL[which(gene_info$Entrezid==gene)]
  }
  
  print(paste((proc.time()[3]-time1)/60,"mins to build dataset to push to welch t test"))
  
  ####welch t test
  time3 = proc.time()[3]
  welch_t_test_summary = list()
  for(gene in genes_of_interest)
  {
    time2 = proc.time()[3]
    welch_t_test_summary[[gene]] = as.data.frame(matrix(NA, (dim(t_TCGA_rnaseq)[2]-length(genes_of_interest)), 1+length(unique_tissue_type)*2+2))
    
    colnames(welch_t_test_summary[[gene]]) = c("ENTREZID"
                                              , paste(unique_tissue_type,"statistic", sep = " /// ")
                                              , paste(unique_tissue_type,"pvalue", sep = " /// ")
                                              , paste("All cancer types","statistic", sep = " /// ")
                                              , paste("All cancer types","pvalue", sep = " /// ")
    )
    welch_t_test_summary[[gene]]$ENTREZID = colnames(t_TCGA_rnaseq)[1:(dim(t_TCGA_rnaseq)[2]-length(genes_of_interest))]
    t_TCGA_rnaseq_v1 = t_TCGA_rnaseq[which(!is.na(t_TCGA_rnaseq[,gene_info$SYMBOL[which(gene_info$Entrezid==gene)]])),]
    
    for(tissue_type in unique_tissue_type)
    {
      t_TCGA_rnaseq_v2 = t_TCGA_rnaseq_v1[which(rownames(t_TCGA_rnaseq_v1)%in%sample_annotation$Name[which(sample_annotation$Type_updated==tissue_type)]),]
      
      gene_amp = which(t_TCGA_rnaseq_v2[,gene_info$SYMBOL[which(gene_info$Entrezid==gene)]]==1)
      gene_normal = which(t_TCGA_rnaseq_v2[,gene_info$SYMBOL[which(gene_info$Entrezid==gene)]]==0)
      
      if((length(gene_amp)>1)&&(length(gene_normal)>1))
      {
        time1 = proc.time()[3]
        no_cores <- 2
        cl <- makeCluster(no_cores, type = "FORK")
        
        x = parLapply(cl, c(1:(dim(t_TCGA_rnaseq_v2)[2]-length(genes_of_interest))),function(i){
          number_of_gene_amp_unique_non_na_samples = sum(!is.na(unique(t_TCGA_rnaseq_v2[gene_amp,i])))
          number_of_gene_normal_unique_non_na_samples = sum(!is.na(unique(t_TCGA_rnaseq_v2[gene_normal,i])))
          
          if((number_of_gene_amp_unique_non_na_samples>1)&&(number_of_gene_normal_unique_non_na_samples>1))
          {
            welch_t_test_result = t.test(t_TCGA_rnaseq_v2[gene_amp,i],t_TCGA_rnaseq_v2[gene_normal,i])
            
            return(c(welch_t_test_result$statistic,
                     welch_t_test_result$p.value
            ))
            
          }else{
            return(c(NA,
                     NA
            ))
          }
        })
        print((proc.time()[3]-time1)/60)
        
        stopCluster(cl)
        
        
        welch_t_test_summary[[gene]][,c(paste(tissue_type,"statistic", sep = " /// ")
                                       , paste(tissue_type,"pvalue", sep = " /// "))] <- as.data.frame(do.call(rbind,lapply(x,matrix,ncol=2,byrow=FALSE)))
        print((proc.time()[3]-time1)/60)
        print(tissue_type)
        print(gene)
        print("Analysis done")
        
      }else{
        print(tissue_type)
        print(gene)
        print("***************Analysis not done*******************")
        
      }
      
    }
    
    t_TCGA_rnaseq_v2 = t_TCGA_rnaseq_v1[which(!is.na(t_TCGA_rnaseq_v1[,gene_info$SYMBOL[which(gene_info$Entrezid==gene)]])),]
    
    gene_amp = which(t_TCGA_rnaseq_v2[,gene_info$SYMBOL[which(gene_info$Entrezid==gene)]]==1)
    gene_normal = which(t_TCGA_rnaseq_v2[,gene_info$SYMBOL[which(gene_info$Entrezid==gene)]]==0)
    
    
    if((length(gene_amp)>1)&&(length(gene_normal)>1))
    {
      time1 = proc.time()[3]
      no_cores <- 2
      cl <- makeCluster(no_cores, type = "FORK")
      
      x = parLapply(cl, c(1:(dim(t_TCGA_rnaseq_v2)[2]-length(genes_of_interest))),function(i){
        number_of_gene_amp_unique_non_na_samples = sum(!is.na(unique(t_TCGA_rnaseq_v2[gene_amp,i])))
        number_of_gene_normal_unique_non_na_samples = sum(!is.na(unique(t_TCGA_rnaseq_v2[gene_normal,i])))
        
        if((number_of_gene_amp_unique_non_na_samples>1)&&(number_of_gene_normal_unique_non_na_samples>1))
        {
          welch_t_test_result = t.test(t_TCGA_rnaseq_v2[gene_amp,i],t_TCGA_rnaseq_v2[gene_normal,i])
          
          return(c(welch_t_test_result$statistic,
                   welch_t_test_result$p.value
          ))
          
        }else{
          return(c(NA,
                   NA
          ))
        }
      })
      print((proc.time()[3]-time1)/60)
      print(gene)
      
      stopCluster(cl)
      
      
      welch_t_test_summary[[gene]][,c(paste("All cancer types","statistic", sep = " /// ")
                                     , paste("All cancer types","pvalue", sep = " /// "))] <- as.data.frame(do.call(rbind,lapply(x,matrix,ncol=2,byrow=FALSE)))
      
    }
    print((proc.time()[3]-time2)/60)
    print(gene)
    print("welch t test done")
    
  }
  
  print((proc.time()[3]-time3)/60)
  print("Whole welch t test analysis for all the genes is done")
  
  time1 = proc.time()[3]
  ######################change the location for saving the RData
  save(welch_t_test_summary, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/welch_t_test_summary_TCGA_for_all_tissue_types_and_combined_",paste(gene_info$SYMBOL[which(gene_info$Entrezid%in%genes_of_interest)], collapse = "_"),".RData", sep = ""))
  
  
  ######################reformat for GSEA
  
  gene_summary = read.table("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/TCGA_rnaseq_gene_summary_per_cancer_types.txt", sep = "\t", header = TRUE)
  
  
  minus_log_pvalue_with_sign_welch_t_test_summary = as.data.frame(matrix(NA, dim(welch_t_test_summary[[genes_of_interest[6]]])[1], 1+length(genes_of_interest)*length(unique_tissue_type)+length(genes_of_interest)))
  unique_tissue_type = sort(unique_tissue_type)
  colnames(minus_log_pvalue_with_sign_welch_t_test_summary) = c("ENTREZID"
                                                                , apply(expand.grid(unique_tissue_type,gene_info$SYMBOL[which(gene_info$Entrezid%in%genes_of_interest)]), 1, paste, collapse=" /// ")
                                                                , apply(expand.grid("All cancer types",gene_info$SYMBOL[which(gene_info$Entrezid%in%genes_of_interest)]), 1, paste, collapse=" /// "))
  
  minus_log_pvalue_with_sign_welch_t_test_summary$ENTREZID = welch_t_test_summary[[genes_of_interest[6]]]$ENTREZID
  for(tissue_type in c(as.character(unique_tissue_type), "All cancer types"))
  {
    for(gene in genes_of_interest)
    {
      minus_log_pvalue_with_sign_welch_t_test_summary[,paste(tissue_type,gene_info$SYMBOL[which(gene_info$Entrezid==gene)],sep = " /// ")] = -log10(ifelse(welch_t_test_summary[[gene]][,paste(tissue_type,"pvalue", sep = " /// ")]==0, 10^(-317),welch_t_test_summary[[gene]][,paste(tissue_type,"pvalue", sep = " /// ")]))*sign(welch_t_test_summary[[gene]][,paste(tissue_type,"statistic", sep = " /// ")])
    }
    
  }
  gene_summary = gene_summary[,c(dim(gene_summary)[2]:(dim(gene_summary)[2]-2),1:(dim(gene_summary)[2]-3))]
  
  minus_log_pvalue_with_sign_welch_t_test_summary_v1 = sqldf("select b.quality_indicator, b.hodgeslehmann_estimate_all_genes as read_count_hl_estimate,a.*  from minus_log_pvalue_with_sign_welch_t_test_summary a
                                                             left join gene_summary b
                                                             on a.ENTREZID = b.ENTREZID")
  minus_log_pvalue_with_sign_welch_t_test_summary_v1 = minus_log_pvalue_with_sign_welch_t_test_summary_v1[which(minus_log_pvalue_with_sign_welch_t_test_summary_v1$quality_indicator==1),]
  minus_log_pvalue_with_sign_welch_t_test_summary_v1 = minus_log_pvalue_with_sign_welch_t_test_summary_v1[,c(3,1:2,4:dim(minus_log_pvalue_with_sign_welch_t_test_summary_v1)[2])]
  
  minus_log_pvalue_with_sign_welch_t_test_summary_v1 = sqldf("select b.SYMBOL, b.GENENAME, b.chromosome_no as CHR, b.BP_mapping,a.*  from minus_log_pvalue_with_sign_welch_t_test_summary_v1 a
                                                             left join gene_info b
                                                             on a.ENTREZID = b.Entrezid")
  minus_log_pvalue_with_sign_welch_t_test_summary_v1 = minus_log_pvalue_with_sign_welch_t_test_summary_v1[,c(5,1:4,6:dim(minus_log_pvalue_with_sign_welch_t_test_summary_v1)[2])]
  
  common_trunk_positive_samples = c(1:dim(minus_log_pvalue_with_sign_welch_t_test_summary_v1)[1])
  common_trunk_negative_samples = c(1:dim(minus_log_pvalue_with_sign_welch_t_test_summary_v1)[1])
  for(gene in genes_of_interest)
  {
    common_trunk_positive_samples = intersect(common_trunk_positive_samples, which(minus_log_pvalue_with_sign_welch_t_test_summary_v1[,paste("All cancer types", gene_info$SYMBOL[which(gene_info$Entrezid==gene)], sep = " /// ")] > 2))
    common_trunk_negative_samples = intersect(common_trunk_negative_samples, which(minus_log_pvalue_with_sign_welch_t_test_summary_v1[,paste("All cancer types", gene_info$SYMBOL[which(gene_info$Entrezid==gene)], sep = " /// ")] <= -2))
    
  }
  
  minus_log_pvalue_with_sign_welch_t_test_summary_v1$common_trunk_positive = 0
  minus_log_pvalue_with_sign_welch_t_test_summary_v1$common_trunk_negative = 0
  
  minus_log_pvalue_with_sign_welch_t_test_summary_v1$common_trunk_positive[common_trunk_positive_samples] = 1
  minus_log_pvalue_with_sign_welch_t_test_summary_v1$common_trunk_negative[common_trunk_negative_samples] = 1
  
  
  ##########################change the location of writing the file##############################
  write.table(minus_log_pvalue_with_sign_welch_t_test_summary_v1, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Results_summary_of_welch_t_test_analysis_on_TCGA_on_all_combinations_",paste(gene_info$SYMBOL[which(gene_info$Entrezid%in%genes_of_interest)], collapse = "_"),".txt",sep = ""), sep = "\t", row.names = FALSE)
  
  if(!is.null(sub_summary_request))
  {
    minus_log_pvalue_with_sign_welch_t_test_summary_v1 = minus_log_pvalue_with_sign_welch_t_test_summary_v1[which(minus_log_pvalue_with_sign_welch_t_test_summary_v1$ENTREZID%in%sub_summary_request),c(1:7, (dim(minus_log_pvalue_with_sign_welch_t_test_summary_v1)[2]-length(genes_of_interest) - 1):dim(minus_log_pvalue_with_sign_welch_t_test_summary_v1)[2])]
    write.table(minus_log_pvalue_with_sign_welch_t_test_summary_v1, file = paste("/Users/arkajyotibhattacharya/Projects/Replication\ stress\ project/Results/Results_summary_of_welch_t_test_analysis_on_TCGA_on_all_combinations_",paste(gene_info$SYMBOL[which(gene_info$Entrezid%in%genes_of_interest)], collapse = "_"),"sub_summary.txt",sep = ""), sep = "\t", row.names = FALSE)
    
  }
  
  print(paste((proc.time()[3]-time1)/60,"mins to obtain summary for GSEA"))
  print(paste((proc.time()[3]-time_full)/60,"mins to run the whole analysis"))
  
  return(1)
}


name_the_function(genes_of_interest = c("1869","4342","3845","2064","4602","595"), sub_summary_request = c("65265", "55661", "79929", "27304", "51678", "55226", "79691", "197407", "22847"))
