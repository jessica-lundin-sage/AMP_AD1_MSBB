## Differential expression analyses on RNAseq data from MAYO
# STEP 02: count and tpm data from synapse
# J Lundin
# June 16 2026


pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser,dplyr,purrr,readr,lubridate,stringr,tibble,ggplot2)
pacman::p_load(httr, purrr)
synLogin()

#work_dir <- ("C:/Users/jlundin/OneDrive - Sage Bionetworks/RNASeq_Harm/MAYO")
#setwd(work_dir)

# get functions from git
urls <- c("https://raw.githubusercontent.com/jessica-lundin-sage/pers_Lundin/main/RNASeq_DE/AMP_AD1/functions/functions_filter_low_count_genes.R")
purrr::walk(urls, function(u) {
  resp <- httr::GET(u, httr::add_headers(Authorization = paste("token", Sys.getenv("GITHUB_PAT"))))
  httr::stop_for_status(resp)
  tmp <- tempfile(fileext = ".R")
  writeLines(httr::content(resp, "text", encoding = "UTF-8"), tmp)
  source(tmp)
})


### pulling count and tpm data from synapse ----
msbb_counts_temp <- synapser::synGet("syn69368998") 
  msbb_counts <- read.csv(msbb_counts_temp$path, sep="\t", header=T, check.names = FALSE)
msbb_tpm_temp <- synapser::synGet("syn69368999") 
  msbb_tpm <- read.csv(msbb_tpm_temp$path, sep="\t", header=T, check.names = FALSE)

# metadata
MSBB_md_temp <- synapser::synGet("syn76887188") #MSBB_md_all
  MSBB_md <- read.csv(MSBB_md_temp$path, header=T)

  
## filtering on low counts ----
filtered_genes_MSBB <- filter_gene_expression(
  tpm_file   = msbb_tpm ,
  reads_file  = msbb_counts,
  metadata  = MSBB_md,
  synid_outfile = c("syn76814282"),
  study_var = c("MSBB_all")
  )

  

  ## process tpms and counts filtering again but stratified by tissue and sex ----
  library(tidyverse)
  
  sex_var <- c("female","male")
  tissue2_var <- c("FP", "IFG", "PG" ,"STG")
  mat1<-NULL
  mat2<-NULL
  for (i in 1:length(tissue2_var)){
    for (j in 1:length(sex_var)){
      msbb_counts_filtered <- read.csv(synapser::synGet('syn76887322')$path, sep="\t", stringsAsFactors = F, check.names = FALSE)
      msbb_tpm_filtered <- read.csv(synapser::synGet('syn76979799')$path, sep="\t", stringsAsFactors = F, check.names = FALSE)
      md_all <- read.csv(synapser::synGet('syn76887188')$path, stringsAsFactors = F)
      
      
      read_thresh = 6
      prop_samples = 0.2
      tpm_thresh = 0.1
      
      # 1. Load Data
      gene_reads <- msbb_counts_filtered
      gene_tpm   <- msbb_tpm_filtered 
      pheno      <- md_all %>% filter(sex == sex_var[j] & tissue2 == tissue2_var[i])
      
      # 2. Match IDs
      # Get list of sample IDs present in phenotype file
      RNA_link <- pheno$specimenID
      
      # Filter columns: Keep 'Name', 'Description', and any columns matching the phenotype IDs
      keep_cols_tpm   <- intersect(names(gene_tpm), RNA_link)
      keep_cols_tpm   <- c("gene_id", keep_cols_tpm)
      
      keep_cols_reads <- intersect(names(gene_reads), RNA_link)
      keep_cols_reads <- c("gene_id", keep_cols_reads)
      
      gene_tpm   <- gene_tpm[, keep_cols_tpm]
      gene_reads <- gene_reads[, keep_cols_reads, drop = FALSE]
      
      
      # 4. Filter by Reads Threshold
      rownames(gene_reads) <- gene_reads$gene_id
      read_matrix <- as.matrix(gene_reads[, 2:ncol(gene_reads)])
      passed_reads <- rowMeans(read_matrix >= read_thresh) >= prop_samples
      cat(paste("Table for sex and tissue specific test of counts >=6 in >=20% of samples for ", sex_var[j]," ",tissue2_var[i]))
      cat(" ")
      table(passed_reads)
      cat(table(passed_reads))
      remove_for_strat <- names(passed_reads)[!passed_reads]
      
      # 3. Filter by TPM Threshold
      # Logic: rowMeans on a logical matrix is much faster than sapply/transpose
      rownames(gene_tpm) <- gene_tpm$gene_id
      tpm_matrix <- as.matrix(gene_tpm[, 2:ncol(gene_tpm)])
      passed_tpm <- rowMeans(tpm_matrix >= tpm_thresh) >= prop_samples
      cat(paste("Table for sex and tissue specific test of tpm >=0.1 in >=20% of samples for ", sex_var[j]," ",tissue2_var[i]))
      cat(" ")
      table(passed_tpm)
      cat(table(passed_tpm))
      remove_for_strat_tpm <- names(passed_tpm)[!passed_tpm]
      
      final_genes_remove <- unique(names(passed_tpm)[!passed_tpm], 
                                   names(passed_reads)[!passed_reads])
      
      if (length(final_genes_remove) != 0){
        mat2 <- data.frame(gene_id = final_genes_remove, stringsAsFactors = FALSE)
        mat2$sex <- sex_var[j]
        mat2$tissue2 <- tissue2_var[i]
        
        mat1 <- rbind(mat1, mat2)
      }
      
      # write to synapse - this does individual files. Below is one file with a column for sex and tissue - can be filtered downstream
    #  file_path <- paste("Remove gene_ids from strat analysis reads and tpm_", sex_var[j],"_",tissue2_var[i],".txt",sep="")
    #  write.table(final_genes_remove, file = file_path, row.names = FALSE)
    #  synapser::synStore(synapser::File(path = file_path, parent = "syn76814282"))
      
    }
  }
  
  table(mat1$sex, mat1$tissue2)
  
  mat1 <- as.data.frame(mat1)
  
  mat1$gene_id_clean <- sub("\\..*", "", mat1$gene_id)
  #BiocManager::install("EnsDb.Hsapiens.v86")  # run once, not in the pipeline script
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
  
  # 3. Define your list of target Gene IDs
  my_ensembl_ids <- unique(mat1$gene_id_clean)
  
  gene_positions <- genes(EnsDb.Hsapiens.v86,
                          filter = GeneIdFilter(my_ensembl_ids),
                          columns = c("gene_id", "seq_name", "gene_seq_start", "gene_seq_end", "seq_strand"))
  gene_positions <- as.data.frame(gene_positions)
  gene_positions <- gene_positions %>% dplyr::rename(chromosome_name = seqnames)
  
  mat1 <- merge(mat1, gene_positions, by.x="gene_id_clean", by.y="gene_id", all=T)
  
  
  # write to synapse
  file_path <- paste("Remove gene_ids from stratified analysis reads and tpm.csv",sep="")
  write.table(mat1, file = file_path, row.names = FALSE, sep="\t")
  synapser::synStore(synapser::File(path = file_path, parent = "syn76814282"))
  
  
  