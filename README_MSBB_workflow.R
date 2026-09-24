## Reprossing RNASeq data for AMP-AD1.0
# J Lundin
# July 14 2026

## WORKFLOW


DE_MSBB01_synapse_metadata.R
  Loads and cleans metadata
  output: final md_all file (syn76887188) (n=xx)

DE_MSBB02_synapse_counts_tpms_from_synapse.R
  Loads and filters counts and tpms 
    counts and tpms filtered for both counts and tpms thresholds and were output as separate files 
  output: msbb_all_counts_filtered.txt (syn76887322)
  output: msbb_all_tpm_filtered.txt (syn76979799)
  
  Check for additional count and tpm filter stratified by tissue and sex
  output: list of gene_ids to remove based on tissue and sex specific filtering (syn76980381)

Technical variables
  technical_stats_multiqc_star.R  Output: MSBB_multiqc_star_technical_stats.csv
  technical_stats_fastqc.R [run on AWS] (MSBB_fq_stats.rds contains: basic_stats.txt, phred_per_base.txt, base_content.txt) 

MSBB_QC.Rmd
  Run checks on RIN thresholds, inferred sex, fast-qc, multiqc, calcuate technical covariate for RNA metrics (plus technical outliers check) -- filter final md and count files
  output: MSBB_md_counts.rds ("metadata", "counts") (syn76887997)

  Ran CQN normalization on QC file
  output: MSBB_md_counts_cqn.rds ("metadata", "counts", "dge_cqn") (syn76980951)  

  Additional filtering based on RNA alignment and QC metrics PCs
  output: MSBB_md_counts_cqn_FINAL.rds ("metadata", "counts", "dge_cqn") (synxx)  

#here for reference but not longer used
DE_QC_MSBB_26Aug2026.Rmd  
   Outlier detection using PCA (crude)
   Calculate SVs - no longer using bc of PCA of RNA metrics
   Variance partitioning visualization - figure saved
   Correlation of model PCs with covariates


   

