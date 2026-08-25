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

MSBB_QC.Rmd
  Run checks on RIN thresholds, inferred sex, fast-qc, multiqc, calcuate technical covariate for RNA metrics (plus technical outliers check) -- filter final md and count files
  output: MSBB_md_counts.rds ("metadata", "counts") (syn76887997)

  Ran CQN normalization on QC file
  output: MSBB_md_counts_cqn_.rds ("metadata", "counts", "dge_cqn", "dge_cqn_df") (syn76980951)  

On AWS:
DE_QC_b1tob4_24July2026.Rmd  #here for reference but not longer used
   Outlier detection using PCA (crude)
   Calculate SVs - no longer using bc of PCA of RNA metrics
   Variance partitioning visualization - figure saved
   Outlier detection using PCA (residualized with incorporation of RNA metrics and SVs) - use AWS version instead (on residualized data)
   #output: ROSMAP_md_sv_final_for_DE_models.csv [final metadata file, includes final metadata with RNA metrics and SVs] ] ("md_sv") (synxxx)
   
On AWS:
DE_with_built_interaction_no_svs.R 
  output "ROSMAP_DE_final.rds" (syn76557299) contains:
      "metadata" = md_sv, 
      "vobj_expr" = voom_gene_expression, 
      "dge_cqn" = dge_cqn, 
      "ebayes" = dream.cont.ebayes6, 
      "fit_contrasts" = fit_contrasts6,
      "males_DLPFC6" = males_DLPFC6, 
      "females_DLPFC6" = females_DLPFC6,  
      "males_PCC6" = males_PCC6, 
      "females_PCC6" = females_PCC6, 
      "males_CN6" = males_CN6, 
      "females_CN6" = females_CN6
   Brings in ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds ("metadata", "counts", "dge_cqn", "dge_cqn_df") (syn76553509)  

      
On AWS:
DE_Residuals_for_sharing.R
   Models technical variables only on cqn normalized counts using dream weights. formula: ~ PMI + RIN + PC1_metrics + PC2_metrics + PC3_metrics + (1|individualID) + (1|final_batch)
   "ROSMAP_DE_res.rds" (syn76564990) contains: 
      "metadata" = md_sv, 
      "vobj_res" = voom_res (voomwithDreamWeights output), 
      "dge_cqn" = dge_cqn (counts with dge_can$E from CQN normalization), 
      "fit_res" = fit_res (dream() with counts, CQN offset (E), and dream weights), 
      "residual_gene_expression" = residual_gene_expression (residuals from dream())
   
   
   form_ck <- ~ 0 + group  + apoe4Status + age_cat + (1|individualID)
       How group was defined: md_sv$group <- interaction(md_sv$diagnosis, md_sv$sex, md_sv$tissue2, drop = TRUE)
                              group <- interaction(md_sv$diagnosis, md_sv$sex, md_sv$tissue2, drop = TRUE)
   "ROSMAP_DE_res2.rds" (syn76565771) contains: 
      "metadata_res3" = md_sv, 
      "fit_res2.dream" = fit_res2 (dream() with residualized counts, formula with additional vars (form_ck), and tissue X sex X diagnosis contrasts), 
      "fit_res3.ebayes" = fit_res3 (eBayes of fit_res2), 
      "males_DLPFC6_res3" = males_DLPFC6 (topTable of fit_res3 for each contrast), 
      "females_DLPFC6_res3" = females_DLPFC6,  
      "males_PCC6_res3" = males_PCC6, 
      "females_PCC6_res3" = females_PCC6, 
      "males_CN6_res3" = males_CN6, 
      "females_CN6_res3" = females_CN6
 
On AWS:
DE_residualized_plots_by_final_batch.R
    Run PCA on residualized data with overlay markers for batch and diagnosis - plots output
    "ROSMAP_residualized_data_by_batch.png"
    "ROSMAP_residualized_data_by_diagnosis.png"

Technical variables
  technical_stats_multiqc_star.R (syn76227881) ROSMAP_multiqc_star_technical_stats.csv
  technical_stats_fastqc.R (syn76283403) [run on AWS] (ROSMAP_fq_stats.rds contains: basic_stats.txt, phred_per_base.txt, base_content.txt) 

      
      
      
Helpful links:
https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
https://diseaseneurogenomics.github.io/variancePartition/articles/dream.html
https://github.com/DiseaseNeuroGenomics/variancePartition/blob/HEAD/R/dream.R
https://github.com/Sage-Bionetworks/amp-rnaseq
https://github.com/Sage-Bionetworks/ampad-rnaseq-reprocessing/tree/main/code/metadata_preprocessing
https://github.com/Sage-Bionetworks/sageseqr/blob/master/R/functions.R
https://github.com/Sage-Bionetworks/ampad-DiffExp/blob/master/gene_level_analysis/MAYO_geneLevel_TMM.Rmd
https://github.com/Sage-Bionetworks/sageseqr/blob/master/R/functions.R




