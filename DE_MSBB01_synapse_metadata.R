
## Differential expression analyses on RNAseq data from MSBB
# STEP 01: metadata from synapse
# J Lundin
# June 16 2026


pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser,dplyr,purrr,readr,lubridate,stringr,tibble,ggplot2)
synLogin()


#MSBB_meta_ind <- read.csv(synapser::synGet('syn6101474')$path, stringsAsFactors = F, check.names = FALSE) # incorrect file - running to compare with previous results

### pulling metadata from synapse ----
MSBB_meta_ind <- read.csv(synapser::synGet('syn73713767')$path, stringsAsFactors = F, check.names = FALSE) #from metadata harmonization study
MSBB_meta_ind <- merge(MSBB_meta_ind, MSBB_dx_prev, by="individualID") # adding in diagnosis_previous from "DE_MSBB_check_prev_md.R"

MSBB_meta_biosp <- read.csv(synapser::synGet('syn21893059')$path, stringsAsFactors = F, check.names = FALSE)
MSBB_meta_assay <- read.csv(synapser::synGet('syn22447899')$path, stringsAsFactors = F, check.names = FALSE)

MSBB_meta_biosp <- MSBB_meta_biosp %>% filter(assay == "rnaSeq"& tissue != 'blood')
MSBB_meta_assay <- MSBB_meta_assay %>% filter(assay == "rnaSeq")
#length(unique(MAYO_meta_biosp$individualID)) #620

table(MSBB_meta_biosp$tissue)
table(MSBB_meta_assay$specimenID %in% MSBB_meta_biosp$specimenID)

metadata_temp <- merge(MSBB_meta_assay, MSBB_meta_biosp, by="specimenID")
metadata_temp2 <- merge(metadata_temp, MSBB_meta_ind, by = "individualID")
comb <- metadata_temp2

#summary_list <- summarize_metadata(comb)
#print_metadata_summary(summary_list)

table(comb$platform, useNA = "always")
table(comb$sequencingBatch, useNA = "always")
table(comb$libraryBatch, useNA = "always")
table(comb$libraryPrep, useNA = "always")
table(comb$libraryPreparationMethod, useNA = "always")
table(comb$runType, useNA = "always")
table(comb$assay.x, useNA = "always")
table(comb$nucleicAcidSource, useNA = "always")
table(comb$dataContributionGroup, useNA = "always")
summary(comb$RIN)
summary(comb$rRNA.rate)
summary(comb$PMI)


# diagnosis
comb$diagnosis <- 'OTHER'
comb[ (comb$Braak == "Stage IV" | comb$Braak == "Stage V" |comb$Braak == "Stage VI" ) &
                 ( comb$amyCerad ==  "Frequent/Definite/C3"| comb$amyCerad== "Moderate/Probable/C2") &
                 (comb$CDR >=1), ]$diagnosis <- 'AD'
comb[ (comb$Braak == "Stage III" |comb$Braak == "Stage II"|comb$Braak == "Stage I"|comb$Braak == "None") &
                 (comb$amyCerad == "Sparse/Possible/C1"  | comb$amyCerad== "None/No AD/C0"  ) & 
                 (comb$CDR <= 0.5), ]$diagnosis <- 'CT'
table(comb$diagnosis, useNA="always")
table(comb$diagnosis, comb$sex, useNA="always")


# diag2
comb$diag2 <- 'OTHER2'
comb[ (comb$Braak == "Stage IV" | comb$Braak == "Stage V" |comb$Braak == "Stage VI" ) &
        (comb$amyCerad == "Frequent/Definite/C3" | comb$amyCerad == "Moderate/Probable/C2" ), ]$diag2 <- 'AD2'
comb[ (comb$Braak == "Stage III" |comb$Braak == "Stage II"|comb$Braak == "Stage I"|comb$Braak == "None") &
        (comb$amyCerad == "Sparse/Possible/C1" | comb$amyCerad == "None/No AD/C0" ), ]$diag2 <- 'CT2'
table(comb$diag2, useNA="always")
table(comb$diag2, comb$sex, useNA="always")

table(comb$diag2, comb$diagnosis_prev)

table(comb$diag2, comb$sex, comb$tissue)
table(comb$diag2, comb$tissue)
table(comb$tissue)
summary(comb$RIN)
summary(comb$PMI)
table(comb$apoe4Status)

#apoe4, is now #apoeGenotype (22, 23, 24, 33, 34, 44) and #apoe4Status (yes (if any 4) or no) with clinical_harmonized dataset
table(comb$apoeGenotype, comb$apoe4Status)

comb <- comb %>% mutate(tissue2 = recode(tissue, "frontal pole"  = "FP",
                                         "inferior frontal gyrus"    = "IFG",
                                         "parahippocampal gyrus"  = "PG"  ,   
                                         "prefrontal cortex"  = "PC",
                                         "superior temporal gyrus" = "STG"))

comb$age_cat <- NA
comb$age_cat[comb$ageDeath == "90+"] <- "90+"
comb$age_cat[as.numeric(comb$ageDeath[comb$ageDeath != "90+"])<80] <- "<80"
comb$age_cat[as.numeric(comb$ageDeath[comb$ageDeath != "90+"])>=85] <- "ge85lt90"
comb$age_cat[(as.numeric(comb$ageDeath[comb$ageDeath != "90+"])>=80 & as.numeric(comb$ageDeath[comb$ageDeath != "90+"])<85)] <- "ge80lt85"
comb$age_cat[comb$ageDeath == "90+"] <- "90+"
table(comb$age_cat, useNA = "always")

comb$age_cat2 <- NA
comb$age_cat2[comb$ageDeath == "90+"] <- "90+"
comb$age_cat2[as.numeric(comb$ageDeath[comb$ageDeath != "90+"])<90] <- "<90"
comb$age_cat2[comb$ageDeath == "90+"] <- "90+"
table(comb$age_cat2, useNA = "always")

# Sequencing Statistics 
metrics <- read.csv(synapser::synGet('syn76845876')$path, stringsAsFactors = F)

metrics2 <- metrics %>% dplyr::select(specimenID, "picard_UNPAIRED_READS_EXAMINED"  ,          
                                      "picard_READ_PAIRS_EXAMINED" ,                  
                                      "picard_UNMAPPED_READS"  ,                         
                                      "picard_PERCENT_DUPLICATION"   ,                   
                                      "picard_READS_UNMAPPED"     ,                      
                                      "rsem_alignable_percent" ,                         
                                      "rsem_uniquely_aligned_percent"  ,   
                                      "samtools_reads_mapped", 
                                      "samtools_reads_duplicated",
                                      "samtools_error_rate",                             
                                      "samtools_average_length"  ,                
                                      "samtools_average_quality"  ,                       
                                      "samtools_insert_size_average" ,       
                                      "samtools_percentage_of_properly_paired_reads_...",
                                      "samtools_reads_mapped_percent"  ,                  
                                      "samtools_reads_mapped_and_paired_percent"  ,      
                                      "samtools_reads_properly_paired_percent"   ,       
                                      "samtools_reads_duplicated_percent"     ,                 
                                      "samtools_percent_mapped_X",                        
                                      "samtools_percent_mapped_Y"      ,                 
                                    #  "cutadapt_percent_trimmed_R2" ,                     
                                    #  "cutadapt_percent_trimmed_R1",                     
                                    #  "cutadapt_mean_percent_trimmed" ,  
                                    "cutadapt_percent_trimmed",
                                      "Average.input.read.length"  ,                     
                                      "Uniquely.mapped.reads.."   ,                      
                                      "Number.of.splices..GT.AG" ,                       
                                      "Number.of.splices..GC.AG" ,                        
                                      "Number.of.splices..AT.AC" ,                       
                                      "Mismatch.rate.per.base..." ,                   
                                      "Deletion.rate.per.base"   ,                                      
                                      "Insertion.rate.per.base"  ,                        
                                      "X..of.reads.mapped.to.multiple.loci"  ,           
                                      "X..of.reads.mapped.to.too.many.loci"    ,         
                                      "X..of.reads.unmapped..too.many.mismatches"   ,    
                                      "X..of.reads.unmapped..too.short" ,                
                                      "X..of.reads.unmapped..other" ,                    
                                      "X..of.chimeric.reads" ) 

md <- merge(comb, metrics2, by.x="specimenID", all.x=T)


# add sequencing statistics below ASAP - until then use this
file_path <- "MSBB_md_all.csv"
write.csv(md, file = file_path, row.names = FALSE)
file <- synapser::synStore(synapser::File(path = file_path, name="MSBB_md_all.csv", parent = "syn76814282"))



# checking technical vars
#source("C:/Users/jlundin/OneDrive/git_code/RNASeq_DE/AMP_AD1/RNASeq_DE/AMP_AD1/functions/functions_summarize_metadata.R")
