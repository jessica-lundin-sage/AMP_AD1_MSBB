
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
MSBB_dx_prev <- read.csv(synapser::synGet('syn77471281')$path, stringsAsFactors = F, check.names = FALSE) #from metadata harmonization study
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

md <- merge(comb, metrics, by.x="specimenID", all.x=T)


## from previous processing of this data
# https://github.com/Sage-Bionetworks/ampad-rnaseq-reprocessing/blob/1d5d2e88987714546c92fd2a8564806f148816fa/code/metadata_preprocessing/msbb_reprocessing.R#L365
sex_swapped <- c( 'hB_RNA_12901', 'hB_RNA_12934', 'BM_36_296', 'hB_RNA_10622',
                  'hB_RNA_10622_L43C014', 'hB_RNA_10702', 'hB_RNA_10702_E007C014', 'hB_RNA_7765', 'hB_RNA_7765_resequenced',
                  'hB_RNA_7995', 'hB_RNA_8015', 'hB_RNA_8025', 'hB_RNA_8025_E007C014', 'hB_RNA_8385', 'hB_RNA_8305')
expression_outliers <- c('hB_RNA_10577', 'hB_RNA_12964', 'hB_RNA_13144', 'hB_RNA_13397', 'hB_RNA_10567', 'hB_RNA_10577', 'hB_RNA_7855', 'hB_RNA_9005')
# Tossed lower RIN, when tied tossed lower AlignmentSummaryMetrics_PF_READS_ALIGNED
resequenced_tossed <- c( 
  "BM_22_178", "BM_22_270", "hB_RNA_7755", "hB_RNA_7995", "hB_RNA_7765_resequenced", "hB_RNA_8025", 
  "hB_RNA_8025_E007C014", "BM_22_254", "hB_RNA_9139", "hB_RNA_9144", "hB_RNA_9147", "BM_22_11", 
  "hB_RNA_9183", "hB_RNA_9186", "hB_RNA_9191", "hB_RNA_9191_L43C014", "BM_22_84", "BM_22_251", 
  "hB_RNA_9222", "BM_22_31", "hB_RNA_10232", "hB_RNA_10242", "hB_RNA_10252", "hB_RNA_12171", 
  "hB_RNA_12181", "hB_RNA_12191", "hB_RNA_12201", "hB_RNA_12211", "hB_RNA_8075_resequenced",
  "hB_RNA_8085_resequenced", "hB_RNA_8115_resequenced", "hB_RNA_8175_resequenced", "hB_RNA_8265",
  "hB_RNA_8265_resequenced", "hB_RNA_8285", "hB_RNA_8295", "hB_RNA_9207_resequenced", "hB_RNA_9210_resequenced",
  "hB_RNA_9212_resequenced", "hB_RNA_8125_resequenced", "hB_RNA_8155_resequenced", "hB_RNA_8215_resequenced",
  "hB_RNA_9209_resequenced", "hB_RNA_10342", "hB_RNA_10302", "hB_RNA_8055", "hB_RNA_8185", "hB_RNA_8225_L43C014",
  "hB_RNA_9115_B82C014", "hB_RNA_9166_L43C014", "hB_RNA_9178", "hB_RNA_9180_L43C014", "hB_RNA_9187_L43C014",
  "hB_RNA_9189_E007C014", "hB_RNA_9208", "hB_RNA_12252", "hB_RNA_12252_resequenced", "hB_RNA_12262", 
  "hB_RNA_12262_resequenced", "BM_36_360", "hB_RNA_12272", "hB_RNA_12282", "hB_RNA_12282_resequenced", 
  "hB_RNA_12292", "hB_RNA_12292_resequenced", "hB_RNA_12312", "hB_RNA_12312_resequenced", "hB_RNA_12322",
  "hB_RNA_12322_resequenced", "hB_RNA_12332", "hB_RNA_12332_B18C014", "hB_RNA_12342", "hB_RNA_12342_resequenced",
  "hB_RNA_12352", "hB_RNA_12352_resequenced", "hB_RNA_12362", "hB_RNA_12372", "hB_RNA_12372_resequenced", "BM_36_407",
  "hB_RNA_12382_resequenced", "hB_RNA_12392", "hB_RNA_12402", "hB_RNA_9085", "hB_RNA_9105", "hB_RNA_10482_resequenced",
  "hB_RNA_10492_resequenced", "hB_RNA_10522_resequenced", "hB_RNA_10532", "hB_RNA_10542_resequenced", "hB_RNA_10552",
  "hB_RNA_10567", "hB_RNA_10577", "hB_RNA_10583", "hB_RNA_10617", "hB_RNA_10632", "hB_RNA_10642", "hB_RNA_10652", 
  "hB_RNA_10662", "hB_RNA_10672", "hB_RNA_10682", "hB_RNA_10692", "hB_RNA_10712_resequenced", 
  "hB_RNA_10722_resequenced", "hB_RNA_10742_resequenced", "hB_RNA_10762_resequenced", "hB_RNA_10822",
  "hB_RNA_10832", "hB_RNA_10842", "hB_RNA_10852", "hB_RNA_10862", "hB_RNA_10872", "hB_RNA_10882",  
  "hB_RNA_10502_resequenced", "hB_RNA_10802_resequenced", "hB_RNA_10512_L43C014", "hB_RNA_10622", 
  "hB_RNA_10702", "hB_RNA_10782", "hB_RNA_10992", "hB_RNA_11002", "hB_RNA_12302_E007C014", 
  "hB_RNA_5041", "hB_RNA_16245_E008C189", "hB_RNA_16715_E009C189", "hB_RNA_16735_E009C189",
  "hB_RNA_16895_E009C189", "hB_RNA_16905_E009C189", "hB_RNA_16965_E009C189",  "hB_RNA_17125_E009C189",
  "hB_RNA_4398_E007C014", "hB_RNA_4631_E007C014", "hB_RNA_4720", "hB_RNA_4751_L43C014", "hB_RNA_4774",
  "hB_RNA_4791", "hB_RNA_4801_L43C014", "hB_RNA_4862", "hB_RNA_4881_L43C014", "hB_RNA_4891_L43C014", 
  "hB_RNA_4923_L43C014", "hB_RNA_4946_L43C014", "hB_RNA_4951_L43C014", "hB_RNA_4961_L43C014", "hB_RNA_4980", 
  "hB_RNA_5011_L43C014", "hB_RNA_5021", "hB_RNA_5031_L43C014", "hB_RNA_8555", "hB_RNA_8935",
  "hB_RNA_8675_L43C014", "BM_10_638", "BM_10_687", "BM_10_634", "hB_RNA_13266", "hB_RNA_13276",
  "hB_RNA_13294", "BM_10_727", "hB_RNA_13389", "hB_RNA_13406", "BM_10_636", "BM_10_598",
  "BM_10_742", "BM_10_554", "BM_10_606", "hB_RNA_13500", "BM_10_620", "hB_RNA_13518", 
  "hB_RNA_13547", "BM_10_627", "hB_RNA_13631", "BM_10_557", "BM_10_673", "hB_RNA_13058",
  "hB_RNA_13068_resequenced", "hB_RNA_13081", "hB_RNA_13048_resequenced", "hB_RNA_13216_resequenced", 
  "hB_RNA_13032_L43C014", "hB_RNA_13375")
total_toss <- c(sex_swapped,expression_outliers,resequenced_tossed)

md2 <- md[!(md$specimenID %in% total_toss),]


# add sequencing statistics below ASAP - until then use this
file_path <- "MSBB_md_all.csv"
write.csv(md, file = file_path, row.names = FALSE)
file <- synapser::synStore(synapser::File(path = file_path, name="MSBB_md_all.csv", parent = "syn76814282"))



# checking technical vars
#source("C:/Users/jlundin/OneDrive/git_code/RNASeq_DE/AMP_AD1/RNASeq_DE/AMP_AD1/functions/functions_summarize_metadata.R")
