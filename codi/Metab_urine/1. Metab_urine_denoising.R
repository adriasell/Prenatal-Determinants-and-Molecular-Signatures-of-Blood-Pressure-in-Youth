###############################################################################
##################### Metab preprocessing & denoising #########################
###############################################################################

library(lumi)
library(doParallel)
library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(dplyr)
library(omics)
library(ggplot2)
library(corrplot)
library(isva)
library(SmartSVA)
library(PCAtools)
library(smplot2)
library(fastDummies)
library(polycor)

setwd("./TFM")
output_path <- "./results/denoising/denoising_metab/urine"
getwd()
##············· PREPARING VARIABLES FOR MODELLING

# /// Load HelixID by group (N/S)

helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)


# /// Loading external functions
source("./script/denoising_2024/functions/extract_data_by_xlsx_styles.R")
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

# /// Loading data
# ···· Metadata .RData filepath
metadataFile <- "./script/denoising_2024/metadata/HELIX_SVA_common_OmicsMetadata_20231026.RData" 

# ···· Omic data (Rdata file with an ExpressionSet or GenomicRatioSet)
# ··················· INDICAR RDATA NUEVO CON LA WINSORIZACIÓN
omicFile <- paste0(output_path,"/urine_winsorized.RDS")

# ···· Phenotype data 
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

rownames(phenotype) <- phenotype$HelixID


names(phenotype)[grep("cohort.x", names(phenotype))] <- "cohort"
names(phenotype)[grep("h_ethnicity_c.x", names(phenotype))] <- "h_ethnicity_c"
names(phenotype)[grep("e3_sex_Time1", names(phenotype))] <- "e3_sex"


# ···· Modify phenotype, cohort and sex, with dummies groups (1 vs. all columns)
phenotype$h_ethnicity_c <- ifelse(phenotype$h_ethnicity_c %in% c("Asian","Pakistani"), "Asian_pakistani", phenotype$h_ethnicity_c)

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_ethnicity_c",
                        remove_selected_columns = TRUE)

phenotype <- dummy_cols(phenotype,
                        select_columns = "cohort",
                        remove_selected_columns = TRUE)

phenotype <- dummy_cols(phenotype,
                        select_columns = "e3_sex",
                        remove_selected_columns = TRUE)

phenotype <- phenotype[,colnames(phenotype)!=c("h_ethnicity_c_NA")]

#Convert var to factor
eth_dummies <- names(phenotype)[grep("h_ethnicity_", names(phenotype))]
phenotype[eth_dummies] <- lapply(phenotype[eth_dummies], as.factor)

coh_dummies <- names(phenotype)[grep("cohort", names(phenotype))]
phenotype[coh_dummies] <- lapply(phenotype[coh_dummies], as.factor)


# /// Main variables definition
# 1. Omic used: Define omics used in analysis omicLayer variable.
#   Possible values are: 
#       Methyl   : Methylation data
#       Gene     : Genomic data
#       Prot     : Proteomic data
#       miRNAs   : micro-RNAs
#       Serum    : Serum methabolome
#       Urine    : Urine metabolome
omicsLayer <- "Urine"

# 2. ID column name in the phenotype dataset 
#   It relates the phenotype to the Omic data. (by default SampleID)
phenotypeID <- "SampleID"


# 3. Variables to be used from phenotype dataset (by default ALL)
#  Examples:
#       phenoVariables <- c('bmi', 'sex', 'age', 'smoke')  # 
#       phenoVariables <- 'ALL'         # If we need ALL the variables present in phenotype dataset
#       phenoVariables <- ''            # If no value is defined then ALL the variables present in the phenotype will be used


# LAS VARIABLES ASIGNADAS EN LOS DATOS FILTRADOS
# c("hs2_zdia_bp.v3_2017_Time1","hs2_zsys_bp.v3_2017_Time1","hs2_zbmi_who_Time1","hs2_height_c_zscore_Time1","e3_sex_Time1","hs2_visit_age_years_Time1")

phenotypeVariables <- c("ALL")

# 4. Ethnicity - Population
#   Ethnic parameter used to make a subset of the data and also to take in to account
#   the PCs related to the population. 
#   For European (ONLY if we are using EUR population and no other) we use the PCs calculated from the EUR population, 
#   otherwise we use the PCs computed from ALL populations 
#   Define ethnicity to filter
#   Possible values are: 
#       AFR     :  African
#       AMR     :  Ad Mixed American
#       EAS     :  East Asian
#       EUR     :  European
#       SAS     :  South Asian
#       MIXED   :  Mixed population (not classified in the other categories)
#       ALL     : No filter
#       
#  Example:
#       ethnic <- 'ALL'                     # All populations - no filter
#       ethnic <- ''                        # All populations - no filter
#       ethnic <- c('AFR', 'AMR', 'SAS')    # Use African, Mixed American and South Asian populations
# 

ethnic <- c('ALL')


#### STEP 3. Data preparation
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

omics <- getfullOmicsPhenotype(omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, ethnic )    

#Imputation hs_dift_mealblood_imp
imp_value<-median(omics$data$hs_dift_mealblood_imp, na.rm = T)
omics$data$hs_dift_mealblood_imp[is.na(omics$data$hs_dift_mealblood_imp)]<-imp_value

omics_n <- omics$data[,which(pData(omics$data)$HelixID.x %in% helixid_n)]
omics_s <- omics$data[,which(pData(omics$data)$HelixID.x %in% helixid_s)]

# ////////////// Denoising using linear mixed effects model

# Denoising variables:
#   Define variables to be used to denoise data. Suggested vars to be used for denoising 
#   each omic are: 
#   
#       * Metylation
#           - nano_conc_ng_ul
#           - extr_batch
#           - Sample_Plate
#           - BATCH_Conversion - IMPORTANT
#           - BATCH_Infinium
#           - Array - MOST IMPORTANT
#           - Slide - IMPORTANT
#           - SampleType
# 
#       * Gene expression:
#           - round_gexp
#           - r_b IMPORTANT
#           - round_rna_extr
#           - extr_batch
#           - bio_conc_denat_ngul IMPORTANT
#           - bio_rin_denat IMPORTANT
#           - extraction_successful_class
#           - sample_milky_blood_class
#           - SampleType
#           
#       * miRNA
#           - initial_exclusion
#           - batch - IMPORTANT
#           - slide - IMPORTANT 
#           - final_selection
#           - bio_conc_denat_ngul - IMPORTANT
#           - bio_rin_denat - IMPORTANT
#           - extraction_successful_class
#           - round_rna_extr
#           - sample_milky_blood_class
#           - SampleType
#           
#       * Proteome
#           - plate - IMPORTANT
#           - SampleType
#           
#       * Urina
#           - Sample.ID
#           - Centre
#           - Child_ID
#           - SamplingType - IMPORTANT
#           
#       * Serum
#           - run_order1 - IMPORTANT
#           - run_order2 - IMPORTANT
#           - plate_batch - IMPORTANT
#           - SampleType
#   
#   
#   NOTE: 
#   
#   Cohort, sex and ethnicity are also important but perhaps they can be used in the final 
#   model during the data analysis process

#   SUGGESTED CONFIGURATION:
#       Please add or remove variables as needed
#       
#       * Methylation:
#           denoise_vars <- c("Methyl.nano_conc_ng_ul", "Methyl.extr_batch", "Methyl.Sample_Plate", 
#                             "Methyl.BATCH_Conversion", "Methyl.BATCH_Infinium", "Methyl.Array", 
#                             "Methyl.Slide", "Methyl.SampleType")
# 
#       * Gene Expression:
#           denoise_vars <- c("Gene.round_gexp", "Gene.r_b", "Gene.round_rna_extr", 
#                             "Gene.extr_batch", "Gene.bio_conc_denat_ngul", "Gene.bio_rin_denat", 
#                             "Gene.extraction_successful_class", "Gene.sample_milky_blood_class", "Gene.SampleType")
#   
#       * miRNAs:
#           denoise_vars <- c("miRNA.initial_exclusion", "miRNA.batch", "miRNA.slide", 
#                             "miRNA.final_selection", "miRNA.bio_conc_denat_ngul", "miRNA.bio_rin_denat", 
#                             "miRNA.extraction_successful_class", "miRNA.round_rna_extr", 
#                             "miRNA.sample_milky_blood_class", "miRNA.SampleType")
# 
#       * Proteome:
#           denoise_vars <- c("Prot.plate", "blood_sam4")
#  
#       * Serum metabolome:
#           denoise_vars <- c("Serum.run_order1", "Serum.run_order2", "Serum.plate_batch", "Serum.SampleType")
#  
#       * Urine metabolome:
#           denoise_vars <- c("Urine.Sample.ID", "Urine.Centre", "Urine.Child_ID", "Urine.SamplingType")
#       
#   

source("./script/denoising_2024/functions/lm_local.R")
source("./script/denoising_2024/functions/mlmer_local.R")
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")

denoise_vars_n <- c("e3_bw", "hs2_visit_age_years_Time1", "e3_sex","h_ethnicity_c_Asian_pakistani",
                    "h_ethnicity_c_Caucasian",'cohort_BIB', "cohort_EDEN", "cohort_KANC")

denoise_vars_s <- c("e3_bw", "hs2_visit_age_years_Time1", "e3_sex", "cohort_SAB")

denoise_vars <- c("e3_bw", "hs2_visit_age_years_Time1", "e3_sex","h_ethnicity_c_Asian_pakistani",
                  "h_ethnicity_c_Caucasian",'cohort_BIB', "cohort_EDEN", "cohort_KANC", "cohort_MOBA","cohort_SAB")

nas_remove_n <- rownames((pData(omics_n))[,denoise_vars_n])[apply((pData(omics_n))[,denoise_vars_n], 1, function(x) any(is.na(x)))]
nas_remove_s <- rownames((pData(omics_s))[,denoise_vars_s])[apply((pData(omics_s))[,denoise_vars_s], 1, function(x) any(is.na(x)))]
nas_remove <- rownames((pData(omics$data))[,denoise_vars])[apply((pData(omics$data))[,denoise_vars], 1, function(x) any(is.na(x)))]

# Removing samples with NAs in interest variables
omics_n <- omics_n[,!colnames(omics_n) %in% nas_remove_n]
omics_s <- omics_s[,!colnames(omics_s) %in% nas_remove_s]
omics$data <- omics$data[,!colnames(omics$data) %in% nas_remove]


#omics_n_denoised <- denoising_mlmer(omics_n, omicsLayer, denoise_vars_n, ncores = 18, blocksize = 100)
#omics_s_denoised <- denoising_mlmer(omics_s, omicsLayer, denoise_vars_s, ncores = 18, blocksize = 100)
omics_denoised <- denoising_mlmer(omics$data, omicsLayer, denoise_vars, ncores = 18, blocksize = 100)


# ///// Save denoised data
#saveRDS(omics_n_denoised, file = paste0(output_path, "/metab_urine_denoised_n.RDS"))
#saveRDS(omics_s_denoised, file = paste0(output_path, "/metab_urine_denoised_s.RDS"))
saveRDS(omics_denoised, file = paste0(output_path, "/metab_urine_denoised.RDS"))


#-------------------------------------------------------------
#               Check noised VS. denoised data
#-------------------------------------------------------------

# Correlation plot variables:
#   Define variables to be used in correlation plot. Please set which variables should be
#   used as numerical and which as a factor
#   
#       num_vars: Numerical variables or variables to be used as numerical
#       factor_vars: Factor variables or variables to be used as factor
#       
#           Example: 
#               num_vars <- c("e3_sex", "blood_sam4", "Prot.plate")
#               factor_vars <- c("cohort", "h_ethnicity_c")
#               factor_vars <- NA                                       # If we don't have factor variables
# Define the number of first principal components to be used in check, 
# If nComponents is NULL or NA, all components will be used 
# 
#       Example: 
#           nComponents <- 10 # Checks first 10 compoents
#           nComponents <- NA # Checks all the components
#           
omics_n_denoised <- readRDS(paste0(output_path, "/metab_urine_denoised_n.RDS"))
omics_s_denoised <- readRDS(paste0(output_path, "/metab_urine_denoised_s.RDS"))
omics_denoised <- readRDS(paste0(output_path, "/metab_urine_denoised.RDS"))


nComponents <- 10

check_noised_denoised(original = omics_n, denoised = omics_n_denoised, 
                      Vars = denoise_vars_n, prefix = omicsLayer, topPCA = nComponents,
                      output_path=paste0(output_path,"/corr_plot_n"))

check_noised_denoised(original = omics_s, denoised = omics_s_denoised, 
                      Vars = denoise_vars_s, prefix = omicsLayer, topPCA = nComponents,
                      output_path=paste0(output_path,"/corr_plot_s"))

check_noised_denoised(original = omics$data, denoised = omics_denoised, 
                      Vars = denoise_vars, prefix = omicsLayer, topPCA = nComponents,
                      output_path=paste0(output_path,"/corr_plot"))


#Scatter plot
source("./codi/functions/scatter_plot_function.R")

dfs <- c("raw"=omics$data,"raw_n"=omics_n, "raw_s"=omics_s, "denoised"=omics_denoised, "denoised_n"=omics_n_denoised, "denoised_s"=omics_s_denoised)
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2")
vars <- rownames(omics$data)


for (i in 1:length(dfs)){
  name <- names(dfs)[i]
  for (var in vars){
    for (outcome in outcomes){
      scatter_plot(data = dfs[[i]], x = var, y = outcome, omicsLayer = omicsLayer) }}}
