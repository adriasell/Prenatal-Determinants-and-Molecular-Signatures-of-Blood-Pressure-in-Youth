################################# Set R environment 
library(writexl)
library(tidyverse)
library(doParallel)
library(parallel)
library(ggsankey)
library(networkD3)
library(sharp)
library(readxl)
library(dplyr)
library(tidyr)
library(fastDummies)
library(htmlwidgets)
library(webshot)
library(purrr)
library(mediation)

################################# Defining working directory.
getwd()
setwd("./TFM")
output_path <-"./results/ExWAS"

source("./codi/Functions/functions_ExWAS.R")
source("./codi/Functions/function_volcano_plot.R")

###Prepare the dataset----------------------------------------------------------
#Load data
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined

#Load projections
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N`
proj_south <- projections$`Projections S`
proj_all <-  projections$`Projections all`

#phenotype data
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)

rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all$prot),]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north$prot),]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south$prot),]

# Convert the specified columns to factors in all three data frames

cols_to_factor <- c("h_edumc_None", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")
# Apply factor conversion only if the columns exist in each dataframe
for (col in cols_to_factor) {
  if (col %in% names(phenotype_n)) phenotype_n[[col]] <- as.factor(phenotype_n[[col]])
  if (col %in% names(phenotype_s)) phenotype_s[[col]] <- as.factor(phenotype_s[[col]])
  if (col %in% names(phenotype_all)) phenotype_all[[col]] <- as.factor(phenotype_all[[col]])
}

#Codebook
codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")

#Load exposure data
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)

#Scale exposures
to_factor <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
               "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter","e3_asmokyn_p_None",
               "e3_alcpreg_yn_None")

exposome[, !colnames(exposome) %in% to_factor] <- scale(exposome[, !colnames(exposome) %in% to_factor], center = T)
exposome[to_factor] <- lapply(exposome[to_factor], factor)
exposome$HelixID <- rownames(exposome)
exposome <- fastDummies::dummy_cols(exposome,
                                    select_columns = to_factor,
                                    remove_first_dummy = TRUE,      
                                    remove_selected_columns = TRUE) 
rownames(exposome) <- exposome$HelixID

exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]

#Prepare df
lat_vars_n <- prepare_df(projections = proj_north, X_data = X_north, exposure_data = exposome_n, response = 4)
lat_vars_s <- prepare_df(projections =  proj_south, X_data = X_south, exposure_data = exposome_s, response = 4)
lat_vars_all <- prepare_df(projections =  proj_all, X_data = X_combined, exposure_data = exposome_all, response = 4)

all(rownames(lat_vars_n)==rownames(phenotype_n))
all(rownames(lat_vars_s)==rownames(phenotype_s))
all(rownames(lat_vars_all)==rownames(phenotype_all))

#Exposures-BPLC
LC_BP <- readRDS("./results/RGCCA/projections_outcome_adol.rds")
LC_BP_n <- LC_BP$`Proj N`
LC_BP_s <- LC_BP$`Proj S`
LC_BP_all <- LC_BP$`Proj all`

dbp_mixture_all <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_all_hs2_zdia_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zdia_bp.v3_2017_Time2")
sbp_mixture_all <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_all_hs2_zsys_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zsys_bp.v3_2017_Time2")
mixture_all <- rbind(dbp_mixture_all, sbp_mixture_all)

#Mixture BPLC
mixture_BP_LC_all <- read.csv2("./results/ExWAS/Mixture_analysis_BP-LC/data/data_all_LC_BP1.csv")
mixture_BP_LC_n <- read.csv2("./results/ExWAS/Mixture_analysis_BP-LC/data/data_n_LC_BP1.csv")
mixture_BP_LC_s <- read.csv2("./results/ExWAS/Mixture_analysis_BP-LC/data/data_s_LC_BP1.csv") 

rownames(mixture_BP_LC_all) <- mixture_BP_LC_all$X

#Inverse diet mixture
mixture_BP_LC_all$diet.mixture0.index_zbmi8 <- mixture_BP_LC_all$diet.mixture0.index_zbmi8*-1
#ExWAS--------------------------------------------------------------------------
rgcca_validation <- readRDS("./results/RGCCA/model/rgcca_final.rds")
response <- rgcca_validation$call$response
ncomp <- rgcca_validation$call$ncomp

library(mediation)
set.seed(291002)


# 1. BMI → prot → SBP
BMI_ProtLC_all <- lm(prot ~ h_mbmi_None + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                   h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                 data = cbind(lat_vars_all, phenotype_all))

BMI_SBP_prot_all <- lm(hs2_zsys_bp.v3_2017_Time2 ~ h_mbmi_None + prot + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                    h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
              data = cbind(lat_vars_all, phenotype_all))

med_BMI_prot_SBP <- mediate(BMI_ProtLC_all, BMI_SBP_prot_all, treat = "h_mbmi_None", mediator = "prot", boot = TRUE, sims = 1000)


# 2. BMI → urine → SBP

BMI_UrineLC_all <- lm(urine ~ h_mbmi_None + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                       h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                     data = cbind(lat_vars_all, phenotype_all))

BMI_SBP_urine_all <- lm(hs2_zsys_bp.v3_2017_Time2 ~ h_mbmi_None + urine + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                    h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                  data = cbind(lat_vars_all, phenotype_all))

med_BMI_urine_SBP <- mediate(BMI_UrineLC_all, BMI_SBP_urine_all, treat = "h_mbmi_None", mediator = "urine", boot = TRUE, sims = 1000)


# 9. MPF → prot → SBP
MPF_prot_all <- lm(prot ~ lifestyle.mixture0.index_zbmi8 + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                     h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                   data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

MPF_SBP_prot <- lm(hs2_zsys_bp.v3_2017_Time2 ~ lifestyle.mixture0.index_zbmi8 + prot + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                     h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                   data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_MPF_prot_SBP <- mediate(MPF_prot_all, MPF_SBP_prot,
                            treat = "lifestyle.mixture0.index_zbmi8",
                            mediator = "prot", boot = TRUE, sims = 1000)

# 10. MPF → urine → SBP
MPF_urine_all <- lm(urine ~ lifestyle.mixture0.index_zbmi8 + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                      h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                    data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

MPF_SBP_urine <- lm(hs2_zsys_bp.v3_2017_Time2 ~ lifestyle.mixture0.index_zbmi8 + urine + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                      h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                    data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_MPF_urine_SBP <- mediate(MPF_urine_all, MPF_SBP_urine,
                             treat = "lifestyle.mixture0.index_zbmi8",
                             mediator = "urine", boot = TRUE, sims = 1000)

# 11. Diet → prot → SBP
diet_prot_all <- lm(prot ~ diet.mixture0.index_zbmi8 + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                      h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                    data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

diet_SBP_prot <- lm(hs2_zsys_bp.v3_2017_Time2 ~ diet.mixture0.index_zbmi8 + prot + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                      h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                    data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_diet_prot_SBP <- mediate(diet_prot_all, diet_SBP_prot,
                             treat = "diet.mixture0.index_zbmi8",
                             mediator = "prot", boot = TRUE, sims = 1000)

# 12. Diet → urine → SBP
diet_urine_all <- lm(urine ~ diet.mixture0.index_zbmi8 + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                       h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                     data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

diet_SBP_urine <- lm(hs2_zsys_bp.v3_2017_Time2 ~ diet.mixture0.index_zbmi8 + urine + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                       h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                     data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_diet_urine_SBP <- mediate(diet_urine_all, diet_SBP_urine,
                              treat = "diet.mixture0.index_zbmi8",
                              mediator = "urine", boot = TRUE, sims = 1000)

# 1. BMI → prot → DBP
BMI_DBP_prot_all <- lm(hs2_zdia_bp.v3_2017_Time2 ~ h_mbmi_None + prot + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                         h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                       data = cbind(lat_vars_all, phenotype_all))

med_BMI_prot_DBP <- mediate(BMI_ProtLC_all, BMI_DBP_prot_all,
                            treat = "h_mbmi_None",
                            mediator = "prot", boot = TRUE, sims = 1000)


# 2. BMI → urine → DBP
BMI_DBP_urine_all <- lm(hs2_zdia_bp.v3_2017_Time2 ~ h_mbmi_None + urine + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                          h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                        data = cbind(lat_vars_all, phenotype_all))

med_BMI_urine_DBP <- mediate(BMI_UrineLC_all, BMI_DBP_urine_all,
                             treat = "h_mbmi_None",
                             mediator = "urine", boot = TRUE, sims = 1000)


# 9. MPF → prot → DBP
MPF_DBP_prot <- lm(hs2_zdia_bp.v3_2017_Time2 ~ lifestyle.mixture0.index_zbmi8 + prot + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                     h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                   data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_MPF_prot_DBP <- mediate(MPF_prot_all, MPF_DBP_prot,
                            treat = "lifestyle.mixture0.index_zbmi8",
                            mediator = "prot", boot = TRUE, sims = 1000)


# 10. MPF → urine → DBP
MPF_DBP_urine <- lm(hs2_zdia_bp.v3_2017_Time2 ~ lifestyle.mixture0.index_zbmi8 + urine + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                      h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                    data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_MPF_urine_DBP <- mediate(MPF_urine_all, MPF_DBP_urine,
                             treat = "lifestyle.mixture0.index_zbmi8",
                             mediator = "urine", boot = TRUE, sims = 1000)


# 11. Diet → prot → DBP
diet_DBP_prot <- lm(hs2_zdia_bp.v3_2017_Time2 ~ diet.mixture0.index_zbmi8 + prot + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                      h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                    data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_diet_prot_DBP <- mediate(diet_prot_all, diet_DBP_prot,
                             treat = "diet.mixture0.index_zbmi8",
                             mediator = "prot", boot = TRUE, sims = 1000)


# 12. Diet → urine → DBP
diet_DBP_urine <- lm(hs2_zdia_bp.v3_2017_Time2 ~ diet.mixture0.index_zbmi8 + urine + h_cohort_INMA + h_cohort_BIB + h_cohort_EDEN +
                       h_cohort_KANC + h_cohort_MOBA + h_parity_None + h_native_None + h_edumc_None + h_age_None,
                     data = cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all))

med_diet_urine_DBP <- mediate(diet_urine_all, diet_DBP_urine,
                              treat = "diet.mixture0.index_zbmi8",
                              mediator = "urine", boot = TRUE, sims = 1000)

extract_mediation_summary <- function(med_obj, analysis_name) {
  s <- summary(med_obj)
  
  data.frame(
    Analysis = analysis_name,
    
    ACME_Est = s$d0, 
    ACME_CI = paste0("[", round(s$d0.ci[1], 3), ", ", round(s$d0.ci[2], 3), "]"),
    ACME_p = s$d0.p,
    
    ADE_Est = s$z0, 
    ADE_CI = paste0("[", round(s$z0.ci[1], 3), ", ", round(s$z0.ci[2], 3), "]"),
    ADE_p = s$z0.p,
    
    Total_Est = s$tau.coef, 
    Total_CI = paste0("[", round(s$tau.ci[1], 3), ", ", round(s$tau.ci[2], 3), "]"),
    Total_p = s$tau.p,
    
    Prop_Med_Est = s$n0, 
    Prop_Med_CI = paste0("[", round(s$n0.ci[1], 3), ", ", round(s$n0.ci[2], 3), "]"),
    Prop_Med_p = s$n0.p,
    
    stringsAsFactors = FALSE
  )
}

# Collect results
mediation_results <- do.call(rbind, list(
  extract_mediation_summary(med_BMI_prot_SBP,   "BMI → prot → SBP"),
  extract_mediation_summary(med_BMI_urine_SBP,  "BMI → urine → SBP"),
  extract_mediation_summary(med_MPF_prot_SBP,   "MPF → prot → SBP"),
  extract_mediation_summary(med_MPF_urine_SBP,  "MPF → urine → SBP"),
  extract_mediation_summary(med_diet_prot_SBP,  "Diet → prot → SBP"),
  extract_mediation_summary(med_diet_urine_SBP, "Diet → urine → SBP"),
  
  extract_mediation_summary(med_BMI_prot_DBP,   "BMI → prot → DBP"),
  extract_mediation_summary(med_BMI_urine_DBP,  "BMI → urine → DBP"),
  extract_mediation_summary(med_MPF_prot_DBP,   "MPF → prot → DBP"),
  extract_mediation_summary(med_MPF_urine_DBP,  "MPF → urine → DBP"),
  extract_mediation_summary(med_diet_prot_DBP,  "Diet → prot → DBP"),
  extract_mediation_summary(med_diet_urine_DBP, "Diet → urine → DBP")
))

# Apply FDR correction afterwards
mediation_results <- mediation_results %>%
  mutate(
    ACME_p_adj     = p.adjust(ACME_p,     method = "BH"),
    ADE_p_adj      = p.adjust(ADE_p,      method = "BH"),
    Total_p_adj    = p.adjust(Total_p,    method = "BH"),
    Prop_Med_p_adj = p.adjust(Prop_Med_p, method = "BH")
  )

print(mediation_results)
mediation_results <- as.data.frame(
  lapply(mediation_results, function(x) {
    if (is.numeric(x)) {
      round(x, 2)
    } else {
      x
    }
  })
)

write_xlsx(mediation_results, "./results/Mediation/mediation_results.xlsx")










#Comprovations_----------------------------------------------------------------
# List all data frames you are cbinding
dfs <- list(lat_vars_all, phenotype_all, mixture_BP_LC_all, LC_BP_all)

# Check if all have the same rownames
all_same_rownames <- all(sapply(dfs, function(df) identical(rownames(df), rownames(dfs[[1]]))))

if (all_same_rownames) {
  message("All data frames have identical rownames in the same order.")
} else {
  message("Rownames differ — check alignment!")
  
  # Optional: find which ones differ
  for (i in seq_along(dfs)) {
    if (!identical(rownames(dfs[[i]]), rownames(dfs[[1]]))) {
      cat("Mismatch in data frame", i, "\n")
    }
  }
}
