require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
library(dplyr)
library(tidyverse)
library(fastDummies)
library(ggeffects)
library(writexl)
setwd("./TFM")

#Log regression exposures-BP_cat
#Load data
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined
#Load data
pheno <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(pheno)<-pheno$HelixID
pheno<-pheno[order(rownames(pheno)),]
pheno <- pheno[c("hs2_BPcat_v3_2017_bin_Time2", "HelixID")]

pheno$hs2_BPcat_v3_2017_bin_Time2 <- as.factor(pheno$hs2_BPcat_v3_2017_bin_Time2)

#Load projections 
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_north) <- names(projections$`Projections N`)
proj_south <- projections$`Projections S` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_south) <- names(projections$`Projections S`)
proj_all <-  projections$`Projections all` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_all) <- names(projections$`Projections all`)

#phenotype data
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)

rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all), 
                         c("HelixID","e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north),
                       c("HelixID","e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC")]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south),
                       c("HelixID","e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA")]

#dataset containing proj and BP trajectories 
pheno_all <- pheno[rownames(pheno) %in% rownames(proj_all),]
pheno_n <- pheno[rownames(pheno) %in% rownames(proj_north),]
pheno_s <- pheno[rownames(pheno) %in% rownames(proj_south),]

all(rownames(proj_all)==rownames(pheno_all))
all(rownames(proj_north)==rownames(pheno_n))
all(rownames(proj_south)==rownames(pheno_s))

proj_all$HelixID <- rownames(proj_all)
proj_north$HelixID <- rownames(proj_north)
proj_south$HelixID <- rownames(proj_south)

df_all <- merge(proj_all, pheno_all, by = "HelixID") %>% mutate(Population = "All")
df_all <- merge(df_all, phenotype_all, by = "HelixID") %>% mutate(Population = "All")
df_n <- merge(proj_north, pheno_n, by = "HelixID") %>% mutate(Population = "North")
df_n <- merge(df_n, phenotype_n, by = "HelixID") %>% mutate(Population = "North")
df_s <- merge(proj_south, pheno_s, by = "HelixID") %>% mutate(Population = "South")
df_s <- merge(df_s, phenotype_s, by = "HelixID") %>% mutate(Population = "South")



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
exposome <- tibble::rownames_to_column(exposome, var = "HelixID")

exposome <- dummy_cols(exposome,
                       select_columns = to_factor,
                       remove_selected_columns = T,
                       remove_first_dummy = T)
rownames(exposome) <- exposome$HelixID

exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]


df_all <- merge(exposome_all, pheno_all, by = "HelixID") %>% mutate(Population = "All")
df_n <- merge(exposome_n, pheno_n, by = "HelixID") %>% mutate(Population = "North")
df_s <- merge(exposome_s, pheno_s, by = "HelixID") %>% mutate(Population = "South")


covariates_all <- c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")
covariates_all <- paste(covariates_all, collapse = " + ")

res_all <- data.frame()

for (x in colnames(exposome_all)[-1]) { 
  formula <- as.formula(paste("hs2_BPcat_v3_2017_bin_Time2 ~", x, "+", covariates_all))
  log_model <- glm(formula, data = cbind(df_all, phenotype_all), family = binomial())
  
  coef_val <- coef(summary(log_model))[2, "Estimate"]
  se_val <- coef(summary(log_model))[2, "Std. Error"]
  p_val <- coef(summary(log_model))[2, "Pr(>|z|)"]
  
  res_all <- rbind(res_all, data.frame(
    Exposure = x,
    Coefficient = coef_val,
    StdError = se_val,
    p = p_val
  ))
}

covariates_n <- c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC")
covariates_n <- paste(covariates_n, collapse = " + ")

res_n <- data.frame()

for (x in colnames(exposome_n)[-1]) { 
  formula <- as.formula(paste("hs2_BPcat_v3_2017_bin_Time2 ~", x, "+", covariates_n))
  log_model <- glm(formula, data = cbind(df_n, phenotype_n), family = binomial())
  
  coef_val <- coef(summary(log_model))[x, "Estimate"]
  se_val <- coef(summary(log_model))[x, "Std. Error"]
  p_val <- coef(summary(log_model))[x, "Pr(>|z|)"]
  
  res_n <- rbind(res_n, data.frame(
    Exposure = x,
    Coefficient = coef_val,
    StdError = se_val,
    p = p_val
  ))
}

covariates_s <- c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA")
covariates_s <- paste(covariates_s, collapse = " + ")

res_s <- data.frame()
for (x in colnames(exposome_s)[-1]) { 
  formula <- as.formula(paste("hs2_BPcat_v3_2017_bin_Time2 ~", x, "+", covariates_s))
  log_model <- glm(formula, data = cbind(df_s, phenotype_s), family = binomial())
  
  coef_val <- coef(summary(log_model))[2, "Estimate"]
  se_val <- coef(summary(log_model))[2, "Std. Error"]
  p_val <- coef(summary(log_model))[2, "Pr(>|z|)"]
  
  res_s <- rbind(res_s, data.frame(
    Exposure = x,
    Coefficient = coef_val,
    StdError = se_val,
    p = p_val
  ))
}

res_n$p_corrected <- p.adjust(res_n$p, method = "fdr")
res_s$p_corrected <- p.adjust(res_s$p, method = "fdr")
res_all$p_corrected <- p.adjust(res_all$p, method = "fdr")

log_reg_BP_cat_exposures <- list("North"=res_n,
                                 "South"=res_s,
                                 "Pooled"=res_all)

write_xlsx(log_reg_BP_cat_exposures, "./results/Log_regression/BP_cat_adol-exposome.xlsx")
