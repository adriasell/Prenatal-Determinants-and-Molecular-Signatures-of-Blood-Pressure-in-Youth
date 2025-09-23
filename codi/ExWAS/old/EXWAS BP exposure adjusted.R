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
library(patchwork)
#install.packages("./packages_R/ggsankey-main", repos = NULL, type = "source")

################################# Defining working directory.
getwd()
setwd("./TFM")
output_path <-"./results/ExWAS"
set.seed(1899)

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
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all$prot), ]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north$prot),]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south$prot),]

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

# Convert the specified columns to factors in all three data frames
cols_to_factor <- c("h_parity_None", "h_native_None", "h_edumc_None")

# Apply factor conversion only if the columns exist in each dataframe
for (col in cols_to_factor) {
  if (col %in% names(phenotype_n)) phenotype_n[[col]] <- as.factor(phenotype_n[[col]])
  if (col %in% names(phenotype_s)) phenotype_s[[col]] <- as.factor(phenotype_s[[col]])
  if (col %in% names(phenotype_all)) phenotype_all[[col]] <- as.factor(phenotype_all[[col]])
}

## Association between outcomes and the exposures------------------------------

##north
res_exposure_N <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c(
    "h_cohort_BIB", 
    "h_cohort_EDEN", 
    "h_cohort_KANC",
    "h_parity_None",
    "h_native_None",
    "h_edumc_None",
    "h_age_None"), collapse = " + "))
  
  res_exposure_N[[outcome]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n),
                                          expos_name = colnames(exposome_n),
                                          form=form) %>%  mutate(outcome=outcome)}

res_exposure_N <- suppressMessages(purrr::reduce(res_exposure_N, full_join)) %>% arrange(p)

##south
exposome_s <- exposome_s[, colSums(exposome_s != 0) > 0]
res_exposure_S <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c("h_cohort_INMA",
                                      "h_parity_None",
                                      "h_native_None",
                                      "h_edumc_None",
                                      "h_age_None")))
  
  res_exposure_S[[outcome]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s),
                                          expos_name = colnames(exposome_s),
                                          form=form) %>%  mutate(outcome=outcome)}

res_exposure_S <- suppressMessages(purrr::reduce(res_exposure_S, full_join)) %>% arrange(p)

##all
res_exposure_all <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c(
    "h_cohort_BIB", 
    "h_cohort_EDEN", 
    "h_cohort_KANC", 
    "h_cohort_MOBA", 
    "h_cohort_INMA",
    "h_parity_None",
    "h_native_None",
    "h_edumc_None",
    "h_age_None"), collapse = " + "))
  
  res_exposure_all[[outcome]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all),
                                            expos_name = colnames(exposome_all),
                                            form=form) %>%  mutate(outcome=outcome)}

res_exposure_all <- suppressMessages(purrr::reduce(res_exposure_all, full_join)) %>% arrange(p)

#Save results
list_res_exposure_BP <- list("res_exposure_N"=res_exposure_N,
                             "res_exposure_S"=res_exposure_S,
                             "res_exposure_all"=res_exposure_all)

write_xlsx(list_res_exposure_BP, paste0(output_path,"/ExWAS_BP_exposome_adjusted.xlsx"))

