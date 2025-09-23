library(tidyverse)
library(RGCCA)
library(Biobase)
library(dplyr)
library(parallel)
library(writexl)
#install.packages("./packages_R/RGCCA-main", repos = NULL, type = "source") RGCCA modified

getwd()
setwd("./TFM")

ncore.use <- parallel::detectCores() - 1
set.seed(1999)

source("./codi/Functions/functions_RGCCA.R")

# Step 0: prepare the dataset---------------------------------------------------

# /// Load HelixID by group (N/S)
helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

# Omic data
### Prot
prot_all <- eset2df("./results/sensitivity_bmi/prot_denoised.RDS")
prot_n <- prot_all[rownames(prot_all) %in% helixid_n,]
prot_s <- prot_all[rownames(prot_all) %in% helixid_s,]

### Serum
serum_all <- eset2df("./results/sensitivity_bmi/metab_serum_denoised.RDS")
serum_n <- serum_all[rownames(serum_all) %in% helixid_n,]
serum_s <- serum_all[rownames(serum_all) %in% helixid_s,]

### Urine
urine_all <- eset2df("./results/sensitivity_bmi/metab_urine_denoised.RDS")
urine_n <- urine_all[rownames(urine_all) %in% helixid_n,]
urine_s <- urine_all[rownames(urine_all) %in% helixid_s,]

#Comprovations
all(rownames(prot_s)==rownames(serum_s))
all(rownames(serum_s)==rownames(urine_s))
all(rownames(prot_n)==rownames(serum_n))
all(rownames(serum_n)==rownames(urine_n))
all(rownames(prot_all)==rownames(serum_all))
all(rownames(urine_all)==rownames(serum_all))


# Outcome
Y <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(Y) <- Y$HelixID
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2","HelixID")
names(Y)[grep("hs2_zsys_bp_v3_2017", names(Y))] <- "hs2_zsys_bp_v3_2017_Time2"
names(Y)[grep("hs2_zdia_bp_v3_2017", names(Y))] <- "hs2_zdia_bp_v3_2017_Time2"
Y <- Y[outcomes]
Y <- Y[complete.cases(Y),]

#Select complete cases
ids_n <- Reduce(intersect, list(prot_n$HelixID, serum_n$HelixID, urine_n$HelixID, Y$HelixID))
ids_s <- Reduce(intersect, list(prot_s$HelixID, serum_s$HelixID, urine_s$HelixID, Y$HelixID))

#Arrange data
prot_n <-  prot_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
prot_s <-  prot_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
prot_all <-  prot_all %>% arrange(HelixID) %>% dplyr::select(-HelixID)
serum_n <-  serum_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
serum_s <-  serum_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
serum_all <-  serum_all %>% arrange(HelixID) %>% dplyr::select(-HelixID)
urine_n <-  urine_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
urine_s <-  urine_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
urine_all <-  urine_all %>% arrange(HelixID) %>% dplyr::select(-HelixID)
Y <- Y %>% 
  arrange(HelixID) %>% dplyr::select(-HelixID) %>%
  mutate(across(everything(), ~ gsub(",", ".", .))) %>%
  mutate(across(everything(), as.numeric))

# Divide it in train(North) test(South)
X_north <- list(prot = prot_n[rownames(prot_n) %in% ids_n,],
                serum= serum_n[rownames(serum_n) %in% ids_n,],
                urine= urine_n[rownames(urine_n) %in% ids_n,],
                Y=Y[rownames(Y) %in% ids_n,])

X_south <- list(prot = prot_s[rownames(prot_s) %in% ids_s,],
                serum= serum_s[rownames(serum_s) %in% ids_s,],
                urine= urine_s[rownames(urine_s) %in% ids_s,],
                Y=Y[rownames(Y) %in% ids_s,])

X_combined <- list(prot = prot_all[rownames(prot_all) %in% c(ids_s, ids_n),],
                   serum = serum_all[rownames(serum_all) %in% c(ids_s, ids_n),],
                   urine = urine_all[rownames(urine_all) %in% c(ids_s, ids_n),],
                   Y = Y[rownames(Y) %in% c(ids_s, ids_n),])

list_X <- list("X_north"=X_north,
               "X_south"=X_south,
               "X_combined"=X_combined)

write_rds(list_X, file = "./results/sensitivity_bmi/RGCCA/list_X.rds")


#Projections with original model(without BMI) in denoised data adjusting BMI
rgcca_res_original <- readRDS("./results/RGCCA/model/rgcca_final.rds")

proj_north1 <- rgcca_predict(rgcca_res_original, blocks_test = X_north, prediction_model = "lm")$projection
proj_south1 <- rgcca_predict(rgcca_res_original, blocks_test = X_south, prediction_model = "lm")$projection
proj_all1 <- rgcca_predict(rgcca_res_original, blocks_test = X_combined, prediction_model = "lm")$projection

projections1<-list("Projections N" = proj_north1,
                   "Projections S" = proj_south1,
                   "Projections all" = proj_all1)

write_rds(projections1, file = "./results/sensitivity_bmi/RGCCA/projections_bmi.rds")

