library(tidyverse)
library(RGCCA)
library(Biobase)
library(dplyr)
library(parallel)
library(writexl)
#devtools::install_github(repo="https://github.com/rgcca-factory/RGCCA.git", force=T)


getwd()
setwd("./TFM")

ncore.use <- parallel::detectCores() - 1
set.seed(1899)

### Step 0: prepare the dataset #####

# Omic data
prot_eset <- readRDS("./results/denoising/denoising_prot/prot_denoised.RDS")
prot <- data.frame(t(exprs(prot_eset)))
prot <- prot[complete.cases(prot),]
prot$HelixID <- pData(prot_eset)$HelixID.x
rownames(prot) <- prot$HelixID


serum_eset <- readRDS("./results/denoising/denoising_metab/serum/metab_serum_denoised.RDS")
serum <- data.frame(t(exprs(serum_eset)))
serum <- serum[complete.cases(serum),]
serum$HelixID <- pData(serum_eset)$HelixID.x
rownames(serum) <- serum$HelixID


urine_eset <- readRDS("./results/denoising/denoising_metab/urine/metab_urine_denoised.RDS")
urine <- data.frame(t(exprs(urine_eset)))
urine <- urine[complete.cases(urine),]
urine$HelixID <- pData(urine_eset)$HelixID.x
rownames(urine) <- urine$HelixID


# Outcome
Y <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2","HelixID")
rownames(Y) <- Y$HelixID
names(Y)[grep("hs2_zsys_bp_v3_2017", names(Y))] <- "hs2_zsys_bp_v3_2017_Time2"
names(Y)[grep("hs2_zdia_bp_v3_2017", names(Y))] <- "hs2_zdia_bp_v3_2017_Time2"
Y <- Y[outcomes]
Y <- Y[complete.cases(Y),]

#Select complete cases

ids <- Reduce(intersect, list(prot$HelixID, serum$HelixID, urine$HelixID, Y$HelixID))
prot <- prot[prot$HelixID %in% ids,]
serum <- serum[serum$HelixID %in% ids,]
urine <- urine[urine$HelixID %in% ids,]
Y <- Y[Y$HelixID %in% ids,]

prot <-  prot %>% arrange(HelixID) %>% dplyr::select(-HelixID) 
serum <-  serum %>% arrange(HelixID) %>% dplyr::select(-HelixID)
urine <-  urine %>% arrange(HelixID) %>% dplyr::select(-HelixID)
Y <- Y %>% arrange(HelixID) %>% dplyr::select(-HelixID) %>%
  mutate(across(everything(), ~ gsub(",", ".", .))) %>%
  mutate(across(everything(), as.numeric)) 

#Comprovacions
all(rownames(prot)==rownames(serum))
all(rownames(urine)==rownames(serum))
all(rownames(prot)==rownames(Y))

# Divide it in train test
id_train <- sample(c(TRUE, FALSE), nrow(Y), replace=TRUE, prob=c(0.7,0.3))


X_train <- list(prot = prot[id_train,],
                serum=serum[id_train,],
                urine=urine[id_train,],
                Y=Y[id_train,])


X_test <- list(prot = prot[!id_train,],
               serum=serum[!id_train,],
               urine=urine[!id_train,],
               Y=Y[!id_train,])

data <-list("X_train"=X_train,
            "X_test"=X_test)

save(data, file= "./results/RGCCA/data_X.Rdata")
