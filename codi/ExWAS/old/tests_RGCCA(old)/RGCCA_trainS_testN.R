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

### Step 0: prepare the dataset
# Omic data
###Prot
prot_n_eset <- readRDS("./results/denoising/denoising_prot/prot_denoised_n.RDS")
prot_n <- data.frame(t(exprs(prot_n_eset)))
prot_n <- prot_n[complete.cases(prot_n),]
prot_n$HelixID <- pData(prot_n_eset)$HelixID.x
rownames(prot_n) <- prot_n$HelixID
prot_s_eset <- readRDS("./results/denoising/denoising_prot/prot_denoised_s.RDS")
prot_s <- data.frame(t(exprs(prot_s_eset)))
prot_s <- prot_s[complete.cases(prot_s),]
prot_s$HelixID <- pData(prot_s_eset)$HelixID.x
rownames(prot_s) <- prot_s$HelixID
###Serum
serum_n_eset <- readRDS("./results/denoising/denoising_metab/serum/metab_serum_denoised_n.RDS")
serum_n <- data.frame(t(exprs(serum_n_eset)))
serum_n <- serum_n[complete.cases(serum_n),]
serum_n$HelixID <- pData(serum_n_eset)$HelixID.x
rownames(serum_n) <- serum_n$HelixID
serum_s_eset <- readRDS("./results/denoising/denoising_metab/serum/metab_serum_denoised_s.RDS")
serum_s <- data.frame(t(exprs(serum_s_eset)))
serum_s <- serum_s[complete.cases(serum_s),]
serum_s$HelixID <- pData(serum_s_eset)$HelixID.x
rownames(serum_s) <- serum_s$HelixID
###Urine
urine_n_eset <- readRDS("./results/denoising/denoising_metab/urine/metab_urine_denoised_n.RDS")
urine_n <- data.frame(t(exprs(urine_n_eset)))
urine_n <- urine_n[complete.cases(urine_n),]
urine_n$HelixID <- pData(urine_n_eset)$HelixID.x
rownames(urine_n) <- urine_n$HelixID
urine_s_eset <- readRDS("./results/denoising/denoising_metab/urine/metab_urine_denoised_s.RDS")
urine_s <- data.frame(t(exprs(urine_s_eset)))
urine_s <- urine_s[complete.cases(urine_s),]
urine_s$HelixID <- pData(urine_s_eset)$HelixID.x
rownames(urine_s) <- urine_s$HelixID

#Comprovations
all(rownames(prot_s)==rownames(serum_s))
all(rownames(serum_s)==rownames(urine_s))
all(rownames(prot_n)==rownames(serum_n))
all(rownames(serum_n)==rownames(urine_n))

# Outcome
Y <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2","HelixID")
rownames(Y) <- Y$HelixID
names(Y)[grep("hs2_zsys_bp_v3_2017", names(Y))] <- "hs2_zsys_bp_v3_2017_Time2"
names(Y)[grep("hs2_zdia_bp_v3_2017", names(Y))] <- "hs2_zdia_bp_v3_2017_Time2"
Y <- Y[outcomes]
Y <- Y[complete.cases(Y),]

#Select complete cases
ids_n <- Reduce(intersect, list(prot_n$HelixID, serum_n$HelixID, urine_n$HelixID, Y$HelixID))
ids_s <- Reduce(intersect, list(prot_s$HelixID, serum_s$HelixID, urine_s$HelixID, Y$HelixID))

prot_n <-  prot_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
prot_s <-  prot_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
serum_n <-  serum_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
serum_s <-  serum_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
urine_n <-  urine_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
urine_s <-  urine_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
Y <- Y %>% arrange(HelixID) %>% dplyr::select(-HelixID) %>%
  mutate(across(everything(), ~ gsub(",", ".", .))) %>%
  mutate(across(everything(), as.numeric)) 

# Divide it in train(North) test(South)
X_test <- list(prot = prot_n[rownames(prot_n) %in% ids_n,],
                serum= serum_n[rownames(serum_n) %in% ids_n,],
                urine= urine_n[rownames(urine_n) %in% ids_n,],
                Y=Y[rownames(Y) %in% ids_n,])

X_train <- list(prot = prot_s[rownames(prot_s) %in% ids_s,],
               serum= serum_s[rownames(serum_s) %in% ids_s,],
               urine= urine_s[rownames(urine_s) %in% ids_s,],
               Y=Y[rownames(Y) %in% ids_s,])


### Step 1: tuning parameters using cross validation for number of components  #

# Sparsity (OPCIONAL)
print("- Tuning of sparsity")
connection = 1 - diag(4)
min.sparsity <- 1 / sqrt(sapply(X_train, ncol))
sparsity_grid <- expand.grid(prot = seq(min.sparsity[1],1, length.out=25),
                             serum = seq(min.sparsity[2],1, length.out=25),
                             urine = seq(min.sparsity[3],1, length.out=25),
                             Y = 0)

set.seed(1899)
cv_sgcca_sparsity <- rgcca_cv(blocks=X_train,
                              response=4,
                              method="sgcca",
                              par_value= sparsity_grid,
                              ncomp=1,
                              par_type="sparsity",
                              n_cores=detectCores(),
                              connection = connection)

sparsity = cv_sgcca_sparsity$best_params
write_rds(cv_sgcca_sparsity, "./results/RGCCA/cv_sgcca_sparsity_SN.rds")
#cv_sgcca_sparsity<-readRDS("./results/RGCCA/cv_sgcca_sparsity_SN.rds")

# Number of components  
print("- Tuning of number of components per blocks")
set.seed(1899)
cv_sgcca_ncomp  <- rgcca_cv(blocks = X_train,
                            response = 4,
                            method = "sgcca",
                            par_type = "ncomp",
                            sparsity = sparsity,
                            par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, Y = 1:3),
                            n_cores = detectCores(),
                            connection = connection)
plot(cv_sgcca_ncomp)
ncomp = cv_sgcca_ncomp$best_params

## Step 2: RGCCA with optimized parameters ##

rgcca_res <- rgcca(cv_sgcca_ncomp)
summary(rgcca_res)
write_rds(rgcca_res, "./results/RGCCA/rgcca_SN.rds")
#rgcca_res <- readRDS("./results/RGCCA/rgcca_SN.rds")

#Stability
#stab<-rgcca_stability(rgcca_res, n_cores = detectCores())
#plot(stab)

## Step 3: Interpretation of the model and performance ##
# Significance: bootstrap (Time estimated: 2 min)
bootstrap <- rgcca_bootstrap(rgcca_res,n_cores = detectCores())
#write_rds(bootstrap, "./results/RGCCA/bootstrap_rgcca_SN.rds")
#bootstrap <- readRDS("./results/RGCCA/bootstrap_rgcca_SN.RData")

# Plots loading for specific components
plot(bootstrap,comp=1,block=1,n_mark=30,display_order=T,type="weight")

plot(bootstrap,comp=1,block=2,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=2,block=2,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=3,block=2,n_mark=30,display_order=T,type="weight")

plot(bootstrap,comp=1,block=3,n_mark=30,display_order=T,type="weight")

plot(bootstrap,comp=1,block=4,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=2,block=4,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=3,block=4,n_mark=30,display_order=T,type="weight")


### Step 4: Prediction ability on test set 

# With rgcca_predict, you can get standard indicators (RMSE, MAE, R2) using different predictor models
rgcca_predict_res <- rgcca_predict(rgcca_res,blocks_test=X_test,prediction_model="lm")
rgcca_predict_res$metric

# Projection of the variables on the latent components for the test set
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()
name_components <- c("Prot","Serum1","Serum2","Serum3","Urine")
colnames(latent_variables)<- name_components
names_out<-colnames(X_test$Y)

# Estimate R2 of each latent variable by runing a linear model with only this component
results_R2 <- data.frame(Outcome = character(), Layer = character(), R2 = numeric())

for (name in names_out) {
  for (comp in name_components) {
    data_r2 <- cbind(latent_variables, X_test$Y) %>% as.data.frame()
    form <- as.formula(paste("X_test$Y$", name, "~", comp))
    r2 <- summary(lm(form, data = data_r2))$r.squared*100
    results_R2 <- rbind(results_R2, data.frame(Outcome = name, Layer = comp, R2 = r2))
  }
}

results_R2

# Estimate correlation between components
cor(latent_variables) %>% round(2)

list_SN<-list("sparsity"=cv_sgcca_sparsity$best_params,
           "ncomp"=cv_sgcca_ncomp$best_params,
           "R2"=results_R2,
           "metrics"=rgcca_predict_res$metric,
           "corr"=cor(latent_variables) %>% round(2),
           "Latent factor"=rgcca_res$Y,
           "Block weights"=rgcca_res$a,
           "projection"=rgcca_predict_res$projection)

save(list_SN, file = "./results/RGCCA/results_rgcca_SN.RData")
load("./results/RGCCA/results_rgcca_SN.RData")

list_SN[1:5]

# Filtra els valors diferents de 0 i manté els rownames
resultat <- lapply(list_SN$`Block weights`, function(x) {
  if (any(x != 0)) {
    x_filtered <- x[x != 0]
    names(x_filtered) <- attr(x, "dimnames")[[1]][x != 0]
    x_filtered
  } else {
    NULL
  }
})
resultat
