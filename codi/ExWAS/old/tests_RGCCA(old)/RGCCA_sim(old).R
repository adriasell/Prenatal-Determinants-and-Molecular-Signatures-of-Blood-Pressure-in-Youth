library(tidyverse)
library(RGCCA)
library(Biobase)
library(dplyr)
#devtools::install_github(repo="https://github.com/rgcca-factory/RGCCA.git")
getwd()
setwd("./TFM")

ncore.use <- parallel::detectCores() - 1
set.seed(1999)

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

#SOLUCIÓ CUTRE ERROR RGCCA_PRED (ERROR A RGCCA_TRANSF0RM, projection)
colnames(X_test$prot) <- paste0("prot_", colnames(X_test$prot))
colnames(X_test$serum) <- paste0("serum_", colnames(X_test$serum))
colnames(X_test$urine) <- paste0("urine_", colnames(X_test$urine))
colnames(X_test$Y) <- paste0("Y_", colnames(X_test$Y))

### Step 1: tuning parameters using cross validation for number of components  #

# Number of components  
print("- Tuning of number of components per blocks")
cv_sgcca_ncomp  <- rgcca_cv(
  blocks = X_train,
  response = 4,
  method = "rgcca",
  connection = 1 - diag(4),
  par_type = "ncomp",
  par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, Y = 1:3),
  n_cores = ncore.use)

plot(cv_sgcca_ncomp)
ncomp = cv_sgcca_ncomp$best_params
write_rds(cv_sgcca_ncomp, "./cv_sgcca_ncomp.rds")

## Step 2: RGCCA with optimized parameters ##
cv_sgcca_ncomp<-readRDS("./cv_sgcca_ncomp.rds")
rgcca_res <- rgcca(cv_sgcca_ncomp)
summary(rgcca_res)
 
## Step 3: Interpretation of the model and performance ##

# Significance: bootstrap (Time estimated: 2 min)
bootstrap <- rgcca_bootstrap(rgcca_res,n_cores = ncore.use)
View(bootstrap$stats)

# Plots loading for all
plot(bootstrap,comp=1,n_mark=20,display_order=T,type="weight")

# Plots loading for specific components
plot(bootstrap,comp=1,block=1,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=1,block=2,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=1,block=3,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=1,block=4,n_mark=30,display_order=T,type="weight")

### Step 4: Prediction ability on test set 

# With rgcca_predict, you can get standard indicators (RMSE, MAE, R2) using different predictor models
rgcca_predict_res <- rgcca_predict(rgcca_res,blocks_test=X_test,prediction_model="lm")
rgcca_predict_res$metric

# Projection of the variables on the latent components for the test set
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()
name_components <- c("Prot","Serum","Urine")
colnames(latent_variables)<- name_components
names_out<-colnames(X_test$Y)

# Estimate R2 of each latent variable by runing a linear model with only this component
for (name in names_out){
  for (comp in name_components){
      data_r2 <- cbind(latent_variables,X_test$Y) %>% as.data.frame() 
      form = paste0(" X_test$Y$",name, "~", comp)
      res <- summary(lm(as.formula(form),data=data_r2))
      r2 <- res$r.squared
      print(paste0("R2 ", name, " of ",comp," = ", 100*signif(r2,2),"%"))
      }
}
# Estimate correlation between components
cor(latent_variables) %>% round(2)

