library(tidyverse)
library(Biobase)
library(dplyr)
library(parallel)
library(writexl)
library(Deriv)
install.packages("./packages_R/RGCCA-main", repos = NULL, type = "source")
library(RGCCA)
getwd()

ncore.use <- parallel::detectCores() - 1
set.seed(1999)

source("./codi/Functions/functions_RGCCA.R")

# Step 1: tuning parameters using cross validation for number of components-----

print("Data loading")
list_X<-readRDS("./results/RGCCA_methyl/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined
print("Data loaded")

# Sparsity
print("- Tuning of sparsity")
connection = 1-diag(5)

min.sparsity <- 1 / sqrt(sapply(X_north, ncol))
sparsity_grid <- expand.grid(prot = seq(min.sparsity[1],1, length.out=10),
                             serum = seq(min.sparsity[2],1, length.out=10),
                             urine = seq(min.sparsity[3],1, length.out=10),
                             methyl = seq(min.sparsity[4],1, length.out=10),
                             Y = 0)

set.seed(1999)
cv_sgcca_sparsity <- rgcca_cv(blocks=X_north,
                              response=5,
                              method="sgcca",
                              par_value= sparsity_grid,
                              ncomp=1,
                              par_type="sparsity",
                              n_cores=detectCores(),
                              connection = connection)
print("Tuning of sparsity DONE")

write_rds(cv_sgcca_sparsity, "./results/RGCCA_methyl/model/cv_sgcca_sparsity.rds")
print("File saved")
#cv_sgcca_sparsity<-readRDS("./results/RGCCA_methyl/model/cv_sgcca_sparsity.rds")
sparsity = cv_sgcca_sparsity$best_params

# Number of components  
print("- Tuning of number of components per blocks")
set.seed(1999)
cv_sgcca_ncomp  <- rgcca_cv(blocks = X_north,
                            response = 5,
                            method = "sgcca",
                            par_type = "ncomp",
                            sparsity =  cv_sgcca_sparsity$best_params,
                            par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, methyl = 1:5, Y = 1:3),
                            n_cores = detectCores(),
                            connection = connection)
print("- Tuning of number of components per blocks DONE")

write_rds(cv_sgcca_ncomp, "./results/RGCCA_methyl/model/cv_sgcca_ncomp.rds")
print("File saved")

#cv_sgcca_ncomp<-readRDS("./results/RGCCA_methyl/model/cv_sgcca_ncomp.rds")
ncomp = cv_sgcca_ncomp$best_params


# Step 2: RGCCA with optimized parameters---------------------------------------
print("Final model")
rgcca_res <- rgcca(cv_sgcca_ncomp)
write_rds(rgcca_res, "./results/RGCCA_methyl/model/rgcca_validation.rds")
#rgcca_res <- readRDS("./results/RGCCA_methyl/model/rgcca_validation.rds")
print("File saved")

summary(rgcca_res)

#Bootstrap
rgcca_res$a
print("Bootstraping")
set.seed(1999)
bootstrap <- rgcca_bootstrap(rgcca_res,n_boot = 5000, n_cores = detectCores())
write_rds(bootstrap, "./results/RGCCA_methyl/model/bootstrap.rds")
print("File saved")

rgcca_res <- bootstrap$rgcca
write_rds(rgcca_res, "./results/RGCCA_methyl/model/rgcca_final.rds")
print("File saved")

#rgcca_res <- readRDS("./results/RGCCA_methyl/model/rgcca_final.rds")


# Step 3: Prediction ability on test set---------------------------------------- 
print("Predicition ability")
# With rgcca_predict, you can get standard indicators (RMSE, MAE, R2) using different predictor models
rgcca_predict_res <- rgcca_predict(rgcca_res,blocks_test=X_south,prediction_model="lm")


# Projection of the variables on the latent components for the test set
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()
name_components <- paste0(rep(names(X_north)[-5],times=cv_sgcca_ncomp$best_params[-5]),unlist(lapply(cv_sgcca_ncomp$best_params[-5],seq)))
colnames(latent_variables)<- name_components
names_out<-colnames(X_south$Y)

# Estimate R2 of each latent variable by runing a linear model with only this component
results_R2_test <- calculate_R2(latent_variables, X_south$Y, names_out, name_components)

#Projections with validated model-----------------------------------------------
print("Obtaining projections")
proj_north <- rgcca_predict(rgcca_res, blocks_test = X_north, prediction_model = "lm")$projection
proj_south <- rgcca_predict(rgcca_res, blocks_test = X_south, prediction_model = "lm")$projection
proj_all <- rgcca_predict(rgcca_res, blocks_test = X_combined, prediction_model = "lm")$projection

projections<-list("Projections N" = proj_north,
                  "Projections S" = proj_south,
                  "Projections all" = proj_all)

write_rds(projections, file = "./results/RGCCA_methyl/projections.rds")
print("Projection saved")

#proj <- readRDS("./results/RGCCA_methyl/projections.rds")
 
print("Other approaches")
# Approaches--------------------------------------------------------------------
### Train N, test All-----------------------------------------------------------
set.seed(1999)
list_N_all<- apply_RGCCA(rgcca_mod = rgcca_res,
                         X_data = X_north,
                         X_data_pred = X_combined, 
                         connection = connection, 
                         response = 5, 
                         ncomp = ncomp, 
                         sparsity = sparsity, 
                         approach = "trainN_testall",
                         output_dir = "RGCCA_methyl")
print("Done")

### All-------------------------------------------------------------------------
list_all<- apply_RGCCA(X_data = X_combined,
                       X_data_pred = X_combined, 
                       connection = connection, 
                       response = 5, 
                       ncomp = ncomp, 
                       sparsity = sparsity, 
                       approach = "all",
                       output_dir = "RGCCA_methyl")
print("Done")

print("End of other approaches")

