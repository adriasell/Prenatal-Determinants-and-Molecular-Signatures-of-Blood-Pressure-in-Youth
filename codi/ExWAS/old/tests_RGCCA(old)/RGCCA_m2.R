library(tidyverse)
library(RGCCA)
library(Biobase)
library(dplyr)
library(parallel)
#devtools::install_github(repo="https://github.com/rgcca-factory/RGCCA.git")
getwd()
setwd("./TFM")
set.seed(1899)

ncore.use <- parallel::detectCores() - 1

### Step 0: prepare the dataset #####
load("./results/RGCCA/data_X.RData")
X_train <- data$X_train
X_test <- data$X_test


#OPCIÓ 2 SPARSITY (ncomp=1:5, y(2:3))-------------------------------------------

### Step 1: tuning parameters using cross validation for number of components  #
cv_sgcca_sparsity<-readRDS("./results/RGCCA/cv_sgcca_sparsity.rds")

# Number of components  
print("- Tuning of number of components per blocks")

connection = 1 - diag(4)
cv_sgcca_ncomp  <- rgcca_cv(blocks = X_train,
                            response = 4,
                            method = "sgcca",
                            connection = connection,
                            par_type = "ncomp",
                            sparsity = cv_sgcca_sparsity$best_params,
                            par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, Y = 2:3),
                            n_cores = detectCores())

plot(cv_sgcca_ncomp)
ncomp = cv_sgcca_ncomp$best_params


## Step 2: RGCCA with optimized parameters ##

rgcca_res <- rgcca(cv_sgcca_ncomp)
summary(rgcca_res)
write_rds(rgcca_res, "./results/RGCCA/rgcca_m2.rds")
stab<-rgcca_stability(rgcca_res, n_cores = detectCores())

## Step 3: Interpretation of the model and performance ##

# Significance: bootstrap (Time estimated: 2 min)
bootstrap <- rgcca_bootstrap(rgcca_res,n_cores = detectCores())
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

m2_list<-list("sparsity"=cv_sgcca_sparsity$best_params,
              "ncomp"=cv_sgcca_ncomp$best_params,
              "R2"=results_R2,
              "metrics"=rgcca_predict_res$metric,
              "corr"=cor(latent_variables) %>% round(2),
              "Latent factor"=rgcca_res$Y,
              "Contributions"=rgcca_res$a,
              "projection"=rgcca_predict_res$projection)

save(m2_list, file = "./results/RGCCA/rgcca_m2.RData")

m2_list[1:5]

