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

load("./results/RGCCA/data_X.RData")
X_train <- data$X_train
X_test <- data$X_test

#OPCIÓ 1 SPARSITY (ncomp=1:5, y(1:3))-------------------------------------------

### Step 1: tuning parameters using cross validation for number of components  #

# Sparsity (OPCIONAL)
print("- Tuning of sparsity")
connection = 1 - diag(4)

min_sparsity <- 1 / sqrt(sapply(X_train, ncol))
sparsity_grid <- expand.grid(prot = seq(min_sparsity[1],1, length.out=10),
                             serum = seq(min_sparsity[2],1, length.out=10),
                             urine = seq(min_sparsity[3],1, length.out=10),
                             Y = seq(min_sparsity[4],1,length.out=10))

cv_sgcca_sparsity <- rgcca_cv(blocks=X_train,
                              response=4,
                              method="sgcca",
                              par_value= sparsity_grid,
                              ncomp=1,
                              par_type="sparsity",
                              n_cores=detectCores(),
                              connection = connection)

sparsity = cv_sgcca_sparsity$best_params
write_rds(cv_sgcca_sparsity, "./results/RGCCA/cv_sgcca_sparsity.rds")
cv_sgcca_sparsity<-readRDS("./results/RGCCA/cv_sgcca_sparsity.rds")

# Number of components  
print("- Tuning of number of components per blocks")
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
write_rds(rgcca_res, "./results/RGCCA/rgcca_m1.rds")
summary(rgcca_res)
stab<-rgcca_stability(rgcca_res, n_cores = detectCores())
plot(stab)

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

m1_list<-list("sparsity"=cv_sgcca_sparsity$best_params,
     "ncomp"=cv_sgcca_ncomp$best_params,
     "R2"=results_R2,
     "metrics"=rgcca_predict_res$metric,
     "corr"=cor(latent_variables) %>% round(2),
     "Latent factor"=rgcca_res$Y,
     "contributions"=rgcca_res$a,
     "projection"=rgcca_predict_res$projection)

save(m1_list, file = "./results/RGCCA/rgcca_m1.RData")
m1_list[7]
