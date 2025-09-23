library(tidyverse)
library(RGCCA)
library(Biobase)
library(dplyr)
library(parallel)
library(writexl)
library(fastDummies)
#install.packages("./packages_R/RGCCA-main", repos = NULL, type = "source") RGCCA modified

getwd()
#setwd("./TFM")

ncore.use <- parallel::detectCores() - 1
set.seed(1999)

source("./codi/Functions/functions_RGCCA.R")

# Step 0: prepare the dataset---------------------------------------------------

# /// Load HelixID by group (N/S)
helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

# Omic data
### Prot
prot_all <- eset2df("./results/denoising/denoising_prot/prot_denoised.RDS")
prot_n <- prot_all[rownames(prot_all) %in% helixid_n,]
prot_s <- prot_all[rownames(prot_all) %in% helixid_s,]

### Serum
serum_all <- eset2df("./results/denoising/denoising_metab/serum/metab_serum_denoised.RDS")
serum_n <- serum_all[rownames(serum_all) %in% helixid_n,]
serum_s <- serum_all[rownames(serum_all) %in% helixid_s,]

### Urine
urine_all <- eset2df("./results/denoising/denoising_metab/urine/metab_urine_denoised.RDS")
urine_n <- urine_all[rownames(urine_all) %in% helixid_n,]
urine_s <- urine_all[rownames(urine_all) %in% helixid_s,]

### Exposome
exposome_all <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
exposome_all$HelixID <- rownames(exposome_all)
to_dummy <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
              "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")
#Scale exposures
exposome_all[, !colnames(exposome_all) %in% c(to_dummy, "HelixID", "e3_asmokyn_p_None", "e3_alcpreg_yn_None")] <- scale(exposome_all[, !colnames(exposome_all) %in% c(to_dummy, "HelixID", "e3_asmokyn_p_None", "e3_alcpreg_yn_None", "h_greenyn300_preg_None")], center = T)

exposome_all <- dummy_cols(exposome_all,
                       select_columns = to_dummy,
                       remove_selected_columns = T,
                       remove_first_dummy = T)

rownames(exposome_all)<-exposome_all$HelixID
exposome_all <- exposome_all[rownames(exposome_all) %in% rownames(prot_all),]
exposome_all <- exposome_all[rownames(prot_all),]
exposome_n <- exposome_all[rownames(exposome_all) %in% rownames(prot_n),]
exposome_s <- exposome_all[rownames(exposome_all) %in% rownames(prot_s),]


#Comprovations
all(rownames(prot_s)==rownames(serum_s))
all(rownames(serum_s)==rownames(urine_s))
all(rownames(exposome_s)==rownames(serum_s))
all(rownames(prot_n)==rownames(serum_n))
all(rownames(serum_n)==rownames(urine_n))
all(rownames(exposome_n)==rownames(serum_n))
all(rownames(prot_all)==rownames(serum_all))
all(rownames(urine_all)==rownames(serum_all))
all(rownames(exposome_all)==rownames(serum_all))

# Outcome
Y <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(Y) <- Y$HelixID
outcomes <- c("hs2_zdia_bp.v3_2017_Time1", "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time2", "hs2_zsys_bp.v3_2017_Time2","HelixID")
names(Y)[grep("hs2_zsys_bp_v3_2017", names(Y))] <- "hs2_zsys_bp_v3_2017_Time2"
names(Y)[grep("hs2_zdia_bp_v3_2017", names(Y))] <- "hs2_zdia_bp_v3_2017_Time2"
Y <- Y[outcomes]
Y <- Y[complete.cases(Y),]

#Select complete cases
ids_n <- Reduce(intersect, list(prot_n$HelixID, serum_n$HelixID, urine_n$HelixID, exposome_n$HelixID, Y$HelixID))
ids_s <- Reduce(intersect, list(prot_s$HelixID, serum_s$HelixID, urine_s$HelixID, exposome_s$HelixID, Y$HelixID))

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
exposome_n <- exposome_n %>% arrange(HelixID) %>% dplyr::select(-HelixID)
exposome_s <- exposome_s %>% arrange(HelixID) %>% dplyr::select(-HelixID)
exposome_all <- exposome_all %>% arrange(HelixID) %>% dplyr::select(-HelixID)
Y <- Y %>% 
  arrange(HelixID) %>% dplyr::select(-HelixID) %>%
  mutate(across(everything(), ~ gsub(",", ".", .))) %>%
  mutate(across(everything(), as.numeric))

# Divide it in train(North) test(South)
X_north <- list(prot = prot_n[rownames(prot_n) %in% ids_n,],
                serum= serum_n[rownames(serum_n) %in% ids_n,],
                urine= urine_n[rownames(urine_n) %in% ids_n,],
                exposome= exposome_n[rownames(exposome_n) %in% ids_n,],
                Y=Y[rownames(Y) %in% ids_n,])

X_south <- list(prot = prot_s[rownames(prot_s) %in% ids_s,],
               serum= serum_s[rownames(serum_s) %in% ids_s,],
               urine= urine_s[rownames(urine_s) %in% ids_s,],
               exposome= exposome_s[rownames(exposome_s) %in% ids_s,],
               Y=Y[rownames(Y) %in% ids_s,])

X_combined <- list(prot = prot_all[rownames(prot_all) %in% c(ids_s, ids_n),],
                   serum = serum_all[rownames(serum_all) %in% c(ids_s, ids_n),],
                   urine = urine_all[rownames(urine_all) %in% c(ids_s, ids_n),],
                   exposome= exposome_all[rownames(exposome_all) %in% c(ids_s, ids_n),],
                   Y = Y[rownames(Y) %in% c(ids_s, ids_n),])

list_X <- list("X_north"=X_north,
               "X_south"=X_south,
               "X_combined"=X_combined)

write_rds(list_X, file = "./results/RGCCA_exposome/list_X.rds")
# Step 1: tuning parameters using cross validation for number of components-----

# Sparsity
print("- Tuning of sparsity")
connection = 1 - diag(5)
min.sparsity <- 1 / sqrt(sapply(X_north, ncol))
sparsity_grid <- expand.grid(prot = seq(min.sparsity[1],1, length.out=15),
                             serum = seq(min.sparsity[2],1, length.out=15),
                             urine = seq(min.sparsity[3],1, length.out=15),
                             exposome = seq(min.sparsity[4],1, length.out=15),
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

write_rds(cv_sgcca_sparsity, "./results/RGCCA_exposome/model/cv_sgcca_sparsity.rds")
cv_sgcca_sparsity<-readRDS("./results/RGCCA_exposome/model/cv_sgcca_sparsity.rds")
sparsity = cv_sgcca_sparsity$best_params
  
# Number of components  
print("- Tuning of number of components per blocks")
set.seed(1999)
cv_sgcca_ncomp  <- rgcca_cv(blocks = X_north,
                            response = 5,
                            method = "sgcca",
                            par_type = "ncomp",
                            sparsity =  cv_sgcca_sparsity$best_params,
                            par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, exposome = 1:5, Y = 1:3),
                            n_cores = detectCores(),
                            connection = connection)

write_rds(cv_sgcca_ncomp, "./results/RGCCA_exposome/model/cv_sgcca_ncomp.rds")
cv_sgcca_ncomp<-readRDS("./results/RGCCA_exposome/model/cv_sgcca_ncomp.rds")
ncomp = cv_sgcca_ncomp$best_params


# Step 2: RGCCA with optimized parameters---------------------------------------
rgcca_res <- rgcca(cv_sgcca_ncomp)
write_rds(rgcca_res, "./results/RGCCA_exposome/model/rgcca_validation.rds")
rgcca_res <- readRDS("./results/RGCCA_exposome/model/rgcca_validation.rds")
summary(rgcca_res)

#Bootstrap
rgcca_res$a

set.seed(1999)
bootstrap <- rgcca_bootstrap(rgcca_res,n_boot = 5000, n_cores = detectCores())
rgcca_res <- bootstrap$rgcca
write_rds(rgcca_res, "./results/RGCCA_exposome/model/rgcca_final.rds")
rgcca_res <- readRDS("./results/RGCCA_exposome/model/rgcca_final.rds")

#LC plot
source("./codi/Functions/bootstrap_function.R")
png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_prot.png")
bootstrap_figure(bootstrap$stats, block = "prot", color = "#E05600")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_serum.png")
bootstrap_figure(bootstrap$stats, block = "serum", color = "#10450F")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_urine.png")
bootstrap_figure(bootstrap$stats, block = "urine", color = "#064550")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_urine2.png")
bootstrap_figure(bootstrap$stats, block = "urine", color = "#064550", comp = 2 )
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_exposome.png")
bootstrap_figure(bootstrap$stats, block = "exposome", color = "#C1C1C1" )
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_exposome2.png")
bootstrap_figure(bootstrap$stats, block = "exposome", color = "#C1C1C1", comp = 2 )
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_Y.png")
bootstrap_figure(bootstrap$stats, block = "Y")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap/bootstrap_Y2.png")
bootstrap_figure(bootstrap$stats, block = "Y", comp = 2)
dev.off()

# Step 3: Prediction ability on test set---------------------------------------- 
# With rgcca_predict, you can get standard indicators (RMSE, MAE, R2) using different predictor models
rgcca_predict_res <- rgcca_predict(rgcca_res,blocks_test=X_south,prediction_model="lm")
rgcca_pred_n <- rgcca_predict(rgcca_res,blocks_test=X_north,prediction_model="lm")
rgcca_pred_all <- rgcca_predict(rgcca_res,blocks_test=X_combined,prediction_model="lm")

# Projection of the variables on the latent components for the test set
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()
name_components <- c("Prot","Serum","Urine1", "Urine2","Exposome1", "Exposome2")
colnames(latent_variables)<- name_components
names_out<-colnames(X_south$Y)

# Projection of the variables on the latent components for the train set
latent_variables_n <- rgcca_pred_n$projection %>% purrr::reduce(cbind)%>% as.data.frame()
colnames(latent_variables_n)<- name_components

# Estimate R2 of each latent variable by runing a linear model with only this component
results_R2_test <- calculate_R2(latent_variables, X_south$Y, names_out, name_components)
results_R2_train <- calculate_R2(latent_variables_n, X_north$Y, names_out, name_components)

# Estimate correlation between components
cor(latent_variables) %>% round(2)

# Result list
results_list_validation <-  list("sparsity"=cv_sgcca_sparsity$best_params,
                                 "ncomp"=cv_sgcca_ncomp$best_params,
                                 "R2 test (south)"=results_R2_test, 
                                 "R2 train (north)"=results_R2_train, 
                                 "metrics"=rgcca_predict_res$metric,
                                 "corr"=cor(latent_variables) %>% round(2),
                                 "Latent factor"=rgcca_res$Y,
                                 "Block weights for block"=rgcca_res$a,
                                 "projection"=rgcca_predict_res$projection)

write_rds(results_list_validation, file = "./results/RGCCA_exposome/model_results/results_rgcca_validation.rds")



#Projections with validated model-----------------------------------------------
proj_north <- rgcca_predict(rgcca_res, blocks_test = X_north, prediction_model = "lm")$projection
proj_south <- rgcca_predict(rgcca_res, blocks_test = X_south, prediction_model = "lm")$projection
proj_all <- rgcca_predict(rgcca_res, blocks_test = X_combined, prediction_model = "lm")$projection

projections<-list("Projections N" = proj_north,
                  "Projections S" = proj_south,
                  "Projections all" = proj_all)

write_rds(projections, file = "./results/RGCCA_exposome/projections.rds")
proj <- readRDS("./results/RGCCA_exposome/projections.rds")

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
                         output_dir = "RGCCA_exposome")

### All-------------------------------------------------------------------------
list_all<- apply_RGCCA(X_data = X_combined,
                       X_data_pred = X_combined, 
                       connection = connection, 
                       response = 5, 
                       ncomp = ncomp, 
                       sparsity = sparsity, 
                       approach = "all",
                       output_dir = "RGCCA_exposome")

### North-----------------------------------------------------------------------
list_N<- apply_RGCCA(rgcca_mod = rgcca_res,
                     X_data = X_north,
                     X_data_pred = X_north, 
                     connection = connection, 
                     response = 5, 
                     ncomp = ncomp, 
                     sparsity = sparsity, 
                     approach = "North",
                     output_dir = "RGCCA_exposome")


###train(S),test(N)-------------------------------------------------------------
list_SN <- apply_RGCCA(X_data = X_south,
                       X_data_pred = X_north, 
                       connection = connection, 
                       response = 5,            
                       ncomp = ncomp, 
                       sparsity = sparsity,            
                       approach = "trainS_testN",
                       output_dir = "RGCCA_exposome")


#Performance RGCCA model--------------------------------------------------------
###North
performance_n <- performance_RGCCA_cv(rgcca_model = rgcca_res, 
                                      X_test = X_north,
                                      response = 5,
                                      sparsity = sparsity,
                                      ncomp = ncomp,
                                      connection = 1-diag(5))
r2_n <- bind_rows(lapply(names(performance_n), function(name) {performance_n[[name]]$r2 %>%
    mutate(name = name)}), .id = "source")

r2_block_n <- bind_rows(lapply(names(performance_n), function(name) {performance_n[[name]]$r2_block %>%
    mutate(name = name)}), .id = "source")

###South
performance_s <- performance_RGCCA_cv(rgcca_model = rgcca_res, 
                                      X_test = X_south,
                                      response = 5,
                                      sparsity = sparsity,
                                      ncomp = ncomp,
                                      connection = 1-diag(5))

r2_s <- bind_rows(lapply(names(performance_s), function(name) {performance_s[[name]]$r2 %>%
    mutate(name = name)}), .id = "source")

r2_block_s <- bind_rows(lapply(names(performance_s), function(name) {performance_s[[name]]$r2_block %>% 
    mutate(name = name)}), .id = "source")

###All
performance_a <- performance_RGCCA_cv(rgcca_model = rgcca_res, 
                                      X_test = X_combined,
                                      response = 5,
                                      sparsity = sparsity,
                                      ncomp = ncomp,
                                      connection = 1-diag(5))

r2_a <- bind_rows(lapply(names(performance_a), function(name) {performance_a[[name]]$r2 %>% 
    mutate(name = name)}), .id = "source")

r2_block_a <- bind_rows(lapply(names(performance_a), function(name) {performance_a[[name]]$r2_block %>% 
    mutate(name = name)}), .id = "source")

#Funció r2 a %
to_perc <- function(df) {
  df <- df %>% mutate(across(c(value, md), ~ . * 100, .names = "{.col}"))
  return(df)}

r2_n <- to_perc(r2_n)
r2_block_n <- to_perc(r2_block_n)

r2_s <- to_perc(r2_s)
r2_block_s <- to_perc(r2_block_s)

r2_a <- to_perc(r2_a)
r2_block_a <- to_perc(r2_block_a)




#Plots--------------------------------------------------------------------------
colors <- c("res_dia1" = "#8CC5E3", "res_dia2" = "#3594CC", "res_sys1" = "#F0B077", "res_sys2" = "#EA801C")

#R2 block
png("./results/RGCCA_exposome/R2/R2_block_all.png", width = 1500, height = 588)
ggplot(data = r2_block_a, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (pooled population)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.title = element_blank()) +
  scale_x_discrete(labels = name_components) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()

png("./results/RGCCA_exposome/R2/R2_block_n.png", width = 1500, height = 588)
ggplot(data = r2_block_n, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (North)") +
  theme(axis.text.x = element_text(angle = 0),legend.title = element_blank()) +
  scale_x_discrete(labels = name_components) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()

png("./results/RGCCA_exposome/R2/R2_block_s.png", width = 1500, height = 588)
ggplot(data = r2_block_s, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (South)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.title = element_blank()) +
  scale_x_discrete(labels = name_components) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()

#R2 
png("./results/RGCCA_exposome/R2/R2_outcome_all.png", width = 1500, height = 588)
ggplot(data = r2_a, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (pooled population)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()

png("./results/RGCCA_exposome/R2/R2_outcome_n.png", width = 1500, height = 588)
ggplot(data = r2_n, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (North)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30) 
dev.off()

png("./results/RGCCA_exposome/R2/R2_outcome_s.png", width = 1500, height = 588)
ggplot(data = r2_s, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (South)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()








#_______________________________________________________________________________
set.seed(123)
train_idx <- sample(rownames(Y), size = 0.8 * nrow(Y))
X_train <- lapply(X_combined, function(x) x[train_idx, , drop = FALSE])
X_test <- lapply(X_combined, function(x) x[setdiff(rownames(Y), train_idx), , drop = FALSE])

connection = 1 - diag(5)
min.sparsity <- 1 / sqrt(sapply(X_north, ncol))
sparsity_grid <- expand.grid(prot = seq(min.sparsity[1],1, length.out=10),
                             serum = seq(min.sparsity[2],1, length.out=10),
                             urine = seq(min.sparsity[3],1, length.out=10),
                             exposome = seq(min.sparsity[4],1, length.out=10),
                             Y = 0)

cv_sgcca_sparsity <- rgcca_cv(blocks=X_train,
                              response=5,
                              method="sgcca",
                              par_value= sparsity_grid,
                              ncomp=1,
                              par_type="sparsity",
                              n_cores=detectCores(),
                              connection = connection)

cv_sgcca_ncomp  <- rgcca_cv(blocks = X_train,
                            response = 5,
                            method = "sgcca",
                            par_type = "ncomp",
                            sparsity =  cv_sgcca_sparsity$best_params,
                            par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, exposome = 1:5, Y = 1:3),
                            n_cores = detectCores(),
                            connection = connection)

write_rds(cv_sgcca_ncomp, "./results/RGCCA_exposome/model/cv_rgcca_all_v1.rds")
rgcca_res_all <- rgcca(cv_sgcca_ncomp)
bootstrap <- rgcca_bootstrap(rgcca_res_all,n_boot = 5000, n_cores = detectCores())
rgcca_res_all <- bootstrap$rgcca
write_rds(rgcca_res_all, "./results/RGCCA_exposome/model/rgcca_all_v1.rds")
summary(rgcca_res_all)
rgcca_res_all$a

#WAITING
rgcca_predict_res <- rgcca_predict(rgcca_res_all,blocks_test=X_test ,prediction_model="lm")
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()

####!!!!!! depen ncomp
name_components <- paste0(rep(names(X_test)[-5],times=cv_sgcca_ncomp$best_params[-5]),unlist(lapply(cv_sgcca_ncomp$best_params[-5],seq)))


colnames(latent_variables)<- name_components
names_out<-colnames(X_train$Y)
R2_results_allv1 <- calculate_R2(latent_variables, X_test$Y, names_out, name_components)
rgcca_predict_res$metric
R2_results_allv1

#LC plot
source("./codi/Functions/bootstrap_function.R")
png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_prot.png")
bootstrap_figure(bootstrap$stats, block = "prot", color = "#E05600")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_serum.png")
bootstrap_figure(bootstrap$stats, block = "serum", color = "#10450F")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_urine.png")
bootstrap_figure(bootstrap$stats, block = "urine", color = "#064550")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_urine2.png")
bootstrap_figure(bootstrap$stats, block = "urine", color = "#064550", comp = 2 )
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_exposome.png")
bootstrap_figure(bootstrap$stats, block = "exposome", color = "#C1C1C1" )
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_exposome2.png")
bootstrap_figure(bootstrap$stats, block = "exposome", color = "#C1C1C1", comp = 2 )
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_Y.png")
bootstrap_figure(bootstrap$stats, block = "Y")
dev.off()

png("./results/RGCCA_exposome/LC_bootstrap_all_v1/bootstrap_Y2.png")
bootstrap_figure(bootstrap$stats, block = "Y", comp = 2)
dev.off()


performance_a_v1 <- performance_RGCCA_cv(rgcca_model = rgcca_res_all, 
                                      X_test = X_train,
                                      response = 5,
                                      sparsity = cv_sgcca_sparsity$best_params,
                                      ncomp = cv_sgcca_ncomp$best_params,
                                      connection = 1-diag(5))
rgcca_res_all$call$ncomp
r2_a_v1 <- bind_rows(lapply(names(performance_a_v1), function(name) {performance_a_v1[[name]]$r2 %>% 
    mutate(name = name)}), .id = "source")

r2_block_a_v1 <- bind_rows(lapply(names(performance_a_v1), function(name) {performance_a_v1[[name]]$r2_block %>% 
    mutate(name = name)}), .id = "source")

r2_a_v1 <- to_perc(r2_a_v1)
r2_block_a_v1 <- to_perc(r2_block_a_v1)

#Plots--------------------------------------------------------------------------
colors <- c("res_dia1" = "#8CC5E3", "res_dia2" = "#3594CC", "res_sys1" = "#F0B077", "res_sys2" = "#EA801C")

#R2 block
png("./results/RGCCA_exposome/R2/R2_block_all.png", width = 1500, height = 588)
ggplot(data = r2_block_a_v1, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (pooled population)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.title = element_blank()) +
  scale_x_discrete(labels = name_components) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()

#R2 
png("./results/RGCCA_exposome/R2/R2_outcome_all.png", width = 1500, height = 588)
ggplot(data = r2_a_v1, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (pooled population)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  ylim(0,30)
dev.off()

