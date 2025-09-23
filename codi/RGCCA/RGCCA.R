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

write_rds(list_X, file = "./results/RGCCA/list_X.rds")
# Step 1: tuning parameters using cross validation for number of components-----

# Sparsity
print("- Tuning of sparsity")
connection = 1 - diag(4)
min.sparsity <- 1 / sqrt(sapply(X_north, ncol))
sparsity_grid <- expand.grid(prot = seq(min.sparsity[1],1, length.out=25),
                             serum = seq(min.sparsity[2],1, length.out=25),
                             urine = seq(min.sparsity[3],1, length.out=25),
                             Y = 0)

set.seed(1999)
cv_sgcca_sparsity <- rgcca_cv(blocks=X_north,
                              response=4,
                              method="sgcca",
                              par_value= sparsity_grid,
                              ncomp=1,
                              par_type="sparsity",
                              n_cores=detectCores(),
                              connection = connection)

write_rds(cv_sgcca_sparsity, "./results/RGCCA/model/cv_sgcca_sparsity.rds")
cv_sgcca_sparsity<-readRDS("./results/RGCCA/model/cv_sgcca_sparsity.rds")
sparsity = cv_sgcca_sparsity$best_params
  
# Number of components  
print("- Tuning of number of components per blocks")
set.seed(1999)
cv_sgcca_ncomp  <- rgcca_cv(blocks = X_north,
                            response = 4,
                            method = "sgcca",
                            par_type = "ncomp",
                            sparsity =  cv_sgcca_sparsity$best_params,
                            par_value = expand.grid(prot = 1:5, serum =1:5, urine = 1:5, Y = 1:3),
                            n_cores = detectCores(),
                            connection = connection)

write_rds(cv_sgcca_ncomp, "./results/RGCCA/model/cv_sgcca_ncomp.rds")
cv_sgcca_ncomp<-readRDS("./results/RGCCA/model/cv_sgcca_ncomp.rds")
ncomp = cv_sgcca_ncomp$best_params


# Step 2: RGCCA with optimized parameters---------------------------------------
rgcca_res <- rgcca(cv_sgcca_ncomp)
write_rds(rgcca_res, "./results/RGCCA/model/rgcca_validation.rds")
rgcca_res <- readRDS("./results/RGCCA/model/rgcca_validation.rds")
summary(rgcca_res)

#Bootstrap
rgcca_res$a

set.seed(1999)
bootstrap <- rgcca_bootstrap(rgcca_res,n_boot = 1000, n_cores = detectCores())
rgcca_res <- bootstrap$rgcca
write_rds(rgcca_res, "./results/RGCCA/model/rgcca_final.rds")
rgcca_res <- readRDS("./results/RGCCA/model/rgcca_final.rds")

#LC plot
source("./codi/Functions/bootstrap_function.R")

pdf("./results/RGCCA/LC_bootstrap/bootstrap_prot.pdf")
bootstrap_figure(bootstrap$stats, block = "prot", color = "#5E0626")
dev.off()

pdf("./results/RGCCA/LC_bootstrap/bootstrap_serum.pdf")
bootstrap_figure(bootstrap$stats, block = "serum", color = "#A7A7A7")
dev.off()

pdf("./results/RGCCA/LC_bootstrap/bootstrap_urine.pdf")
bootstrap_figure(bootstrap$stats, block = "urine", color = "#3C6B66")
dev.off()

pdf("./results/RGCCA/LC_bootstrap/bootstrap_Y.pdf")
bootstrap_figure(bootstrap$stats, block = "Y")
dev.off()

# Step 3: Prediction ability on test set---------------------------------------- 
# With rgcca_predict, you can get standard indicators (RMSE, MAE, R2) using different predictor models
rgcca_predict_res <- rgcca_predict(rgcca_res,blocks_test=X_south,prediction_model="lm")
rgcca_pred_n <- rgcca_predict(rgcca_res,blocks_test=X_north,prediction_model="lm")
rgcca_pred_all <- rgcca_predict(rgcca_res,blocks_test=X_combined,prediction_model="lm")
rgcca_pred_all$metric
rgcca_predict_res$metric
# Projection of the variables on the latent components for the test set
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()
name_components <- c("Prot","Serum","Urine")
colnames(latent_variables)<- name_components
write.csv2(latent_variables, "./results/RGCCA/model_results/lat_vars.csv")
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

write_rds(results_list_validation, file = "./results/RGCCA/model_results/results_rgcca_validation.rds")



#Projections with validated model-----------------------------------------------
proj_north <- rgcca_predict(rgcca_res, blocks_test = X_north, prediction_model = "lm")$projection
proj_south <- rgcca_predict(rgcca_res, blocks_test = X_south, prediction_model = "lm")$projection
proj_all <- rgcca_predict(rgcca_res, blocks_test = X_combined, prediction_model = "lm")$projection

projections<-list("Projections N" = proj_north,
                  "Projections S" = proj_south,
                  "Projections all" = proj_all)

write_rds(projections, file = "./results/RGCCA/projections.rds")
proj <- readRDS("./results/RGCCA/projections.rds")


#EXTRACT outcome projections
#Validation
Y_validation <- scale(
  as.matrix(X_north$Y),
  center = attr(rgcca_res$blocks$Y, "scaled:center"),
  scale = attr(rgcca_res$blocks$Y, "scaled:scale"))

all(Y_validation %*% rgcca_res$a$Y==rgcca_res$Y$Y)

Y_s_scaled <-  scale(
  as.matrix(X_south$Y),
  center = attr(rgcca_res$blocks$Y, "scaled:center"),
  scale = attr(rgcca_res$blocks$Y, "scaled:scale"))

Y_all_scaled <-scale(
  as.matrix(X_combined$Y),
  center = attr(rgcca_res$blocks$Y, "scaled:center"),
  scale = attr(rgcca_res$blocks$Y, "scaled:scale"))

proj_outcome <- list("Proj N"=Y_validation %*% rgcca_res$a$Y,
                     "Proj S"=Y_s_scaled %*% rgcca_res$a$Y,
                      "Proj all"=Y_all_scaled %*% rgcca_res$a$Y)

write_rds(proj_outcome, file = "./results/RGCCA/projections_outcome.rds")

#EXTRACT outcome projections (only adol)
#Validation
Y_validation <- scale(
  as.matrix(X_north$Y[3:4]),
  center = attr(rgcca_res$blocks$Y, "scaled:center")[3:4],
  scale = attr(rgcca_res$blocks$Y, "scaled:scale")[3:4])

Y_s_scaled <-  scale(
  as.matrix(X_south$Y[3:4]),
  center = attr(rgcca_res$blocks$Y, "scaled:center")[3:4],
  scale = attr(rgcca_res$blocks$Y, "scaled:scale")[3:4])

Y_all_scaled <-scale(
  as.matrix(X_combined$Y[3:4]),
  center = attr(rgcca_res$blocks$Y, "scaled:center")[3:4],
  scale = attr(rgcca_res$blocks$Y, "scaled:scale")[3:4])

proj_outcome_adol <- list("Proj N"=Y_validation %*% rgcca_res$a$Y[3:4, drop=F],
                     "Proj S"=Y_s_scaled %*% rgcca_res$a$Y[3:4, drop=F],
                     "Proj all"=Y_all_scaled %*% rgcca_res$a$Y[3:4, drop=F])

write_rds(proj_outcome_adol, file = "./results/RGCCA/projections_outcome_adol.rds")

# Approaches--------------------------------------------------------------------
### Train N, test All-----------------------------------------------------------
set.seed(1999)
list_N_all<- apply_RGCCA(rgcca_mod = rgcca_res,
                         X_data = X_north,
                         X_data_pred = X_combined, 
                         connection = connection, 
                         response = 4, 
                         ncomp = ncomp, 
                         sparsity = sparsity, 
                         approach = "trainN_testall")

### All-------------------------------------------------------------------------
list_all<- apply_RGCCA(X_data = X_combined,
                       X_data_pred = X_combined, 
                       connection = connection, 
                       response = 4, 
                       ncomp = ncomp, 
                       sparsity = sparsity, 
                       approach = "all")

### North-----------------------------------------------------------------------
list_N<- apply_RGCCA(rgcca_mod = rgcca_res,
                     X_data = X_north,
                     X_data_pred = X_north, 
                     connection = connection, 
                     response = 4, 
                     ncomp = ncomp, 
                     sparsity = sparsity, 
                     approach = "North")


###train(S),test(N)-------------------------------------------------------------
list_SN <- apply_RGCCA(X_data = X_south,
                       X_data_pred = X_north, 
                       connection = connection, 
                       response = 4,            
                       ncomp = ncomp, 
                       sparsity = sparsity,            
                       approach = "trainS_testN")


#Performance RGCCA model--------------------------------------------------------
###North
performance_n <- performance_RGCCA_cv(rgcca_model = rgcca_res, 
                                      X_test = X_north,
                                      response = 4,
                                      sparsity = sparsity,
                                      ncomp = ncomp,
                                      connection = 1-diag(4))
r2_n <- bind_rows(lapply(names(performance_n), function(name) {performance_n[[name]]$r2 %>%
    mutate(name = name)}), .id = "source") %>% group_by(name)

r2_block_n <- bind_rows(lapply(names(performance_n), function(name) {performance_n[[name]]$r2_block %>%
    mutate(name = name)}), .id = "source") %>% group_by(Indicator, name)

###South
performance_s <- performance_RGCCA_cv(rgcca_model = rgcca_res, 
                                      X_test = X_south,
                                      response = 4,
                                      sparsity = sparsity,
                                      ncomp = ncomp,
                                      connection = 1-diag(4))

r2_s <- bind_rows(lapply(names(performance_s), function(name) {performance_s[[name]]$r2 %>%
    mutate(name = name)}), .id = "source") %>% group_by(name)

r2_block_s <- bind_rows(lapply(names(performance_s), function(name) {performance_s[[name]]$r2_block %>% 
    mutate(name = name)}), .id = "source") %>% group_by(Indicator, name)

###All
performance_a <- performance_RGCCA_cv(rgcca_model = rgcca_res, 
                                      X_test = X_combined,
                                      response = 4,
                                      sparsity = sparsity,
                                      ncomp = ncomp,
                                      connection = 1-diag(4))

r2_a <- bind_rows(lapply(names(performance_a), function(name) {performance_a[[name]]$r2 %>% 
    mutate(name = name)}), .id = "source") %>% group_by(name)

r2_block_a <- bind_rows(lapply(names(performance_a), function(name) {performance_a[[name]]$r2_block %>% 
    mutate(name = name)}), .id = "source") %>% group_by(Indicator, name)

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

# R2 per block (pooled population)
pdf("./results/RGCCA/R2/R2_block_all.pdf", width = 6, height = 5)
ggplot(data = r2_block_a, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (pooled population)") +
  theme(axis.text.x = element_text(angle = 0, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank()) +
  scale_x_discrete(labels = c("Prot-LC", "Serum-LC", "Urine-LC")) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  geom_text(data = r2_block_a, 
            aes(x = Indicator, y = md, label = round(md, 2)), 
            position = position_dodge(width = 0.75), 
            size = 3, 
            vjust = -0.5)
dev.off()

# R2 per block (North)
pdf("./results/RGCCA/R2/R2_block_n.pdf", width = 6, height = 5)
ggplot(data = r2_block_n, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (North)") +
  theme(axis.text.x = element_text(angle = 0, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_x_discrete(labels = c("Prot-LC", "Serum-LC", "Urine-LC")) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  geom_text(data = r2_block_n, 
            aes(x = Indicator, y = md, label = round(md, 2)), 
            position = position_dodge(width = 0.75), 
            size = 3, 
            vjust = -0.5)
dev.off()

# R2 per block (South)
pdf("./results/RGCCA/R2/R2_block_s.pdf", width = 6, height = 5)
ggplot(data = r2_block_s, aes(x = Indicator, y = value, fill = name)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Latent component", y = "R2 (%)", title = "R2 per block (South)") +
  theme(axis.text.x = element_text(angle = 0, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none") + 
  scale_x_discrete(labels = c("Prot-LC", "Serum-LC", "Urine-LC")) +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  geom_text(data = r2_block_s, 
            aes(x = Indicator, y = md, label = round(md, 2)), 
            position = position_dodge(width = 0.75), 
            size = 3, 
            vjust = -0.5)
dev.off()

# R2 per outcome (pooled population)
pdf("./results/RGCCA/R2/R2_outcome_all1.pdf", width = 3, height = 5)
ggplot(data = r2_a, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.position = "none") +  
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (pooled population)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  geom_text(data = r2_a, 
            aes(x = name, y = md, label = round(md, 2)), 
            position = position_dodge(width = 0.75), 
            size = 3, 
            vjust = -0.5)
dev.off()

# R2 per outcome (North)
pdf("./results/RGCCA/R2/R2_outcome_n.pdf", width = 3, height = 5)
ggplot(data = r2_n, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.position = "none") +  
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (North)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  geom_text(data = r2_n, 
            aes(x = name, y = md, label = round(md, 2)), 
            position = position_dodge(width = 0.75), 
            size = 3, 
            vjust = -0.5)
dev.off()

# R2 per outcome (South)
pdf("./results/RGCCA/R2/R2_outcome_s.pdf", width = 3, height = 5)
ggplot(data = r2_s, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.position = "none") +  
  labs(x = "", y = "R2 (%)", title = "R2 per outcome (South)") +
  scale_fill_manual(values = colors,
                    labels = c("res_dia1" = "Diastolic BP SD score in childhood", 
                               "res_dia2" = "Diastolic BP SD score in adolescence", 
                               "res_sys1" = "Systolic BP SD score in childhood", 
                               "res_sys2" = "Systolic BP SD score in adolescence")) +
  geom_text(data = r2_s, 
            aes(x = name, y = md, label = round(md, 2)), 
            position = position_dodge(width = 0.75), 
            size = 3, 
            vjust = -0.5)
dev.off()

