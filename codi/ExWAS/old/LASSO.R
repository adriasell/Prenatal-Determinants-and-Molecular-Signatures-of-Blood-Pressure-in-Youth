#LASSO (sharp)------------------------------------------------------------------
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
to_dummy <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
              "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")
#Scale exposures
exposome[, !colnames(exposome) %in% to_factor] <- scale(exposome[, !colnames(exposome) %in% to_factor], center = T)
ids<-rownames(exposome)

exposome <- dummy_cols(exposome,
                       select_columns = to_dummy,
                       remove_selected_columns = T,
                       remove_first_dummy = T)
rownames(exposome)<-ids
exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]

#NORTH
res_exp_n<-list()
for (m in 1:length(proj_north)){
  penalty <- c(rep(1,ncol(exposome_n)),rep(0,ncol(phenotype_n)))
  
  res <- VariableSelection(xdata = cbind(exposome_n, phenotype_n),
                           ydata = proj_north[[m]],
                           family = "gaussian", 
                           Lambda = LambdaSequence(lmax =  1e-1, lmin =  1e-3, cardinal = 100),
                           seed = c(12345),
                           penalty.factor = penalty,
                           tau=0.8, 
                           K = 500, 
                           n_cat=3 , 
                           pi_list=seq(0.51, 0.99, by = 0.01), 
                           beep=NULL, 
                           verbose=TRUE)
  
  res_exp_n[[m]] <- res
  names(res_exp_n)[m]<-names(proj_north)[m]
}   

#SOUTH
res_exp_s<-list()
for (m in 1:length(proj_south)){
  penalty <- c(rep(1,ncol(exposome_s)),rep(0,ncol(phenotype_s)))
  
  res <- VariableSelection(xdata = cbind(exposome_s, phenotype_s),
                           ydata = proj_south[[m]],
                           family = "gaussian", 
                           Lambda = LambdaSequence(lmax =  1e-1, lmin =  1e-3, cardinal = 100),
                           seed = c(12345),
                           penalty.factor = penalty,
                           tau=0.8, 
                           K = 500, 
                           n_cat=3 , 
                           pi_list=seq(0.51, 0.99, by = 0.01), 
                           beep=NULL, 
                           verbose=TRUE)
  
  res_exp_s[[m]] <- res
  names(res_exp_s)[m]<-names(proj_south)[m]
}   


#ALL
res_exp_all<-list()
for (m in 1:length(proj_all)){
  penalty <- c(rep(1,ncol(exposome_all)),rep(0,ncol(phenotype_all)))
  
  res <- VariableSelection(xdata = cbind(exposome_all, phenotype_all),
                           ydata = proj_all[[m]],
                           family = "gaussian", 
                           Lambda = LambdaSequence(lmax =  1e-1, lmin = 1e-3, cardinal = 100),
                           seed = c(12345), 
                           penalty.factor = penalty,
                           tau=0.8, 
                           K = 500, 
                           n_cat=3 , 
                           pi_list=seq(0.51, 0.99, by = 0.01), 
                           beep=NULL, 
                           verbose=TRUE)
  
  res_exp_all[[m]] <- res
  names(res_exp_all)[m]<-names(proj_all)[m]
}   

#RESULTS 

#north
for (m in 1:length(res_exp_n)){
  name<- names(res_exp_n[m])
  stab_m <- res_exp_n[[m]]
  class(stab_m) <- "variable_selection"
  
  selected_m <- SelectedVariables(stab_m)
  selected_names <- names(selected_m[selected_m==1])
  selprop_m <- SelectionProportions(stab_m)
  selected_selprop <- selprop_m[selected_m==1]
  pi_list=seq(0.51, 0.99, by = 0.01)
  argmax_id <- ArgmaxId(stab_m)[2]
  
  tmp <- t(stab_m$Beta[argmax_id, colnames(stab_m$selprop), ])
  beta_m <- apply(tmp, 2, FUN = function(x) {mean(x[x != 0])})
  selected_beta <- beta_m[selected_m==1]
  
  write.csv2(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/ExWAS/LASSO/Results/n_exp_",name,".csv")
  )
  svg(paste0("./results/ExWAS/LASSO/Calibration_Plots/Cal_Plot_N_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
  
}


#South
for (m in 1:length(res_exp_s)){
  name<- names(res_exp_s[m])
  stab_m <- res_exp_s[[m]]
  class(stab_m) <- "variable_selection"
  
  selected_m <- SelectedVariables(stab_m)
  selected_names <- names(selected_m[selected_m==1])
  selprop_m <- SelectionProportions(stab_m)
  selected_selprop <- selprop_m[selected_m==1]
  pi_list=seq(0.51, 0.99, by = 0.01)
  argmax_id <- ArgmaxId(stab_m)[2]
  
  tmp <- t(stab_m$Beta[argmax_id, colnames(stab_m$selprop), ])
  beta_m <- apply(tmp, 2, FUN = function(x) {mean(x[x != 0])})
  selected_beta <- beta_m[selected_m==1]
  
  write.csv2(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/ExWAS/LASSO/Results/s_exp_",name,".csv")
  )
  svg(paste0("./results/ExWAS/LASSO/Calibration_Plots/Cal_Plot_S_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
  
}

#all
for (m in 1:length(res_exp_all)){
  name<- names(res_exp_all[m])
  stab_m <- res_exp_all[[m]]
  class(stab_m) <- "variable_selection"
  
  selected_m <- SelectedVariables(stab_m)
  selected_names <- names(selected_m[selected_m==1])
  selprop_m <- SelectionProportions(stab_m)
  selected_selprop <- selprop_m[selected_m==1]
  pi_list=seq(0.51, 0.99, by = 0.01)
  argmax_id <- ArgmaxId(stab_m)[2]
  
  tmp <- t(stab_m$Beta[argmax_id, colnames(stab_m$selprop), ])
  beta_m <- apply(tmp, 2, FUN = function(x) {mean(x[x != 0])})
  selected_beta <- beta_m[selected_m==1]
  
  write.csv2(
    cbind("Names"=selected_names,
          "1_Selprop"=1-selected_selprop,
          "Beta"=selected_beta), row.names = FALSE,
    file = paste0("./results/ExWAS/LASSO/Results/all_exp_",name,".csv")
  )
  svg(paste0("./results/ExWAS/LASSO/Calibration_Plots/Cal_Plot_all_", name, ".svg"))
  CalibrationPlot(stab_m)
  dev.off()
}

##FOREST LASSO PLOT COMP-EXPOSURE-----------------------------------------------
all_exp_prot <- read.csv2("results/ExWAS/LASSO/Results/all_exp_prot.csv") %>% mutate(source = "all", comp = "prot")
all_exp_serum <- read.csv2("results/ExWAS/LASSO/Results/all_exp_serum.csv") %>% mutate(source = "all", comp = "serum")
all_exp_urine <- read.csv2("results/ExWAS/LASSO/Results/all_exp_urine.csv") %>% mutate(source = "all", comp = "urine")
n_exp_prot <- read.csv2("results/ExWAS/LASSO/Results/n_exp_prot.csv") %>% mutate(source = "N", comp = "prot")
n_exp_serum <- read.csv2("results/ExWAS/LASSO/Results/n_exp_serum.csv") %>% mutate(source = "N", comp = "serum")
n_exp_urine <- read.csv2("results/ExWAS/LASSO/Results/n_exp_urine.csv") %>% mutate(source = "N", comp = "urine")
s_exp_prot <- read.csv2("results/ExWAS/LASSO/Results/s_exp_prot.csv") %>% mutate(source = "S", comp = "prot")
s_exp_serum <- read.csv2("results/ExWAS/LASSO/Results/s_exp_serum.csv") %>% mutate(source = "S", comp = "serum")
s_exp_urine <- read.csv2("results/ExWAS/LASSO/Results/s_exp_urine.csv") %>% mutate(source = "S", comp = "urine")

df_lasso <- rbind(
  all_exp_prot, 
  all_exp_serum, 
  all_exp_urine, 
  n_exp_prot, 
  n_exp_serum, 
  n_exp_urine, 
  s_exp_prot, 
  s_exp_serum, 
  s_exp_urine)

df_lasso$Beta <- as.numeric(df_lasso$Beta)
df_lasso <- inner_join(df_lasso, codebook[, 8:9], by = c("Names" = "Variable_name_TRANS"))
write.csv2(df_lasso, "./results/ExWAS/LASSO_comp_exposome.csv")

#PLOT 95%????? COM FERHO??
#PROT
png("./results/ExWAS/Forest_plot/LASSO_prot.png", width = 1300, height = 900)
ggplot(df_lasso[df_lasso$comp=="prot",],
       aes(y = reorder(Label, Beta), x = Beta, xmin = Beta, xmax = Beta, shape = source, colour = source)) +
  scale_color_manual(values = colors) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "LASSO Beta exposure-prot component") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#SERUM
png("./results/ExWAS/Forest_plot/LASSO_serum.png", width = 1300, height = 900)
ggplot(df_lasso[df_lasso$comp=="serum",],
       aes(y = reorder(Label, Beta), x = Beta, xmin = Beta, xmax = Beta, shape = source, colour = source)) +
  scale_color_manual(values = colors) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "LASSO Beta exposure-serum component") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#URINE
png("./results/ExWAS/Forest_plot/LASSO_urine.png", width = 1300, height = 900)
ggplot(df_lasso[df_lasso$comp=="urine",],
       aes(y = reorder(Label, Beta), x = Beta, xmin = Beta, xmax = Beta, shape = source, colour = source)) +
  scale_color_manual(values = colors) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "LASSO Beta exposure-urine component") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()