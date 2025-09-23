################################# Set R environment 
library(writexl)
library(tidyverse)
library(doParallel)
library(parallel)
library(ggsankey)
library(networkD3)
library(sharp)
library(readxl)
library(dplyr)
library(tidyr)
library(fastDummies)
library(htmlwidgets)
library(webshot)
library(patchwork)
library(stargazer)
library(broom)
#install.packages("./packages_R/ggsankey-main", repos = NULL, type = "source")

################################# Defining working directory.
getwd()
setwd("./TFM")
output_path <-"./results/ExWAS"
set.seed(1899)

source("./codi/Functions/functions_ExWAS.R")
source("./codi/Functions/function_volcano_plot.R")

###Prepare the dataset----------------------------------------------------------
#Load data
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined

#Load projections
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N`
proj_south <- projections$`Projections S`
proj_all <-  projections$`Projections all`

#phenotype data
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)

rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all$prot),]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north$prot),]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south$prot),]

#Codebook
codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")

#Load exposure data
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)

#Scale exposures
to_factor <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
               "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter","e3_asmokyn_p_None",
               "e3_alcpreg_yn_None")

exposome[, !colnames(exposome) %in% to_factor] <- scale(exposome[, !colnames(exposome) %in% to_factor], center = T)
exposome[to_factor] <- lapply(exposome[to_factor], factor)

exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]

#Prepare df
lat_vars_n <- prepare_df(projections = proj_north, X_data = X_north, exposure_data = exposome_n, response = 4)
lat_vars_s <- prepare_df(projections =  proj_south, X_data = X_south, exposure_data = exposome_s, response = 4)
lat_vars_all <- prepare_df(projections =  proj_all, X_data = X_combined, exposure_data = exposome_all, response = 4)

all(rownames(lat_vars_n)==rownames(phenotype_n))
all(rownames(lat_vars_s)==rownames(phenotype_s))
all(rownames(lat_vars_all)==rownames(phenotype_all))

# Convert the specified columns to factors in all three data frames

cols_to_factor <- c("h_parity_None", "h_native_None", "h_edumc_None", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")
# Apply factor conversion only if the columns exist in each dataframe
for (col in cols_to_factor) {
  if (col %in% names(phenotype_n)) phenotype_n[[col]] <- as.factor(phenotype_n[[col]])
  if (col %in% names(phenotype_s)) phenotype_s[[col]] <- as.factor(phenotype_s[[col]])
  if (col %in% names(phenotype_all)) phenotype_all[[col]] <- as.factor(phenotype_all[[col]])
}

#ExWAS--------------------------------------------------------------------------
rgcca_validation <- readRDS("./results/RGCCA/model/rgcca_final.rds")
response <- rgcca_validation$call$response
ncomp <- rgcca_validation$call$ncomp

## Association between components and the exposures------------------------------
name_components <- paste0(rep(names(X_combined)[-response],times=ncomp[-response]))

##north
res_exposure_N <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_exposure_N[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n),
                                     expos_name = colnames(exposome_n),
                                     form=form) %>%  mutate(comp=comp)}

res_exposure_N <- suppressMessages(purrr::reduce(res_exposure_N, full_join)) %>% arrange(p)

##south
exposome_s <- exposome_s[, colSums(exposome_s != 0) > 0]
res_exposure_S <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_INMA",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_exposure_S[[comp]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s),
                                       expos_name = colnames(exposome_s),
                                       form=form) %>%  mutate(comp=comp)}

res_exposure_S <- suppressMessages(purrr::reduce(res_exposure_S, full_join)) %>% arrange(p)

##all
res_exposure_all <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC", 
                                   "h_cohort_MOBA", 
                                   "h_cohort_INMA",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_exposure_all[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all),
                                          expos_name = colnames(exposome_all),
                                          form=form) %>%  mutate(comp=comp)}

res_exposure_all <- suppressMessages(purrr::reduce(res_exposure_all, full_join)) %>% arrange(p)

#Save results
list_res_exposure <- list("res_exposure_N"=res_exposure_N,
                         "res_exposure_S"=res_exposure_S,
                         "res_exposure_all"=res_exposure_all)

write_xlsx(list_res_exposure, paste0(output_path,"/ExWAS_comp_exposome.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_comp_exposome.xlsx"))
list_res_exposure <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_comp_exposome.xlsx"))
names(list_res_exposure)<-sheets

## Association between components and outcomes----------------------------------
##north
res_outcome_N <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c("h_cohort_BIB", 
                                      "h_cohort_EDEN", 
                                      "h_cohort_KANC",
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_outcome_N[[outcome]]  <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n),
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_N <- suppressMessages(purrr::reduce(res_outcome_N, full_join)) %>% arrange(p)

##south
res_outcome_S <- list()
for (outcome in colnames(X_south$Y)){
  form <- paste(outcome, "~ h_cohort_INMA + h_parity_None + h_native_None + h_edumc_None + h_age_None")
  
  res_outcome_S[[outcome]]  <-ExWAS_mixed(data=cbind(lat_vars_s, phenotype_s),
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_S <- suppressMessages(purrr::reduce(res_outcome_S, full_join)) %>% arrange(p)

##all
res_outcome_all  <- list()
for (outcome in colnames(X_combined$Y)){
  form <- paste(outcome,  "~", paste(c("h_cohort_BIB", 
                                       "h_cohort_EDEN", 
                                       "h_cohort_KANC", 
                                       "h_cohort_MOBA", 
                                       "h_cohort_INMA",
                                       "h_parity_None", 
                                       "h_native_None", 
                                       "h_edumc_None", 
                                       "h_age_None"), collapse = " + "))
  
  res_outcome_all[[outcome]]  <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all),
                                            expos_name = name_components,
                                            form=form) %>% mutate(outcome=outcome)}

res_outcome_all <- suppressMessages(purrr::reduce(res_outcome_all, full_join)) %>% arrange(p)

#Save results
list_res_outcome <- list("res_outcome_N"=res_outcome_N,
                         "res_outcome_S"=res_outcome_S,
                         "res_outcome_all"=res_outcome_all)

write_xlsx(list_res_outcome, paste0(output_path,"/ExWAS_comp_outcome.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_comp_outcome.xlsx"))
list_res_outcome <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_comp_outcome.xlsx"))
names(list_res_outcome)<-sheets

## Association between outcomes and the exposures------------------------------

##north
res_exposure_N <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c("h_cohort_BIB", 
                                      "h_cohort_EDEN", 
                                      "h_cohort_KANC",
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_exposure_N[[outcome]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n),
                                       expos_name = colnames(exposome_n),
                                       form=form) %>%  mutate(outcome=outcome)}

res_exposure_N <- suppressMessages(purrr::reduce(res_exposure_N, full_join)) %>% arrange(p)

##south
exposome_s <- exposome_s[, colSums(exposome_s != 0) > 0]
res_exposure_S <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c("h_cohort_INMA",
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_exposure_S[[outcome]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s),
                                       expos_name = colnames(exposome_s),
                                       form=form) %>%  mutate(outcome=outcome)}

res_exposure_S <- suppressMessages(purrr::reduce(res_exposure_S, full_join)) %>% arrange(p)

##all
res_exposure_all <- list()
for (outcome in colnames(X_north$Y)){
  form <- paste(outcome, "~", paste(c("h_cohort_BIB", 
                                      "h_cohort_EDEN", 
                                      "h_cohort_KANC",
                                      "h_cohort_INMA",
                                      "h_cohort_MOBA",
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_exposure_all[[outcome]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all),
                                         expos_name = colnames(exposome_all),
                                         form=form) %>%  mutate(outcome=outcome)}

res_exposure_all <- suppressMessages(purrr::reduce(res_exposure_all, full_join)) %>% arrange(p)

#Save results
list_res_exposure_BP <- list("res_exposure_N"=res_exposure_N,
                          "res_exposure_S"=res_exposure_S,
                          "res_exposure_all"=res_exposure_all)

write_xlsx(list_res_exposure_BP, paste0(output_path,"/ExWAS_BP_exposome.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_BP_exposome.xlsx"))
list_res_exposure_BP <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_BP_exposome.xlsx"))
names(list_res_exposure_BP)<-sheets

## Association between outcomes-LC and the exposures----------------------------
LC_BP <- readRDS("./results/RGCCA/projections_outcome_adol.rds")
LC_BP_n <- LC_BP$`Proj N`
LC_BP_s <- LC_BP$`Proj S`
LC_BP_all <- LC_BP$`Proj all`

res_exp_N <- ExWAS_mixed(data=cbind(lat_vars_n[8:ncol(lat_vars_n)], phenotype_n, LC_BP_n),
                         expos_name = colnames(exposome_n),
                         form <- paste("LC_BP_n", "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + ")))


res_exp_S <- ExWAS_mixed(data=cbind(lat_vars_s[8:ncol(lat_vars_s)], phenotype_s, LC_BP_s),
                         expos_name = colnames(exposome_s),
                         form <- paste("LC_BP_s", "~", paste(c("h_cohort_INMA",
                                                               "h_parity_None", 
                                                               "h_native_None", 
                                                               "h_edumc_None", 
                                                               "h_age_None"), collapse = " + ")))


res_exp_all <- ExWAS_mixed(data=cbind(lat_vars_all[8:ncol(lat_vars_all)], phenotype_all, LC_BP_all),
                           expos_name = colnames(exposome_all),
                           form <- paste("LC_BP_all", "~", paste(c("h_cohort_BIB", 
                                                                   "h_cohort_EDEN", 
                                                                   "h_cohort_KANC",
                                                                   "h_cohort_INMA",
                                                                   "h_cohort_MOBA",
                                                                   "h_parity_None", 
                                                                   "h_native_None", 
                                                                   "h_edumc_None", 
                                                                   "h_age_None"), collapse = " + ")))

#Save results
list_res_BP_LC_exposure <- list("res_exposure_N"=res_exp_N,
                                "res_exposure_S"=res_exp_S,
                                "res_exposure_all"=res_exp_all)

write_xlsx(list_res_BP_LC_exposure, paste0(output_path,"/ExWAS_BP-LC_exposome.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_BP-LC_exposome.xlsx"))
list_res_BP_LC_exposure <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_BP-LC_exposome.xlsx"))
names(list_res_BP_LC_exposure)<-sheets

## Association between outcomes-LC and signatures---------------------------------
LC_BP <- readRDS("./results/RGCCA/projections_outcome_adol.rds")
LC_BP_n <- LC_BP$`Proj N`
LC_BP_s <- LC_BP$`Proj S`
LC_BP_all <- LC_BP$`Proj all`

res_BP_LC_N <- ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n, LC_BP_n),
                         expos_name = name_components,
                         form <- paste("LC_BP_n", "~", paste(c("h_cohort_BIB", 
                                                               "h_cohort_EDEN", 
                                                               "h_cohort_KANC",
                                                               "h_parity_None", 
                                                               "h_native_None", 
                                                               "h_edumc_None", 
                                                               "h_age_None"), collapse = " + ")))

res_BP_LC_S <- ExWAS_mixed(data=cbind(lat_vars_s, phenotype_s, LC_BP_s),
                         expos_name = name_components,
                         form <- paste("LC_BP_s", "~", paste(c("h_cohort_INMA", 
                                                               "h_parity_None", 
                                                               "h_native_None", 
                                                               "h_edumc_None", 
                                                               "h_age_None"), collapse = " + ")))


res_BP_LC_all <- ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all, LC_BP_all),
                           expos_name = name_components,
                           form <- paste("LC_BP_all", "~", paste(c("h_cohort_BIB", 
                                                                   "h_cohort_EDEN", 
                                                                   "h_cohort_KANC",
                                                                   "h_cohort_MOBA",
                                                                   "h_cohort_INMA",
                                                                   "h_parity_None", 
                                                                   "h_native_None", 
                                                                   "h_edumc_None", 
                                                                   "h_age_None"), collapse = " + ")))

#Save results
list_res_BP_LC_signatures <- list("res_exposure_N"=res_BP_LC_N,
                                "res_exposure_S"=res_BP_LC_S,
                                "res_exposure_all"=res_BP_LC_all)

write_xlsx(list_res_BP_LC_signatures, paste0(output_path,"/ExWAS_BP-LC_signatures.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_BP-LC_signatures.xlsx"))
list_res_BP_LC_signatures <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_BP-LC_signatures.xlsx"))
names(list_res_BP_LC_signatures)<-sheets


#Sankeyplot---------------------------------------------------------------------
#Preparar codebook
zdia <- codebook[codebook$Variable_name_TRANS == "hs_zdia_bp",]
zsys <- codebook[codebook$Variable_name_TRANS == "hs_zsys_bp",]

zdia_t1 <- transform(zdia, Variable_name_TRANS = "hs2_zdia_bp.v3_2017_Time1", 
                     description = paste0(zdia$description, " Time1"),
                     Label.short..e.g..for.figures.= "Diastolic BP SD score in childhood")
zdia_t2 <- transform(zdia, Variable_name_TRANS = "hs2_zdia_bp.v3_2017_Time2",
                     description = paste0(zdia$description, " Time2"),
                     Label.short..e.g..for.figures.="Diastolic BP SD score in adolescence")
zsys_t1 <- transform(zsys, Variable_name_TRANS = "hs2_zsys_bp.v3_2017_Time1",
                     description = paste0(zsys$description, " Time1"),
                     Label.short..e.g..for.figures.="Systolic BP SD score in childhood")
zsys_t2 <- transform(zsys, Variable_name_TRANS = "hs2_zsys_bp.v3_2017_Time2",
                     description = paste0(zsys$description, " Time2"),
                     Label.short..e.g..for.figures.="Systolic BP SD score in adolescence")

codebook <- rbind(codebook, zdia_t1, zdia_t2, zsys_t1, zsys_t2)

#PLOTS--------------------------------------------------------------------------
##SANKEY OUTCOME-COMP------------------------------------------------------------

df <- list_res_outcome$res_outcome_all

df <- df %>%
  mutate(outcome = case_when(
    outcome == "hs2_zdia_bp.v3_2017_Time1" ~ "Diastolic BP SD score in childhood",
    outcome == "hs2_zdia_bp.v3_2017_Time2" ~ "Diastolic BP SD score in adolescence",
    outcome == "hs2_zsys_bp.v3_2017_Time1" ~ "Systolic BP SD score in childhood",
    outcome == "hs2_zsys_bp.v3_2017_Time2" ~ "Systolic BP SD score in adolescence",
    TRUE ~ outcome))

df <- df[order(df$variable),]
df$variable <- c(rep(c("Prot-LC"),4), rep(c("Serum-LC"),4), rep(c("Urine-LC"),4))

df <- df %>% filter(p_corrected < 0.05)


nodes <- data.frame(name = unique(c(df$variable, df$outcome)), group = c("#5E0626", "#A7A7A7", "#3C6B66", "#F0B077", "#8CC5E3", "#EA801C", "#3594CC"))
nodes <- nodes[c(1,2,3,5,7,4,6),]

links <- df %>%
  mutate(
    source = match(variable, nodes$name) - 1,
    target = match(outcome, nodes$name) - 1,
    value = abs(beta),
    color = ifelse(beta < 0, "olivedrab", "coral")
  ) %>%
  dplyr::select(source, target, value, color)

# Sankey plot
sankey <- sankeyNetwork(Links = links,
              Nodes = nodes, 
              Source = "source", 
              Target = "target", 
              Value = "value", 
              NodeID = "name", 
              units = "Beta", 
              fontSize = 20, 
              fontFamily = "serif", 
              nodeWidth = 20,
              sinksRight = F,
              iterations = 0, colourScale <- JS(
  "d3.scaleOrdinal()
  .domain(['#5E0626', '#A7A7A7', '#3C6B66', '#8CC5E3', '#3594CC', '#F0B077', '#EA801C', 'olivedrab', 'coral'])
  .range(['#5E0626', '#A7A7A7', '#3C6B66', '#8CC5E3', '#3594CC', '#F0B077', '#EA801C', 'olivedrab', 'coral'])"
), NodeGroup = "group", LinkGroup = "color")

saveWidget(sankey, file = "./results/ExWAS/Sankey_LC_outcome.html", selfcontained = TRUE)
webshot("./results/ExWAS/Sankey_LC_outcome.html", 
        file = "./results/ExWAS/Sankey_LC_outcome.tiff",
        selector = "body", delay = 1)

##FOREST ExWAS PLOT COMP-EXPOSURE-----------------------------------------------
df_N <- list_res_exposure$res_exposure_N %>% mutate(source = "North")
df_S <- list_res_exposure$res_exposure_S %>% mutate(source = "South")
df_all <- list_res_exposure$res_exposure_all %>% mutate(source = "Pooled")

df_combined <- rbind(df_N, df_S, df_all)
df_combined <- df_combined %>%  group_by(variable, comp, modality) %>%  filter(any(p < 0.05))
df_combined$source <- factor(df_combined$source, 
                             levels = c("North", "South", "Pooled"), 
                             labels = c("North", "South", "Pooled"))

codebook <- codebook %>% rename(Label = `Label.short..e.g..for.figures.`)
df_combined <- inner_join(df_combined, codebook[, 8:9], by = c("variable" = "Variable_name_TRANS"))

df_combined <- df_combined %>%
  mutate(Label = case_when(modality == 1 ~ paste(Label, "- Yes or No"),
                           modality == 2 ~ paste(Label, "intake - Medium vs low"),
                           modality == 3 ~ paste(Label, "intake - High vs low"),
                           TRUE ~ Label))
df_combined$Label[df_combined$Label=="Traffic__100m"]<-"Traffic (100m)"

write.csv2(df_combined, "./results/ExWAS/ExWAS_comp_exposome_fil.csv")

#PROT
colors <- c("North" = "#3A0417",
            "South" = "#A03358",  
            "Pooled" = "#5E0626") 

pdf("./results/ExWAS/Forest_plot/ExWAS_prot.pdf")
prot<-ggplot(df_combined[df_combined$comp=="prot",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, colour = source)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Prot-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#SERUM
colors <- c("North" = "#707070",
            "South" = "#D3D3D3",  
            "Pooled" = "#A7A7A7") 

pdf("./results/ExWAS/Forest_plot/ExWAS_serum.pdf")
serum <- ggplot(df_combined[df_combined$comp=="serum",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Serum-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#URINE
colors <- c("North" = "#1F3F3C",
            "South" = "#72A09C",  
            "Pooled" = "#3C6B66") 

pdf("./results/ExWAS/Forest_plot/ExWAS_urine.pdf")
urine <- ggplot(df_combined[df_combined$comp=="urine",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Urine-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

prot | serum | urine

##FOREST ExWAS PLOT EXPOSURE-OUTCOME---------------------------------------------
df_N <- list_res_exposure_BP$res_exposure_N %>% mutate(source = "North")
df_S <- list_res_exposure_BP$res_exposure_S %>% mutate(source = "South")
df_all <- list_res_exposure_BP$res_exposure_all %>% mutate(source = "Pooled")

df_combined <- rbind(df_N, df_S, df_all)
df_combined <- df_combined[!(df_combined$variable %in% c("h_edumc_None", "h_age_None")), ]

df_combined <- df_combined %>%  group_by(variable, modality, outcome) %>%  filter(any(p < 0.05))
df_combined$source <- factor(df_combined$source, 
                             levels = c("North", "South", "Pooled"), 
                             labels = c("North", "South", "Pooled"))

codebook <- codebook %>% rename(Label = `Label.short..e.g..for.figures.`)
df_combined <- inner_join(df_combined, codebook[, 8:9], by = c("variable" = "Variable_name_TRANS"))

df_combined <- df_combined %>%
  mutate(Label = case_when(modality == 1 ~ paste(Label, "- Yes or No"),
                           modality == 2 ~ paste(Label, "intake - Medium vs low"),
                           modality == 3 ~ paste(Label, "intake - High vs low"),
                           TRUE ~ Label))
df_combined$Label[df_combined$Label=="Traffic__100m"]<-"Traffic (100m)"
df_combined$Label[df_combined$Label=="Facility_rich"]<-"Facility richness"
df_combined$Label[df_combined$Label=="Facility_dens"]<-"Facility density"


write.csv2(df_combined, "./results/ExWAS/ExWAS_BP_exposome_fil.csv")

#sbpt1
colors <- c("North"="#D88E4B", "South"="#F8D4B0", "Pooled"="#F0B077")


pdf("./results/ExWAS/Forest_plot/ExWAS_SBPt1.pdf")
sbp1<-ggplot(df_combined[df_combined$outcome=="hs2_zsys_bp.v3_2017_Time1",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, colour = source)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "SBP (childhood)") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#sbp_t2
colors <- c("North"="#BF620D", "South"="#F4AE6D", "Pooled"="#EA801C")

pdf("./results/ExWAS/Forest_plot/ExWAS_SBPt2.pdf")
sbp2<-ggplot(df_combined[df_combined$outcome=="hs2_zsys_bp.v3_2017_Time2",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "SBP (adolescence)") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#dbp_t1
colors <- c("North"="#62A2CC", "South"="#B4DAF0", "Pooled"="#8CC5E3")

pdf("./results/ExWAS/Forest_plot/ExWAS_DBPt1.pdf")
dbp1<-ggplot(df_combined[df_combined$outcome=="hs2_zdia_bp.v3_2017_Time1",],
       aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "DBP (childhood)") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

#dbp_t2
colors <- c("North"="#1A6C9B", "South"="#72B9E5", "Pooled"="#3594CC")

pdf("./results/ExWAS/Forest_plot/ExWAS_DBPt2.pdf")
dbp2<-ggplot(df_combined[df_combined$outcome=="hs2_zdia_bp.v3_2017_Time2",],
             aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "DBP (adolescence)") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()

(sbp1 | sbp2) / (dbp1 | dbp2)


#FOREST PLOT EXWAS LC-BP OUTCOMES
df_N <- list_res_BP_LC_exposure$res_exposure_N %>% mutate(source = "North")
df_S <- list_res_BP_LC_exposure$res_exposure_S %>% mutate(source = "South")
df_all <- list_res_BP_LC_exposure$res_exposure_all %>% mutate(source = "Pooled")

df_combined <- rbind(df_N, df_S, df_all)
df_combined <- df_combined %>%  group_by(variable, modality) %>%  filter(any(p < 0.05))
df_combined$source <- factor(df_combined$source, 
                             levels = c("North", "South", "Pooled"), 
                             labels = c("North", "South", "Pooled"))

codebook <- codebook %>% rename(Label = `Label.short..e.g..for.figures.`)
df_combined <- inner_join(df_combined, codebook[, 8:9], by = c("variable" = "Variable_name_TRANS"))

df_combined <- df_combined %>%
  mutate(Label = case_when(modality == 1 ~ paste(Label, "- Yes or No"),
                           modality == 2 ~ paste(Label, "intake - Medium vs low"),
                           modality == 3 ~ paste(Label, "intake - High vs low"),
                           TRUE ~ Label))
df_combined$Label[df_combined$Label=="Traffic__100m"]<-"Traffic (100m)"

write.csv2(df_combined, "./results/ExWAS/ExWAS_BP-LC_exposome_fil.csv")

pdf("./results/ExWAS/Forest_plot/ExWAS_BP-LC_exposure.pdf")
ggplot(df_combined,
             aes(y = reorder(Label, beta), x = beta, xmin = `CI 2.5`, xmax = `CI 97.5`, shape = source, , colour = source)) +
  scale_color_manual(values = colors) + 
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "BP-LC") +
  guides(shape = guide_legend("Population"), color = guide_legend("Population")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))
dev.off()








##########OLD#############################################################################
#EXTRA ANALYSIS-----------------------------------------------------------------
#ExWAS childhood diet-LC________________________________________________________
exposome1 <- read.csv2("./db/exposome/data_HELIXsubcohort_18_09_23_FILTERED_F.csv")
rownames(exposome1) <- exposome1$HelixID
exposome1 <- exposome1[c("h_bf_None",
                         "h_bfdur_Ter",  
                         "hs_bakery_prod_Ter",  
                         "hs_beverages_Ter",  
                         "hs_break_cer_Ter",  
                         "hs_caff_drink_Ter",  
                         "hs_dairy_Ter",  
                         "hs_fastfood_Ter",  
                         "hs_org_food_Ter",  
                         "hs_proc_meat_Ter",  
                         "hs_readymade_Ter",  
                         "hs_total_bread_Ter",  
                         "hs_total_cereal_Ter",  
                         "hs_total_fish_Ter",  
                         "hs_total_fruits_Ter",  
                         "hs_total_lipids_Ter",  
                         "hs_total_meat_Ter",  
                         "hs_total_potatoes_Ter",  
                         "hs_total_sweets_Ter",  
                         "hs_total_veg_Ter",  
                         "hs_total_yog_Ter",  
                         "h_cereal_post_Log",  
                         "h_dairy_post_Log",  
                         "h_fastfood_post_Log",  
                         "h_fish_post_Log",  
                         "h_fruit_post_Log",  
                         "h_legume_post_Log",  
                         "h_meat_post_Log",  
                         "h_nonalc_post_Log",  
                         "h_sugar_post_Log",  
                         "h_veg_post_Log")]

diet_factors<-c("h_bf_None","h_bfdur_Ter", "hs_bakery_prod_Ter", "hs_beverages_Ter", "hs_break_cer_Ter",
                "hs_caff_drink_Ter", "hs_dairy_Ter", "hs_fastfood_Ter", "hs_org_food_Ter",
                "hs_proc_meat_Ter", "hs_readymade_Ter", "hs_total_bread_Ter", "hs_total_cereal_Ter", "hs_total_fish_Ter",
                "hs_total_fruits_Ter", "hs_total_lipids_Ter", "hs_total_meat_Ter", "hs_total_potatoes_Ter", 
                "hs_total_sweets_Ter", "hs_total_veg_Ter", "hs_total_yog_Ter")

exposome1[, !colnames(exposome1) %in% diet_factors] <- scale(exposome1[,!colnames(exposome1) %in% diet_factors])
exposome1[, colnames(exposome1) %in% diet_factors] <- lapply(exposome1[,diet_factors], as.factor)


exposome1_n <- exposome1[rownames(exposome1) %in% rownames(X_north$prot),]
exposome1_s <- exposome1[rownames(exposome1) %in% rownames(X_south$prot),]
exposome1_all <- exposome1[rownames(exposome1) %in% rownames(X_combined$prot),]

lat_vars_n1 <- prepare_df(projections = proj_north["urine"], X_data = X_north[3:4], exposure_data = exposome1_n, response = 2)
lat_vars_s1 <- prepare_df(projections = proj_south["urine"], X_data = X_south[3:4], exposure_data = exposome1_s, response = 2)
lat_vars_all1 <- prepare_df(projections = proj_all["urine"], X_data = X_combined[3:4], exposure_data = exposome1_all, response = 2)


res_diet_urine_n <- ExWAS_mixed(data=cbind(lat_vars_n1, phenotype_n),
                                expos_name = colnames(exposome1),
                                form=paste("urine"," ~1")) %>% mutate(Group = "North")

res_diet_urine_s <- ExWAS_mixed(data=cbind(lat_vars_s1, phenotype_s),
                                expos_name = colnames(exposome1),
                                form=paste("urine"," ~1")) %>% mutate(Group = "South")

res_diet_urine_all <- ExWAS_mixed(data=cbind(lat_vars_all1, phenotype_all),
                                  expos_name = colnames(exposome1),
                                  form=paste("urine"," ~1")) %>% mutate(Group = "Pooled")
df_res_diet_urine <- rbind(res_diet_urine_n, res_diet_urine_s, res_diet_urine_all) %>% arrange(p)
write_xlsx(df_res_diet_urine, "./results/ExWAS/ExWAS_UrineLC_dietvars_childhood.xlsx")

#LM pregnancy hg ~ fish_ter_____________________________________________________
#Load data
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

rownames(phenotype) <- phenotype$HelixID
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]
phenotype_all <- phenotype[rownames(phenotype) %in% rownames(X_combined$prot),]

all(rownames(exposome_all)==rownames(phenotype_all))
df <- cbind(exposome_all, phenotype_all)
df <- df[, !duplicated(names(df))]

#Linear regression
df$h_fish_preg_Ter <- as.factor(df$h_fish_preg_Ter)
lm_model <- lm(hs_hg_m_Log2 ~ h_fish_preg_Ter + h_cohort + e3_sex.x + hs2_visit_age_years_Time1, data = df)
summary(lm_model)
#summary(lm(hs_hg_m_Log2 ~ h_fish_preg_Ter, data = df))

#Boxplot

ggplot(df, aes(x = factor(h_fish_preg_Ter), y = hs_hg_m_Log2, fill = factor(h_fish_preg_Ter))) +
  geom_boxplot() +
  labs(title = "Mercury Levels by Fish Intake",
       x = "Fish Intake During Pregnancy",
       y = "Log2 Mercury Levels") +
  theme_bw() +
  theme(legend.position = "none") 

#Sankey_________________________________________________________________________
#Sankey plot North

sankey_res_N <- plotSankey(res_outcome = list_res_outcome$res_outcome_N, 
                           res_exposure = list_res_exposure$res_exposure_N,
                           name_outcomes = colnames(X_north$Y), 
                           variable_labels = codebook, 
                           path = "./results/Sankey_plot/North",
                           p_val = 0.05)

#Sankey plot south
sankey_res_S <-plotSankey(res_outcome = list_res_outcome$res_outcome_S,
                          res_exposure = list_res_exposure$res_exposure_S,
                          name_outcomes = colnames(X_south$Y), 
                          variable_labels = codebook, 
                          path = "./results/Sankey_plot/South",
                          p_val = 0.05)


#Sankey plot all
sankey_res_all <-plotSankey(res_outcome = list_res_outcome$res_outcome_all, 
                            res_exposure = list_res_exposure$res_exposure_all,
                            name_outcomes = colnames(X_combined$Y), 
                            variable_labels = codebook, 
                            path = "./results/Sankey_plot/all",
                            p_val = 0.05)


## Association between outcomes and mixtures------------------------------------
dbp_mixture_all <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_all_hs2_zdia_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zdia_bp.v3_2017_Time2")
sbp_mixture_all <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_all_hs2_zsys_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zsys_bp.v3_2017_Time2")
dbp_mixture_n <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_n_hs2_zdia_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zdia_bp.v3_2017_Time2")
sbp_mixture_n <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_n_hs2_zsys_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zsys_bp.v3_2017_Time2")
dbp_mixture_s <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_s_hs2_zdia_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zdia_bp.v3_2017_Time2")
sbp_mixture_s <- read_csv2("./results/ExWAS/Mixture_analysis_BP/data/data_s_hs2_zsys_bp.v3_2017_Time21.csv") %>% mutate(Outcome="hs2_zsys_bp.v3_2017_Time2")

mixtures <- c("air_pollution.mixture0.index_zbmi8", "built_env.mixture0.index_zbmi8", "traffic.mixture0.index_zbmi8", "naturalspaces.mixture0.index_zbmi8", "metals.mixture0.index_zbmi8","OCs.mixture0.index_zbmi8", "PFAS.mixture0.index_zbmi8", "lifestyle.mixture0.index_zbmi8","demographic.mixture0.index_zbmi8", "diet.mixture0.index_zbmi8", "meteo.mixture0.index_zbmi8")

mixture_all <- rbind(dbp_mixture_all, sbp_mixture_all)
mixture_n <- rbind(dbp_mixture_n, sbp_mixture_n)
mixture_s <- rbind(dbp_mixture_s, sbp_mixture_s)


##north
all(mixture_n[mixture_n$Outcome=="hs2_zdia_bp.v3_2017_Time2",]$...1==rownames(phenotype_n))
all(mixture_n[mixture_n$Outcome=="hs2_zsys_bp.v3_2017_Time2",]$...1==rownames(phenotype_n))
all(mixture_n[mixture_n$Outcome=="hs2_zdia_bp.v3_2017_Time2",]$...1==rownames(lat_vars_n))
all(mixture_n[mixture_n$Outcome=="hs2_zsys_bp.v3_2017_Time2",]$...1==rownames(lat_vars_n))
all(mixture_s[mixture_s$Outcome=="hs2_zdia_bp.v3_2017_Time2",]$...1==rownames(phenotype_s))
all(mixture_s[mixture_s$Outcome=="hs2_zsys_bp.v3_2017_Time2",]$...1==rownames(phenotype_s))
all(mixture_s[mixture_s$Outcome=="hs2_zdia_bp.v3_2017_Time2",]$...1==rownames(lat_vars_s))
all(mixture_s[mixture_s$Outcome=="hs2_zsys_bp.v3_2017_Time2",]$...1==rownames(lat_vars_s))
all(mixture_all[mixture_all$Outcome=="hs2_zdia_bp.v3_2017_Time2",]$...1==rownames(phenotype_all))
all(mixture_all[mixture_all$Outcome=="hs2_zsys_bp.v3_2017_Time2",]$...1==rownames(phenotype_all))
all(mixture_all[mixture_all$Outcome=="hs2_zdia_bp.v3_2017_Time2",]$...1==rownames(lat_vars_all))
all(mixture_all[mixture_all$Outcome=="hs2_zsys_bp.v3_2017_Time2",]$...1==rownames(lat_vars_all))

res_mixture_N <- list()
for (outcome in colnames(X_north$Y)[3:4]){
  form <- paste(outcome, "~", paste(c("h_cohort_BIB", 
                                      "h_cohort_EDEN", 
                                      "h_cohort_KANC",
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_mixture_N[[outcome]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n, mixture_n[mixture_n$Outcome==outcome,]),
                                         expos_name = mixtures,
                                         form=form) %>%  mutate(outcome=outcome)}

res_mixture_N <- suppressMessages(purrr::reduce(res_mixture_N, full_join)) %>% arrange(p)

##south
res_mixture_S <- list()
for (outcome in colnames(X_north$Y)[3:4]){
  form <- paste(outcome, "~", paste(c("h_cohort_INMA", 
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_mixture_S[[outcome]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s, mixture_s[mixture_s$Outcome==outcome,]),
                                         expos_name = mixtures,
                                         form=form) %>%  mutate(outcome=outcome)}

res_mixture_S <- suppressMessages(purrr::reduce(res_mixture_S, full_join)) %>% arrange(p)

##all
res_mixture_all <- list()
for (outcome in colnames(X_north$Y)[3:4]){
  form <- paste(outcome, "~", paste(c("h_cohort_BIB", 
                                      "h_cohort_EDEN", 
                                      "h_cohort_KANC",
                                      "h_cohort_INMA",
                                      "h_parity_None", 
                                      "h_native_None", 
                                      "h_edumc_None", 
                                      "h_age_None"), collapse = " + "))
  
  res_mixture_all[[outcome]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all, mixture_all[mixture_all$Outcome==outcome,]),
                                           expos_name = mixtures,
                                           form=form) %>%  mutate(outcome=outcome)}

res_mixture_all <- suppressMessages(purrr::reduce(res_mixture_all, full_join)) %>% arrange(p)

#Save results
list_res_BP_mixture <- list("res_mixture_N"=res_mixture_N,
                            "res_mixture_S"=res_mixture_S,
                            "res_mixture_all"=res_mixture_all)

write_xlsx(list_res_BP_mixture, paste0(output_path,"/ExWAS_BP_mixtures.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_BP_mixtures.xlsx"))
list_res_BP_mixture <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_BP_mixtures.xlsx"))
names(list_res_BP_mixture)<-sheets

## Association between signatures and mixtures----------------------------------
LC_BP <- readRDS("./results/RGCCA/projections_outcome_adol.rds")
LC_BP_n <- LC_BP$`Proj N`
LC_BP_s <- LC_BP$`Proj S`
LC_BP_all <- LC_BP$`Proj all`

res_mixture_N <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_mixture_N[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n, mixture_n[mixture_n$Outcome=="hs2_zdia_bp.v3_2017_Time2",]),
                                      expos_name = mixtures,
                                      form=form) %>%  mutate(comp=comp)}

res_mixture_N <- suppressMessages(purrr::reduce(res_mixture_N, full_join)) %>% arrange(p)

##south
res_mixture_S <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_INMA", 
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_mixture_S[[comp]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s, mixture_s[mixture_s$Outcome=="hs2_zdia_bp.v3_2017_Time2",]),
                                      expos_name = mixtures,
                                      form=form) %>%  mutate(comp=comp)}

res_mixture_S <- suppressMessages(purrr::reduce(res_mixture_S, full_join)) %>% arrange(p)

##all
res_mixture_all <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC",
                                   "h_cohort_INMA",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_mixture_all[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all, mixture_all[mixture_all$Outcome=="hs2_zdia_bp.v3_2017_Time2",]),
                                        expos_name = mixtures,
                                        form=form) %>%  mutate(comp=comp)}

res_mixture_all <- suppressMessages(purrr::reduce(res_mixture_all, full_join)) %>% arrange(p)

#Save results
list_res_LC_mixture <- list("res_mixture_N"=res_mixture_N,
                            "res_mixture_S"=res_mixture_S,
                            "res_mixture_all"=res_mixture_all)

write_xlsx(list_res_LC_mixture, paste0(output_path,"/ExWAS_LC_mixtures.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_LC_mixtures.xlsx"))
list_res_LC_mixture <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_LC_mixtures.xlsx"))
names(list_res_LC_mixture)<-sheets

#DIFFERENT MIXTURES (BP-LC as outcome)
## Association between BP-LC and mixtures---------------------------------------
mixture_BP_LC_all <- read.csv2("./results/ExWAS/Mixture_analysis_BP-LC/data/data_all_LC_BP1.csv")
mixture_BP_LC_n <- read.csv2("./results/ExWAS/Mixture_analysis_BP-LC/data/data_n_LC_BP1.csv")
mixture_BP_LC_s <- read.csv2("./results/ExWAS/Mixture_analysis_BP-LC/data/data_s_LC_BP1.csv")



res_exp_N <- ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n, LC_BP_n, mixture_BP_LC_n),
                         expos_name = mixtures,
                         form <- paste("LC_BP_n", "~", paste(c("h_cohort_BIB", 
                                                               "h_cohort_EDEN", 
                                                               "h_cohort_KANC",
                                                               "h_parity_None", 
                                                               "h_native_None", 
                                                               "h_edumc_None", 
                                                               "h_age_None"), collapse = " + ")))

res_exp_S <- ExWAS_mixed(data=cbind(lat_vars_s, phenotype_s, LC_BP_s, mixture_BP_LC_s),
                         expos_name = mixtures,
                         form <- paste("LC_BP_s", "~", paste(c("h_cohort_INMA",
                                                               "h_parity_None", 
                                                               "h_native_None", 
                                                               "h_edumc_None", 
                                                               "h_age_None"), collapse = " + ")))


res_exp_all <- ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all, LC_BP_all, mixture_BP_LC_all),
                           expos_name = mixtures,
                           form <- paste("LC_BP_all","~", paste(c("h_cohort_BIB", 
                                                                  "h_cohort_EDEN", 
                                                                  "h_cohort_KANC",
                                                                  "h_cohort_INMA",
                                                                  "h_parity_None", 
                                                                  "h_native_None", 
                                                                  "h_edumc_None", 
                                                                  "h_age_None"), collapse = " + ")))

#Save results
list_res_BP_LC_mixture <- list("res_exposure_N"=res_exp_N,
                               "res_exposure_S"=res_exp_S,
                               "res_exposure_all"=res_exp_all)

write_xlsx(list_res_BP_LC_mixture, paste0(output_path,"/ExWAS_BP-LC_mixture.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_BP-LC_mixture.xlsx"))
list_res_BP_LC_mixture <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_BP-LC_mixture.xlsx"))
names(list_res_BP_LC_mixture)<-sheets

## Association between signatures and mixtures(BP_LC based)---------------------
LC_BP <- readRDS("./results/RGCCA/projections_outcome_adol.rds")
LC_BP_n <- LC_BP$`Proj N`
LC_BP_s <- LC_BP$`Proj S`
LC_BP_all <- LC_BP$`Proj all`

res_mixture_N <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_mixture_N[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n, mixture_BP_LC_n),
                                      expos_name = mixtures,
                                      form=form) %>%  mutate(comp=comp)}

res_mixture_N <- suppressMessages(purrr::reduce(res_mixture_N, full_join)) %>% arrange(p)

##south
res_mixture_S <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_INMA",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_mixture_S[[comp]] <-ExWAS_mixed(data= cbind(lat_vars_s, phenotype_s, mixture_BP_LC_s),
                                      expos_name = mixtures,
                                      form=form) %>%  mutate(comp=comp)}

res_mixture_S <- suppressMessages(purrr::reduce(res_mixture_S, full_join)) %>% arrange(p)

##all
res_mixture_all <- list()
for (comp in name_components){
  form <- paste(comp, "~", paste(c("h_cohort_BIB", 
                                   "h_cohort_EDEN", 
                                   "h_cohort_KANC",
                                   "h_cohort_INMA",
                                   "h_parity_None", 
                                   "h_native_None", 
                                   "h_edumc_None", 
                                   "h_age_None"), collapse = " + "))
  
  res_mixture_all[[comp]] <-ExWAS_mixed(data=cbind(lat_vars_all, phenotype_all, mixture_BP_LC_all),
                                        expos_name = mixtures,
                                        form=form) %>%  mutate(comp=comp)}

res_mixture_all <- suppressMessages(purrr::reduce(res_mixture_all, full_join)) %>% arrange(p)

#Save results
list_res_LC_mixture <- list("res_mixture_N"=res_mixture_N,
                            "res_mixture_S"=res_mixture_S,
                            "res_mixture_all"=res_mixture_all)

write_xlsx(list_res_LC_mixture, paste0(output_path,"/ExWAS_LC_mixtures_vBP_LC.xlsx"))
sheets <- excel_sheets(paste0(output_path,"/ExWAS_LC_mixtures_vBP_LC.xlsx"))
list_res_LC_mixture <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_LC_mixtures_vBP_LC.xlsx"))
names(list_res_LC_mixture)<-sheets

#MODELS OUTCOME SIGNATURES------------------------------------------------------
#LM SIGNATURES ~ BP with interactions
sbpt1_mod <- lm(hs2_zsys_bp.v3_2017_Time1 ~ prot * serum * urine, data = lat_vars_all)
dbpt1_mod <- lm(hs2_zdia_bp.v3_2017_Time1 ~ prot * serum * urine, data = lat_vars_all)
sbpt2_mod <- lm(hs2_zsys_bp.v3_2017_Time2 ~ prot * serum * urine, data = lat_vars_all)
dbpt2_mod <- lm(hs2_zdia_bp.v3_2017_Time2 ~ prot * serum * urine, data = lat_vars_all)

stargazer(sbpt1_mod, dbpt1_mod, sbpt2_mod, dbpt2_mod,
          type = "text", 
          title = "Blood Pressure Models",
          dep.var.labels = c("SBP T1", "DBP T1", "SBP T2", "DBP T2"),
          covariate.labels = c("Protein", "Serum", "Urine",
                               "Protein x Serum", "Protein x Urine",
                               "Serum x Urine", "Protein x Serum x Urine"),
          out = paste0(output_path,"/bp_models_interactions.rtf"))

#Quadratic 
sbpt1_mod_quad <- lm(hs2_zsys_bp.v3_2017_Time1 ~ 
                       prot + I(prot^2) + 
                       serum + I(serum^2) + 
                       urine + I(urine^2),
                     data = lat_vars_all)

dbpt1_mod_quad <- lm(hs2_zdia_bp.v3_2017_Time1 ~ 
                       prot + I(prot^2) + 
                       serum + I(serum^2) + 
                       urine + I(urine^2),
                     data = lat_vars_all)

sbpt2_mod_quad <- lm(hs2_zsys_bp.v3_2017_Time2 ~ 
                       prot + I(prot^2) + 
                       serum + I(serum^2) + 
                       urine + I(urine^2),
                     data = lat_vars_all)

dbpt2_mod_quad <- lm(hs2_zdia_bp.v3_2017_Time2 ~ 
                       prot + I(prot^2) + 
                       serum + I(serum^2) + 
                       urine + I(urine^2),
                     data = lat_vars_all)

stargazer(sbpt1_mod_quad, dbpt1_mod_quad, sbpt2_mod_quad, dbpt2_mod_quad,
          type = "text", 
          title = "Blood Pressure Models - Quadratic Terms (No Interactions)",
          dep.var.labels = c("SBP T1", "DBP T1", "SBP T2", "DBP T2"),
          covariate.labels = c("Protein", "Protein²", 
                               "Serum", "Serum²", 
                               "Urine", "Urine²"),
          out = paste0(output_path, "/bp_models_quadratic.rtf"))


# Quadratic models WITH INTERACTIONS
sbpt1_mod_quad_int <- lm(hs2_zsys_bp.v3_2017_Time1 ~ 
                           (prot + I(prot^2)) * 
                           (serum + I(serum^2)) * 
                           (urine + I(urine^2)),
                         data = lat_vars_all)

dbpt1_mod_quad_int <- lm(hs2_zdia_bp.v3_2017_Time1 ~ 
                           (prot + I(prot^2)) * 
                           (serum + I(serum^2)) * 
                           (urine + I(urine^2)),
                         data = lat_vars_all)

sbpt2_mod_quad_int <- lm(hs2_zsys_bp.v3_2017_Time2 ~ 
                           (prot + I(prot^2)) * 
                           (serum + I(serum^2)) * 
                           (urine + I(urine^2)),
                         data = lat_vars_all)

dbpt2_mod_quad_int <- lm(hs2_zdia_bp.v3_2017_Time2 ~ 
                           (prot + I(prot^2)) * 
                           (serum + I(serum^2)) * 
                           (urine + I(urine^2)),
                         data = lat_vars_all)

# Get a tidy dataframe of model fit statistics
sbpt1_fit <- glance(sbpt1_mod_quad_int)
dbpt1_fit <- glance(dbpt1_mod_quad_int)
sbpt2_fit <- glance(sbpt2_mod_quad_int)
dbpt2_fit <- glance(dbpt2_mod_quad_int)

# Combine them into one table for easy comparison
all_model_fits <- bind_rows(
  "SBP Time1" = sbpt1_fit,
  "DBP Time1" = dbpt1_fit,
  "SBP Time2" = sbpt2_fit,
  "DBP Time2" = dbpt2_fit,
  .id = "Model"
)
print(all_model_fits, width = Inf)

#modelo con interacciones (no prot:serum)+ urine cuadratico
sbpt1_mod2 <- lm(hs2_zsys_bp.v3_2017_Time1 ~ prot + serum + urine + I(urine^2) +
                   prot:urine + serum:urine + prot:serum:urine, 
                 data = lat_vars_all)

dbpt1_mod2 <- lm(hs2_zdia_bp.v3_2017_Time1 ~ prot + serum + urine + I(urine^2) +
                   prot:urine + serum:urine + prot:serum:urine, 
                 data = lat_vars_all)

sbpt2_mod2 <- lm(hs2_zsys_bp.v3_2017_Time2 ~ prot + serum + urine + I(urine^2) +
                   prot:urine + serum:urine + prot:serum:urine, 
                 data = lat_vars_all)

dbpt2_mod2 <- lm(hs2_zdia_bp.v3_2017_Time2 ~ prot + serum + urine + I(urine^2) +
                   prot:urine + serum:urine + prot:serum:urine, 
                 data = lat_vars_all)
stargazer(sbpt1_mod2, dbpt1_mod2, sbpt2_mod2, dbpt2_mod2,
          type = "text", 
          title = "Blood Pressure Models (sin prot*serum, con urine²)",
          dep.var.labels = c("SBP T1", "DBP T1", "SBP T2", "DBP T2"),
          covariate.labels = c("Protein", "Serum", "Urine", "Urine²",
                               "Protein x Urine", "Serum x Urine",
                               "Protein x Serum x Urine"),
          out = paste0(output_path,"/bp_models_interactions2.rtf"))




#PLOTS BP-SIGNATURES------------------------------------------------------------
library(ggplot2)
library(patchwork)

# Funció per fer els 3 subgràfics d'una variable BP contra prot, serum, urine
make_bp_plots <- function(df, bp_var, title) {
  p1 <- ggplot(df, aes(x = prot, y = .data[[bp_var]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(x = "Protein", y = "BP", title = paste(title, "- prot")) +
    theme_minimal()
  
  p2 <- ggplot(df, aes(x = serum, y = .data[[bp_var]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(x = "Serum", y = "BP", title = paste(title, "- serum")) +
    theme_minimal()
  
  p3 <- ggplot(df, aes(x = urine, y = .data[[bp_var]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "green4") +
    labs(x = "Urine", y = "BP", title = paste(title, "- urine")) +
    theme_minimal()
  
  p1 + p2 + p3  # els junta en una fila
}

# Ara generem els 4 gràfics compostos
bp1 <- make_bp_plots(lat_vars_all, "hs2_zsys_bp.v3_2017_Time1", "SBP Time 1")
bp2 <- make_bp_plots(lat_vars_all, "hs2_zdia_bp.v3_2017_Time1", "DBP Time 1")
bp3 <- make_bp_plots(lat_vars_all, "hs2_zsys_bp.v3_2017_Time2", "SBP Time 2")
bp4 <- make_bp_plots(lat_vars_all, "hs2_zdia_bp.v3_2017_Time2", "DBP Time 2")

# Mostrar tots junts en un layout de 2x2
(bp1 / bp2) | (bp3 / bp4)
