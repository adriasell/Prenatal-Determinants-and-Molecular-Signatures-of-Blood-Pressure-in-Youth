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
#install.packages("./packages_R/ggsankey-main", repos = NULL, type = "source")

################################# Defining working directory.
getwd()
setwd("./TFM")
output_path <-"./results/sensitivity_bmi/ExWAS"
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
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all$prot), 
                         c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA", "hs2_zbmi_who_Time1")]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north$prot),
                       c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "hs2_zbmi_who_Time1")]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south$prot),
                       c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA", "hs2_zbmi_who_Time1")]


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

#ExWAS--------------------------------------------------------------------------
rgcca_validation <- readRDS("./results/RGCCA/model/rgcca_final.rds")
response <- rgcca_validation$call$response
ncomp <- rgcca_validation$call$ncomp
name_components <- paste0(rep(names(X_combined)[-response],times=ncomp[-response]))

## Association between components and outcomes----------------------------------
##north
res_outcome_N <- list()
for (outcome in colnames(X_north$Y)){
  form = paste(outcome," ~ 1 + hs2_zbmi_who_Time1")
  res_outcome_N[[outcome]]  <-ExWAS_mixed(data=cbind(lat_vars_n, phenotype_n),
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_N <- suppressMessages(purrr::reduce(res_outcome_N, full_join)) %>% arrange(p)

##south
res_outcome_S <- list()
for (outcome in colnames(X_south$Y)){
  form = paste(outcome," ~ 1 + hs2_zbmi_who_Time1")
  res_outcome_S[[outcome]]  <-ExWAS_mixed(data=cbind(lat_vars_s, phenotype_s),
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_S <- suppressMessages(purrr::reduce(res_outcome_S, full_join)) %>% arrange(p)

##all
res_outcome_all  <- list()
for (outcome in colnames(X_combined$Y)){
  form = paste(outcome," ~ 1 + hs2_zbmi_who_Time1")
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


nodes <- data.frame(name = unique(c(df$variable, df$outcome)), group = c("#5E0626", "#3C6B66", "#8CC5E3", "#3594CC", "#F0B077"))

links <- df %>%
  mutate(
    source = match(variable, nodes$name) - 1,
    target = match(outcome, nodes$name) - 1,
    value = abs(beta),
    color = ifelse(beta < 0, "olivedrab", "coral")
  ) %>%
  select(source, target, value, color)

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

saveWidget(sankey, file = "./results/sensitivity_BMI/ExWAS/Sankey_LC_outcome.html", selfcontained = TRUE)
webshot("./results/sensitivity_BMI/ExWAS/Sankey_LC_outcome.html", 
        file = "./results/sensitivity_BMI/ExWAS/Sankey_LC_outcome.tiff",
        selector = "body", delay = 1)

