################################# Set R environment 
library(writexl)
library(tidyverse)
library(doParallel)
library(parallel)
library(ggsankey)
library(networkD3)
library(fastDummies)
#install.packages("./packages_R/ggsankey-main", repos = NULL, type = "source")

################################# Defining working directory.
getwd()
setwd("./TFM")
output_path <-"./results/ExWAS"
set.seed(1899)

source("./codi/Functions/functions_ExWAS.R")

###Prepare the dataset----------------------------------------------------------
#Load data
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined

#Codebook
codebook <- read.csv2("./db/exposome/pregnancy_codebook_CURATED.csv")
to_factor <- codebook$Variable_name_TRANS[codebook$Type2 %in% c("Level2")]
dummies <- codebook$Variable_name_TRANS[codebook$Type2 %in% c("Level3","Level4")]

#Load projections
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N`
proj_south <- projections$`Projections S`
proj_all <-  projections$`Projections all`

#Load exposure data
exposome <- read.csv2("./db/exposome/data_HELIXsubcohort_18_09_23_FILTERED_F.csv", row.names = 1)
exposome <- exposome[, -setdiff(1:9, 4)]
exposome <- exposome[colnames(exposome) %in% c(codebook$Variable_name_TRANS,"HelixID")]
vars <- colnames(exposome)

exposome<- dummy_cols(exposome, 
                      select_columns = dummies, 
                      remove_selected_columns = T)

new_dummies <- names(exposome)[!names(exposome) %in% vars]
cols <- c(new_dummies, to_factor)
exposome[cols] <- lapply(exposome[cols], as.factor)

rownames(exposome) <- exposome$HelixID
exposome$HelixID <- NULL
exposome[exposome == Inf] <- 0
colnames(exposome) <- gsub("[- ]", "_", colnames(exposome))
colnames(exposome) <- gsub("<", "less", colnames(exposome))
colnames(exposome) <- gsub(">", "more", colnames(exposome))

exposome_n <- exposome[rownames(exposome) %in% rownames(X_north$prot),]
exposome_s <- exposome[rownames(exposome) %in% rownames(X_south$prot),]
exposome_all <- exposome[rownames(exposome) %in% rownames(X_combined$prot),]

#Prepare df
lat_vars_n <- prepare_df(projections = proj_north, X_data = X_north, exposure_data = exposome_n, response = 4)
lat_vars_s <- prepare_df(projections =  proj_south, X_data = X_south, exposure_data = exposome_s, response = 4)
lat_vars_all <- prepare_df(projections =  proj_all, X_data = X_combined, exposure_data = exposome_all, response = 4)

#ExWAS--------------------------------------------------------------------------
rgcca_validation <- readRDS("./results/RGCCA/model/rgcca_validation.rds")
response <- rgcca_validation$call$response
ncomp <- rgcca_validation$call$ncomp

## Association between components and the exposures------------------------------
name_components <- paste0(rep(names(X_combined)[-response],times=ncomp[-response]))

##north
exposome_n <- exposome_n[, colSums(exposome_n != 0) > 0]
res_exposure_N <- list()
for (comp in name_components){
  form=paste(comp," ~1")
  res_exposure_N[[comp]] <-ExWAS_mixed(data=lat_vars_n,
                                       expos_name = colnames(exposome_n),
                                       form=form) %>%  mutate(comp=comp)}

res_exposure_N <- suppressMessages(purrr::reduce(res_exposure_N, full_join)) %>% arrange(p)

##south
exposome_s <- exposome_s[, colSums(exposome_s != 0) > 0]
res_exposure_S <- list()
for (comp in name_components){
  form=paste(comp," ~1")
  res_exposure_S[[comp]] <-ExWAS_mixed(data=lat_vars_s,
                                       expos_name = colnames(exposome_s),
                                       form=form) %>%  mutate(comp=comp)}

res_exposure_S <- suppressMessages(purrr::reduce(res_exposure_S, full_join)) %>% arrange(p)

##all
res_exposure_all <- list()
for (comp in name_components){
  form=paste(comp," ~1")
  res_exposure_all[[comp]] <-ExWAS_mixed(data=lat_vars_all,
                                         expos_name = colnames(exposome_all),
                                         form=form) %>%  mutate(comp=comp)}

res_exposure_all <- suppressMessages(purrr::reduce(res_exposure_all, full_join)) %>% arrange(p)

#Save results
list_res_exposure <- list("res_exposure_N"=res_exposure_N,
                          "res_exposure_S"=res_exposure_S,
                          "res_exposure_all"=res_exposure_all)

#write_rds(list_res_exposure, paste0(output_path,"/ExWAS_comp_exposome_unfiltered.rds"))
list_res_exposure<-readRDS(paste0(output_path,"/ExWAS_comp_exposome_unfiltered.rds"))

## Association between components and outcomes----------------------------------
##north
res_outcome_N <- list()
for (outcome in colnames(X_north$Y)){
  form = paste(outcome," ~ 1")
  res_outcome_N[[outcome]]  <-ExWAS_mixed(data=lat_vars_n,
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_N <- suppressMessages(purrr::reduce(res_outcome_N, full_join)) %>% arrange(p)

##south
res_outcome_S <- list()
for (outcome in colnames(X_south$Y)){
  form = paste(outcome," ~ 1")
  res_outcome_S[[outcome]]  <-ExWAS_mixed(data=lat_vars_s,
                                          expos_name = name_components,
                                          form=form) %>% mutate(outcome=outcome)}

res_outcome_S <- suppressMessages(purrr::reduce(res_outcome_S, full_join)) %>% arrange(p)

##all
res_outcome_all  <- list()
for (outcome in colnames(X_combined$Y)){
  form = paste(outcome," ~ 1")
  res_outcome_all[[outcome]]  <-ExWAS_mixed(data=lat_vars_all,
                                            expos_name = name_components,
                                            form=form) %>% mutate(outcome=outcome)}

res_outcome_all <- suppressMessages(purrr::reduce(res_outcome_all, full_join)) %>% arrange(p)

#Save results
list_res_outcome <- list("res_outcome_N"=res_outcome_N,
                         "res_outcome_S"=res_outcome_S,
                         "res_outcome_all"=res_outcome_all)

#write_rds(list_res_outcome, paste0(output_path,"/ExWAS_comp_outcome_unfiltered.rds"))
list_res_outcome<-readRDS(paste0(output_path,"/ExWAS_comp_outcome_unfiltered.rds"))

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

p_vals <- c(0.05,0.1)
#Sankey plot North

for (p in p_vals){
  sankey_res_N <- plotSankey(res_outcome = list_res_outcome$res_outcome_N, 
                           res_exposure = list_res_exposure$res_exposure_N,
                           name_outcomes = colnames(X_north$Y), 
                           variable_labels = codebook, 
                           path = "./results/Sankey_plot/North_unfil",
                           p_val = p)

#Sankey plot south
  sankey_res_S <-plotSankey(res_outcome = list_res_outcome$res_outcome_S,
                          res_exposure = list_res_exposure$res_exposure_S,
                          name_outcomes = colnames(X_south$Y), 
                          variable_labels = codebook, 
                          path = "./results/Sankey_plot/South_unfil",
                          p_val = p)


#Sankey plot all
  sankey_res_all <-plotSankey(res_outcome = list_res_outcome$res_outcome_all, 
                            res_exposure = list_res_exposure$res_exposure_all,
                            name_outcomes = colnames(X_combined$Y), 
                            variable_labels = codebook, 
                            path = "./results/Sankey_plot/all_unfil",
                            p_val = p)
}

