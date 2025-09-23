# # # # Loading packages.
list.of.packages <- c(
  "RMTL","data.table","rrpack","MCMCpack","GIGrvg","utils","MASS","nortest","MBSP","foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra","haven","corrplot","pheatmap", "FactoMineR", "factoextra", "fastDummies")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
for(package.i in list.of.packages){ #loading packages
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE ))}

options(max.print=100000)
setwd("./TFM")
getwd()

#### Loading data.--------------------------------------------------------------
helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)

#Proj
projections <- readRDS("./results/RGCCA/projections.rds")
proj_all <-  projections$`Projections all` %>% purrr::reduce(cbind) %>% as.data.frame()
proj_north <-  projections$`Projections N` %>% purrr::reduce(cbind) %>% as.data.frame()
proj_south <-  projections$`Projections S` %>% purrr::reduce(cbind) %>% as.data.frame()
colnames(proj_all) <- colnames(proj_south) <- colnames(proj_north) <- c("Prot", "Serum", "Urine")

#Diet vars
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
diet_vars <- c("h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter","h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")
exposome <- exposome[diet_vars] %>%  mutate(across(everything(), as.factor))

exposome_all <- exposome[rownames(exposome) %in% rownames(proj_all),]
exposome_n <- exposome[rownames(exposome) %in% rownames(proj_north),]
exposome_s <- exposome[rownames(exposome) %in% rownames(proj_south),]

#BP
list_X <- readRDS("./results/RGCCA/list_X.rds")
BP <- list_X$X_combined$Y

#Covs
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)
rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)

phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all), 
                         c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")]

phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north),
                       c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC")]

phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south),
                       c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA")]


# # # # MCA tests http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/
res.mca_all <- MCA(exposome_all, graph = F)
res.mca_n <- MCA(exposome_n, graph = F)
res.mca_s <- MCA(exposome_s, graph = F)

#Eigen values
eig.val_all <- get_eigenvalue(res.mca_all)
eig.val_n <- get_eigenvalue(res.mca_n)
eig.val_s <- get_eigenvalue(res.mca_s)

#Get results for vars
var_all <- get_mca_var(res.mca_all)
var_n <- get_mca_var(res.mca_n)
var_s <- get_mca_var(res.mca_s)

# Extract 1:3 dimensions
MCA_all <- as.data.frame(res.mca_all$ind$coord)
MCA_n <- as.data.frame(res.mca_n$ind$coord)
MCA_s <- as.data.frame(res.mca_s$ind$coord)
colnames(MCA_all) <- colnames(MCA_n) <- colnames(MCA_s) <- paste0("hs_RA_PC", 1:ncol(MCA_all))

MCA_all$HelixID <- rownames(MCA_all)
BP$HelixID <- rownames(BP)

all(rownames(phenotype_all)==rownames(proj_all))
all(rownames(phenotype_all)==rownames(MCA_all))
all(rownames(phenotype_n)==rownames(proj_north))
all(rownames(phenotype_n)==rownames(MCA_n))
all(rownames(phenotype_s)==rownames(proj_south))
all(rownames(phenotype_s)==rownames(MCA_s))

res <- list()

for (LC in c("Prot", "Serum", "Urine")){
  for (PC in colnames(MCA_all)[1:5]){
    formula_all = as.formula(paste0(LC, "~", PC,"+ e3_sex_Time1 + hs2_visit_age_years_Time1 + h_cohort_BIB + h_cohort_EDEN + h_cohort_KANC + h_cohort_MOBA + h_cohort_INMA"))
    formula_n = as.formula(paste0(LC, "~", PC,"+ e3_sex_Time1 + hs2_visit_age_years_Time1 + h_cohort_BIB + h_cohort_EDEN + h_cohort_KANC"))
    formula_s = as.formula(paste0(LC, "~", PC,"+ e3_sex_Time1 + hs2_visit_age_years_Time1 + h_cohort_INMA"))
    
    res[[LC]][[PC]] <- list(All = summary(lm(formula_all, data = cbind(phenotype_all, proj_all, MCA_all))),
                            North = summary(lm(formula_n, data = cbind(phenotype_n, proj_north, MCA_n ))),
                            South = summary(lm(formula_s, data = cbind(phenotype_s, proj_south, MCA_s))))
  }
}



















#PLOTS--------------------------------------------------------------------------
fviz_screeplot(res.mca_all, addlabels = TRUE, ylim = c(0, 45))
fviz_screeplot(res.mca_n, addlabels = TRUE, ylim = c(0, 45))
fviz_screeplot(res.mca_s, addlabels = TRUE, ylim = c(0, 45))

fviz_mca_var(res.mca_all, repel = TRUE)
fviz_mca_var(res.mca_n, repel = TRUE)
fviz_mca_var(res.mca_s, repel = TRUE)

fviz_mca_var(res.mca_all, repel = TRUE, choice = "mca.cor")
fviz_mca_var(res.mca_n, repel = TRUE, choice = "mca.cor")
fviz_mca_var(res.mca_s, repel = TRUE, choice = "mca.cor")

fviz_mca_var(res.mca_all, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())
fviz_cos2(res.mca_all, choice = "var", axes = 1:3)

fviz_mca_var(res.mca_n, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())
fviz_cos2(res.mca_n, choice = "var", axes = 1:3)

fviz_mca_var(res.mca_s, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())
fviz_cos2(res.mca_s, choice = "var", axes = 1:3)

Ys <- merge(BP, MCA_all, by= "HelixID")
toplot <- cor(data.matrix(Ys))
corrplot(toplot)