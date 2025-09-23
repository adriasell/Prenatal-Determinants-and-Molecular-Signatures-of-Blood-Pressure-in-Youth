library(ggplot2)
library(tidyr)
library(dplyr)
library(corrplot)
library(pheatmap)
library(reshape2)
library(gridExtra)
library(patchwork)

#Table 2 Metrics________________________________________________________________
rgcca_val <- readRDS("./results/RGCCA/model_results/results_rgcca_validation.rds")
rgcca_N_all <- readRDS("./results/RGCCA/model_results/results_rgcca_trainN_testall.rds")

metric_N <- rgcca_val$metrics$train %>% mutate(Population = "N")
metric_S <- rgcca_val$metrics$test %>% mutate(Population = "S")
metric_all <- rgcca_N_all$`RGCCA pred (X_pred)`$metric$test %>% mutate(Population = "Pooled")

metrics_combined <- rbind(metric_N, metric_S, metric_all)
metrics_combined$Population <- factor(metrics_combined$Population,
                                      levels = c("N", "S", "Pooled"), 
                                      labels = c("North", "South", "Pooled"))

metrics_combined <- metrics_combined[order(rownames(metrics_combined)), ] %>% t() %>% data.frame()
metrics_combined1 <- rbind(metrics_combined["Population",], as.data.frame(lapply(metrics_combined[-5, ], function(x) round(as.numeric(x), 3)))) %>%
  as.data.frame()

rownames(metrics_combined1)[2:5] <- c("Diastolic BP SD score in childhood", "Systolic BP SD score in childhood",
                                     "Diastolic BP SD score in adolescence", "Systolic BP SD score in adolescence")


write.csv2(metrics_combined1,"./results/Paper_figures/metrics.csv")

#Supplementary figure____________________________________________________________
#Hierarchical clustering and correlation LC & biomarkers selected


set.seed(1999)

#Load data----------------------------------------------------------------------
#Latent vars
latent_vars <- read.csv2("./results/RGCCA/model_results/lat_vars.csv", row.names = 1)
colnames(latent_vars) <- c("Prot-LC", "Serum-LC", "Urine-LC")

#Biomarkers
list_X <- readRDS("./results/RGCCA/list_X.rds")
X_combined <- list_X$X_combined
df_biomarkers <- cbind(X_combined$prot[c("IL6", "IL1beta", "HGF", "BAFF", "TNFalfa", "IL8")],
                       X_combined$serum["log.PC.aa.C38.3"],
                       X_combined$urine[c("p.cresol.sulfate", "X3.Indoxylsulfate")]) %>% scale()

#Projections
proj <- readRDS("./results/RGCCA/projections.rds")
projections <-proj$`Projections all` %>% purrr::reduce(cbind)
colnames(projections) <- c("Prot-LC", "Serum-LC", "Urine-LC")

#Plots--------------------------------------------------------------------------
#Corr plot LC
pdf("./results/Paper_figures/Extra/correlation_LC.pdf")
corrplot.mixed(cor(latent_vars), upper = "ellipse", lower = "number",
               tl.pos = "lt", tl.col = "black", tl.offset=1, tl.srt = 40, 
               tl.cex = 0.5, number.cex = 0.5)
dev.off()

#Corr plot selected biomarkers
pdf("./results/Paper_figures/Extra/correlation_biomarkers.pdf")
corrplot.mixed(cor(scale(df_biomarkers)), upper = "ellipse", lower = "number",
               tl.pos = "lt", tl.col = "black", tl.offset=1, tl.srt = 40, 
               tl.cex = 0.5, number.cex = 0.5)
dev.off()

#Hierarchical clustering LC
pdf("./results/Paper_figures/Extra/hclust_LC.pdf")
pheatmap(t(projections), show_colnames = F)
dev.off()

#Hierarchical clustering biomarkers+
pdf("./results/Paper_figures/Extra/hclust_biomarkers.pdf")
pheatmap(t(df_biomarkers),show_colnames = F)
dev.off()


#BP distribution----------------------------------------------------------------
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

colors <- c(
  "res_dia1" = "#8CC5E3",  
  "res_dia2" = "#3594CC",  
  "res_sys1" = "#F0B077",  
  "res_sys2" = "#EA801C"  
)

p1 <- ggplot(phenotype, aes(x = hs2_zdia_bp.v3_2017_Time1)) +
  geom_histogram(color = "black", fill = colors["res_dia1"]) +
  ggtitle("z-Diastolic BP - Childhood") +
  theme_minimal()

p2 <- ggplot(phenotype, aes(x = hs2_zdia_bp.v3_2017_Time2)) +
  geom_histogram(color = "black", fill = colors["res_dia2"]) +
  ggtitle("z-Diastolic BP - Adolescence") +
  theme_minimal()

p3 <- ggplot(phenotype, aes(x = hs2_zsys_bp.v3_2017_Time1)) +
  geom_histogram(color = "black", fill = colors["res_sys1"]) +
  ggtitle("z-Systolic BP - Childhood") +
  theme_minimal()

p4 <- ggplot(phenotype, aes(x = hs2_zsys_bp.v3_2017_Time2)) +
  geom_histogram(color = "black", fill = colors["res_sys2"]) +
  ggtitle("z-Systolic BP - Adolescence") +
  theme_minimal()

# Arrange in 2x2 grid
grid.arrange(p1, p2, p3, p4, ncol = 2)

#Correlation heat map

corr_data <- phenotype[, c(
  "hs2_zdia_bp.v3_2017_Time1",
  "hs2_zdia_bp.v3_2017_Time2",
  "hs2_zsys_bp.v3_2017_Time1",
  "hs2_zsys_bp.v3_2017_Time2"
)]

colnames(corr_data) <- c("zDBP - Childhood", "zDBP - Adolescence", "zSBP - Childhood", "zSBP - Adolescence")

corr_matrix <- cor(corr_data, use = "pairwise.complete.obs")
corr_melt <- melt(corr_matrix)

ggplot(corr_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
    limit = c(-1, 1), name = "Pearson\nr"
  ) +
  theme_minimal() +
  labs(title = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


#Covariates distribution--------------------------------------------------------

# Barplot of Serum Run Order
library(lumi)
library(doParallel)
library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(dplyr)
library(omics)
library(ggplot2)
library(corrplot)
library(isva)
library(SmartSVA)
library(PCAtools)
library(smplot2)
library(fastDummies)
library(polycor)

setwd("./TFM")
output_path <- "./results/denoising/denoising_metab/serum"
getwd()
##············· PREPARING VARIABLES FOR MODELLING

# /// Load HelixID by group (N/S)

helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)


# /// Loading external functions
source("./script/denoising_2024/functions/extract_data_by_xlsx_styles.R")
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

# /// Loading data
# ···· Metadata .RData filepath
metadataFile <- "./script/denoising_2024/metadata/HELIX_SVA_common_OmicsMetadata_20231026.RData" 

# ···· Omic data (Rdata file with an ExpressionSet or GenomicRatioSet)
# ··················· INDICAR RDATA NUEVO CON LA WINSORIZACIÓN
omicFile <- paste0(output_path,"/serum_winsorized.RDS")

# ···· Phenotype data 
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

rownames(phenotype) <- phenotype$HelixID

names(phenotype)[grep("cohort.x", names(phenotype))] <- "cohort"
names(phenotype)[grep("h_ethnicity_c.x", names(phenotype))] <- "h_ethnicity_c"
names(phenotype)[grep("e3_sex_Time1", names(phenotype))] <- "e3_sex"


# ···· Modify phenotype, cohort and sex, with dummies groups (1 vs. all columns)
phenotype$h_ethnicity_c <- ifelse(phenotype$h_ethnicity_c %in% c("Asian","Pakistani"), "Asian_pakistani", phenotype$h_ethnicity_c)

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_ethnicity_c",
                        remove_selected_columns = TRUE)

phenotype <- dummy_cols(phenotype,
                        select_columns = "cohort",
                        remove_selected_columns = T)

phenotype <- dummy_cols(phenotype,
                        select_columns = "e3_sex",
                        remove_selected_columns = TRUE)

phenotype <- phenotype[,colnames(phenotype)!=c("h_ethnicity_c_NA")]

#Convert var to factor
eth_dummies <- names(phenotype)[grep("h_ethnicity_", names(phenotype))]
phenotype[eth_dummies] <- lapply(phenotype[eth_dummies], as.factor)

coh_dummies <- names(phenotype)[grep("cohort", names(phenotype))]
phenotype[coh_dummies] <- lapply(phenotype[coh_dummies], as.factor)

omicsLayer <- "Serum"
phenotypeID <- "SampleID"
phenotypeVariables <- c("ALL")
ethnic <- c('ALL')


#### STEP 3. Data preparation
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

omics <- getfullOmicsPhenotype(omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, ethnic )
ids1<-colnames(omics$data)
serumrunorder <- pData(omics$data)["Serum.run_order1"]
serumrunorder$Serum.run_order1 <- as.factor(serumrunorder$Serum.run_order1)

library(lumi)
library(doParallel)
library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(dplyr)
library(omics)
library(ggplot2)
library(corrplot)
library(isva)
library(SmartSVA)
library(PCAtools)
library(smplot2)
library(fastDummies)
library(polycor)

setwd("./TFM")
output_path <- "./results/denoising/denoising_prot"
getwd()
##············· PREPARING VARIABLES FOR MODELLING

# /// Load HelixID by group (N/S)

helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)


# /// Loading external functions
source("./script/denoising_2024/functions/extract_data_by_xlsx_styles.R")
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

# /// Loading data
# ···· Metadata .RData filepath
metadataFile <- "./script/denoising_2024/metadata/HELIX_SVA_common_OmicsMetadata_20231026.RData" 

# ···· Omic data (Rdata file with an ExpressionSet or GenomicRatioSet)
# ··················· INDICAR RDATA NUEVO CON LA WINSORIZACIÓN
omicFile <- paste0(output_path,"/prot_winsorized.RDS")

# ···· Phenotype data 
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

rownames(phenotype) <- phenotype$HelixID

names(phenotype)[grep("cohort.x", names(phenotype))] <- "cohort"
names(phenotype)[grep("h_ethnicity_c.x", names(phenotype))] <- "h_ethnicity_c"
names(phenotype)[grep("e3_sex_Time1", names(phenotype))] <- "e3_sex"


# ···· Modify phenotype, cohort and sex, with dummies groups (1 vs. all columns)
phenotype$h_ethnicity_c <- ifelse(phenotype$h_ethnicity_c %in% c("Asian","Pakistani"), "Asian_pakistani", phenotype$h_ethnicity_c)

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_ethnicity_c",
                        remove_selected_columns = T)

phenotype <- dummy_cols(phenotype,
                        select_columns = "cohort",
                        remove_selected_columns = T)

phenotype <- dummy_cols(phenotype,
                        select_columns = "e3_sex",
                        remove_selected_columns = TRUE)

phenotype <- phenotype[,colnames(phenotype)!=c("h_ethnicity_c_NA")]

#Convert var to factor
eth_dummies <- names(phenotype)[grep("h_ethnicity_", names(phenotype))]
phenotype[eth_dummies] <- lapply(phenotype[eth_dummies], as.factor)

coh_dummies <- names(phenotype)[grep("cohort", names(phenotype))]
phenotype[coh_dummies] <- lapply(phenotype[coh_dummies], as.factor)

omicsLayer <- "Prot"
phenotypeID <- "SampleID"
phenotypeVariables <- c("ALL")
ethnic <- c('ALL')

#### STEP 3. Data preparation
source("./script/denoising_2024/functions/generic_functions_denoising_v2.R")
source("./script/denoising_2024/functions/mlmer_local.R")

omics <- getfullOmicsPhenotype(omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, ethnic )    
ids1<-colnames(omics$data)

#Imputation hs_dift_mealblood_imp
imp_value<-median(omics$data$hs_dift_mealblood_imp, na.rm = T)
omics$data$hs_dift_mealblood_imp[is.na(omics$data$hs_dift_mealblood_imp)]<-imp_value

# ───────────────────────────────
# Extract confounders
denoise_vars <- c("Prot.plate", "e3_bw", "hs2_visit_age_years_Time1",
                  "hs_dift_mealblood_imp", "e3_sex", "h_ethnicity_c", "h_cohort")

df <- pData(omics$data)[, denoise_vars]

df_long <- df %>%
  select(e3_bw, hs2_visit_age_years_Time1, hs_dift_mealblood_imp) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

df_long1 <- df_long %>%
  filter(Variable %in% c("e3_bw", "hs2_visit_age_years_Time1"))

df_long2 <- df_long %>%
  filter(Variable == "hs_dift_mealblood_imp")

df_cat <- df %>%
  select(where(is.factor)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

df_long1$Variable <- factor(df_long1$Variable,
                            levels = c("e3_bw", "hs2_visit_age_years_Time1"),
                            labels = c("Birthweight", "Age"))

df_cat$Variable <- factor(df_cat$Variable,
                          levels = c("e3_sex", "h_cohort", "h_ethnicity_c", "Prot.plate"),
                          labels = c("Sex", "Cohort", "Ethnicity", "Protein plate"))



# p1: Serum.run_order1 barplot with black border and uniform fill
p1 <- ggplot(serumrunorder, aes(x = Serum.run_order1)) +
  geom_bar(color = "black", fill = "steelblue") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  labs(title = "Serum run order")

# p2: histograms for e3_bw and hs2_visit_age_years_Time1 with black border, same fill
p2 <- ggplot(df_long1, aes(x = Value)) +
  geom_histogram(bins = 30, color = "black", fill = "steelblue", alpha = 0.8) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# p3: histogram for hs_dift_mealblood_imp with black border, same fill
p3 <- ggplot(df_long2, aes(x = Value)) +
  geom_histogram(bins = 30, color = "black", fill = "steelblue", alpha = 0.8) +
  facet_wrap(~ Variable, scales = "free", labeller = labeller(Variable = function(x) "")) +  # remove facet title here
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  labs(title = "Time to last meal")

# p4: categorical barplots with black border, uniform fill
p4 <- ggplot(df_cat, aes(x = Value)) +
  geom_bar(color = "black", fill = "steelblue") +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# Combine all plots with consistent layout
combined_plot <- (p1 + p3) / p2 / p4 +
  plot_layout(heights = c(1, 1.5, 1.5))

# Display the combined plot
combined_plot







