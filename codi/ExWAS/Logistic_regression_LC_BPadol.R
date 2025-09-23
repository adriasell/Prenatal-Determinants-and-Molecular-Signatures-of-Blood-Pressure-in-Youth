require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
library(dplyr)
library(tidyverse)
library(fastDummies)
library(ggeffects)
library(writexl)
setwd("./TFM")

#Load data
pheno <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(pheno)<-pheno$HelixID
pheno<-pheno[order(rownames(pheno)),]
pheno <- pheno[c("hs2_BPcat_v3_2017_bin_Time2", "HelixID")]

pheno$hs2_BPcat_v3_2017_bin_Time2 <- as.factor(pheno$hs2_BPcat_v3_2017_bin_Time2)

#Load projections 
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_north) <- names(projections$`Projections N`)
proj_south <- projections$`Projections S` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_south) <- names(projections$`Projections S`)
proj_all <-  projections$`Projections all` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_all) <- names(projections$`Projections all`)

#phenotype data
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)

rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)
phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all), 
                         c("HelixID","e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA")]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north),
                       c("HelixID","e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC")]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south),
                       c("HelixID","e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA")]

#dataset containing proj and BP trajectories 
pheno_all <- pheno[rownames(pheno) %in% rownames(proj_all),]
pheno_n <- pheno[rownames(pheno) %in% rownames(proj_north),]
pheno_s <- pheno[rownames(pheno) %in% rownames(proj_south),]

all(rownames(proj_all)==rownames(pheno_all))
all(rownames(proj_north)==rownames(pheno_n))
all(rownames(proj_south)==rownames(pheno_s))

proj_all$HelixID <- rownames(proj_all)
proj_north$HelixID <- rownames(proj_north)
proj_south$HelixID <- rownames(proj_south)

df_all <- merge(proj_all, pheno_all, by = "HelixID") %>% mutate(Population = "All")
df_all <- merge(df_all, phenotype_all, by = "HelixID") %>% mutate(Population = "All")
df_n <- merge(proj_north, pheno_n, by = "HelixID") %>% mutate(Population = "North")
df_n <- merge(df_n, phenotype_n, by = "HelixID") %>% mutate(Population = "North")
df_s <- merge(proj_south, pheno_s, by = "HelixID") %>% mutate(Population = "South")
df_s <- merge(df_s, phenotype_s, by = "HelixID") %>% mutate(Population = "South")



#Cases for each population
table(df_all$hs2_BPcat_v3_2017_bin_Time2)
table(df_n$hs2_BPcat_v3_2017_bin_Time2)
table(df_s$hs2_BPcat_v3_2017_bin_Time2)


# Process predictions function
process_predictions <- function(predictions) {
  predictions %>%
    mutate(Normal = 1 - Altered) %>%
    pivot_longer(cols = c(Normal, Altered), names_to = "BP_trajectory", values_to = "probability") %>%
    mutate(conf.low = ifelse(BP_trajectory == "Normal", 1 - conf.low, conf.low),
           conf.high = ifelse(BP_trajectory == "Normal", 1 - conf.high, conf.high),
           BP_trajectory = factor(BP_trajectory, levels = c("Normal", "Altered")))
}

# Plot predictions function
plot_predictions <- function(pred_df, exposure_label, pop_label) {
  ggplot(pred_df, aes(x = x, y = probability, fill = BP_trajectory)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
    labs(
      x = paste0(exposure_label, " LC"),
      y = "Probabilities",
      title = paste("Predicted Probabilities of BP Alteration in Adolescence\nPopulation:", pop_label),
      fill = "BP trajectory"
    ) +
    scale_fill_manual(values = c("Normal" = "#089912", "Altered" = "#950C00")) +
    theme_minimal()
}

plots <- list()
populations <- list(
  n = df_n,
  s = df_s,
  all = df_all
)

covariates_per_pop <- list(
  n = paste(c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC"), collapse = " + "),
  s = paste(c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_INMA"), collapse = " + "),
  all = paste(c("e3_sex_Time1", "hs2_visit_age_years_Time1", "h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA", "h_cohort_INMA"), collapse = " + ")
)

signatures <- c("prot", "serum", "urine")

outcome <- "hs2_BPcat_v3_2017_bin_Time2"



# Initialize results dataframes per population
res_n <- data.frame()
res_s <- data.frame()
res_all <- data.frame()

# For easier assignment
res_list <- list(n = res_n, s = res_s, all = res_all)

# Loop over populations
for (pop_name in names(populations)) {
  df <- populations[[pop_name]]                          # Select data
  covars <- covariates_per_pop[[pop_name]]               # Covariates
  
  # Fit combined model once per population
  formula_str <- paste(outcome, "~", paste(signatures, collapse = " + "))
  model <- glm(as.formula(formula_str), data = df, family = binomial())
  
  # Loop over signatures
  for (sig in signatures) {
    pred <- ggemmeans(model, terms = paste0(sig, " [all]")) %>%
      rename(Altered = predicted) %>%
      data.frame()
    
    pred_long <- process_predictions(pred)               # Process predictions
    
    # Generate plot
    plot_obj <- plot_predictions(pred_long, exposure_label = sig, pop_label = toupper(pop_name))
    plots[[pop_name]][[sig]] <- plot_obj                  # Save plot
    
    # Save plot as PDF
    ggsave(filename = paste0("./results/Log_regression/logreg_", pop_name, "_", sig, ".pdf"),
           plot = plot_obj, width = 6, height = 4)
    
    # Extract coefficient summary for this exposure only
    coef_summary <- summary(model)$coefficients
    # Coef row name for the signature
    coef_name <- sig
    
    if (coef_name %in% rownames(coef_summary)) {
      coef_val <- coef_summary[coef_name, "Estimate"]
      se_val <- coef_summary[coef_name, "Std. Error"]
      p_val <- coef_summary[coef_name, "Pr(>|z|)"]
      
      # Create temp result row
      temp_res <- data.frame(
        Population = toupper(pop_name),
        Exposure = sig,
        Coefficient = coef_val,
        StdError = se_val,
        p = p_val
      )
      
      # Append to results dataframe in res_list
      res_list[[pop_name]] <- rbind(res_list[[pop_name]], temp_res)
    }
  }
}

# Assign back results dataframes
res_n <- res_list$n
res_s <- res_list$s
res_all <- res_list$all

res_n$p_corrected <- p.adjust(res_n$p, method = "fdr")
res_s$p_corrected <- p.adjust(res_s$p, method = "fdr")
res_all$p_corrected <- p.adjust(res_all$p, method = "fdr")


log_reg_BP_cat_signatures <- list("North"=res_n,
                                 "South"=res_s,
                                 "Pooled"=res_all)

write_xlsx(log_reg_BP_cat_signatures, "./results/Log_regression/BP_cat_adol-LC.xlsx")





