###############################################################################
######################### Prenatal Exposome Analysis ##########################
###############################################################################
library(dplyr)
library(fastDummies)
library(corrplot)

setwd("./TFM")
exposome <- read.csv2("./db/exposome/data_HELIXsubcohort_18_09_23_FILTERED_F.csv")
getwd()
vars <- c(
  # Air pollution   # Traffic load
  "h_no2_ratio_preg_Log",
  "h_pm10_ratio_preg_None",
  "h_pm25_ratio_preg_None",
  "h_trafload_preg_pow1over3",
  "h_trafnear_preg_pow1over3",
  
  # Built environment + noise
  "h_builtdens300_preg_Sqrt",
  "h_fdensity300_preg_Log",
  "h_frichness300_preg_None",
  "h_landuseshan300_preg_None",
  "h_popdens_preg_Sqrt",
  "h_lden_preg_None",
  "h_ln_preg_None",
  
  # Natural spaces
  "h_walkability_mean_preg_None",
  "h_ndvi100_preg_None",

  # Sociodemographic
  "h_age_None",
  "e3_gac_None",
  "e3_asmokyn_p_None",
  "e3_alcpreg_yn_None",
  "h_mbmi_None",
  "hs_wgtgain_None",
  "h_edumc_None",
  
  # Metals
  "hs_cd_m_Log2",
  "hs_hg_m_Log2",
  "hs_pb_m_Log2",
  
  # OCs
  "hs_dde_madj_Log2",
  "hs_hcb_madj_Log2",
  "hs_pcb138_madj_Log2",
  "hs_pcb153_madj_Log2",
  "hs_pcb180_madj_Log2",
  
  # PFAS
  "hs_pfhxs_m_Log2",
  "hs_pfna_m_Log2",
  "hs_pfoa_m_Log2",
  "hs_pfos_m_Log2",
  
  # Diet
  "h_fish_preg_Ter",
  "h_fruit_preg_Ter",
  "h_legume_preg_Ter",
  "h_veg_preg_Ter",
  "h_dairy_preg_Ter",
  "h_meat_preg_Ter",
  
  #Meteorological vars
  "h_humidity_preg_None",
  "h_pressure_preg_None",
  "h_temperature_preg_None",
  "hs_uvdvf_mt_hs_h_None"
  )

exposome <- exposome[c(vars,"X.7")]
exposome$e3_asmokyn_p_None <- as.integer(ifelse(exposome$e3_asmokyn_p_None == "yes", 1, 0))

rownames(exposome) <- exposome$X.7
exposome$X.7 <- NULL
write.csv2(exposome, "./results/ExWAS/exposome_filtered.csv", row.names = T)
str(exposome)


#
codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")


codebook_fil<-codebook[codebook$Variable_name_TRANS %in% vars,] %>%arrange(Group)
write_xlsx(codebook_fil[6:12], "./results/ExWAS/codebook_fil.xlsx")


#Corrplot exposome vars
int_vars <- exposome %>% select(where(is.integer)) %>% names()
df_factored <- exposome %>% mutate(across(all_of(int_vars), as.factor))
df_dummies <- dummy_cols(df_factored, select_columns = int_vars, remove_selected_columns = TRUE, remove_first_dummy = T)
name_map <- setNames(codebook$Label.short..e.g..for.figures., codebook$Variable_name_TRANS)

ordered_vars <- c(
  # Air pollution 
  "h_no2_ratio_preg_Log",
  "h_pm10_ratio_preg_None",
  "h_pm25_ratio_preg_None",
  
  # Built environment 
  "h_builtdens300_preg_Sqrt",
  "h_fdensity300_preg_Log",
  "h_frichness300_preg_None",
  "h_landuseshan300_preg_None",
  "h_popdens_preg_Sqrt",
  "h_lden_preg_None",
  "h_ln_preg_None",
  
  # Traffic load
  "h_trafload_preg_pow1over3",
  "h_trafnear_preg_pow1over3",
  
  # Green spaces and walkability
  "h_ndvi100_preg_None",
  "h_walkability_mean_preg_None",
  
  # Metals
  "hs_cd_m_Log2",
  "hs_hg_m_Log2",
  "hs_pb_m_Log2",
  
  # OCs
  "hs_dde_madj_Log2",
  "hs_hcb_madj_Log2",
  "hs_pcb138_madj_Log2",
  "hs_pcb153_madj_Log2",
  "hs_pcb180_madj_Log2",
  
  # PFAS
  "hs_pfhxs_m_Log2",
  "hs_pfna_m_Log2",
  "hs_pfoa_m_Log2",
  "hs_pfos_m_Log2",
  
  # Lifestyle factors
  "e3_asmokyn_p_None",
  "e3_alcpreg_yn_None",
  "h_mbmi_None",
  "hs_wgtgain_None",
  
  # Demographic factors
  "e3_gac_None",
  "h_age_None",
  "h_edumc_None", 
  
  # Diet
  "h_fish_preg_Ter",  
  "h_fruit_preg_Ter", 
  "h_legume_preg_Ter",
  "h_veg_preg_Ter", 
  "h_dairy_preg_Ter", 
  "h_meat_preg_Ter",
  
  # Meteorological factors
  "h_humidity_preg_None",
  "h_pressure_preg_None",
  "h_temperature_preg_None",
  "hs_uvdvf_mt_hs_h_None"
)

ordered_vars_present <- ordered_vars[ordered_vars %in% colnames(df_dummies)]
df_ordered <- df_dummies[, ordered_vars_present]
colnames(df_ordered) <- name_map[colnames(df_ordered)]
corr_matrix <- cor(df_ordered, use = "pairwise.complete.obs")
ordered_labels <- name_map[ordered_vars_present]
corr_matrix <- corr_matrix[ordered_labels, ordered_labels]

corrplot.mixed(
  corr_matrix, upper = "ellipse", lower = "number",
  tl.pos = "lt", tl.col = "black", tl.offset=1, tl.srt = 60, 
  tl.cex = 0.5, number.cex = 0.5)



