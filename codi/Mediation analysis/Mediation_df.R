library(dplyr)
library(writexl)
library(purrr)
library(readxl)
library(mediation)

getwd()
setwd("./TFM")

#overlaps btw exposome-BP and exposome-signatures-BP----------------------------
ExWAS_BP_exposome <- read.csv2("./results/ExWAS/ExWAS_BP_exposome_fil.csv")
ExWAS_BP_exposome$signif_fdr <- ifelse(ExWAS_BP_exposome$p_corrected < 0.05, "*", "")
ExWAS_LC_exposome <- read.csv2("./results/ExWAS/ExWAS_comp_exposome_fil.csv")
ExWAS_LC_exposome <- ExWAS_LC_exposome[ExWAS_LC_exposome$p<0.05,]
ExWAS_LC_exposome$signif_fdr <- ifelse(ExWAS_LC_exposome$p_corrected < 0.05, "*", "")

sheet_names <- excel_sheets("./results/ExWAS/ExWAS_comp_outcome.xlsx")
ExWAS_BP_LC <- lapply(sheet_names, function(sheet) {read_xlsx("./results/ExWAS/ExWAS_comp_outcome.xlsx", sheet = sheet)})
names(ExWAS_BP_LC) <- sheet_names

ExWAS_BP_LC$res_outcome_N <- ExWAS_BP_LC$res_outcome_N %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "North")

ExWAS_BP_LC$res_outcome_S <- ExWAS_BP_LC$res_outcome_S %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "South")

ExWAS_BP_LC$res_outcome_all <- ExWAS_BP_LC$res_outcome_all %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "Pooled")

ExWAS_BP_combined <- bind_rows(
  ExWAS_BP_LC$res_outcome_N,
  ExWAS_BP_LC$res_outcome_S,
  ExWAS_BP_LC$res_outcome_all
)

df_SBP <- list("Exposome-Signatures"=  ExWAS_LC_exposome %>% arrange(Label),
               "Exposome-BP"= ExWAS_BP_exposome[ExWAS_BP_exposome$outcome=="hs2_zsys_bp.v3_2017_Time2",] %>% arrange(Label),
               "Signatures-BP"=ExWAS_BP_combined[ExWAS_BP_combined$outcome=="hs2_zsys_bp.v3_2017_Time2",]%>% arrange(variable))

df_DBP <- list("Exposome-Signatures"=  ExWAS_LC_exposome %>% arrange(Label),
               "Exposome-BP"= ExWAS_BP_exposome[ExWAS_BP_exposome$outcome=="hs2_zdia_bp.v3_2017_Time2",] %>% arrange(Label),
               "Signatures-BP"=ExWAS_BP_combined[ExWAS_BP_combined$outcome=="hs2_zdia_bp.v3_2017_Time2",]%>% arrange(variable))


write_xlsx(df_SBP, "./results/Mediation/df_SBP.xlsx")
write_xlsx(df_DBP, "./results/Mediation/df_DBP.xlsx")

common_SBP <- merge(df_SBP$`Exposome-Signatures`,
                df_SBP$`Exposome-BP`,
                by = c("variable", "modality", "source"))

colnames(common_SBP) <- gsub("\\.x$", "_expLC", colnames(common_SBP))
colnames(common_SBP) <- gsub("\\.y$", "_direct", colnames(common_SBP))


common_DBP <- merge(df_DBP$`Exposome-Signatures`,
                df_DBP$`Exposome-BP`,
                by = c("variable", "modality", "source"))

colnames(common_DBP) <- gsub("\\.x$", "_expLC", colnames(common_DBP))
colnames(common_DBP) <- gsub("\\.y$", "_direct", colnames(common_DBP))


cols <- c("Label_expLC", "variable", "modality", "source", "comp", "outcome", "beta..CI.95.._expLC", "p_expLC", "p_corrected_expLC", "signif_fdr_expLC", "beta..CI.95.._direct", "p_direct",  "p_corrected_direct", "signif_fdr_direct")
df_SBP_common <- list("Common Exposures"= common_SBP[,cols],
                      "Signatures-BP"=ExWAS_BP_combined[ExWAS_BP_combined$outcome=="hs2_zsys_bp.v3_2017_Time2",]%>% arrange(variable))
df_DBP_common <- list("Common Exposures"= common_DBP[,cols],
                      "Signatures-BP"=ExWAS_BP_combined[ExWAS_BP_combined$outcome=="hs2_zdia_bp.v3_2017_Time2",]%>% arrange(variable))
write_xlsx(df_SBP_common, "./results/Mediation/df_SBP_overlapped.xlsx")
write_xlsx(df_DBP_common, "./results/Mediation/df_DBP_overlapped.xlsx")


#overlaps btw exposome-BP_LC and exposome-signatures-BP_LC----------------------
ExWAS_BP_exposome <- read.csv2("./results/ExWAS/ExWAS_BP-LC_exposome_fil.csv")
ExWAS_BP_exposome$signif_fdr <- ifelse(ExWAS_BP_exposome$p_corrected < 0.05, "*", "")
ExWAS_LC_exposome <- read.csv2("./results/ExWAS/ExWAS_comp_exposome_fil.csv")
ExWAS_LC_exposome <- ExWAS_LC_exposome[ExWAS_LC_exposome$p<0.05,]
ExWAS_LC_exposome$signif_fdr <- ifelse(ExWAS_LC_exposome$p_corrected < 0.05, "*", "")

sheet_names <- excel_sheets("./results/ExWAS/ExWAS_BP-LC_signatures.xlsx")
ExWAS_BP_LC <- lapply(sheet_names, function(sheet) {read_xlsx("./results/ExWAS/ExWAS_BP-LC_signatures.xlsx", sheet = sheet)})
names(ExWAS_BP_LC) <- sheet_names

ExWAS_BP_LC$res_exposure_N <- ExWAS_BP_LC$res_exposure_N %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "North")

ExWAS_BP_LC$res_exposure_S <- ExWAS_BP_LC$res_exposure_S %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "South")

ExWAS_BP_LC$res_exposure_all <- ExWAS_BP_LC$res_exposure_all %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "Pooled")

ExWAS_BP_LC_combined <- bind_rows(
  ExWAS_BP_LC$res_exposure_N,
  ExWAS_BP_LC$res_exposure_S,
  ExWAS_BP_LC$res_exposure_all)

df_LC_BP <- list("Exposome-Signatures"=  ExWAS_LC_exposome %>% arrange(comp) %>% arrange(source),
               "Exposome-BP"= ExWAS_BP_exposome %>% arrange(source),
               "Signatures-BP"=ExWAS_BP_LC_combined %>% arrange(variable))

write_xlsx(df_LC_BP, "./results/Mediation/df_LC_BP.xlsx")

common_BP_LC <- merge(df_LC_BP$`Exposome-Signatures`,
                    df_LC_BP$`Exposome-BP`,
                    by = c("variable", "modality", "source"))

colnames(common_BP_LC) <- gsub("\\.x$", "_expLC", colnames(common_BP_LC))
colnames(common_BP_LC) <- gsub("\\.y$", "_direct", colnames(common_BP_LC))

cols <- c("Label_expLC", "variable", "modality", "source", "comp", "beta..CI.95.._expLC", "p_expLC", "p_corrected_expLC", "signif_fdr_expLC", "beta..CI.95.._direct", "p_direct",  "p_corrected_direct", "signif_fdr_direct")
df_common_BP_LC <- list("Common Exposures"= common_BP_LC[,cols] %>% arrange(comp) %>% arrange(source),
                      "Signatures-BP"=ExWAS_BP_LC_combined %>% arrange(variable))

write_xlsx(df_common_BP_LC, "./results/Mediation/df_BP-LC_overlapped.xlsx")

#overlaps btw mixtures-BP and mixtures-signatures-BP----------------------------
sheet_names <- excel_sheets("./results/ExWAS/ExWAS_BP_mixtures.xlsx")
ExWAS_BP_mixtures <- lapply(sheet_names, function(sheet) {read_xlsx("./results/ExWAS/ExWAS_BP_mixtures.xlsx", sheet = sheet)})
names(ExWAS_BP_mixtures) <- sheet_names

ExWAS_BP_mixtures$res_mixture_N <- ExWAS_BP_mixtures$res_mixture_N %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "North")

ExWAS_BP_mixtures$res_mixture_S <- ExWAS_BP_mixtures$res_mixture_S %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "South")

ExWAS_BP_mixtures$res_mixture_all <- ExWAS_BP_mixtures$res_mixture_all %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "Pooled")

ExWAS_BP_mixtures <- bind_rows(
  ExWAS_BP_mixtures$res_mixture_N,
  ExWAS_BP_mixtures$res_mixture_S,
  ExWAS_BP_mixtures$res_mixture_all)

ExWAS_BP_mixtures <- ExWAS_BP_mixtures %>% filter(p<0.05)

sheet_names <- excel_sheets("./results/ExWAS/ExWAS_LC_mixtures.xlsx")
ExWAS_LC_mixtures <- lapply(sheet_names, function(sheet) {read_xlsx("./results/ExWAS/ExWAS_LC_mixtures.xlsx", sheet = sheet)})
names(ExWAS_LC_mixtures) <- sheet_names

ExWAS_LC_mixtures$res_mixture_N <- ExWAS_LC_mixtures$res_mixture_N %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "North")

ExWAS_LC_mixtures$res_mixture_S <- ExWAS_LC_mixtures$res_mixture_S %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "South")

ExWAS_LC_mixtures$res_mixture_all <- ExWAS_LC_mixtures$res_mixture_all %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "Pooled")

ExWAS_LC_mixtures <- bind_rows(
  ExWAS_LC_mixtures$res_mixture_N,
  ExWAS_LC_mixtures$res_mixture_S,
  ExWAS_LC_mixtures$res_mixture_all
)
ExWAS_LC_mixtures <- ExWAS_LC_mixtures %>% filter(p<0.05)

df_mixtures_BP <- list("Mixtures-Signatures"=  ExWAS_BP_mixtures %>% arrange(comp) %>% arrange(Population),
                 "Mixtures-BP"= ExWAS_LC_mixtures %>% arrange(Population),
                 "Signatures-BP"=ExWAS_BP_combined %>% arrange(variable))

write_xlsx(df_mixtures_BP, "./results/Mediation/df_mixtures_BP.xlsx")

common_mixtures_BP <- merge(df_mixtures_BP$`Mixtures-Signatures`,
                            df_mixtures_BP$`Mixtures-BP`,
                            by = c("variable", "Population"))

colnames(common_mixtures_BP) <- gsub("\\.x$", "_mixture-LC", colnames(common_mixtures_BP))
colnames(common_mixtures_BP) <- gsub("\\.y$", "_direct", colnames(common_mixtures_BP))


cols <- c(
  "variable", "Population", "comp", "outcome",
  "beta (CI 95%)_mixture-LC", "p_mixture-LC", "p_corrected_mixture-LC", "signif_fdr_mixture-LC",
  "beta (CI 95%)_direct", "p_direct", "p_corrected_direct", "signif_fdr_direct"
)

df_common_mixtures_BP <- list("Common Exposures"= common_mixtures_BP[,cols] %>% arrange(comp) %>% arrange(Population),
                              "Signatures-BP"=ExWAS_BP_combined %>% arrange(variable))


# Create separate data frames based on the outcome column
df_DBP_mixtures <- list(
  "Common Exposures" = df_common_mixtures_BP[["Common Exposures"]] %>%
    filter(outcome == "hs2_zdia_bp.v3_2017_Time2"),
  
  "Signatures-BP" = df_common_mixtures_BP[["Signatures-BP"]] %>%
    filter(outcome == "hs2_zdia_bp.v3_2017_Time2")
)

write_xlsx(df_DBP_mixtures, "./results/Mediation/df_DBP-mixtures_overlapped.xlsx")

# SBP outcomes with p < 0.05
df_SBP_mixtures <- list(
  "Common Exposures" = df_common_mixtures_BP[["Common Exposures"]] %>%
    filter(outcome == "hs2_zsys_bp.v3_2017_Time2"),
  
  "Signatures-BP" = df_common_mixtures_BP[["Signatures-BP"]] %>%
    filter(outcome == "hs2_zsys_bp.v3_2017_Time2") 
)

write_xlsx(df_SBP_mixtures, "./results/Mediation/df_SBP-mixtures_overlapped.xlsx")


#overlaps btw mixtures-BP_LC and mixtures-signatures-BP_LC----------------------
sheets <- excel_sheets(paste0(output_path,"/ExWAS_LC_mixtures_vBP_LC.xlsx"))
ExWAS_LC_mixtures <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_LC_mixtures_vBP_LC.xlsx"))
names(ExWAS_LC_mixtures)<-sheets

ExWAS_LC_mixtures$res_mixture_N <- ExWAS_LC_mixtures$res_mixture_N %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "North")

ExWAS_LC_mixtures$res_mixture_S <- ExWAS_LC_mixtures$res_mixture_S %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "South")

ExWAS_LC_mixtures$res_mixture_all <- ExWAS_LC_mixtures$res_mixture_all %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "Pooled")

ExWAS_LC_mixtures <- bind_rows(
  ExWAS_LC_mixtures$res_mixture_N,
  ExWAS_LC_mixtures$res_mixture_S,
  ExWAS_LC_mixtures$res_mixture_all)

ExWAS_LC_mixtures <- ExWAS_LC_mixtures %>% filter(p<0.05)

sheets <- excel_sheets(paste0(output_path,"/ExWAS_BP-LC_mixture.xlsx"))
ExWAS_BP_LC_mixture <- lapply(sheets, read_excel, path = paste0(output_path,"/ExWAS_BP-LC_mixture.xlsx"))
names(ExWAS_BP_LC_mixture)<-sheets

ExWAS_BP_LC_mixture$res_exposure_N <- ExWAS_BP_LC_mixture$res_exposure_N %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "North")

ExWAS_BP_LC_mixture$res_exposure_S <- ExWAS_BP_LC_mixture$res_exposure_S %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "South")

ExWAS_BP_LC_mixture$res_exposure_all <- ExWAS_BP_LC_mixture$res_exposure_all %>%
  mutate(signif_fdr = ifelse(p_corrected < 0.05, "*", ""),
         Population = "Pooled")

ExWAS_BP_LC_mixture <- bind_rows(
  ExWAS_BP_LC_mixture$res_exposure_N,
  ExWAS_BP_LC_mixture$res_exposure_S,
  ExWAS_BP_LC_mixture$res_exposure_all)

ExWAS_BP_LC_mixture <- ExWAS_BP_LC_mixture %>% filter(p<0.05)

df_LC_BP_mixtures <- list("Mixtures-Signatures"= ExWAS_LC_mixtures %>% arrange(comp) %>% arrange(Population),
                 "Mixtures-BP_LC"= ExWAS_BP_LC_mixture %>% arrange(Population),
                 "Signatures-BP_LC"=ExWAS_BP_LC_combined %>% arrange(variable))

write_xlsx(df_LC_BP_mixtures, "./results/Mediation/df_LC_BP_mixtures.xlsx")

common_BP_LC_mixture <- merge(df_LC_BP_mixtures$`Mixtures-Signatures`,
                              df_LC_BP_mixtures$`Mixtures-BP`,
                      by = c("variable", "Population"))

colnames(common_BP_LC_mixture) <- gsub("\\.x$", "_expLC", colnames(common_BP_LC_mixture))
colnames(common_BP_LC_mixture) <- gsub("\\.y$", "_direct", colnames(common_BP_LC_mixture))

cols <- c("variable", "Population", "comp", "beta (CI 95%)_expLC", "p_expLC", 
          "p_corrected_expLC", "signif_fdr_expLC", "beta (CI 95%)_direct", 
          "p_direct", "p_corrected_direct", "signif_fdr_direct")

df_common_BP_LC <- list("Common Exposures"= common_BP_LC_mixture[,cols] %>% arrange(comp) %>% arrange(Population),
                        "Signatures-BP"=ExWAS_BP_LC_combined %>% arrange(variable))

write_xlsx(df_common_BP_LC, "./results/Mediation/df_BP-LC_mixture_overlapped.xlsx")



