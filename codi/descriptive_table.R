#Load packages------------------------------------------------------------------
library(dplyr)
library(table1)
library(writexl)
library(tibble)
library(nortest)
#Set wd-------------------------------------------------------------------------
setwd("./TFM")

#Load data----------------------------------------------------------------------
load("./script/denoising_2024/metadata/HELIX_SVA_common_OmicsMetadata_20231026.RData")

list_X <- readRDS("./results/RGCCA/list_X.rds")
X_north <- list_X$X_north
X_south <- list_X$X_south
X_combined <- list_X$X_combined

exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
exposome$HelixID <- rownames(exposome)

phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
phenotype<- merge(phenotype, omics_metadata_final, by = "HelixID")
phenotype<- merge(phenotype, exposome, by = "HelixID")

phenotype <- phenotype[phenotype$HelixID %in% rownames(X_combined$prot),]
phenotype <- phenotype %>% distinct(HelixID, .keep_all = TRUE)

#Data transformation------------------------------------------------------------
rownames(phenotype) <- phenotype$HelixID
#Group
phenotype$Group <- as.factor(ifelse(phenotype$HelixID %in% rownames(X_north$prot), "N", "S"))

#Cohort
phenotype$h_cohort <- factor(phenotype$h_cohort, 
                             levels = c("BIB", "EDEN", "KANC", "MOBA", "INMA", "RHEA"),
                             labels = c("BIB", "EDEN", "KANC", "MOBA", "INMA", "RHEA"))

#Education lvl
phenotype$h_edumc_None.x <- factor(phenotype$h_edumc_None.x,
                                 levels = c(1, 2, 3),
                                 labels = c("Low (Primary school)",
                                            "Medium (Secondary school)",
                                            "High (University degree or higher)"))

#Ethnicity
phenotype$h_ethnicity_c <- factor(ifelse(phenotype$h_ethnicity_c %in% c("Asian", "Pakistani"), "Asian or pakistani",
                                           ifelse(phenotype$h_ethnicity_c == "Caucasian", "Caucasian", "Other")),
                                    levels = c( "Caucasian", "Asian or pakistani", "Other"))

#family native from the country 
phenotype$h_native_None <- factor(phenotype$h_native_None,
                                  levels = c(0, 1, 2),
                                  labels = c("Any", "One parent", "Both parents"))

#Sex
phenotype$e3_sex_Time1 <- factor(phenotype$e3_sex_Time1,
                                 levels = c("male", "female"),
                                 labels = c("Male", "Female"))

#Desc table------------------------------------------------------------------
#Set df
y_t1 <- phenotype[c("hs2_visit_age_years_Time1", "hs2_height_c_Time1", "hs2_c_bmi_None_Time1", 
                    "hs2_zsys_bp.v3_2017_Time1", "hs2_zdia_bp.v3_2017_Time1", "hs2_BPcat_v3_2017_bin_Time1",  
                    "e3_sex_Time1", "h_cohort", "h_edumc_None.x", "h_mbmi_None.x", "h_age_None.x", 
                    "e3_gac_None", "e3_bw", "h_native_None", "h_ethnicity_c", "Group")] %>% mutate (Time= "T1")
y_t2 <- phenotype[c("hs2_visit_age_years_Time2", "hs2_height_c_Time2", "hs2_c_bmi_None_Time2", 
                    "hs2_zsys_bp.v3_2017_Time2", "hs2_zdia_bp.v3_2017_Time2", "hs2_BPcat_v3_2017_bin_Time2", "Group", "h_cohort")] %>% mutate (Time= "T2")

colnames(y_t1)[1:6]<-c("hs2_visit_age_years", "hs2_height", "hs2_c_bmi", "hs2_zsys_bp.v3_2017", "hs2_zdia_bp.v3_2017", "hs2_BPcat_v3_2017")
colnames(y_t2)[1:6]<-c("hs2_visit_age_years", "hs2_height", "hs2_c_bmi", "hs2_zsys_bp.v3_2017", "hs2_zdia_bp.v3_2017", "hs2_BPcat_v3_2017")
all(rownames(y_t1)==rownames(y_t2))
df <- bind_rows(y_t1, y_t2) 
df$Time <- as.factor(df$Time)

#Labels
label(df$hs2_zsys_bp.v3_2017) <- "Systolic BP SD score"
label(df$hs2_zdia_bp.v3_2017) <- "Diastolic BP SD score"
label(df$hs2_BPcat_v3_2017) <- "BP category"
label(df$hs2_visit_age_years) <- "Child age (years)"
label(df$hs2_height) <- "Height (m)"
label(df$hs2_c_bmi) <- "BMI (kg/m²)"
label(df$e3_sex_Time1) <- "Sex"
label(df$h_cohort) <- "Cohort"
label(df$h_edumc_None.x) <- "Maternal education"
label(df$h_mbmi_None.x) <- "Maternal BMI (kg/m²)"
label(df$h_age_None.x) <- "Maternal age (years)"
label(df$e3_gac_None) <- "Gestational age (weeks)"
label(df$e3_bw) <- "Birth weight (g)"
label(df$h_native_None) <- "Is the family native from the country of recruitment?"
label(df$h_ethnicity_c) <- "Ethnicity"

#Create table
tb1<-as.data.frame(table1( ~ hs2_visit_age_years + hs2_height + hs2_c_bmi + hs2_zsys_bp.v3_2017 + 
                             hs2_zdia_bp.v3_2017 + hs2_BPcat_v3_2017 + e3_sex_Time1 + h_cohort +
                             h_edumc_None.x + h_mbmi_None.x + h_age_None.x + e3_gac_None + e3_bw +
                             h_native_None + h_ethnicity_c | Group + Time , data = df))

tb1 <- as.data.frame(lapply(tb1, as.character), stringsAsFactors = FALSE)

#Add p-vals
tb1 <- tb1 %>%
  add_column(p.value = NA, .after = 3) %>% 
  add_column(p.value = NA, .after = 6) %>% 
  add_column(p.value_overall = NA, .after = 7) %>% 
  add_column(p.value = NA, .after = 10)

#age
tb1[2,4] <- t.test(df$hs2_visit_age_years[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_visit_age_years[df$Group == "N" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[2,7] <- t.test(df$hs2_visit_age_years[df$Group == "S" & df$Time == "T1"], 
                   df$hs2_visit_age_years[df$Group == "S" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[2,8] <- t.test(df$hs2_visit_age_years[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_visit_age_years[df$Group == "S" & df$Time == "T1"])$p.value

tb1[2,11] <- t.test(df$hs2_visit_age_years[df$Time == "T1"], 
                    df$hs2_visit_age_years[df$Time == "T2"], 
                    paired = TRUE)$p.value

#Height
tb1[5,4] <- t.test(df$hs2_height[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_height[df$Group == "N" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[5,7] <- t.test(df$hs2_height[df$Group == "S" & df$Time == "T1"], 
                   df$hs2_height[df$Group == "S" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[5,8] <- t.test(df$hs2_height[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_height[df$Group == "S" & df$Time == "T1"])$p.value

tb1[5,11] <- t.test(df$hs2_height[df$Time == "T1"], 
                    df$hs2_height[df$Time == "T2"], 
                    paired = TRUE)$p.value

#bmi
tb1[8,4] <- t.test(df$hs2_c_bmi[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_c_bmi[df$Group == "N" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[8,7] <- t.test(df$hs2_c_bmi[df$Group == "S" & df$Time == "T1"], 
                   df$hs2_c_bmi[df$Group == "S" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[8,8] <- t.test(df$hs2_c_bmi[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_c_bmi[df$Group == "S" & df$Time == "T1"])$p.value

tb1[8,11] <- t.test(df$hs2_c_bmi[df$Time == "T1"], 
                    df$hs2_c_bmi[df$Time == "T2"], 
                    paired = TRUE)$p.value


#zsys
tb1[11,4] <- t.test(df$hs2_zsys_bp.v3_2017[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_zsys_bp.v3_2017[df$Group == "N" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[11,7] <- t.test(df$hs2_zsys_bp.v3_2017[df$Group == "S" & df$Time == "T1"], 
                   df$hs2_zsys_bp.v3_2017[df$Group == "S" & df$Time == "T2"], 
                   paired = TRUE)$p.value

tb1[11,8] <- t.test(df$hs2_zsys_bp.v3_2017[df$Group == "N" & df$Time == "T1"], 
                   df$hs2_zsys_bp.v3_2017[df$Group == "S" & df$Time == "T1"])$p.value

tb1[11,11] <- t.test(df$hs2_zsys_bp.v3_2017[df$Time == "T1"], 
                    df$hs2_zsys_bp.v3_2017[df$Time == "T2"], 
                    paired = TRUE)$p.value

#zdia
tb1[14,4] <- t.test(df$hs2_zdia_bp.v3_2017[df$Group == "N" & df$Time == "T1"], 
                    df$hs2_zdia_bp.v3_2017[df$Group == "N" & df$Time == "T2"], 
                    paired = TRUE)$p.value

tb1[14,7] <- t.test(df$hs2_zdia_bp.v3_2017[df$Group == "S" & df$Time == "T1"], 
                    df$hs2_zdia_bp.v3_2017[df$Group == "S" & df$Time == "T2"], 
                    paired = TRUE)$p.value

tb1[14,8] <- t.test(df$hs2_zdia_bp.v3_2017[df$Group == "N" & df$Time == "T1"], 
                    df$hs2_zdia_bp.v3_2017[df$Group == "S" & df$Time == "T1"])$p.value

tb1[14,11] <- t.test(df$hs2_zdia_bp.v3_2017[df$Time == "T1"], 
                     df$hs2_zdia_bp.v3_2017[df$Time == "T2"], 
                     paired = TRUE)$p.value

#BPCAT
tb1[17,4] <- chisq.test(df$hs2_BPcat_v3_2017[df$Group == "N"], df$Time[df$Group == "N"])$p.value
tb1[17,7] <- chisq.test(df$hs2_BPcat_v3_2017[df$Group == "S"], df$Time[df$Group == "S"])$p.value
tb1[17,8] <- chisq.test(df$hs2_BPcat_v3_2017[df$Time == "T1"], df$Group[df$Time == "T1"])$p.value
tb1[17,11] <- chisq.test(df$hs2_BPcat_v3_2017, df$Time)$p.value


tb1[20,8] <- chisq.test(df$e3_sex_Time1[df$Time == "T1"], df$Group[df$Time == "T1"])$p.value
tb1[31,8] <- chisq.test(df$h_edumc_None.x[df$Time == "T1"], df$Group[df$Time == "T1"])$p.value

tb1[36,8] <- t.test(df$h_mbmi_None.x[df$Group == "N" & df$Time == "T1"], 
                    df$h_mbmi_None.x[df$Group == "S" & df$Time == "T1"])$p.value

tb1[40,8] <- t.test(df$h_age_None.x[df$Group == "N" & df$Time == "T1"], 
                    df$h_age_None.x[df$Group == "S" & df$Time == "T1"])$p.value

tb1[44,8] <- t.test(df$e3_gac_None[df$Group == "N" & df$Time == "T1"], 
                    df$e3_gac_None[df$Group == "S" & df$Time == "T1"])$p.value

tb1[48,8] <- t.test(df$e3_bw[df$Group == "N" & df$Time == "T1"], 
                    df$e3_bw[df$Group == "S" & df$Time == "T1"])$p.value

tb1[57,8] <-  chisq.test(df$h_ethnicity_c[df$Time == "T1"], df$Group[df$Time == "T1"])$p.value

tb1[tb1 == "0 (0%)" | tb1 == "NA (NA)"] <- NA
tb1[c("p.value", "p.value.1", "p.value_overall", "p.value.2")] <- lapply(tb1[c("p.value", "p.value.1", "p.value_overall", "p.value.2")], function(x) ifelse(x < 0.01, "<0.01", round(x,2)))

#Simplify table
cols<-c(2, 3, 5, 6, 9, 10)
tb1[2, cols] <- tb1[3, cols]
tb1[5, cols] <- tb1[6, cols]
tb1[8, cols] <- tb1[9, cols]
tb1[11, cols] <- tb1[12, cols]
tb1[14, cols] <- tb1[15, cols]
tb1[36, cols] <- tb1[37, cols]
tb1[40, cols] <- tb1[41, cols]
tb1[44, cols] <- tb1[45, cols]
tb1[48, cols] <- tb1[49, cols]

#Delete cohort data at adol
tb1[25:30,c(3,6,10)]<-NA

# Delete Missing rows
tb1 <- tb1[!grepl("Missing|Mean|Median", tb1$X.), ]

colnames(tb1) <- c("", rep(c("Childhood", "Adolescence", "p-value"),2), "p-value overall", c("Childhood", "Adolescence", "p-value"))
write_xlsx(tb1, "./results/Paper_figures/desc_table.xlsx")


#Normality test
nor_df <- phenotype[c(
  "hs2_visit_age_years_Time1",
  "hs2_visit_age_years_Time2",
  "hs2_height_c_Time1",
  "hs2_height_c_Time2",
  "hs2_c_bmi_None_Time1",
  "hs2_c_bmi_None_Time2",
  "hs2_zsys_bp.v3_2017_Time1",
  "hs2_zsys_bp.v3_2017_Time2",
  "hs2_zdia_bp.v3_2017_Time1",
  "hs2_zdia_bp.v3_2017_Time2",
  "h_mbmi_None.x",
  "h_age_None.x",
  "e3_gac_None",
  "e3_bw")]

res_lil <- apply(nor_df, 2, function(x) lillie.test(x))

