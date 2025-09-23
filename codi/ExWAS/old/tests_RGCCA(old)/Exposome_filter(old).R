###############################################################################
######################### Prenatal Exposome Analysis ##########################
###############################################################################
library(fastDummies)

setwd("./TFM")
getwd()

#Load data
PREG_childhood <- read.csv2("./results/ExWAS/PREG_childhood_onlineEXWAS_28_10_24.csv")
PREG_adol <- read.csv2("./results/ExWAS/PREG_adolesc_onlineEXWAS_28_10_24.csv")
var_names <- read.csv2("./db/exposome/FIXED_NAMES_DIET_VARS_preg_y_pos.csv")
codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")
exposome <- read.csv2("./db/exposome/data_HELIXsubcohort_18_09_23_FILTERED_F.csv")

#Select relevant outcomes
Health_out <- c("SBP Z-Score - Childhood", "DBP Z-Score - Childhood", "SBP Z-Score - Adolescence", "DBP Z-Score - Adolescence")
PREG_childhood_fil<- PREG_childhood[PREG_childhood$Health.Outcome %in% Health_out,]
PREG_adol_fil<- PREG_adol[PREG_adol$Health.Outcome %in% Health_out,]

#Select significant outcomes
PREG_childhood_sig <- PREG_childhood_fil[PREG_childhood_fil$P.Value<0.1,]
PREG_adol_sig <- PREG_adol_fil[PREG_adol_fil$P.Value<0.1,]
exposure_list_child <- (PREG_childhood_sig$Exposure)
exposure_list_adol <- (PREG_adol_sig$Exposure)

#var_names
vars_child<-var_names[var_names$long_name_figures %in% exposure_list_child,]$db_name
vars_adol<-var_names[var_names$long_name_figures %in% exposure_list_adol,]$db_name
vars_child[!(vars_child %in% colnames(exposome))]
vars_child <- c(vars_child, "h_fish_preg_Ter")
vars_adol[!(vars_adol %in% colnames(exposome))]
vars_adol <- c(vars_adol, "hs_sibpos", "hs_globalsmok_m_None")
exposome_c <- exposome[colnames(exposome) %in% c(vars_child, vars_adol, "X.7")]
str(exposome_c)

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


#Dummies
exposome_c <- dummy_cols(exposome_c, 
                         select_columns = "h_fish_preg_Ter", 
                         remove_selected_columns = TRUE)


cols<-colnames(exposome_c)[grepl("hs_sibpos|h_fish_preg_Ter", colnames(exposome_c))]

#To factor
to_factor <- c("hs_globalsmok_m_None", "e3_psmokanyt_None", "h_mwork_None",cols)
exposome_c[to_factor] <- lapply(exposome_c[to_factor], factor)


rownames(exposome_c) <- exposome$X.7
exposome_c <- exposome_c[colnames(exposome_c) %in% codebook$Variable_name_TRANS[codebook$Period=="Pregnancy"]]

write.csv2(exposome_c, "./results/ExWAS/exposome_filtered.csv", row.names = T)
