#Division btw north (BiB, EDEN, KANC, MoBa) and south group (INMA/SAB, RHEA)
setwd("./TFM")
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
x<-readRDS(omicFile)
getwd()
north<-c("BIB","EDEN","KANC","MOBA")
south<-c("SAB","RHEA")

phenotype$cohort.x
n_helixid<-phenotype$HelixID[phenotype$cohort.x %in% north]
s_helixid<-phenotype$HelixID[phenotype$cohort.x %in% south]

write.csv(n_helixid,"./helixid_n.csv")
write.csv(s_helixid,"./helixid_s.csv")




table(phenotype$ancestry_prediction[phenotype$cohort.x %in% north])
table(phenotype$ancestry_prediction[phenotype$cohort.x %in% south])









omics_metadata_final_N <- omics_metadata_final[omics_metadata_final$cohort %in% north,]
phenotype_N <- phenotype[phenotype$cohort.x %in% north,]
omicFile_N <- omicFile[ , which(pData(omicFile)$cohort %in% north)]

omics_metadata_final_S <- omics_metadata_final[omics_metadata_final$cohort %in% south,]
phenotype_S <- phenotype[phenotype$cohort.x %in% south,]
omicFile_S <- omicFile[ , which(pData(omicFile)$cohort %in% south)]

saveRDS(omics_metadata_final_N, "./omics_metadata_final_N.rds")
saveRDS(omics_metadata_final_S, "./omics_metadata_final_S.rds")
saveRDS(phenotype_N, "./phenotype_N.rds")
saveRDS(phenotype_S, "./phenotype_S.rds")
saveRDS(omicFile_N, "./omicFile_N.rds")
saveRDS(omicFile_S, "./omicFile_S.rds")