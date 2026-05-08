##############################################################################
##
##  Description: 
##  This function gets: 
##      - Specific-omic metadata
##      - What samples are in pane / subcohort or in both datasets.
##  
##  Input:       
##      @subcohortData: string path to subcohort data file
##      @panelData: string path to panel data file
##      @technicalVars: Specific-omic technical variables (variables to extract)
##      @Prefix: Prefix to assign to the variables
##      
##  Output:
##      
##      dataset with the technical variables + sampletype (panel/subcohort/both) 
##      where the prefix has been added to the name to disinguish them from the
##      other omics.
##
##
## @Author: Dolors Pelegri-Sisó - Bioinformatics Core Facility - ISGlobal
## 
##############################################################################


library(minfi)
library(data.table)



loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

getOmicsSpecificMetadata <- function( subcohortData, panelData, technicalVars, Prefix) {
    
    # browser()
    
    metadata_subc <- loadRData(subcohortData)
    metadata_panel <- loadRData(panelData)
    
    # Get unique dataframe from urine panel and subcohort
    meta_panel <- pData(metadata_panel)
    meta_subc <- pData(metadata_subc)
    
    rm( list = c("metadata_panel", "metadata_subc")); gc();
    
    # Add variable to distinguish between panel and Subcohort samples
    meta_panel$SampleType <- "panel"
    meta_subc$SampleType <- "subcohort"
    
    # Get only one copy of the sample if samples is duplicated in subcohort and panel data and update protSubCohort variable
    commonSamples <- intersect(meta_panel$SampleID, meta_subc$SampleID)
    currentOmic <- rbind( meta_panel[-which(meta_panel$SampleID %in% commonSamples),], meta_subc )
    currentOmic[which( currentOmic$SampleID %in% commonSamples ), "SampleType"] <- "both"
    currentOmic$SampleType = as.factor(currentOmic$SampleType)
    
    # browser()
    
    # Add prot prefix to omics variables
    colnames(currentOmic)[which(colnames(currentOmic) %in%  c(technicalVars, "SampleType")  & !colnames(currentOmic) %like% paste0(Prefix,".")) ] <-  paste( Prefix, 
                                                                                                                                                             c(technicalVars[ which(!technicalVars %like% paste0(Prefix,".")) ], "SampleType"), 
                                                                                                                                                             sep = ".") 
    return( currentOmic[ , c("HelixID", "SampleID", "cohort", colnames(currentOmic)[which(colnames(currentOmic) %like% paste0(Prefix,".") )] )] )
    
}