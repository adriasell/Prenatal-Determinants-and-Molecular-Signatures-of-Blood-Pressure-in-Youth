library(stringr)

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}



# Create data with Phenotypes - Technical Variables and Omics Data
getfullOmicsPhenotype <- function( omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, population   ) 
{
    
    # Load omic Data
    if( tolower(tools::file_ext(omicFile)) == 'rds') {
        omicData <- readRDS(omicFile)
    } else if( tolower(tools::file_ext(omicFile)) == "rdata") {
        omicData <- get(load(omicFile))
    } else{
        stop("File extension unknown, allowed formats: RDS and .RData.")
    }
    
    if( !inherits(omicData, "ExpressionSet") & !inherits(omicData, "GenomicRatioSet") ) {
        stop("The omics data is not of the expected type, a GenomicRatioSet object or an ExpressionSet is expected")
    } 
    
    # Transform str phenotypes to factor phenotypes
    phenotype <- as.data.frame(unclass(phenotype),stringsAsFactors=TRUE)
    
    # Load metadata with technical variables
    omicMetadata <- get(load(metadataFile)[which(load(metadataFile) %like% 'omics_metadata')])
    technicVars <- colnames(omicMetadata)[ which(colnames(omicMetadata) %like% omicsLayer) ]
    
    # Removing subjects withous variables from expressionSet or SummarizedExperiment
    omicData <- omicData[ ,sampleNames(omicData) %in% omicMetadata$SampleID] 
    
    # Filter by population (if not 'ALL' populations are selected)
    if( !'ALL' %in% population) {
        # Filter samples by population
        Sample.filter <- omicMetadata[, "SampleID"][which(omicMetadata$FINAL_ancestry %in% population)]
        if( length(Sample.filter)>0 ){
            omicData <- omicData[ ,sampleNames(omicData) %in% Sample.filter ] 
            omicMetadata <- omicMetadata[which(omicMetadata$SampleID %in% Sample.filter), ]
        } else {
            stop("No samples with ",population," ancestry  ")
        }
    }
    
    # Remove metadata from samples not present in GenomicRatioSet or ExpressionSet
    omicMetadata.f <- merge( data.frame( SampleID = pData(omicData)[, "SampleID"]) , omicMetadata, by = "SampleID")
    omicMetadata.f.ord <- omicMetadata.f[ order( match(omicMetadata.f$SampleID, pData(omicData)$SampleID)), ] # Sort metadata SampleID == methylome SampleID
    
    # Remove variables not related with Current Omics
    toRemove <- sample_collect_cdb[ which(sample_collect_cdb$RelatedOmics != omicsLayer ), ]$var
    if ( omicsLayer != 'Urine') {
        toRemove <- c(toRemove, "uricomments")
    } 
    
    omicMetadata.f.ord <- omicMetadata.f.ord[ , -which(colnames(omicMetadata.f.ord) %in% toRemove)]    
    
    # Get common Vars in omics and phenotypes for merging purposes
    omicMetadataVars <- colnames(omicMetadata.f.ord)
    phenoVars <- colnames(phenotype)
    
    ## #############################################################################################
    ##IMPORTANT!!!: LI ESTIC DONANT MÉS PES A LES DADES FENOTIPIQUES QUE PASSA
    ## L'USUARI!!! HA DE SER AIXÍ O HAURIEN DE PREVALER LES DADES QUE ESTAN A LES METADADES
    ## D'HELIX PERQUÈ SON LES OFICIALS ????!!!!!
    ## 
    ## PREGUNTAR A MARIONA
    ## #############################################################################################
    toRemove <- intersect(omicMetadataVars, phenoVars)[ -which( intersect(omicMetadataVars, phenoVars) == "SampleID")]
    toKeep <- colnames(omicMetadata.f.ord)[  - which(colnames(omicMetadata.f.ord) %in% toRemove)]
    omicMetadata.f.ord <- omicMetadata.f.ord[, toKeep]
    
    if( phenotypeVariables == 'ALL' || 'ALL' %in% phenotypeVariables){
        omicMetadata.f.ord <- merge( omicMetadata.f.ord, phenotype, by="SampleID")
        final_phenoVars <- colnames(phenotype)[-which(colnames(phenotype) == "SampleID")]
    } else {
        omicMetadata.f.ord <- merge( omicMetadata.f.ord, phenotype[,c("SampleID", phenotypeVariables)], by="SampleID")
        final_phenoVars <- phenotypeVariables
    }
    
    omicData <- omicData[ ,sampleNames(omicData) %in% omicMetadata.f.ord$SampleID]
    rownames(omicMetadata.f.ord) <- omicMetadata.f.ord$SampleID
    
    if( inherits(omicData, "ExpressionSet") ) {
        phenoData(omicData) <- new("AnnotatedDataFrame", omicMetadata.f.ord)        
    } else if (inherits(omicData, "GenomicRatioSet")) {
        omicData <- minfi:::.pDataAdd(omicData, omicMetadata.f.ord)
    } 
    
    sampleNames(omicData) <- omicMetadata.f.ord$SampleID

    return(list( data = omicData, 
                 sampes = ncol(omicData),
                 nfeatures = nrow(omicData),
                 nvariables = ncol(pData(omicData)),
                 phenoVars = final_phenoVars,
                 technicVars = technicVars))
}


# Create formula to get SVA for each method: 
#       Interest: With all variables - only protect the phenotype of interest
#       Full: Protect all covariables - protect the phenotype of interest and all the covariates - apply SVA only in technical variables
#       Partial: Protect some covariables - protect the phenotype of interest and all come covariates - apply SVA in technical variables and some covariates

get_SVAFormulas <- function(svaMode, interest_var, covariates, to_protect_covariates, oType) {
    
    browser()
    if(oType == 'ExpressionSet') {
        leftForm <- "t(exprs(omicSVA$data))"
    } else {
        # leftForm <- "assay(omicSVA$data)"
        leftForm <- "t(getM(omicSVA$data))" # Mandatory o deixar-ho com un opcional??
    }
    
    if('Interest' %in% svaMode) {
        if( interest_var !='' && !is.na(interest_var) ) {
            rightForm_inter <- interest_var
        }
    }
    
    if('Full' %in% svaMode) {
        if(length(covariates)>0 && covariates[1]!='' && !is.na(covariates[1]) && interest_var !='' && !is.na(interest_var)) {
            rightForm_full <- paste0(interest_var, covariates, collapse = " + ")
        }
    }
    
    if('Partial' %in% svaMode) {
        if(length(to_protect_covariates)>0 && to_protect_covariates[1]!='' && !is.na(to_protect_covariates[1]) ) {
            rightForm_partial <- paste0(interest_var, to_protect_covariates, collapse = " + ")
        }
    }
    
    return( list( f_Interest = ifelse( exists("rightForm_inter"), paste0(leftForm, " ~ ", rightForm_inter), NA ) , 
                  f_Full = ifelse( exists("rightForm_full"), paste0(leftForm, " ~ ", rightForm_full), NA ) , 
                  f_Partial = ifelse( exists( "rightForm_partial" ), paste0(leftForm, " ~ ", rightForm_partial), NA )
    )
    )
}



## ////////////////////////////////////////////////////////////////////////////////////////////////////////////
## ///                              COMPUTE SURROGATE VARIABLES (SVA)                                       ///
## ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Compute SVA variables
get_SVA <- function(svaMode, formulas, omicSVA, omicsLayer, interest_var, prefix="", respath = NULL, maxAssocPlots = 10, PCs = FALSE) {
    
    SurrogateVars <- sapply(svaMode, function(mode){
        
        if(prefix!="") {
            bfilename = paste0( "SVs_", prefix, "_", omicsLayer, "_", mode)    
        } else {
            bfilename = paste0( "SVs_", omicsLayer, "_", mode)
        }
        
        useFormula <- names(formulas)[which(str_replace(names(formulas), "f_", "") == mode)]
        if( !is.na(formulas[[useFormula]]) ) {
            
            if(class(omicSVA$data)=="GenomicRatioSet") {
                assay(omicSVA$data)[assay(omicSVA$data)==0] <- 0.0001 #esto es importantisimo para trabajar luego con Mvalues (sino, aparecen -Inf)
            }
            
            Y.r <- t( resid( lm(  as.formula(formulas[[useFormula]]) , data = pData(omicSVA$data) ) ))
            
            n.sv <- EstDimRMT( Y.r, TRUE )$dim + 1 # Adding one extra dimension to compensate potential loss of 1 degree of freedom
            mod <- model.matrix( as.formula(gsub('^.*(~)', '\\1',formulas[[useFormula]] )), data = pData(omicSVA$data) ) # Creating model
            
            data_object <- gsub( "\\)$" ,"" ,  gsub( "^t\\("  ,"" , gsub(".~.*","",  formulas[[useFormula]])))
            
            sv.obj <- smartsva.cpp( eval(parse(text = data_object)), mod, mod0 = NULL, n.sv = n.sv )   # creating SVs
            
            # Create RObject with results
            save(sv.obj, file = paste0( bfilename ,".RData") )    
            
            # Get data to be tested
            test_sv <- sv.obj$sv
            test_sv <- as.data.frame(test_sv)
            rownames(test_sv) <- rownames(pData(omicSVA$data))
            colnames(test_sv) <- paste0("SV",1:ncol(test_sv))
            
            # Add phenotype data
            test_sv$varInterest <- pData(omicSVA$data)[ , interest_var]
            test_sv$cohort <- pData(omicSVA$data)$cohort
            
            if(is.null(respath) || respath==""){
                respath == "./"
            } else {
                if(str_sub(respath, -1) != '/')
                    respath <- paste0(respath, "/")
            }
            
            ggplot(test_sv, aes(SV1, SV2) )+
                geom_point(aes(colour = cohort))
            ggsave(paste0(respath, bfilename, "_SV1_2_cohort.pdf"), width = 4, height = 4)
            
            ggplot(test_sv, aes(SV2, SV3) )+
                geom_point(aes(colour = cohort))
            ggsave(paste0(respath, bfilename, "_SV2_3_cohort.pdf"), width = 4, height = 4)
            
            test_sv$e3_sex <- pData(omicSVA$data)$e3_sex
            
            ggplot(test_sv, aes(SV1, SV2) )+
                geom_point(aes(colour = e3_sex))
            ggsave(paste0(respath, bfilename, "_SV1_2_sex.pdf"), width = 4, height = 4)
            
            ggplot(test_sv, aes(SV2, SV3) )+
                geom_point(aes(colour = e3_sex))
            ggsave(paste0(respath, bfilename, "_SV2_3_sex.pdf"), width = 4, height = 4)
            
            pve <- variance_explained_sv(eval(parse(text = data_object)), sv.obj$pprob.gam, sv.obj$pprob.b, n.sv, bfilename, interest_var)
            
            if( length(omicSVA$phenoVars) > 0 ){
                associat <- sapply( omicSVA$phenoVars, association, model = mode, df = pData(omicSVA$data), svs=sv.obj$sv, curvar = interest_var, maxsvs = maxAssocPlots, PCs = PCs, simplify = F)
                
                if( length(associat) > 0 ) {
                    write.csv2( do.call(rbind,associat), paste0( interest_var, "/Summary_association_", bfilename, ".csv"))
                }
            }
            
            return( SVs = sv.obj)
            
        }    
    }, simplify = F)
    
    # Create RObject with all the surrogate variables obtained with all methods
    save(SurrogateVars, file = paste0( respath, interest_var, "/SurrogateVars.RData") )    
    
    return(SurrogateVars)
    
}



## ///////////////////////////////////////////////////////////////////////////////////////////
## ///                              RESIDUALIZE DATA                                       ///
## ///////////////////////////////////////////////////////////////////////////////////////////

# Returns residualized data 
get_SVsResiduals <- function( resSVA, n_SVs, dd, oType) {
    
    if(is.na(n_SVs)) {
        n_SVs <- resSVA$n.sv
    }
    
    if(oType == 'ExpressionSet') {
        
        residualized <- sapply( resSVA, function(data, n.sv, dd) {
            
            model <- model.matrix( ~ data$sv[,1:n_SVs] )
            
            res <- residuals( lmFit( exprs(dd), model ), exprs(dd) )
            
            ## -- ALL DATAT
            residual_SVA <- ExpressionSet(
                assayData = res,
                phenoData = phenoData( dd),
                featureData = featureData( dd )
            )
            return(residual_SVA)
            
        }, n.sv = n_SVs, dd = dd)
        
    } else {
        
        residualized <- sapply( resSVA, function(data, n.sv, dd) {
            
            model <- model.matrix( ~ data$sv[,1:n_SVs] )
            assay(dd)[assay(dd)==0] <- 0.0001 
            res <- residuals( lmFit( getM(dd), model ), getM(dd) )
            
            ## -- ALL DATAT
            assay(dd) <- res
            return(dd)
            
        }, n.sv = n_SVs, dd = dd)
        
    }
    
    
    # Create RObject with all the surrogate variables obtained with all methods
    save(residualized, file = paste0( "Residualized.RData") )    
    
    return(residualized)
    
}
