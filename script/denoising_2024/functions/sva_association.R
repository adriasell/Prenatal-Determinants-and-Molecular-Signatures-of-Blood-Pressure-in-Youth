#### libraries
library(minfi)
library(sva)
library(isva)
library(limma)
library(SmartSVA)
library(Biobase)
library(xlsx)
library(dplyr)

#..# devtools::install_github('smin95/smplot')
library(smplot)
library(ggplot2)






#### % VARIANCE EXPLAINED BY SVs
####        based on : https://support.bioconductor.org/p/88553/
variance_explained_sv <- function(dat, pprob.gam, pprob.b, n.sv, model, curvar, PCs = FALSE) 
{
    vartype = "SV"
    if( PCs ) { vartype="PC" }
    
    pprob <- pprob.gam * (1-pprob.b)
    dats <- dat * pprob
    dats <- dats - rowMeans(dats)
    uu <- eigen(t(dats) %*% dats)
    
    #  could express the fraction of each eigenvalue over the total sum. 
    #  Can be interpreted as the "Percent of explained variance" for each SV
    uu_val <- uu$values / sum(uu$values)
    
    # Plot
    
    x_val <- 1:ncol(dats)
    expl_var_plot <- as.data.frame(cbind(x_val,uu_val))
    colnames(expl_var_plot) <- c(vartype, "PVE")
    
    
    p <- ggplot(expl_var_plot[1:n.sv,], aes(SV,PVE)) +
        geom_point(size=2, pch=19, color="blue")  +
        ggtitle(paste0(model, " : ", curvar)) +
        scale_color_gradient() +
        geom_text(aes(label = rownames(expl_var_plot[1:n.sv,])),
                  check_overlap = TRUE,
                  hjust=0.5,
                  vjust=2,
                  size=2) +
        xlab(vartype) +
        ylab("% explained variance (estimate)")
    
    pdf(paste0( curvar,"/", model, "_Percent_variance_explained_SV.pdf" ))
    print(p)
    dev.off()
    return(uu_val)
    
}





#### 
#### -- FUNCTION FOR ASSOCIATION TEST --
#### 

association <- function(var, df, svs, model, curvar, maxsvs=NA, PCs=FALSE ){
    
    
    vartype = "SV"
    if( PCs ) { vartype="PC" }
    
    
    # Check if variable is factor and n_factors>=2
    if( inherits(df[,var], "factor") && nlevels( droplevels(df[,var]) ) < 2) {
        return(NA)    
    }
    
    if( ( inherits(df[,var], "integer") || inherits(df[,var], "numeric")) && length(unique(df[,var])) < 2) {
        return(NA)    
    }
    
    imgPath <- file.path(getwd(), curvar, "association_plots", vartype)
    dir.create(imgPath, showWarnings = FALSE, recursive = TRUE)
    
    if( is.na(maxsvs) || maxsvs > dim(svs)[2] ) {
        nsvs <- dim(svs)[2]
    } else {
        nsvs <- maxsvs
    }
    
    res <- sapply( 1:nsvs, function(i) {
        
        df2 <- data.frame( as.numeric(svs[,i]), df[,var])
        colnames(df2) <- c(paste0(vartype,i), var)
        
        frm <- as.formula(paste0(vartype,i, " ~ ",var))
        lm.res <- lm(df2[,1] ~ df2[,2]);
        lm.res <- lm(frm, data = df2)
        
        
        #### Plot only association with SV1 : SV5
        if( i %in% seq_along(1:maxsvs) ) {
            if( inherits(df[,var], "factor")) {
                if( length(unique(df[,var])) > 20 ) {
                    warning("Warning: Factor variable?, maybe this is a continuous var=")
                }
                
                if(length(unique(df[,var])) >= 2 ){
                    p <- ggplot(data = df2, mapping = aes( df2[,2], df2[,1], colour = df2[,2])) +
                        geom_boxplot() +
                        xlab(var) +
                        ylab(paste0(vartype,i)) +
                        guides(fill=guide_legend(title=var))
                } else {
                    warning("Warning: Only one factor present ")
                }
                
            } else {
                eq <- substitute(italic(paste0(vartype,i)) == a + b %.% italic(paste0(var))*","~~italic(r)^2~"="~r2,
                                 list(a = format(unname(coef(lm.res)[1]), digits = 2),
                                      b = format(unname(coef(lm.res)[2]), digits = 2),
                                      r2 = format(summary(lm.res)$r.squared, digits = 3)))
                
                p <- ggplot(data = df2, mapping = aes(x = df2[,2], y = df2[,1])) +
                    geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
                    sm_corr_theme() +
                    sm_statCorr() +
                    scale_x_continuous(var) +
                    scale_y_continuous(paste0(vartype,i)) +
                    sm_statCorr()
            }
            
            pdf(paste0(imgPath, "/", model,"_", vartype, i, "_",var,".pdf"))
            print(p)
            dev.off()
        }
        
        print(paste0("Processing ", var, " - ", vartype, i))
        
        
        sum_lm <- summary(lm.res)
        
        if(inherits(coefficients(sum_lm)[2:nrow(coefficients(sum_lm)),], "numeric")) {
            res <- data.frame(t( coefficients(sum_lm)[2:nrow(coefficients(sum_lm)),]) , 
                              R2 = sum_lm$r.squared, 
                              Global_pvalue = pf(summary((lm.res))$fstatistic[1],summary((lm.res))$fstatistic[2],summary((lm.res))$fstatistic[3],lower.tail=F),
                              model)
            res$coefficient <- attributes(sum_lm$coefficients)$dimnames[[1]][2]
            res$SV <- paste0(vartype,i)
            res$association <- paste0(vartype,i," - ",var)
            
        } else {
            if( length(unique(df[,var])) >= 2 ){
                
                res <-  data.frame( coefficients(sum_lm)[2:nrow(coefficients(sum_lm)),] , 
                                    R2 = sum_lm$r.squared, 
                                    Global_pvalue = pf(summary((lm.res))$fstatistic[1],summary((lm.res))$fstatistic[2],summary((lm.res))$fstatistic[3],lower.tail=F), 
                                    model) 
                res$coefficient <- attributes(sum_lm$coefficients)$dimnames[[1]][2:nrow(coefficients(sum_lm))]
                res$SV <- paste0(vartype,i)
                res$association <- paste0(vartype,i," - ",var)
            } else {
                warning("Warning: Only one factor present ")
            }
        }
        return(res) 
        
    } ,simplify = F)
    
    do.call(rbind,res)
    
}
