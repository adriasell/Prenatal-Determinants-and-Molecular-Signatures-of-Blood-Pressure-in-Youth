#### libraries
library(limma)
library(Biobase)
library(dplyr)




variance_explained_sv <- function(dat, pprob.gam, pprob.b, n.sv, model, curvar) 
{
    
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
    colnames(expl_var_plot) <- c("SV", "PVE")
    
    
    p <- ggplot(expl_var_plot[1:n.sv,], aes(SV,PVE)) +
        geom_point(size=2, pch=19, color="blue")  +
        ggtitle(paste0(model, " : ", curvar)) +
        scale_color_gradient() +
        geom_text(aes(label = rownames(expl_var_plot[1:n.sv,])),
                  check_overlap = TRUE,
                  hjust=0.5,
                  vjust=2,
                  size=2) +
        xlab("SV") +
        ylab("% explained variance (estimate)")
    
    dir.create(file.path(getwd(), curvar), showWarnings = FALSE)
    pdf(paste0( curvar,"/", model, "_Percent_variance_explained_SV.pdf" ))
        print(p)
    dev.off()
    return(uu_val)
    
}