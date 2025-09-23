################################################################
#################### Serum  preproc ############################
################################################################

library(lumi)
library(doParallel)
library(minfi)
library(data.table)
library(readxl)
library(xlsx)
library(dplyr)
library(omics)
library(ggplot2)
library(corrplot)
library(isva)
library(SmartSVA)
library(PCAtools)
library(smplot2)
library(DescTools)
library(fastDummies)
library(polycor)

getwd()
setwd("./TFM")
output_path <- "./results/denoising/denoising_metab/serum"

##············· DATA PREPROCESSING

# ///  Data winsorization 
# Loading omics data
data<-readRDS("./db/metab/BP_metab_serum_subcohort_N533_sep23.Rds")

# Winsorization functions
getBeta2 <- function( ms, offset  = 0.001) {
  bs <- exprs( ms )
  bs[ bs == 0 ] <- offset
  bs[ bs == 1 ] <- 1 - offset
  return(bs)
}

process_outliers <- function(bs){
  bs_without_outliers <- bs |> t() |>
    as.data.frame()
  bs_without_outliers[,] <- mclapply(bs_without_outliers,Winsorize) |>
    t()
  return(bs_without_outliers)
}

# Obtaining winsorization results
winsor_m <- data|>
  getBeta2(offset  = 0.001) |>
  process_outliers() |>
  t()

# Inserting winsorizated values into the array
exprs(data) <- winsor_m


# Saving filtered winsorizated data
saveRDS(data, file = paste0(output_path,"/serum_winsorized.RDS"))

