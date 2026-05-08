## #########################################################
## 
## Get background color from excel to extract data
## 
## - If background color == color then get values from cell and return a 
## dataframe with cell values
## 
## 
## #########################################################

library(xlsx)




# Get Styles from Cells
cellColor <- function(style) 
{
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
    rgb <- paste(rgb, collapse = "")
    return(rgb)
}


extract_main_technical_variables <- function( fileVariables, color, ncols, colnames) 
{
    
    # Get important variables (in yellow)
    wb     <- xlsx::loadWorkbook(fileVariables)
    covariates <- xlsx::getSheets(wb)[[3]]
    
    rows  <- getRows(covariates) # get all rows
    cells <- getCells(rows)     # get all cells
    
    styles <- sapply(cells, getCellStyle) # Get cell styles
    
    BackColor <- list(main = color)
    m     <- match(sapply(styles, cellColor), BackColor)
    labs  <- names(BackColor)[m]
    cells <- cells[which(labs=="main")]
    cells <- cells[-1] # Remove first cell with background color (yellow) === main technical variables
    
    technicalVars <- as.data.frame(matrix(sapply(cells, getCellValue), 
                                          ncol = ncols, byrow = T)) 
    colnames(technicalVars) <- colnames
    
    return(technicalVars)
    
}