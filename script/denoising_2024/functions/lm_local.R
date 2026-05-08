# Function from omics package
# Removed unnecessary calculus to obtain only residuals.

lm_local <- function(formula, data=NULL, vars, lrt=TRUE, save.residuals=FALSE,
                        save.ranks=TRUE, ncores = 2L) {
  if (missing(data)) {
    data <- NULL
  }
  
  Y <- as.data.frame(get(omics:::response.name(formula, data=data), envir=environment(formula)))
  if (is.null(data)) {
    data <- environment(formula)
  }
  
  options(na.action='na.pass')
  mm <- model.matrix(formula, data=data)
  mm <- as.matrix(mm[,1])
  idx <- which(colnames(mm) %in% vars)
  
  
  
  if (length(idx) == 0 & !save.residuals) {
    stop("No variables selected")
  }
  
  vars <- colnames(mm)[idx]
  colnames(mm) <- sprintf("V%d", 1:ncol(mm))
  
  re.labs <- all.vars(formula)[-1]
  formula <- as.formula(paste0("y ~ ",
                                paste0(all.vars(formula)[-1], collapse=" + "), " + ", colnames(mm), " -1"))

  model.data <- data.frame(mm, mget(re.labs, as.environment(data)))
  
  if (save.residuals) {
    residuals <- matrix(NA, nrow(Y), ncol(Y), dimnames=dimnames(Y))
  }
  res <- lapply( 1:ncol(Y), function(i) {
    model.data$y <- Y[,i]
    model <- try(
      lm(formula, data=model.data, na.action=na.exclude),
      silent=TRUE )
    if (inherits(model, "try-error")) {
      return(NA)
    }
    
    return( list(residual = resid(model, type="response"),
                 name = colnames(Y)[i] ))
  })
  
  residuals = do.call(rbind,  lapply(res, "[[", "residual"))
  rownames(residuals) = do.call(rbind,  lapply(res, "[[", "name"))
  
  return( residuals)
}
