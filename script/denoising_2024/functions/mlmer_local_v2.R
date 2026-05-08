mlmer_local <- function(
    formula,
    data,
    vars,
    train_ids = helixid_n,
    test_ids  = helixid_s,   
    save.residuals = TRUE
) {
  
  helixid_n <- c(read.csv("./helixid_n.csv", row.names = 1)$x)
  helixid_s <- c(read.csv("./helixid_s.csv", row.names = 1)$x)
  
  
  # ---------------------------
  # 1 — Train/test split
  # ---------------------------
  if (!is.null(train_ids)) {
    df_train <- data[data$helixid %in% train_ids, ]
  } else df_train <- data
  
  if (!is.null(test_ids)) {
    df_test <- data[data$helixid %in% test_ids, ]
  } else df_test <- NULL
  
  # ---------------------------
  # 2 — Extract Y
  # ---------------------------
  Yname <- omics:::response.name(formula, data = df_train)
  Y <- get(Yname, envir = environment(formula))
  
  Y_train <- Y[rownames(df_train), , drop=FALSE]
  Y_test  <- if(!is.null(df_test)) Y[rownames(df_test), , drop=FALSE] else NULL
  
  # ---------------------------
  # 3 — Fixed + Random terms
  # ---------------------------
  tmp <- formula
  tmp[[2]] <- NULL
  
  lf <- lme4::lFormula(
    tmp, df_train, REML=FALSE, na.action=na.pass,
    control=lme4::lmerControl(
      check.nobs.vs.nlev="ignore",
      check.nobs.vs.nRE="ignore",
      check.rankX="ignore",
      check.scaleX="ignore",
      check.formula.LHS="ignore"
    )
  )
  
  labs <- labels(terms(lme4::nobars(formula), data=df_train))
  if (missing(vars)) vars <- labs
  
  mm <- lf$X
  idx <- which(attr(lf$X, "assign") %in% match(vars, labs))
  vars <- colnames(mm)[idx]
  colnames(mm) <- sprintf("V%d", 1:ncol(mm))
  
  formula2 <- as.formula(sprintf(
    "y ~ %s - 1",
    paste0(c(sprintf("(%s)", lme4::findbars(formula)), colnames(mm)), collapse=" + ")
  ))
  
  re.labs <- as.character(attr(terms(lf$fr), "predvars.random")[-1])
  model.data <- data.frame(mm, mget(re.labs, as.environment(lf$fr)))
  
  # ---------------------------
  # 4 — Initialize residual matrix
  # ---------------------------
  all_rows <- c(rownames(df_train), if(!is.null(df_test)) rownames(df_test) else NULL)
  residuals_matrix <- matrix(
    NA,
    nrow = ncol(Y),
    ncol = length(all_rows),
    dimnames = list(colnames(Y), all_rows)
  )
  
  # ---------------------------
  # 5 — Fit models per column
  # ---------------------------
  for (i in seq_len(ncol(Y_train))) {
    
    model.data$y <- Y_train[, i]
    
    model <- try(
      lme4::lmer(formula2,
                 data = model.data,
                 REML = FALSE,
                 na.action = na.exclude),
      silent = TRUE
    )
    
    yname <- colnames(Y_train)[i]
    if (inherits(model, "try-error")) next
    
    # ---- TRAIN residuals ----
    r_train <- resid(model)
    residuals_matrix[yname, rownames(df_train)] <- r_train
    
    # ---- TEST residuals / predictions ----
    if (!is.null(df_test)) {
      pred_test <- try(
        predict(model, newdata = df_test, allow.new.levels = TRUE),
        silent = TRUE
      )
      if (!inherits(pred_test, "try-error")) {
        # If Y_test exists → compute true residuals
        if (!is.null(Y_test)) {
          r_test <- Y_test[, i] - pred_test
          residuals_matrix[yname, rownames(df_test)] <- r_test
        } else {
          # No Y_test → store predictions
          residuals_matrix[yname, rownames(df_test)] <- pred_test
        }
      }
    }
  }
  
  return(residuals_matrix)
}

