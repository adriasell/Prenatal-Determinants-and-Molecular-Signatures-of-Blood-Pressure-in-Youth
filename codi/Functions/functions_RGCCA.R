library(tidyverse)
library(RGCCA)
library(Biobase)
library(dplyr)
library(parallel)
library(writexl)

#apply_RGCCA--------------------------------------------------------------------
apply_RGCCA <- function(rgcca_mod = NULL, 
                        X_data, 
                        X_data_pred, 
                        connection, 
                        response, 
                        ncomp, 
                        sparsity, 
                        approach,
                        output_dir="RGCCA"){
  
  #RGCCA model
  set.seed(1999)
  if(is.null(rgcca_mod)){
  rgcca_model <- rgcca(blocks = X_data, 
                       connection = connection, 
                       method = "sgcca", 
                       response = response, 
                       ncomp = ncomp, 
                       sparsity = sparsity)
  #Save model
  filename <- paste0("./results/",output_dir,"/model/rgcca_", approach, ".rds")
  write_rds(rgcca_model, filename)
  } else{
    rgcca_model <- rgcca_mod
  }
  
  
  #Bootstrap
  bootstrap_res <- rgcca_bootstrap(rgcca_model, n_boot = 1000, n_cores = detectCores())
  
  RGCCA_mod <- list("sparsity" = rgcca_model$call$sparsity,
                    "ncomp" = rgcca_model$call$ncomp,
                    "Latent factor" = rgcca_model$Y,
                    "Block weights for block" = rgcca_model$a,
                    "Bootstrap results" = bootstrap_res)
  
####Prediction with X_data######################################################
  #Prediction
  rgcca_predict_res <- rgcca_predict(rgcca_model,blocks_test = X_data, prediction_model = "lm")
  latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind) %>% as.data.frame()
  
  #Calculate R2
  if (any(ncomp>1)){
    name_components <- paste0(rep(names(X_data)[-response],times=ncomp[-response]),unlist(lapply(ncomp[-response],seq)))
  } else {
    name_components <- paste0(rep(names(X_data)[-response],times=ncomp[-response]))
  }
  
  colnames(latent_variables) <- name_components
  names_out <- colnames(X_data$Y)
  results_R2 <- calculate_R2(latent_variables, X_data$Y, names_out, name_components)
  
  #Calculate correlation
  corr_matrix <- cor(latent_variables) %>% round(2)
  
  rgcca_xdat <- list("R2" = results_R2,
                     "metrics" = rgcca_predict_res$metric,
                     "corr" = corr_matrix,
                     "projection" = rgcca_predict_res$projection)
  
  #List of metrics
  list_result <- list("RGCCA model"= RGCCA_mod,
                      "RGCCA pred (X_data)"= rgcca_xdat)

    ####Prediction with X_data_pred################################################
  if (!identical(X_data, X_data_pred)){
    rgcca_predict_res_1 <- rgcca_predict(rgcca_model,blocks_test = X_data_pred, prediction_model = "lm")
    latent_variables_1 <- rgcca_predict_res_1$projection %>% purrr::reduce(cbind) %>% as.data.frame()
  
      #Calculate R2
    if (any(ncomp>1)){
      name_components <- paste0(rep(names(X_data)[-response],times=ncomp[-response]),unlist(lapply(ncomp[-response],seq)))
    } else {
      name_components <- paste0(rep(names(X_data)[-response],times=ncomp[-response]))
    }
    colnames(latent_variables_1) <- name_components
    names_out <- colnames(X_data_pred$Y)
    results_R2_1<- calculate_R2(latent_variables_1, X_data_pred$Y, names_out, name_components)
    
    #Calculate correlation
    corr_matrix_1 <- cor(latent_variables_1) %>% round(2)
    
    rgcca_xpred <- list("R2" = results_R2_1,
                        "metrics" = rgcca_predict_res_1$metric,
                        "corr" = corr_matrix_1,
                        "projection" = rgcca_predict_res_1$projection)
    #List of metrics
    list_result <- c(list_result, list("RGCCA pred (X_pred)"=rgcca_xpred))
  }
  

  
  filename1 <-  paste0("./results/", output_dir,"/model_results/results_rgcca_", approach, ".rds")
  write_rds(list_result, file = filename1)
  
  return(list_result)
}

#Function read data-------------------------------------------------------------
eset2df <- function(eset_path) {
  eset <- readRDS(eset_path)
  if(class(eset)=="GenomicRatioSet"){
    df <- data.frame(t(getBeta(eset)))
  } else{
    df <- data.frame(t(exprs(eset)))
  }
  df <- df[complete.cases(df), ]
  df$HelixID <- pData(eset)$HelixID.x
  rownames(df) <- df$HelixID
  return(df)
}

#Function calculate R2----------------------------------------------------------
calculate_R2 <- function(latent_vars, Y_data, outcomes, components) {
  results_R2 <- data.frame(Outcome = character(), Layer = character(), R2 = numeric())
  for (name in outcomes) {
    for (comp in components) {
      data_r2 <- cbind(latent_vars, Y_data) %>% as.data.frame()
      form <- as.formula(paste("Y_data$", name, "~", comp))
      r2 <- summary(lm(form, data = data_r2))$r.squared * 100
      results_R2 <- rbind(results_R2, data.frame(Outcome = name, Layer = comp, R2 = r2))
    }
  }
  return(results_R2)
}

#Function performance RGCCA-----------------------------------------------------
performance_RGCCA_cv <- function(rgcca_model,
                                 X_test=X,
                                 response=4,
                                 sparsity ,
                                 ncomp ,
                                 connection,
                                 metrics=c("Rsquared","RMSE","MAE"),
                                 Nfold=5,
                                 n_run=5,
                                 ncore.use=ncore.use){
  set.seed(1999)
  
  n_block = length(X_test)
  # Names of components (un peu barbare mais j'ai pas su simplifier)
  
  if (any(ncomp>1)){
    name_components <- paste0(rep(names(X_test)[-response],times=ncomp[-response]),unlist(lapply(ncomp[-response],seq)))
  } else {
    name_components <- paste0(rep(names(X_test)[-response],times=ncomp[-response]))
  }
  outcomes<-names(X_test$Y)
  
  quality_cv_test = matrix(nrow = Nfold,
                           ncol=length(metrics)+2*length(name_components),
                           data=NA,
                           dimnames=list(1:(Nfold),c(metrics,paste0("cor_",name_components),paste0("R2_",name_components))))
  
  quality_cv_test_dia1 <- quality_cv_test_sys1 <- quality_cv_test_dia2 <- quality_cv_test_sys2 <- quality_cv_test
  
  # Function for 1 k-fold validation (to rep n_run time, parrallelized)
  cross_validation_rep_i <- function(rep){
    id_cv <- sample(1:Nfold,nrow(X_test[[1]]),replace=TRUE)
    
    for (i in 1:Nfold){
      n <- i
      print(paste0("rep=",rep,", i=",i))
      # Generate train and test datasets
      X_test_i = list()
      for (k in 1:n_block){
        if (k!=response){
          X_test_i[[names(X_test[k])]] <- X_test[[k]][id_cv==i,]
        }else{
          X_test_i[[names(X_test[k])]] <- X_test[[k]][id_cv==i,]
        }
        
      }
      
      # Run RGCCA
      results_i <- rgcca_res
      
      # Estimate prediction error on (K-1) block = train, and the Kth block = test
      pred_quality <- rgcca_predict(results_i,blocks_test=X_test_i,prediction_model="lm")
      # Estimate correlation between latent variables and the outcome
      latent_variables_test <- pred_quality$projection %>% purrr::reduce(cbind)
      colnames(latent_variables_test)  <- name_components
      cor_values<-cor(latent_variables_test,X_test_i[[response]])%>% as.data.frame()%>% signif(2) %>% t()
      
      for (outcome in outcomes){
        for (metric in metrics){
          if (outcome=="hs2_zdia_bp.v3_2017_Time1"){
            quality_cv_test_dia1[n,metric] <- (pred_quality$metric$test[metric,outcome])
            quality_cv_test_dia1[n,paste0("cor_",name_components)] <- cor_values[outcome,]
          } 
          if (outcome=="hs2_zsys_bp.v3_2017_Time1"){
            quality_cv_test_sys1[n,metric] <- (pred_quality$metric$test[metric,outcome])
            quality_cv_test_sys1[n,paste0("cor_",name_components)] <- cor_values[outcome,]
            
          }
          if (outcome=="hs2_zdia_bp.v3_2017_Time2"){
            quality_cv_test_dia2[n,metric] <- (pred_quality$metric$test[metric,outcome])
            quality_cv_test_dia2[n,paste0("cor_",name_components)] <- cor_values[outcome,]
            
          }
          if (outcome=="hs2_zsys_bp.v3_2017_Time2"){
            quality_cv_test_sys2[n,metric] <- (pred_quality$metric$test[metric,outcome])
            quality_cv_test_sys2[n,paste0("cor_",name_components)] <- cor_values[outcome,]
          }}}
      
      
      # Estimate R2 of each latent variable by runing a linear model with only this component
      data_r2 <- cbind(latent_variables_test,X_test_i[[response]]) %>% as.data.frame()
      
      for (comp in name_components){
        for (outcome in outcomes){
          form = paste0(outcome," ~", comp)
          res <- summary(lm(as.formula(form),data=data_r2))
          
          if (outcome=="hs2_zdia_bp.v3_2017_Time1"){
            quality_cv_test_dia1[n,paste0("R2_",comp)] <- res$r.squared
          } 
          
          if (outcome=="hs2_zsys_bp.v3_2017_Time1"){
            quality_cv_test_sys1[n,paste0("R2_",comp)] <- res$r.squared
          }
          
          if (outcome=="hs2_zdia_bp.v3_2017_Time2"){
            quality_cv_test_dia2[n,paste0("R2_",comp)] <- res$r.squared
          }
          
          if (outcome=="hs2_zsys_bp.v3_2017_Time2"){
            quality_cv_test_sys2[n,paste0("R2_",comp)] <- res$r.squared
          }}}
      
    }
    list_res <- list("Dia1"=quality_cv_test_dia1,
                     "Sys1"=quality_cv_test_sys1,
                     "Dia2"=quality_cv_test_dia2,
                     "Sys2"=quality_cv_test_sys2)
    
    return(list_res)
  }
  res<-list()
  for (i in 1:n_run) { 
    res[[i]] <- cross_validation_rep_i(rep = i) 
  }
  res1 <- lapply(names(res[[1]]), function(outcome) {lapply(res, function(x) x[[outcome]])})
  names(res1) <- names(res[[1]])
  
  res_dia1 <- purrr::reduce(res1[[1]],rbind)
  res_sys1 <- purrr::reduce(res1[[2]],rbind)
  res_dia2 <- purrr::reduce(res1[[3]],rbind)
  res_sys2 <- purrr::reduce(res1[[4]],rbind)
  
  ## Make plot
  # Long version of the results
  list_res <- list("res_dia1"=res_dia1, "res_sys1"=res_sys1, "res_dia2"=res_dia2, "res_sys2"=res_sys2)
  
  result_list <- lapply(names(list_res), function(name) {
    res_long <- list_res[[name]] %>%
      as.data.frame() %>%
      pivot_longer(cols = colnames(list_res[[name]]), names_to = "Indicator", values_to = "value") %>%
      mutate(
        set = factor(Indicator),
        group = case_when(
          Indicator %in% metrics ~ Indicator,
          Indicator %in% paste0("cor_", name_components) ~ "Correlation with the outcome",
          Indicator %in% paste0("R2_", name_components) ~ "R2 of each component",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(Indicator != "Correlation with the outcome")  # # Line at zero for the correlations
    
    # Medians of each group   
    medians <- res_long %>%
      group_by(Indicator, group) %>%
      dplyr::summarise(md = median(value), .groups = "drop") %>%
      dplyr::mutate(md = round(md, 3))
    
    # Multiply by -1 for negative correlation (all components positively associated with outcome)   
    res_long <- res_long %>%
      left_join(medians, by = c("Indicator", "group")) %>%
      mutate(value = ifelse(Indicator %in% paste0("cor_", name_components), value * sign(md), value))
    
    # Extract R2 values   
    r2 <- res_long %>% filter(Indicator == "Rsquared")
    r2_block <- res_long %>% filter(Indicator %in% paste0("R2_",name_components))
    
    # Final plots   
    final_plot <- ggplot() +
      geom_boxplot(data = res_long, aes(x = Indicator, y = value)) +
      #geom_hline(data = line_zero,aes(yintercept = yint),lty=2) + 
      ggforce::facet_row(vars(group), scales = 'free', space = 'free') +
      ggrepel::geom_text_repel(data = medians, aes(x = Indicator, y = abs(md), label = abs(md)), direction = "y") +
      xlab("Cross validation") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    return(list(plot = final_plot,
                median_quality = medians,
                quality_all = res_long,
                r2 = r2,
                r2_block = r2_block))
  })
  
  names(result_list) <- names(list_res)
  
  return(result_list)
}
