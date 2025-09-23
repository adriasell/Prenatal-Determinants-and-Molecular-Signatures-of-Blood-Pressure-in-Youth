library(RGCCA)
library(writexl)
library(tidyverse)
library(ggsankey)
library(networkD3)
library(webshot)
library(flextable)

#Function prepare_df (with latent vars, outcome and exposure data)--------------
prepare_df <- function(projections, X_data, exposure_data, response, ncomp=NULL){
  set.seed(1899)
  latent_vars <- projections %>% purrr::reduce(cbind)%>% as.data.frame()
  colnames(latent_vars)<- c(paste0(rep(names(X_data)[-response])))
  
  if (!is.null(ncomp)){
    colnames(latent_vars)<- c(paste0(rep(names(X_data)[-response],times = ncomp[-response])))
    colnames(latent_vars) <- paste0(colnames(latent_vars),unlist(lapply(ncomp[-response],seq)))
  }
  
  all(rownames(latent_vars)==rownames(X_data$Y))
  #Add outcome to latent vars
  latent_variables_with_outcome <-latent_vars %>%  cbind(X_data$Y)
  colnames(latent_variables_with_outcome)<- c(colnames(latent_vars),colnames(X_data$Y))

  #Add exposure to latent vars
  latent_variables_with_outcome <- latent_variables_with_outcome[order(rownames(latent_variables_with_outcome)), ]
  exposure_data <- exposure_data[order(rownames(exposure_data)), ]
  all(rownames(exposure_data)==rownames(latent_variables_with_outcome))
  latent_variables_with_exposures <- cbind(latent_variables_with_outcome, exposure_data)

  return(latent_variables_with_exposures)
}

#ExWAS(Ines)--------------------------------------------------------------------------

ExWAS_mixed <- function(data,expos_name,form,correction="BH",M0=NULL) {
  # 4 - ExWAS core function : run one regression model by exposure and 
  #     return the table of results (coeff, CI, raw p values)
  # Inputs
  # - data = imputed datasets (exposures, covariates, outcomes)
  # - expos_name = list of exposures to consider
  # - form = string "outcome ~ covariates"
  # Output
  # - table with beta, CI, raw p value and corrected p value for each regression
  
  ExWAS_parralel_mixed = function(expo) {
    j=1
    # Formula : outcome ~ expo + covariates 
    form_exwas <- paste(form, "+", expo)
    n = length(expo)
    # Regression with lm on each imputed data (with()) + Rubin's rule to
    # pool the results (pool()) + summary to get the coefficients
    a <-summary(lm(as.formula(form_exwas),data=data))$coefficients %>% as.data.frame()
    a$term <- rownames(a)
    ## Fill of the result table
    data_expo = data[,expo]
    
    # Categorical variables with more than 2 modality: k-1 beta and p value
    if (!is.numeric(data_expo) & length(table(data[, expo])) > 2) {
      # Checking that the variable had category with n=0 individuals
      if (any(table(data[,expo])==0)){
        print(paste0("Category with no data for variable ",expo," (removed)"))
        res <- data.frame(variable=expo,modality="",`beta (CI 95%)`=NA,
                          beta=NA,sd=NA,`CI 2.5`=NA, `CI 97.5`=NA)
        # We begin from the second modality (first = ref)
      } else {
        res <- matrix(NA, ncol = 8, nrow = (length(table(data_expo))-1)) %>%
          as.data.frame()
        colnames(res) <-  c("variable","modality", "beta (CI 95%)","p","beta", "sd", "CI 2.5","CI 97.5")
        cat_names <- rownames(table(data[, expo]))[-1]
        # And iterate until the last modality
        for (modality in cat_names) {
          name_modality = paste(expo,modality,sep="")
          beta_expo <- a[a$term == name_modality , "Estimate"] #beta
          sd_expo <- a[a$term == name_modality , "Std. Error"] #std.error
          p_value <- signif(a[a$term == name_modality , "Pr(>|t|)"],2) # p value
          IC_low <- as.numeric(beta_expo)-1.96*as.numeric(sd_expo) #CI 2.5%
          IC_high <- as.numeric(beta_expo)+1.96*as.numeric(sd_expo) #CI 75%
          beta_IC <- paste0(signif(as.numeric(beta_expo),2)," [",
                            signif(as.numeric(IC_low),2),"; ",
                            signif(as.numeric(IC_high),2),"]") # beta (CI)
          modality_expo = modality
          res[j,] <- res[j,]  %>%
            mutate(variable = expo,
                   modality = modality_expo,
                   `beta (CI 95%)`= beta_IC,
                   p = as.numeric(p_value),
                   beta =  signif(beta_expo,2),
                   sd =  signif(sd_expo,2),
                   `CI 2.5`= signif(IC_low,2),
                   `CI 97.5`= signif(IC_high,2))
          
          j <- j + 1
        }
      }
      # Continous variables or categorical with 2 modality: only one beta and p value
    } else {
      res <- matrix(NA, ncol = 8, nrow = 1) %>% as.data.frame()
      colnames(res) <-  c("variable","modality", "beta (CI 95%)","p","beta", "sd", "CI 2.5","CI 97.5")
      
      # Categorical variable with 2 modality
      if (!is.numeric(data_expo)) {
        modality_expo = rownames(table(data[, expo]))[2]
        expo_raw = expo
        expo_name = paste0(expo, " - ",modality_expo)
        expo <- paste(expo, modality_expo, sep = "")
        # Continuous variable
      } else { 
        expo_raw = expo
        modality_expo = ""
      }
      beta_expo <- a[a$term == expo , "Estimate"] #beta
      sd_expo <- a[a$term == expo , "Std. Error"] #std.error
      IC_low <- as.numeric(beta_expo)-1.96*as.numeric(sd_expo) #CI 2.5%
      IC_high <- as.numeric(beta_expo)+1.96*as.numeric(sd_expo) #CI 75%
      p_value <- signif(a[a$term == expo , "Pr(>|t|)"],2) # p value
      beta_IC <- paste0(signif(as.numeric(beta_expo),2)," [",
                        signif(as.numeric(IC_low),2),"; ",
                        signif(as.numeric(IC_high),2),"]") # beta (CI)
      res[j,] <- res[j,]  %>%
        mutate(variable = expo_raw,
               modality = modality_expo,
               `beta (CI 95%)`= beta_IC,
               p = as.numeric(p_value),
               beta =  signif(beta_expo,2),
               sd =  signif(sd_expo,2),
               `CI 2.5`= signif(IC_low,2),
               `CI 97.5`= signif(IC_high,2))
      
      j <- j + 1
    }
    return(as.data.frame(res))
  }
  
  # RUN LINEAR REGRESSION FUNCTION ON EACH EXPO
  n_expo = length(expos_name)
  res <- foreach::foreach(i = 1:n_expo,.packages = c("dplyr")) %dopar%{
    return(ExWAS_parralel_mixed(expo=expos_name[i]))}
  res <- suppressMessages(purrr::reduce(res,full_join))
  
  res <- res |> mutate(across(c(beta,p,sd,`CI 2.5`,`CI 97.5`), ~as.numeric(.x)))
  
  # Correction on multiple testing
  if (correction=="BH"){
    res$p_corrected <- p.adjust(res$p,method="BH")
  } else if (correction=="Li"){
    if (is.null(M0)){print("Correction has not been performed as M0 was not given")}else{
      res <- res |> mutate(
        p = as.numeric(p),
        p_corrected = ifelse(p*M0>1,1,p*M0))
    }
  }
  
  return(res)
}



cor_mixed <- function(data,method="spearman"){
  ## Purpose: 2 -  use the function mixedCor from psych to calculate a 
  #                correlation matrix on mixed type variables (continuous and categorical).
  #                Note: mixedCor consider categorical variable as ordered factors.
  ## Inputs: - data: dataframe (not mice) with only the variables to consider 
  ##                 (mixed type allowed)
  ## Output: - correlation matrix
  
  ## STEP 1 : TYPE OF VARIABLES
  continuous_var = which(sapply(data, class) == "numeric")
  names(continuous_var)=NULL
  # Categorical var
  categorical_var = which(sapply(data, class) == "factor")
  #  - 2 levels only
  binary_var = categorical_var[sapply(data[,categorical_var],nlevels)==2] 
  binary_var = binary_var[!is.na(binary_var)]
  names(binary_var)=NULL
  #  - More than 2 levels (but less than 8)
  poly_var = categorical_var[(sapply(data[,categorical_var],nlevels)>2 & sapply(data[,categorical_var],nlevels)<8)] %>% na.exclude()
  names(poly_var)=NULL
  
  ## STEP 2 : CORRELATION MATRIX USING MIXEDCOR FUNCTION (FROM PSYCH)
  # data converted in numeric (necessary)
  data[,] = lapply(data[,],as.numeric)
  # Correlation matrix
  cor = data  %>% 
    mixedCor(c=continuous_var,p=poly_var,d=binary_var,use="pairwise.complete.obs",method=method)%>% pluck('rho')
  return(cor)
}

#### SANKEY PLOTS ####
plotSankey <- function(res_outcome,res_exposure,name_outcomes, p_val=0.05,
                       variable_labels, label_components = NULL, path=path, color=T){
  # res_outcome = ExWAS between outcome (last layer) and omics (middle layers)
  # res_exposure = ExWAS between omics (middle layer) and the exposures (first layer)
  # name_outcome = name of the outcome (last layer)
  # variable_labels = dataframe with the labels associated with each exposure (first layer) and outcome (last layer), and their group
  # label_components = dataframe with the labels associated with each component
  # path = where to save the data
  # Data 
  res_outcome_reshaped <- res_outcome %>%
    # dplyr::mutate(outcome=name_outcome) %>%
    filter(outcome%in% name_outcomes)%>%
    dplyr::rename(beta_comp_outcome=beta,Comp=variable,p_comp_outcome=p,p_comp_outcome_adj=p_corrected) %>%
    dplyr::select(Comp,outcome,beta_comp_outcome,p_comp_outcome, p_comp_outcome_adj) %>%
    dplyr::left_join(variable_labels, by = join_by("outcome"=="Variable_name_TRANS")) %>%
    rename(label_outcomes = Label.short..e.g..for.figures.) %>%
    select(-c("Subgroup","Group","Period"))
  res_exposure_reshaped <- res_exposure %>%
    dplyr::rename(beta_expos_comp=beta,Exposure=variable,p_expos_comp=p,Comp=comp, p_comp_expos_adj=p_corrected) %>%
    dplyr::left_join(variable_labels, by = join_by("Exposure"=="Variable_name_TRANS")) %>%
    dplyr::select(Exposure,modality,Comp,beta_expos_comp,p_expos_comp,Label.short..e.g..for.figures.,Group,Period,p_comp_expos_adj) %>%
    rename(label_exposures = Label.short..e.g..for.figures.) 
  
  res_total_unfil <- full_join(res_exposure_reshaped,res_outcome_reshaped, relationship = "many-to-many")
  write.csv2(res_total_unfil, paste0(path,"/sankeyresults_unfil.csv"))
 
  res_total <- full_join(res_exposure_reshaped,res_outcome_reshaped, relationship = "many-to-many")
  
  res_total <- res_total %>%  dplyr::filter(p_expos_comp < p_val & p_comp_outcome < p_val)

  if (!is.null(label_components)){
    res_total <- left_join(res_total,label_components) %>%
      mutate(Comp = coalesce(label_comp,Comp),.keep="unused")
  }
  res_total <- res_total %>%
    dplyr::mutate(
      Group2 = Comp) %>%
    group_by(Exposure)%>%
    reframe(Group2 = rep(dplyr::first(Group2),n()), Comp = Comp, beta_expos_comp = beta_expos_comp,
            p_expos_comp = p_expos_comp, p_comp_expos_adj=p_comp_expos_adj, label_exposures = label_exposures, Group =Group, 
            Period = Period, outcome = outcome, beta_comp_outcome = beta_comp_outcome,
            p_comp_outcome = p_comp_outcome,  p_comp_outcome_adj = p_comp_outcome_adj, label_outcomes = label_outcomes)
  
  ## Nodes = names of nodes
  Nodes <- data.frame(
    name = c(res_total$label_exposures,
             res_total$Comp,
             res_total$label_outcomes),
    node.group = c(paste0(res_total$Group,res_total$Group2),
                   res_total$Comp,res_total$outcome))%>%
    arrange(node.group) %>%
    distinct()
  ## Links: data with (x,node) = (type of 1st node,name of node)
  ##                  (next_x,next_node) = same for the node following
  ##                  value = link value
  Links <- ggsankey::make_long(res_total,
                               label_exposures, 
                               Comp, 
                               label_outcomes, 
                               value = "beta_expos_comp")
  Links$value[seq(2, nrow(Links), by = 3)] <- res_total$beta_comp_outcome
  Links$value[seq(3, nrow(Links), by = 3)] <- NA_real_
  Links <- Links %>%
  mutate(
    source = match(node, Nodes$name) - 1,
    target = match(next_node, Nodes$name) - 1, 
    Value = 1 # abs(value)
  ) %>%
  mutate(
    link.group = case_when(
      value < 0 & !grepl("BP SD score", next_node) ~ "negative",
      value > 0 & !grepl("BP SD score", next_node) ~ "positive",
      value == 0 ~ "categorical",
      value < 0 & grepl("BP SD score", next_node) ~ "negative_out",
      value > 0 & grepl("BP SD score", next_node) ~ "positive_out",
      is.na(value) ~ "NA"
    )
  ) %>%
  filter(!is.na(value)) %>%
  as.data.frame() %>%
  group_by(x, node, next_x, next_node, value, source, target, Value, link.group) %>%
  dplyr::summarize(Value = n()) # first(Value)) %>%
  # group_by(next_x)%>%
  # dplyr::reframe(x=x,node=node,next_x=next_x,source=source,
  #                  target=target,link.group=link.group)#,
  #                  Value=(Value)/sum(Value)) 
  # 
  if (color){
    my_color <- 'd3.scaleOrdinal()
    .domain(["negative", "positive", "negative_out", "positive_out", "categorical", "NA"])
    .range(["#70A9D1", "#D4B96C", "olivedrab", "coral","#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080", "#808080", "#808080","#808080","#808080", "#808080","#808080","#808080", "#808080","#808080","#808080", "#808080","#808080", "#808080", "#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080"])'
    } else {
    my_color <- 'd3.scaleOrdinal() .domain(["negative", "positive","categorical","NA"]) .range(["white","white","white","white","white","white","white","white","white","white","white","white","white","white"])'
    
    }

  # create the plot:
  p <- networkD3::sankeyNetwork(
    Links = (Links), Source = "source", Target = "target", Value = "Value",
    Nodes = Nodes, NodeID = "name", 
    LinkGroup = "link.group", NodeGroup = "node.group",colourScale = my_color, 
    fontSize = 20, fontFamily = "serif",
    nodeWidth = 20,
    sinksRight = F, iterations = 0) 
  

    print(p)
  write.csv2(res_total, paste0(path,"/sankeyresults_",p_val,".csv"))
  
  flextable(res_total) %>%
    flextable::theme_vanilla() %>% 
    flextable::bg(bg = "white", part = "all") %>% 
    flextable::save_as_image(paste0(path, "/sankeyresults_",p_val,".png"))
  
  name_sankey <- paste0("/sankeyplot",p_val,".html")
  networkD3::saveNetwork(p, paste0(path,name_sankey))
  
  # 
  # df <- Links %>%
  #   left_join(Nodes,join_by("node"=="name")) %>%
  #   mutate(node = reorder(node,as.numeric(as.factor(node.group))))
  # 
  # ggplot(df, aes(x = x, 
  #                next_x = next_x, 
  #                node = node, 
  #                next_node = next_node,
  #                fill = link.group,
  #                label = node,
  #                color = node.group#,
  #               # value=value
  #               
  #              )) +
  #   geom_sankey(color="white")+#,width = 0.2,linewidth=10) +
  #   #geom_sankey_text(aes(group = node.group, label = node,color=node.group), size = 3, color = "black",hjust = 0, position = position_nudge(x = 0.1))+
  #   geom_sankey_label(fill="white",size=5) +
  #   scale_fill_manual(values=c("black","coral","olivedrab"))+
  #   theme_sankey(base_line_size = 10,base_rect_size = 0.1)
  # 
  # ggsave(paste0(path,"/sankeyplot_ggsankey.png"),
  #        height=7,width=20)
  
  return(res_total)
}

  