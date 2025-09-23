library(ggplot2)
library(ggrepel)
library(dplyr)


volcano_plot <- function(res_list=stab_m, prop_label= "Selection prop", output_path = ".", top_n = 20, filename) {
  
  #create df
  umbral_m <- pi_list[ArgmaxId(stab_m)[2]]
  mat_volcan <- data.frame(prop = as.numeric(selprop_m), beta_mean = as.numeric(beta_m))
  rownames(mat_volcan) <- names(selprop_m)
  mat_volcan <- mat_volcan %>%
    mutate(Associations = case_when(beta_mean > 0 & prop >= umbral_m ~ "Risk-factor",
                                    beta_mean < 0 & prop >= umbral_m ~ "Protective-factor",
                                    TRUE ~ "No assoc"))
  
  mat_volcan$Associations <- factor(mat_volcan$Associations, levels = c("No assoc", "Protective-factor", "Risk-factor"))
  
  # Selecció dels top N
  top_symbols <- bind_rows(mat_volcan %>% filter(Associations == 'Risk-factor') 
                           %>% arrange(desc(abs(beta_mean))) %>% head(top_n),
                           mat_volcan %>% filter(Associations == 'Protective-factor') 
                           %>% arrange(abs(beta_mean)) %>% head(top_n))
  #Simbols
  match_cpgs <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))[,c("Name", "UCSC_RefGene_Name")]
  match_cpgs$UCSC_RefGene_Name <- gsub("\\;.*$", "", match_cpgs$UCSC_RefGene_Name)
  match_cpgs$cg_Name <- paste(match_cpgs$Name, match_cpgs$UCSC_RefGene_Name, sep = "_")
  top_symbols$symbol <- match_cpgs[rownames(top_symbols),"cg_Name"]
  
  if(grepl("zdia t1", name)){
    lab <- "Beta (DBP Z-Score childhood)"}
  
  if(grepl("zsys t1", name)){
    lab <- "Beta (SBP Z-Score childhood)"}
  
  if(grepl("zdia t2", name)){
    lab <- "Beta (DBP Z-Score adolescence)"}
  
  if(grepl("zsys t2", name)){
    lab <- "Beta (SBP Z-Score adolescence)"}
  

  #plot
  plot <- ggplot(mat_volcan, aes(beta_mean, prop)) +
    geom_point(aes(color = Associations), size = 2) +
    xlab(lab) +
    ylab(prop_label) +
    scale_color_manual(values = c("gray50", "dodgerblue3", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    geom_label_repel(data = top_symbols,
                     aes(beta_mean, prop, label = symbol),
                     size = 3) +
    theme_minimal()
  
  ggsave(paste0("./results/EWAS/volcanoplot/",filename), plot, width = 10, height = 8, dpi = 150, units = "in", device = "svg")
  

  return(plot)
}
