##################################################
##################################################
##########Project: HELIX study - childhood (mxitures and BMI)
##########Creator: NuRia Güil-Oumrait
##########Date: 15/03/24
##########Last update: 5/9/24
##########Aim: groupBWQS (prenatal mixt and adolescence BMI) ADJUSTED MODEL. EXAMPLE FOR BETHANY KNOX. 
##################################################
###Install/Upload these packages 
library("rstan")     
library('ggplot2')
library('BWQS')
library('mvtnorm')
library('tidyverse')
library("fastDummies")
library("foreach")
library("doParallel")
library("writexl")
library("dplyr")
#install_github("https://github.com/ElenaColicino/bwqs.git")
##In my example, I have 10 mixture groups. If you have less, delete the lines of the elements with the numbers you won't need. 
#For instance, if you have 8 mixture groups delete line #33, #34, #47, #48, #59, #60, #85, #86, and in line #94, delete:" + (XC9*WC9)*beta[9]+ (XC10*WC10)*beta[10]) ". Delente lines 111, 112, 126, 127

setwd("./TFM")
getwd()
set.seed(1999)

#Load projections
projections <- readRDS("./results/RGCCA/projections.rds")

proj_north <- projections$`Projections N` %>% purrr::reduce(cbind) 
colnames(proj_north) <- names(projections$`Projections N`)
proj_south <- projections$`Projections S` %>% purrr::reduce(cbind)
colnames(proj_south) <- names(projections$`Projections S`)
proj_all <-  projections$`Projections all` %>% purrr::reduce(cbind)
colnames(proj_all) <- names(projections$`Projections all`)

#Load exposome
exposome <- read.csv2("./results/ExWAS/exposome_filtered.csv", row.names = 1)
#inverse gestational age
exposome$e3_gac_None <- exposome$e3_gac_None*-1

ids<-rownames(exposome)
to_dummy <- c("h_edumc_None", "h_fish_preg_Ter",  "h_fruit_preg_Ter", "h_legume_preg_Ter",
              "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")

#Scale exposures
exposome[, !colnames(exposome) %in% c(to_dummy,"e3_asmokyn_p_None", "e3_alcpreg_yn_None")] <- scale(exposome[, !colnames(exposome) %in% c(to_dummy,"e3_asmokyn_p_None", "e3_alcpreg_yn_None", "h_greenyn300_preg_None")], center = T)

#exposome <- dummy_cols(exposome,
#                       select_columns = to_dummy,
#                       remove_selected_columns = T,
#                       remove_first_dummy = T)

rownames(exposome)<-ids

exposome_all <- exposome[rownames(exposome) %in% rownames(proj_all),]
exposome_n <- exposome[rownames(exposome) %in% rownames(proj_north),]
exposome_s <- exposome[rownames(exposome) %in% rownames(proj_south),]

#edumc_dummies <- colnames(exposome)[grepl("h_edumc_None", colnames(exposome))]
#diet_dummies <- colnames(exposome)[grepl(paste(to_dummy[to_dummy!="h_edumc_None"], collapse = "|"), colnames(exposome))]

#Load phenotype
phenotype <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")

phenotype <- dummy_cols(phenotype,
                        select_columns = "h_cohort",
                        remove_selected_columns = T,
                        remove_first_dummy = F)

rownames(phenotype) <- phenotype$HelixID
phenotype$e3_sex_Time1 <- as.numeric(phenotype$e3_sex_Time1)

phenotype_all<-phenotype[rownames(phenotype) %in% rownames(proj_all),]
phenotype_n<-phenotype[rownames(phenotype) %in% rownames(proj_north),]
phenotype_s<-phenotype[rownames(phenotype) %in% rownames(proj_south),]

all(rownames(phenotype_all)==rownames(exposome_all))
all(rownames(phenotype_all)==rownames(proj_all))
all(rownames(phenotype_n)==rownames(exposome_n))
all(rownames(phenotype_n)==rownames(proj_north))
all(rownames(phenotype_s)==rownames(exposome_s))
all(rownames(phenotype_s)==rownames(proj_south))

df_all<-cbind(proj_all, exposome_all, phenotype_all)
df_n<-cbind(proj_north, exposome_n, phenotype_n)
df_s<-cbind(proj_south, exposome_s, phenotype_s)


model_bwqs_gaussian_group_lasso <- "data {                                                         

int<lower=0> N;              // Sample size

int<lower=0> C1;             // number of element of first mix
int<lower=0> C2;             // number of element of second mix
int<lower=0> C3;             // number of element of third mix
int<lower=0> C4;             // number of element of fourth mix
int<lower=0> C5;             // number of element of fifth mix
int<lower=0> C6;             // number of element of six mix
int<lower=0> C7;             // number of element of seven mix
int<lower=0> C8;             // number of element of eight mix
int<lower=0> C9;             // number of element of ninth mix
int<lower=0> C10;             // number of element of ninth mix



int<lower=0> K;              // number of covariates
int<lower=0> G;              //  number of Groups

matrix[N,C1] XC1;	           // matrix of first mix
matrix[N,C2] XC2;	           // matrix of second mix
matrix[N,C3] XC3;	           // matrix of third mix
matrix[N,C4] XC4;	           // matrix of fourth mix
matrix[N,C5] XC5;	           // matrix of fifth mix
matrix[N,C6] XC6;	           // matrix of six mix
matrix[N,C7] XC7;	           // matrix of seven mix
matrix[N,C8] XC8;	           // matrix of eight mix
matrix[N,C9] XC9;	           // matrix of ninth mix
matrix[N,C10] XC10;	           // matrix of ninth mix


vector[C1] DalpC1;           // vector of the Dirichlet coefficients for first mix
vector[C2] DalpC2;           // vector of the Dirichlet coefficients for second mix
vector[C3] DalpC3;           // vector of the Dirichlet coefficients for third mix
vector[C4] DalpC4;           // vector of the Dirichlet coefficients for fourth mix
vector[C5] DalpC5;           // vector of the Dirichlet coefficients for fifth mix
vector[C6] DalpC6;           // vector of the Dirichlet coefficients for sixth mix
vector[C7] DalpC7;           // vector of the Dirichlet coefficients for seven mix
vector[C8] DalpC8;           // vector of the Dirichlet coefficients for eight mix
vector[C9] DalpC9;           // vector of the Dirichlet coefficients for ninth mix
vector[C10] DalpC10;           // vector of the Dirichlet coefficients for ninth mix


matrix[N,K] KV;	             // matrix of covariates
real y[N];                   // outcome gaussian variable
}

parameters {

real <lower=0> sigma;
real mu;               // intercepts

vector[G] beta;     // indiv coeffs by group
vector[K] delta;         // covariates coefficients

vector<lower=0>[G] tau_squared;     
real<lower=0> lambda_squared;         // penalization factor

simplex[C1] WC1;          // weights of first mix
simplex[C2] WC2;          // weights of second mix
simplex[C3] WC3;          // weights of second mix
simplex[C4] WC4;          // weights of second mix
simplex[C5] WC5;          // weights of second mix
simplex[C6] WC6;          // weights of second mix
simplex[C7] WC7;          // weights of second mix
simplex[C8] WC8;          // weights of second mix
simplex[C9] WC9;          // weights of second mix
simplex[C10] WC10;          // weights of second mix




}
transformed parameters {

vector[N] Xb;

Xb = (mu + (XC1*WC1)*beta[1] + (XC2*WC2)*beta[2] + (XC3*WC3)*beta[3] + (XC4*WC4)*beta[4] + (XC5*WC5)*beta[5]
                 + (XC6*WC6)*beta[6]+ (XC7*WC7)*beta[7]+ (XC8*WC8)*beta[8]+ (XC9*WC9)*beta[9]+(XC10*WC10)*beta[10]) + KV*delta;    
}
model {

mu ~ normal(0, 100);

sigma ~ inv_gamma(0.01,0.01);

lambda_squared ~ gamma(2,0.5);

tau_squared[1] ~ gamma((C1+1)/2.0, (lambda_squared)/2);
tau_squared[2] ~ gamma((C2+1)/2.0, (lambda_squared)/2);
tau_squared[3] ~ gamma((C3+1)/2.0, (lambda_squared)/2);
tau_squared[4] ~ gamma((C4+1)/2.0, (lambda_squared)/2);
tau_squared[5] ~ gamma((C5+1)/2.0, (lambda_squared)/2);
tau_squared[6] ~ gamma((C6+1)/2.0, (lambda_squared)/2);
tau_squared[7] ~ gamma((C7+1)/2.0, (lambda_squared)/2);
tau_squared[8] ~ gamma((C8+1)/2.0, (lambda_squared)/2);
tau_squared[9] ~ gamma((C9+1)/2.0, (lambda_squared)/2);
tau_squared[10] ~ gamma((C10+1)/2.0, (lambda_squared)/2);

beta ~ multi_normal(rep_vector(0,G),diag_matrix(tau_squared));
delta ~ multi_normal(rep_vector(0,K),diag_matrix(rep_vector(100, K)));

WC1 ~ dirichlet(DalpC1);
WC2 ~ dirichlet(DalpC2);
WC3 ~ dirichlet(DalpC3);
WC4 ~ dirichlet(DalpC4);
WC5 ~ dirichlet(DalpC5);
WC6 ~ dirichlet(DalpC6);
WC7 ~ dirichlet(DalpC7);
WC8 ~ dirichlet(DalpC8);
WC9 ~ dirichlet(DalpC9);
WC10 ~ dirichlet(DalpC10);


y ~  normal(Xb, sigma);

}

"
m_lasso <- rstan::stan_model(model_code =  model_bwqs_gaussian_group_lasso)

#Make a new dataset with your exposures variables (each line with each mixture group, followed by your covariates, and in the last line your outcome)
vars<-c(
  # Air pollution 
  "h_no2_ratio_preg_Log",
  "h_pm10_ratio_preg_None",
  "h_pm25_ratio_preg_None",
  
  # Built environment 
  "h_builtdens300_preg_Sqrt",
  "h_fdensity300_preg_Log",
  "h_frichness300_preg_None",
  "h_landuseshan300_preg_None",
  "h_popdens_preg_Sqrt",
  "h_lden_preg_None",
  "h_ln_preg_None",
  
  # Traffic load
  "h_trafload_preg_pow1over3",
  "h_trafnear_preg_pow1over3",
  
  # Green spaces and walkability
  "h_ndvi100_preg_None",
  "h_walkability_mean_preg_None",
  
  # Metals
  "hs_cd_m_Log2",
  "hs_hg_m_Log2",
  "hs_pb_m_Log2",
  
  # OCs
  "hs_dde_madj_Log2",
  "hs_hcb_madj_Log2",
  "hs_pcb138_madj_Log2",
  "hs_pcb153_madj_Log2",
  "hs_pcb180_madj_Log2",
  
  # PFAS
  "hs_pfhxs_m_Log2",
  "hs_pfna_m_Log2",
  "hs_pfoa_m_Log2",
  "hs_pfos_m_Log2",
  
  # Lifestyle factors
  "e3_asmokyn_p_None",
  "e3_alcpreg_yn_None",
  "h_mbmi_None",
  "hs_wgtgain_None",
  "e3_gac_None",
  
  #Diet
  "h_fish_preg_Ter",  
  "h_fruit_preg_Ter", 
  "h_legume_preg_Ter",
  "h_veg_preg_Ter", 
  "h_dairy_preg_Ter", 
  "h_meat_preg_Ter",
  
  #Environment
  "h_humidity_preg_None",
  "h_pressure_preg_None",
  "h_temperature_preg_None",
  "hs_uvdvf_mt_hs_h_None",
  
  
  #COVARS (subcohorts)
  "h_parity_None",
  "h_native_None",
  "h_edumc_None",
  "h_age_None",
  "h_cohort_BIB", 
  "h_cohort_EDEN",
  "h_cohort_KANC", 
  "h_cohort_MOBA", 
  "h_cohort_INMA",
  
  
  #Outcome
  "hs2_zdia_bp.v3_2017_Time2",
  "hs2_zsys_bp.v3_2017_Time2"
)

data_all <- df_all[,vars]
data_n <- df_n[,vars]
data_s <- df_s[,vars]

data_n <- data_n[,!colnames(data_n) %in% c("h_cohort_INMA","h_cohort_MOBA")]
data_s <- data_s[,!colnames(data_s) %in% c("h_cohort_BIB", "h_cohort_EDEN", "h_cohort_KANC", "h_cohort_MOBA")]

data_list <- list("data_all"=data_all,
                  "data_n"=data_n,
                  "data_s"=data_s)

comps <- c("hs2_zdia_bp.v3_2017_Time2","hs2_zsys_bp.v3_2017_Time2" )

q = 4  # number of quantiles
G = 10 # Number of mixture groups

res_list <- c()

for (x in names(data_list)){
  data <- data_list[[x]]
  
  for (comp in comps){
    ##Add your chemicals belongign to each mixture and delete the parts you don't need (if you have less than 10 groups)
    if (x=="data_all"){
      formula = as.formula(paste0(comp, "~ h_cohort_BIB + h_cohort_EDEN + h_cohort_KANC + h_cohort_MOBA + h_cohort_INMA + h_parity_None + h_native_None + h_edumc_None + h_age_None")) # write only the name of the covariates
    }
    if (x=="data_n"){
      formula = as.formula(paste0(comp, "~ h_cohort_BIB + h_cohort_EDEN + h_cohort_KANC + h_parity_None + h_native_None + h_edumc_None + h_age_None")) # write only the name of the covariates
    }
    if (x=="data_s"){
      formula = as.formula(paste0(comp, "~ 1 + h_cohort_INMA + h_parity_None + h_native_None + h_edumc_None + h_age_None")) # write only the name of the covariates
    }
    y_name  <- all.vars(formula)[1]
    KV_name <- all.vars(formula)[-1]
    
    mix_name_1 <- c("h_no2_ratio_preg_Log", "h_pm10_ratio_preg_None", "h_pm25_ratio_preg_None")
    mix_name_2 <- c("h_builtdens300_preg_Sqrt", "h_fdensity300_preg_Log", "h_frichness300_preg_None", "h_landuseshan300_preg_None", "h_popdens_preg_Sqrt", "h_lden_preg_None", "h_ln_preg_None")
    mix_name_3 <- c("h_trafload_preg_pow1over3", "h_trafnear_preg_pow1over3")
    mix_name_4 <- c("h_ndvi100_preg_None", "h_walkability_mean_preg_None") 
    mix_name_5 <- c("hs_cd_m_Log2", "hs_hg_m_Log2", "hs_pb_m_Log2") 
    mix_name_6 <- c("hs_dde_madj_Log2", "hs_hcb_madj_Log2", "hs_pcb138_madj_Log2", "hs_pcb153_madj_Log2", "hs_pcb180_madj_Log2")
    mix_name_7 <- c("hs_pfhxs_m_Log2", "hs_pfna_m_Log2", "hs_pfoa_m_Log2", "hs_pfos_m_Log2") 
    mix_name_8 <- c("e3_asmokyn_p_None", "e3_alcpreg_yn_None", "h_mbmi_None", "hs_wgtgain_None", "e3_gac_None")
    mix_name_9 <- c("h_fish_preg_Ter", "h_fruit_preg_Ter", "h_legume_preg_Ter", "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")
    mix_name_10  <- c("h_humidity_preg_None", "h_pressure_preg_None", "h_temperature_preg_None", "hs_uvdvf_mt_hs_h_None")
    
    X1 = BWQS::quantile_split(data = data, mix_name = mix_name_1, q)[, mix_name_1]
    X2 = BWQS::quantile_split(data = data, mix_name = mix_name_2, q)[, mix_name_2]
    X3 = BWQS::quantile_split(data = data, mix_name = mix_name_3, q)[, mix_name_3]
    X4 = BWQS::quantile_split(data = data, mix_name = mix_name_4, q)[, mix_name_4]
    X5 = BWQS::quantile_split(data = data, mix_name = mix_name_5, q)[, mix_name_5]
    X6 = BWQS::quantile_split(data = data, mix_name = mix_name_6, q)[, mix_name_6]
    X7 = BWQS::quantile_split(data = data, mix_name = mix_name_7, q)[, mix_name_7]
    X8 = BWQS::quantile_split(data = data, mix_name = mix_name_8, q)[, mix_name_8]
    X9 = BWQS::quantile_split(data = data, mix_name = mix_name_9, q)[, mix_name_9]
    X10 = BWQS::quantile_split(data = data, mix_name = mix_name_10, q)[, mix_name_10]
    
    #Same, delete the parts you don't need if you have less than 10 mixture groups
    #Probably you will need to change G and set the number of groups you have.
    data_reg <- list(
      N   = nrow(data),
      
      C1  = length(mix_name_1),
      C2  = length(mix_name_2),
      C3  = length(mix_name_3),
      C4  = length(mix_name_4),
      C5  = length(mix_name_5),
      C6  = length(mix_name_6),
      C7  = length(mix_name_7),
      C8  = length(mix_name_8),
      C9  = length(mix_name_9),
      C10  = length(mix_name_10),
      
      XC1 = cbind(X1),
      XC2 = cbind(X2),
      XC3 = cbind(X3),
      XC4 = cbind(X4),
      XC5 = cbind(X5),
      XC6 = cbind(X6),
      XC7 = cbind(X7),
      XC8 = cbind(X8),
      XC9 = cbind(X9),
      XC10 = cbind(X10),
      
      DalpC1 = rep(1, length(mix_name_1)),
      DalpC2 = rep(1, length(mix_name_2)),
      DalpC3 = rep(1, length(mix_name_3)),
      DalpC4 = rep(1, length(mix_name_4)),
      DalpC5 = rep(1, length(mix_name_5)),
      DalpC6 = rep(1, length(mix_name_6)),
      DalpC7 = rep(1, length(mix_name_7)),
      DalpC8 = rep(1, length(mix_name_8)),
      DalpC9 = rep(1, length(mix_name_9)),
      DalpC10 = rep(1, length(mix_name_10)),
      
      KV = as.matrix(data[, KV_name]),
      K   = length(KV_name),
      G = 10,
      y = as.vector(data[, y_name])
    )
    
    #TRY TO RUN THIS AND CHANGE the number of iterations from 10000 to 100 just to see if the code is running. Then, you can repeat the model with 10000.
    fit_lasso_mixt0_zbmi8 <- rstan::sampling(
      m_lasso,
      data = data_reg,
      chains = 4,
      iter = 10000,
      thin = 1,
      refresh = 0,
      verbose = T,
      control = list(max_treedepth = 20, adapt_delta = 0.999999999999999))
    
    result_name <- paste0(x, "_", comp)
    saveRDS(fit_lasso_mixt0_zbmi8, file = paste0("./results/ExWAS/Mixture_analysis_BP_adol/fit_",result_name,".rds"))
    res_list[[result_name]] <- fit_lasso_mixt0_zbmi8
    
    #Results
    s_pchem_zbmi8 <- summary(fit_lasso_mixt0_zbmi8, probs = c(0.025, 0.975))
    results_chemp_zbmi8 <- as.data.frame(s_pchem_zbmi8$summary)
    
    #Select rows corresponding to contributions (THIS MAY CHANGE DEPENDING ON THE N, N covars, and THE MIXTURES GROUPS!!!!!!!!!!!!)
    if (x=="data_all"){
      results_chemp_zbmi8 <-results_chemp_zbmi8[c(rownames(results_chemp_zbmi8)[3:12], rownames(results_chemp_zbmi8)[33:73]),]###############3
    }
    if (x=="data_n"){
      results_chemp_zbmi8 <-results_chemp_zbmi8[c(rownames(results_chemp_zbmi8)[3:12], rownames(results_chemp_zbmi8)[31:71]),]###############3
    }
    if (x=="data_s"){
      results_chemp_zbmi8 <-results_chemp_zbmi8[c(rownames(results_chemp_zbmi8)[3:12], rownames(results_chemp_zbmi8)[29:69]),]###############3
    }
    
    rownames(results_chemp_zbmi8) = c("Air pollution", "Built environment", "Traffic", "Green spaces and walkability", 
                                      "Metals", "OCs", "PFAS", "Maternal personal factors", "Diet", "Meteorological factors",
                                      mix_name_1, mix_name_2, mix_name_3, mix_name_4, mix_name_5, 
                                      mix_name_6, mix_name_7, mix_name_8, mix_name_9, mix_name_10)
    
    data$air_pollution.mixture0.index_zbmi8 = results_chemp_zbmi8["Air pollution","mean"] * as.matrix(data_reg$XC1) %*% matrix(results_chemp_zbmi8[mix_name_1,"mean"], nrow = length(mix_name_1))
    data$built_env.mixture0.index_zbmi8 = results_chemp_zbmi8["Built environment","mean"] * as.matrix(data_reg$XC2) %*% matrix(results_chemp_zbmi8[mix_name_2,"mean"], nrow = length(mix_name_2))
    data$traffic.mixture0.index_zbmi8 = results_chemp_zbmi8["Traffic","mean"] * as.matrix(data_reg$XC3) %*% matrix(results_chemp_zbmi8[mix_name_3,"mean"], nrow = length(mix_name_3))
    data$naturalspaces.mixture0.index_zbmi8 = results_chemp_zbmi8["Green spaces and walkability","mean"] * as.matrix(data_reg$XC4) %*% matrix(results_chemp_zbmi8[mix_name_4,"mean"], nrow = length(mix_name_4))
    data$metals.mixture0.index_zbmi8 = results_chemp_zbmi8["Metals","mean"] * as.matrix(data_reg$XC5) %*% matrix(results_chemp_zbmi8[mix_name_5,"mean"], nrow = length(mix_name_5))
    data$OCs.mixture0.index_zbmi8 = results_chemp_zbmi8["OCs","mean"] * as.matrix(data_reg$XC6) %*% matrix(results_chemp_zbmi8[mix_name_6,"mean"], nrow = length(mix_name_6))
    data$PFAS.mixture0.index_zbmi8 = results_chemp_zbmi8["PFAS","mean"] * as.matrix(data_reg$XC7) %*% matrix(results_chemp_zbmi8[mix_name_7,"mean"], nrow = length(mix_name_7))
    data$lifestyle.mixture0.index_zbmi8 = results_chemp_zbmi8["Maternal personal factors","mean"] * as.matrix(data_reg$XC8) %*% matrix(results_chemp_zbmi8[mix_name_8,"mean"], nrow = length(mix_name_8))
    data$diet.mixture0.index_zbmi8 = results_chemp_zbmi8["Diet","mean"] * as.matrix(data_reg$XC9) %*% matrix(results_chemp_zbmi8[mix_name_9,"mean"], nrow = length(mix_name_9))
    data$meteo.mixture0.index_zbmi8 = results_chemp_zbmi8["Meteorological factors","mean"] * as.matrix(data_reg$XC10) %*% matrix(results_chemp_zbmi8[mix_name_10,"mean"], nrow = length(mix_name_10))
    
    
    data$air_pollution.mixture0.index_zbmi8 <- as.numeric(data$air_pollution.mixture0.index_zbmi8)
    data$built_env.mixture0.index_zbmi8 <- as.numeric(data$built_env.mixture0.index_zbmi8)
    data$traffic.mixture0.index_zbmi8 <- as.numeric(data$traffic.mixture0.index_zbmi8)
    data$naturalspaces.mixture0.index_zbmi8 <- as.numeric(data$naturalspaces.mixture0.index_zbmi8)
    data$metals.mixture0.index_zbmi8 <- as.numeric(data$metals.mixture0.index_zbmi8)
    data$OCs.mixture0.index_zbmi8 <- as.numeric(data$OCs.mixture0.index_zbmi8)
    data$PFAS.mixture0.index_zbmi8 <- as.numeric(data$PFAS.mixture0.index_zbmi8)
    data$lifestyle.mixture0.index_zbmi8 <- as.numeric(data$lifestyle.mixture0.index_zbmi8)
    data$diet.mixture0.index_zbmi8 <- as.numeric(data$diet.mixture0.index_zbmi8)
    data$meteo.mixture0.index_zbmi8 <- as.numeric(data$meteo.mixture0.index_zbmi8)
    
    data$overall_mixture <- data$air_pollution.mixture0.index_zbmi8 + data$built_env.mixture0.index_zbmi8 + data$traffic.mixture0.index_zbmi8 + data$naturalspaces.mixture0.index_zbmi8 +
      data$metals.mixture0.index_zbmi8 + data$OCs.mixture0.index_zbmi8 + data$PFAS.mixture0.index_zbmi8 + data$lifestyle.mixture0.index_zbmi8 +  data$diet.mixture0.index_zbmi8 + data$meteo.mixture0.index_zbmi8 
    
    data$R = data$overall_mixture - data[[comp]]
    
    write.csv2(results_chemp_zbmi8, paste0("./results/ExWAS/Mixture_analysis_BP_adol/res_",result_name,"1.csv"))
    write.csv2(data, paste0("./results/ExWAS/Mixture_analysis_BP_adol/data/",result_name,"1.csv"))
    write_rds(data_reg, paste0("./results/ExWAS/Mixture_analysis_BP_adol/data/data_reg_",result_name,"1.RDS"))
    
  }
}

#Forest plot--------------------------------------------------------------------

res_all_diat2 <- read.csv2("./results/ExWAS/Mixture_analysis_BP_adol/res_data_all_hs2_zdia_bp.v3_2017_Time21.csv")%>% mutate(Population = "Pooled")
res_all_syst2 <- read.csv2("./results/ExWAS/Mixture_analysis_BP_adol/res_data_all_hs2_zsys_bp.v3_2017_Time21.csv")%>% mutate(Population = "Pooled")

res_n_diat2 <- read.csv2("./results/ExWAS/Mixture_analysis_BP_adol/res_data_n_hs2_zdia_bp.v3_2017_Time21.csv")%>% mutate(Population = "North")
res_n_syst2 <- read.csv2("./results/ExWAS/Mixture_analysis_BP_adol/res_data_n_hs2_zsys_bp.v3_2017_Time21.csv")%>% mutate(Population = "North")

res_s_diat2 <- read.csv2("./results/ExWAS/Mixture_analysis_BP_adol/res_data_s_hs2_zdia_bp.v3_2017_Time21.csv")%>% mutate(Population = "South")
res_s_syst2 <- read.csv2("./results/ExWAS/Mixture_analysis_BP_adol/res_data_s_hs2_zsys_bp.v3_2017_Time21.csv")%>% mutate(Population = "South")

df_diat2 <- rbind(res_all_diat2, res_n_diat2, res_s_diat2)
df_syst2 <- rbind(res_all_syst2, res_n_syst2, res_s_syst2)

df_diat2$Population <- factor(df_diat2$Population, 
                              levels =  c("North", "South", "Pooled"),
                              labels = c("North", "South", "Pooled"))

df_syst2$Population <- factor(df_syst2$Population, 
                              levels =  c("North", "South", "Pooled"),
                              labels = c("North", "South", "Pooled"))

df_diat2$X[df_diat2$X == "Lifestyle factors"] <- "Maternal personal factors"
df_syst2$X[df_syst2$X == "Lifestyle factors"] <- "Maternal personal factors"

codebook <- read.csv2("./db/exposome/CODEBOOK_ANALYSIS_AUGUSTO.csv")
colnames(codebook)[9] <- c("Label")



df_diat2 <- df_diat2 %>%
  left_join(codebook %>% dplyr::select(Variable_name_TRANS, Label), by = c("X" = "Variable_name_TRANS"))


df_syst2 <- df_syst2 %>%
  left_join(codebook %>% dplyr::select(Variable_name_TRANS, Label), by = c("X" = "Variable_name_TRANS"))


df_diat2$Label[df_diat2$Label == "" |is.na(df_diat2$Label)] <- df_diat2$X[df_diat2$Label == "" |is.na(df_diat2$Label)]

df_syst2$Label[df_syst2$Label == "" |is.na(df_syst2$Label)] <- df_syst2$X[df_syst2$Label == "" |is.na(df_syst2$Label)]


write_xlsx(list("DBP t2"=df_diat2,
                "SBP t2"=df_syst2),
           "./results/ExWAS/Mixture_analysis_BP_adol/bwqs_results.xlsx")

#SBP t1
colors <- c("North"="#D88E4B", "South"="#F8D4B0", "Pooled"="#F0B077")

sbpt1<- ggplot(df_syst1[df_syst1$X %in% c("Air pollution", "Built environment", "Traffic", "Green spaces and walkability", "Metals", "OCs", "PFAS", "Maternal personal factors", "Diet", "Meteorological factors"),],
               aes(y = X, x = mean, xmin = `X2.5.`, xmax = `X97.5.`, shape = Population, colour = Population)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Mixture analysis Beta exposure-SBP in childhood") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))

#SBP t2
colors <- c("North"="#BF620D", "South"="#F4AE6D", "Pooled"="#EA801C")

sbpt2 <- ggplot(df_syst2[df_syst2$X %in% c("Air pollution", "Built environment", "Traffic", "Green spaces and walkability", "Metals", "OCs", "PFAS", "Maternal personal factors", "Diet", "Meteorological factors"),],
                aes(y = X, x = mean, xmin = `X2.5.`, xmax = `X97.5.`, shape = Population, colour = Population)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Mixture analysis Beta exposure-SBP in adolescence") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))

#DBP t1
colors <- c("North"="#62A2CC", "South"="#B4DAF0", "Pooled"="#8CC5E3")

dbpt1 <- ggplot(df_diat1[df_diat1$X %in% c("Air pollution", "Built environment", "Traffic", "Green spaces and walkability", "Metals", "OCs", "PFAS", "Lifestyle factors", "Diet", "Meteorological factors"),],
                aes(y = X, x = mean, xmin = `X2.5.`, xmax = `X97.5.`, shape = Population, colour = Population)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Mixture analysis Beta exposure-DBP in childhood") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))

#DBP t2
colors <- c("North"="#1A6C9B", "South"="#72B9E5", "Pooled"="#3594CC")

dbpt2 <- ggplot(df_diat2[df_diat2$X %in% c("Air pollution", "Built environment", "Traffic", "Green spaces and walkability", "Metals", "OCs", "PFAS", "Maternal personal factors", "Diet", "Meteorological factors"),],
                aes(y = X, x = mean, xmin = `X2.5.`, xmax = `X97.5.`, shape = Population, colour = Population)) +
  geom_pointrange(size = 0.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors) + 
  labs(y = "Exposure", x = "Beta (CI 95%)", title = "Mixture analysis Beta exposure-DBP in adolescence") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.1) +  
  theme_minimal() +
  theme(text = element_text(size = 14))



library(patchwork)
(sbpt1 | sbpt2) / (dbpt1 | dbpt2)

mix_name_1 <- c("h_no2_ratio_preg_Log", "h_pm10_ratio_preg_None", "h_pm25_ratio_preg_None")
mix_name_2 <- c("h_builtdens300_preg_Sqrt", "h_fdensity300_preg_Log", "h_frichness300_preg_None", "h_landuseshan300_preg_None", "h_popdens_preg_Sqrt", "h_lden_preg_None", "h_ln_preg_None")
mix_name_3 <- c("h_trafload_preg_pow1over3", "h_trafnear_preg_pow1over3")
mix_name_4 <- c("h_ndvi100_preg_None", "h_walkability_mean_preg_None") 
mix_name_5 <- c("hs_cd_m_Log2", "hs_hg_m_Log2", "hs_pb_m_Log2") 
mix_name_6 <- c("hs_dde_madj_Log2", "hs_hcb_madj_Log2", "hs_pcb138_madj_Log2", "hs_pcb153_madj_Log2", "hs_pcb180_madj_Log2")
mix_name_7 <- c("hs_pfhxs_m_Log2", "hs_pfna_m_Log2", "hs_pfoa_m_Log2", "hs_pfos_m_Log2") 
mix_name_8 <- c("e3_asmokyn_p_None", "e3_alcpreg_yn_None", "h_mbmi_None", "hs_wgtgain_None", "e3_gac_None")
mix_name_9 <- c("h_fish_preg_Ter", "h_fruit_preg_Ter", "h_legume_preg_Ter", "h_veg_preg_Ter", "h_dairy_preg_Ter", "h_meat_preg_Ter")
mix_name_10  <- c("h_humidity_preg_None", "h_pressure_preg_None", "h_temperature_preg_None", "hs_uvdvf_mt_hs_h_None")

#Metals DBP t1
colors <- c("North"="#62A2CC", "South"="#B4DAF0", "Pooled"="#8CC5E3")
p1 <- ggplot(df_diat1[df_diat1$X %in% mix_name_5, ], aes(
  x = Label,
  y = mean,
  fill = Population
)) +
  geom_bar(stat = "identity", position = "dodge") +   geom_errorbar(
    aes(ymin = `X2.5.`, ymax = `X97.5.`),
    position = position_dodge(width = 0.9),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Weights relative to mixture (95% CI)", title = "Metals (DBP in childhood)") +
  geom_hline(
    yintercept = 1 / length(mix_name_5),
    linetype = "dashed",
    color = "black",
    size = 0.1
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14))

#MPF SBP t2
colors <- c("North"="#BF620D", "South"="#F4AE6D", "Pooled"="#EA801C")

p2 <- ggplot(df_syst2[df_syst2$X %in% mix_name_8, ], aes(
  x = Label,
  y = mean,
  fill = Population
)) +
  geom_bar(stat = "identity", position = "dodge") +   geom_errorbar(
    aes(ymin = `X2.5.`, ymax = `X97.5.`),
    position = position_dodge(width = 0.9),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Weights relative to mixture (95% CI)", title = "Maternal personal factors (SBP in adolescence)") +
  geom_hline(
    yintercept = 1 / length(mix_name_8),
    linetype = "dashed",
    color = "black",
    size = 0.1
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14))

#Diet SBP t2
colors <- c("North"="#BF620D", "South"="#F4AE6D", "Pooled"="#EA801C")

p3 <- ggplot(df_syst2[df_syst2$X %in% mix_name_9, ], aes(
  x = Label,
  y = mean,
  fill = Population
)) +
  geom_bar(stat = "identity", position = "dodge") +   geom_errorbar(
    aes(ymin = `X2.5.`, ymax = `X97.5.`),
    position = position_dodge(width = 0.9),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Weights relative to mixture (95% CI)", title = "Diet (SBP in adolescence)") +
  geom_hline(
    yintercept = 1 / length(mix_name_9),
    linetype = "dashed",
    color = "black",
    size = 0.1
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14))
(p1  | p3) / p2
