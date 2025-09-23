require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
library(dplyr)
library(tidyverse)
library(fastDummies)
library(ggeffects)

setwd("./TFM")

#Load data
pheno <- readRDS("./db/pheno/final/bp_wide_validN5332023-10-16.rds")
rownames(pheno)<-pheno$HelixID
pheno<-pheno[order(rownames(pheno)),]
pheno <- pheno[c("hs2_BPcat_v3_2017_bin_Time1", "hs2_BPcat_v3_2017_bin_Time2", "HelixID")]
pheno$hs2_TRAJBPcat_v3_2017_Time <- ifelse(pheno$hs2_BPcat_v3_2017_bin_Time1=="Normal" & pheno$hs2_BPcat_v3_2017_bin_Time2=="Normal", 1,
                                           ifelse(pheno$hs2_BPcat_v3_2017_bin_Time1=="Altered"  & pheno$hs2_BPcat_v3_2017_bin_Time2=="Normal", 2,
                                                  ifelse(pheno$hs2_BPcat_v3_2017_bin_Time1=="Normal" & pheno$hs2_BPcat_v3_2017_bin_Time2=="Altered", 3, 4)))

pheno$hs2_TRAJBPcat_v3_2017_Time <- as.factor(pheno$hs2_TRAJBPcat_v3_2017_Time)

#Load projections 
projections <- readRDS("./results/RGCCA/projections.rds")
proj_north <- projections$`Projections N` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_north) <- names(projections$`Projections N`)
proj_south <- projections$`Projections S` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_south) <- names(projections$`Projections S`)
proj_all <-  projections$`Projections all` %>% purrr::reduce(cbind) %>% data.frame()
colnames(proj_all) <- names(projections$`Projections all`)
 
#dataset containing proj and BP trajectories 
pheno_all <- pheno[rownames(pheno) %in% rownames(proj_all),]
pheno_n <- pheno[rownames(pheno) %in% rownames(proj_north),]
pheno_s <- pheno[rownames(pheno) %in% rownames(proj_south),]

all(rownames(proj_all)==rownames(pheno_all))
all(rownames(proj_north)==rownames(pheno_n))
all(rownames(proj_south)==rownames(pheno_s))

proj_all$HelixID <- rownames(proj_all)
proj_north$HelixID <- rownames(proj_north)
proj_south$HelixID <- rownames(proj_south)

df_all <- merge(proj_all, pheno_all, by = "HelixID") %>% mutate(Population = "All")
df_n <- merge(proj_north, pheno_n, by = "HelixID") %>% mutate(Population = "North")
df_s <- merge(proj_south, pheno_s, by = "HelixID") %>% mutate(Population = "South")

df_combined <- rbind(df_all, df_n, df_s) %>% select(-c("HelixID"))

#Multinomial regression_________________________________________________________
multinom_all <- multinom(hs2_TRAJBPcat_v3_2017_Time ~ prot + serum + urine, data = df_all, model = T, Hess=T) 
predictions.multinom_all_prot <- ggeffects::ggemmeans(multinom_all, terms = "prot [all]") %>% rename("BP trajectory" = response.level) %>% data.frame()
predictions.multinom_all_serum <- ggeffects::ggemmeans(multinom_all, terms = "serum [all]") %>% rename("BP trajectory" = response.level) %>% data.frame()
predictions.multinom_all_urine <- ggeffects::ggemmeans(multinom_all, terms = "urine [all]") %>% rename("BP trajectory" = response.level) %>% data.frame()

multinom_n <- multinom(hs2_TRAJBPcat_v3_2017_Time ~ prot + serum + urine, data = df_n, model = T, Hess=T)
predictions.multinom_n_prot <- ggeffects::ggemmeans(multinom_n, terms = "prot [all]") %>% rename("BP trajectory" = response.level)
predictions.multinom_n_serum <- ggeffects::ggemmeans(multinom_n, terms = "serum [all]") %>% rename("BP trajectory" = response.level)
predictions.multinom_n_urine <- ggeffects::ggemmeans(multinom_n, terms = "urine [all]") %>% rename("BP trajectory" = response.level)

multinom_s <- multinom(hs2_TRAJBPcat_v3_2017_Time ~ prot + serum + urine, data = df_s, model = T, Hess=T)
predictions.multinom_s_prot <- ggeffects::ggemmeans(multinom_s, terms = "prot [all]") %>% rename("BP trajectory" = response.level)
predictions.multinom_s_serum <- ggeffects::ggemmeans(multinom_s, terms = "serum [all]") %>% rename("BP trajectory" = response.level)
predictions.multinom_s_urine <- ggeffects::ggemmeans(multinom_s, terms = "urine [all]") %>% rename("BP trajectory" = response.level)

#p_vals
z <- summary(multinom_all)$coefficients/summary(multinom_all)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

#Relative odd ratio
predictions.multinom_all_prot$BP.trajectory <- as.numeric(as.character(predictions.multinom_all_prot$BP.trajectory))
predictions.multinom_all_serum$BP.trajectory <- as.numeric(as.character(predictions.multinom_all_serum$BP.trajectory))
predictions.multinom_all_urine$BP.trajectory <- as.numeric(as.character(predictions.multinom_all_urine$BP.trajectory))
predictions.multinom_all_prot$odd_ratio<- NA
predictions.multinom_all_serum$odd_ratio<- NA
predictions.multinom_all_urine$odd_ratio <- NA

predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==1,]$odd_ratio <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==1,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot, mean)[1,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==2,]$odd_ratio <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==2,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot, mean)[2,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==3,]$odd_ratio <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==3,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot, mean)[3,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==4,]$odd_ratio <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==4,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot, mean)[4,2]

predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==1,]$conf.high <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==1,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot, mean)[1,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==2,]$conf.high <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==2,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot, mean)[2,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==3,]$conf.high <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==3,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot, mean)[3,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==4,]$conf.high <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==4,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot, mean)[4,2]

predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==1,]$conf.low <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==1,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot, mean)[1,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==2,]$conf.low <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==2,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot, mean)[2,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==3,]$conf.low <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==3,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot, mean)[3,2]
predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==4,]$conf.low <- predictions.multinom_all_prot[predictions.multinom_all_prot$BP.trajectory==4,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot, mean)[4,2]

predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==1,]$odd_ratio <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==1,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum, mean)[1,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==2,]$odd_ratio <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==2,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum, mean)[2,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==3,]$odd_ratio <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==3,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum, mean)[3,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==4,]$odd_ratio <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==4,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum, mean)[4,2]

predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==1,]$conf.high <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==1,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum, mean)[1,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==2,]$conf.high <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==2,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum, mean)[2,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==3,]$conf.high <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==3,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum, mean)[3,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==4,]$conf.high <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==4,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum, mean)[4,2]

predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==1,]$conf.low <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==1,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum, mean)[1,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==2,]$conf.low <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==2,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum, mean)[2,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==3,]$conf.low <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==3,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum, mean)[3,2]
predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==4,]$conf.low <- predictions.multinom_all_serum[predictions.multinom_all_serum$BP.trajectory==4,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum, mean)[4,2]

predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==1,]$odd_ratio <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==1,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine, mean)[1,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==2,]$odd_ratio <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==2,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine, mean)[2,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==3,]$odd_ratio <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==3,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine, mean)[3,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==4,]$odd_ratio <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==4,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine, mean)[4,2]

predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==1,]$conf.high <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==1,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine, mean)[1,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==2,]$conf.high <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==2,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine, mean)[2,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==3,]$conf.high <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==3,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine, mean)[3,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==4,]$conf.high <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==4,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine, mean)[4,2]

predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==1,]$conf.low <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==1,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine, mean)[1,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==2,]$conf.low <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==2,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine, mean)[2,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==3,]$conf.low <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==3,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine, mean)[3,2]
predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==4,]$conf.low <- predictions.multinom_all_urine[predictions.multinom_all_urine$BP.trajectory==4,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine, mean)[4,2]

#Plots 
predictions.multinom_all_prot$BP.trajectory <- as.factor(predictions.multinom_all_prot$BP.trajectory)
predictions.multinom_all_serum$BP.trajectory <- as.factor(predictions.multinom_all_serum$BP.trajectory)
predictions.multinom_all_urine$BP.trajectory <- as.factor(predictions.multinom_all_urine$BP.trajectory)

labs <- c("1" = "Normal to normal",
          "2" = "Altered to normal",
          "3" = "Normal to altered",
          "4" = "Altered to altered")

colors <- c("1"="#089912",
            "2"="#C1C1C1",
            "3"="#FAB7A7",
            "4"="#950C00")

ggplot(predictions.multinom_all_prot, aes(x = x, y = odd_ratio, colour =BP.trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill =BP.trajectory), alpha = 0.2) + 
  labs(x = "Prot LC", y = "Odd ratio",
       title = "Predicted Probabilities of BP trajectories") +
  scale_color_manual(values = colors, labels = labs) + 
  scale_fill_manual(values = colors, labels = labs) +
  theme_minimal()

ggplot(predictions.multinom_all_serum, aes(x = x, y = odd_ratio, colour =BP.trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill =BP.trajectory), alpha = 0.2) + 
  labs(x = "Serum LC", y = "Odd ratio",
       title = "Predicted Probabilities of BP trajectories") +
  scale_color_manual(values = colors, labels = labs) + 
  scale_fill_manual(values = colors, labels = labs) +
  theme_minimal()

ggplot(predictions.multinom_all_urine, aes(x = x, y = odd_ratio, colour =BP.trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill =BP.trajectory), alpha = 0.2) + 
  labs(x = "Urine LC", y = "Odd ratio",
       title = "Predicted Probabilities of BP trajectories") +
  scale_color_manual(values = colors, labels = labs) + 
  scale_fill_manual(values = colors, labels = labs) +
  theme_minimal()

#Cases for each population
table(df_all$hs2_TRAJBPcat_v3_2017_Time)
table(df_n$hs2_TRAJBPcat_v3_2017_Time)
table(df_s$hs2_TRAJBPcat_v3_2017_Time)

#Univar multinomial reg_________________________________________________________
multinom_all_prot <- multinom(hs2_TRAJBPcat_v3_2017_Time ~ prot, data = df_all, model = T, Hess=T) 
multinom_all_serum <- multinom(hs2_TRAJBPcat_v3_2017_Time ~ serum, data = df_all, model = T, Hess=T) 
multinom_all_urine <- multinom(hs2_TRAJBPcat_v3_2017_Time ~ urine, data = df_all, model = T, Hess=T) 

predictions.multinom_all_prot1 <- ggeffects::ggemmeans(multinom_all_prot, terms = "prot [all]") %>% rename("BP trajectory" = response.level) %>% data.frame()
predictions.multinom_all_serum1 <- ggeffects::ggemmeans(multinom_all_serum, terms = "serum [all]") %>% rename("BP trajectory" = response.level) %>% data.frame()
predictions.multinom_all_urine1 <- ggeffects::ggemmeans(multinom_all_urine, terms = "urine [all]") %>% rename("BP trajectory" = response.level) %>% data.frame()

z <- summary(multinom_all_prot)$coefficients/summary(multinom_all_prot)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

z <- summary(multinom_all_serum)$coefficients/summary(multinom_all_serum)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

z <- summary(multinom_all_urine)$coefficients/summary(multinom_all_urine)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

#Relative odd ratio
predictions.multinom_all_prot1$BP.trajectory <- as.numeric(as.character(predictions.multinom_all_prot1$BP.trajectory))
predictions.multinom_all_serum1$BP.trajectory <- as.numeric(as.character(predictions.multinom_all_serum1$BP.trajectory))
predictions.multinom_all_urine1$BP.trajectory <- as.numeric(as.character(predictions.multinom_all_urine1$BP.trajectory))
predictions.multinom_all_prot1$odd_ratio<- NA
predictions.multinom_all_serum1$odd_ratio<- NA
predictions.multinom_all_urine1$odd_ratio <- NA

predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==1,]$odd_ratio <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==1,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[1,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==2,]$odd_ratio <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==2,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[2,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==3,]$odd_ratio <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==3,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[3,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==4,]$odd_ratio <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==4,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[4,2]

predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==1,]$conf.high <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==1,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[1,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==2,]$conf.high <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==2,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[2,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==3,]$conf.high <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==3,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[3,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==4,]$conf.high <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==4,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[4,2]

predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==1,]$conf.low <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==1,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[1,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==2,]$conf.low <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==2,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[2,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==3,]$conf.low <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==3,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[3,2]
predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==4,]$conf.low <- predictions.multinom_all_prot1[predictions.multinom_all_prot1$BP.trajectory==4,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_prot1, mean)[4,2]

predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==1,]$odd_ratio <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==1,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[1,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==2,]$odd_ratio <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==2,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[2,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==3,]$odd_ratio <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==3,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[3,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==4,]$odd_ratio <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==4,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[4,2]

predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==1,]$conf.high <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==1,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[1,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==2,]$conf.high <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==2,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[2,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==3,]$conf.high <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==3,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[3,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==4,]$conf.high <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==4,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[4,2]

predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==1,]$conf.low <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==1,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[1,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==2,]$conf.low <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==2,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[2,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==3,]$conf.low <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==3,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[3,2]
predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==4,]$conf.low <- predictions.multinom_all_serum1[predictions.multinom_all_serum1$BP.trajectory==4,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_serum1, mean)[4,2]

predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==1,]$odd_ratio <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==1,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[1,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==2,]$odd_ratio <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==2,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[2,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==3,]$odd_ratio <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==3,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[3,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==4,]$odd_ratio <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==4,]$predicted/aggregate(predicted ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[4,2]

predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==1,]$conf.high <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==1,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[1,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==2,]$conf.high <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==2,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[2,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==3,]$conf.high <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==3,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[3,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==4,]$conf.high <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==4,]$conf.high/aggregate(conf.high ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[4,2]

predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==1,]$conf.low <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==1,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[1,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==2,]$conf.low <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==2,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[2,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==3,]$conf.low <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==3,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[3,2]
predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==4,]$conf.low <- predictions.multinom_all_urine1[predictions.multinom_all_urine1$BP.trajectory==4,]$conf.low/aggregate(conf.low ~BP.trajectory, data = predictions.multinom_all_urine1, mean)[4,2]

#Plots 
predictions.multinom_all_prot1$BP.trajectory <- as.factor(predictions.multinom_all_prot1$BP.trajectory)
predictions.multinom_all_serum1$BP.trajectory <- as.factor(predictions.multinom_all_serum1$BP.trajectory)
predictions.multinom_all_urine1$BP.trajectory <- as.factor(predictions.multinom_all_urine1$BP.trajectory)

ggplot(predictions.multinom_all_prot1, aes(x = x, y = odd_ratio, colour =BP.trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill =BP.trajectory), alpha = 0.2) + 
  labs(x = "Prot LC", y = "Odd ratio",
       title = "Predicted Probabilities of BP trajectories") +
  scale_color_manual(values = colors, labels = labs) + 
  scale_fill_manual(values = colors, labels = labs) +
  theme_minimal()

ggplot(predictions.multinom_all_serum1, aes(x = x, y = odd_ratio, colour =BP.trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill =BP.trajectory), alpha = 0.2) + 
  labs(x = "Serum LC", y = "Odd ratio",
       title = "Predicted Probabilities of BP trajectories") +
  scale_color_manual(values = colors, labels = labs) + 
  scale_fill_manual(values = colors, labels = labs) +
  theme_minimal()

ggplot(predictions.multinom_all_urine1, aes(x = x, y = odd_ratio, colour =BP.trajectory)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill =BP.trajectory), alpha = 0.2) + 
  labs(x = "Urine LC", y = "Odd ratio",
       title = "Predicted Probabilities of BP trajectories") +
  scale_color_manual(values = colors, labels = labs) + 
  scale_fill_manual(values = colors, labels = labs) +
  theme_minimal()
