library(MASS)
library(Matrix)
library(dlnm)
library(splines)
library(mgcv)
library(stringr)
library(sNPLS)
library(plyr)

################################# EvaluateModel.cv.spls() 

EvaluateModel.cv.spls <- function(x,y,max.K=10,nfolds=5,plot_RMSE_eta=T) {
  set.seed(1899)
  m <- spls::cv.spls(x=x,y=y,K=1:max.K,eta=c(seq(1e-6,1-1e-6,(1-1e-6-1e-6)/1000)),fold=nfolds,plot.it=F,scale.x=F,scale.y=F)   
  RMSE_H0 <- sqrt(mean((y-mean(y))^2))
  K.min <- m$K.opt   
  eta.min <- m$eta.opt   
  if (RMSE_H0<=min(m$mspemat)) {     
    K.min <- 0   
  }   
  if (K.min>0) {     
    model <- spls::spls(x=x,y=y,K=K.min,eta=eta.min)     
    coef <- as.numeric(coef(model))   
  }else {     
    coef <- rep(0,ncol(x))   
  }   
  names(coef) <- paste0("C",1:length(coef))   
  results <- cbind(data.frame(model="SPLS",type="MIN"),t(coef))
  results_with_parameters<-list("results"=results, "K"=K.min, "eta"=eta.min)

  if (plot_RMSE_eta){
    #Graph RMSE over Lambda(eta)
    df <- reshape2::melt(m$mspemat)
    df$Var1_num <- as.numeric(gsub("eta= ", "", df$Var1))
    df$Var2_num <- as.numeric(gsub("K = ", "", df$Var2))
    LambdaVal <- unique(round(unique(df$Var1_num),1))
    
    fun_color_range <- colorRampPalette(c("#08737f", "#f95559"))
    my_colors <- fun_color_range(length(unique(df$Var1_num)))
    graph <- ggplot(df, aes(x = Var2_num, y = value, fill = as.factor(Var1_num))) +
      geom_boxplot() +
      labs(x = "Averaged lambda over folds", y = "Mean squared error", fill = 'Lambda value') +
      scale_x_continuous(breaks = seq(0,max.K,1),labels = LambdaVal) +
      scale_fill_manual(values = my_colors) +
      theme(legend.position = "none")
    
    print(graph)
  }
  return(results_with_parameters) 
}

##### Apply sPLS

applySPLS<-function(data.Y,data.X) {
  set.seed(1899)
  beta <- NA
  # beta <- data.Y$beta
  # names(beta)<-data.Y$true.pred
  init <- colnames(data.X)
  RES<-as.data.frame(init)
  names(RES) <- "var"
  # sPLS
  SPLS <- EvaluateModel.cv.spls(x=data.X,y=data.Y)
  resSPLS<- SPLS$results
  colnames(resSPLS)[-(1:2)] <- colnames(data.X) 
  RES <- merge(RES,cbind(var=colnames(resSPLS),t(resSPLS)),by='var',all.x=TRUE) 
  colnames(RES)[ncol(RES)+1-nrow(resSPLS):1] <- paste(resSPLS[,1],"_", resSPLS[,2],sep="") 
  
  #summary
  RES2<-RES
  RES2$SPLS_MIN.TP<-as.factor(ifelse(is.na(RES2$SPLS_MIN),"0",ifelse(RES2$SPLS_MIN%in%"0","0","1")))
  RES3<-list("results"=RES2,"eta"=data.frame(SPLS$eta), "K"=data.frame(SPLS$K))

  return(RES3)
}
