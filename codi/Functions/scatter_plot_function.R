library(ggplot2)

unite <- function(data_, x_, y_){
  outcome <- as.data.frame(data_)[[y_]]
  var <- as.data.frame(data_)[[x_]]
  return(data.frame("x"=var,"y"=outcome))
}

scatter_plot <- function(data, x, y, omicsLayer) {
  #Save names
  x <- as.character(x)
  y <- as.character(y)

  #Create df with data
  data<-unite(data, x, y)
  colnames(data)<-c(x,y)
  
  #Scatter plot
  plot <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point() +
    labs(title = paste0("Scatter plot ", name), x = paste0(x), y = paste0(y))
  
  #Save the plot
  path <- paste0(getwd(), "/results/scatter_plot/", omicsLayer, "/", name, "_", x, "_", y, ".png")
  ggsave(path, plot)
}
