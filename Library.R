#################################################################
# Code for article: Regular pattern formation regulates population dynamics: logistic growth in cellular automata
# Author : Jingyao Sun
# Date: 2019-10-22
#################################################################
# Library
#################################################################
# To get extended matrix to simulate periodic boundary condition
matrix.ext <- function(matrix,x_ext=2,y_ext=2){
  d_x <- dim(matrix)[1]
  d_y <- dim(matrix)[1]
  left  <- rbind(diag(1,d_x)[(d_x-x_ext+1):d_x,],diag(1,d_x),diag(1,d_x)[1:(1+x_ext-1),])
  right <- cbind(diag(1,d_y)[,(d_y-y_ext+1):d_y],diag(1,d_y),diag(1,d_y)[,1:(1+y_ext-1)])
  matrix_ext <- left %*% matrix %*% right
  return(matrix_ext)
}
# Get total cell number for given neighborhood
get.total.num <- function(a = 5,b = 3){
  # a is axis in x direction; b is axis in y direction
  test <- matrix(1,2*a+1,2*b+1)
  for(i in 1:(2*a+1)){
    for(j in 1:(2*b+1)){
      if(  ((i-a-1)^2)/(a^2) + ((j-b-1)^2)/(b^2) > 1 
      ){
        test[i,j] <- 0
      }
    }
  }
  return(sum(test))
}
# Get number of occupied cells in a neiborhood
get.oval.r <- function(matrix_ext,x_ext,y_ext,
                       x,y,a=5,b=3,
                       skewx = 0,
                       skewy = 0,which=1){
  x1 <- x + x_ext + skewx
  y1 <- y + y_ext + skewy
  #c  <- sqrt(a^2-b^2)
  matrix_cut <- matrix_ext[(x1-a):(x1+a),(y1-b):(y1+b)]
  for(i in 1:(2*a+1)){
    for(j in 1:(2*b+1)){
      if(((i-a-1)^2)/(a^2) + ((j-b-1)^2)/(b^2) > 1 ){ 
        #sqrt((i-(a+1+c))^2+(j-b-1)^2) + 
        # sqrt((i-(a+1-c))^2+(j-b-1)^2) 
        #  > 2*a ){
        matrix_cut[i,j] <- NA
      }
    }
  }
  result <- sum(matrix_cut==which,na.rm=T)
  return(result)
}
# Main function for states updating
next_step <- function(states,
                      alpha = 70,
                      beta = 42,
                      deltat = 0.5,
                      growth = 2,
                      num_all = 113,
                      a = 6,
                      b = 6,
                      skewx=0,
                      skewy=0){
  d_x <- dim(states)[1]
  d_y <- dim(states)[2]
  skewm <- max(skewx,skewy)
  # The extended matrix is adopted to simulate periodic boundary condition
  states_ext  <- matrix.ext(states,a+b+skewm,a+b+skewm)
  states_next <- matrix(0,d_x,d_y)
  speed <- growth * deltat
  for(i in 1:d_x){
    for(j in 1:d_y){
      skewxs= sample(size=1,x=c(-skewx,skewx),prob=c(0.5,0.5)) # Runoff blockage 
      num0 <- get.oval.r(matrix_ext = states_ext,x_ext = a+b+skewm,y_ext = a+b+skewm,
                         x = i,y = j,a=a,b=b,skewx=skewxs,skewy = skewy,which = 0)
      num1 <- num_all - num0
      if((states[i,j] == 0) & ((1- (alpha*num1)/(beta*num_all))>0)){
        rate <- speed*num1*(1- (alpha*num1)/(beta*num_all))/num0
        if(rate<0){rate <- 0}
        if(rate>1){rate <- 1}
        states_next[i,j] <- sample(size=1,x=c(0,1), prob=c(1-rate,rate),replace=T)
      }else if((states[i,j] ==1) & ((1- (alpha*num1)/(beta*num_all))<0)){
        rate <- (-1)*speed*num1*(1- (alpha*num1)/(beta*num_all))/num1
        if(rate<0){rate <- 0}
        if(rate>1){rate <- 1}
        states_next[i,j] <- sample(size=1,x=c(0,1), prob=c(rate,1-rate),replace=T)
      }else{
        states_next[i,j] <- states[i,j]
      }
    }
  }
  return(states_next)
}
# Plotting
library("ggplot2")
automata.plot <- function(a,main="automata"){
  library(ggplot2)
  data.trans <- function(data){
    l <- dim(data)[1]
    y <- rep(1:l,each=l)
    x <- rep(1:l,l)
    state <- c(data)
    return(data.frame(x,y,state))
  }
  test  <- data.trans(a)
  test$state <- as.factor(test$state)
  ggplot(data = test,aes(x=x,y=y,fill=state))+
    geom_tile(colour='grey')+
    labs(title=element_blank())+
    theme_bw() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid =element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")+
    scale_fill_manual(breaks = c("0", "1"),
                      values=c("#FFFFFF", "#666666"))+
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))
}
