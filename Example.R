#################################################################
# Code for article: Regular pattern formation regulates population dynamics: logistic growth in cellular automata
# Author : Jingyao Sun
# Date: 2019-10-22
#################################################################
# Example 1
#################################################################
source("Library.R")
matrix_ini <- matrix(rbinom(2500, size=1, prob=0.01),50,50) # Creat initial matrix
matrix_temp <- matrix_ini 
num_all <- get.total.num(6,6)  # Get total cell number of neiborhood
for(i in 1:200){matrix_temp <- next_step(matrix_temp,
                                         alpha = 70,
                                         beta = 15,
                                         deltat = 0.5,
                                         growth = 2,
                                         num_all = num_all,
                                         a = 6,
                                         b = 6,
                                         skewx=0,
                                         skewy=0)}
automata.plot(matrix_temp) # Spots pattern 

matrix_ini <- matrix(rbinom(2500, size=1, prob=0.01),50,50) # Creat initial matrix
matrix_temp <- matrix_ini 
num_all <- get.total.num(6,6)  # Get total cell number of neiborhood
for(i in 1:200){matrix_temp <- next_step(matrix_temp,
                                         alpha = 70,
                                         beta = 35,
                                         deltat = 0.5,
                                         growth = 2,
                                         num_all = num_all,
                                         a = 6,
                                         b = 6,
                                         skewx=0,
                                         skewy=0)}
automata.plot(matrix_temp) # labyrinths pattern 

matrix_ini <- matrix(rbinom(2500, size=1, prob=0.01),50,50) # Creat initial matrix
matrix_temp <- matrix_ini 
num_all <- get.total.num(6,6)  # Get total cell number of neiborhood
for(i in 1:200){matrix_temp <- next_step(matrix_temp,
                                         alpha = 70,
                                         beta = 55,
                                         deltat = 0.5,
                                         growth = 2,
                                         num_all = num_all,
                                         a = 6,
                                         b = 6,
                                         skewx=0,
                                         skewy=0)}
automata.plot(matrix_temp) # Gaps pattern 

matrix_ini <- matrix(rbinom(2500, size=1, prob=0.01),50,50) # Creat initial matrix
matrix_temp <- matrix_ini 
num_all <- get.total.num(6,6)  # Get total cell number of neiborhood
for(i in 1:200){matrix_temp <- next_step(matrix_temp,
                                         alpha = 70,
                                         beta = 35,
                                         deltat = 0.5,
                                         growth = 2,
                                         num_all = num_all,
                                         a = 6,
                                         b = 6,
                                         skewx=2,
                                         skewy=1)}
automata.plot(matrix_temp) # Bands pattern 
