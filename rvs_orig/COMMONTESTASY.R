
setwd("C:/Users/Scott/Desktop/RVS-master/")
source("TESTHELPER.R")

X <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/X.txt", sep="\t", header = F)
Y <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Y.txt", sep="\t", header = F)
P <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/P.txt", sep="\t", header = F)
M <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/M.txt", sep="\t", header = T)
RVS=TRUE

P <- P[1,]
X <- X[1,]



n = length(X[,1])

for (i in 1:100){

  print(RVS_asy(Y,X[i,],P[i,]))
  
}