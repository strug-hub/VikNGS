library("numDeriv")

#'  Function to calcualte the robust variance of \eqn{E(Gij|Dij)} for case
#'
#'   Use formule \eqn{Var(X)=E(X^2)-E(X)^2} to calculate the variance of genotypes. It is called from \code{\link{RVS_asy}},\code{\link{RVS_btrap}}, and \code{\link{calc_ScoreV_RVS}}.
#' @param P  the genotype frequency, calcualted from EM algorithm
#' @return the Robust variance of E(Gij|Dij)
calc_robust_Var = function(P){
  Sq = 4*P[3] + P[2]
  Sm = 2*P[3] + P[2]
  S = Sq - Sm^2
  return (S)
}

#For Binomial

binom.b=function(eta) {
  return(log(1+exp(eta)) )
}

calc_mean=function(Y,X,Z,rd.group) {
  
  #First, remove NA
  data1=na.omit(data.frame(cbind(Y,X,Z,rd.group)))  #X is the SNP data
  rd.group = data1$rd.group
  Y = data1$Y
  
  Z = Z[!is.na(X),]
  
  X = data1$X  
  
  
  #Fit a glm to get the MLEs for beta0, beta2 and beta3 - The following code fits the model Y=beta0+beta2*A+beta3*B. We basically need to get the MLEs of the coefficients of the covariates under H0: beta1=0.
  
  fit <- summary(glm(Y~., data=Z, family='gaussian'))
  
  
  #output of the 'fit$coefficients' is the MLEs of betas.
  
  Z1=cbind(1,Z)     #Add the "1" column to Z for beta0
  ncov=ncol(Z1)
  
  fitted.val=0       # fitted.val=hat(beta0)+hat(beta2)*A+hat(beta3)*B
  for (i in 1:ncov){	
    fitted.val=fitted.val+fit$coefficients[i,1]*Z1[,i]
  }
  
  #Fitted value= Z1*beta. Now, find the mean function. 
  
  mean.values=NULL
  
  n=length(fitted.val)
  
  for (i in 1:n){
    
    fitted=fitted.val[i]
    
    mean.values[i]=fitted
    #mean.values[i]=grad(binom.b,fitted)
    #grad(norm.b) takes the gradient of norm.b and plug in "fitted" value. 
    
  }
  return(list(Y=Y,X=X,rd.group=rd.group,mean.values=mean.values))
}

int = 10

X <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/X.txt", sep="\t", header = F)
Y <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Y.txt", sep="\t", header = F)
P <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/P.txt", sep="\t", header = F)
M <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/M.txt", sep="\t", header = T)
Z <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Z.txt", sep="\t", header = F)

rd.group <- M$hrg
X = t(X[int,])
P = P[int,]
colnames(Y) = "Y"
colnames(Z) = c("Z1", "Z2", "Z3")

values=calc_mean(Y,X,Z,rd.group)  

mean.values=values$mean.values
Y=values$Y
X=values$X
rd.group=values$rd.group

X.hrd = X[rd.group == 1]
X.lrd = X[rd.group == 0]

Y.hrd = Y[rd.group == 1]
Y.lrd = Y[rd.group == 0]

mean.values.hrd = mean.values[rd.group == 1]
mean.values.lrd = mean.values[rd.group == 0]

nhrd = length(Y.hrd)
nlrd = length( Y.lrd)

score=sum( (Y-mean.values)*X)

if(RVS %in% c('TRUE','True','true','T','t')) {
  v = sum( (Y.lrd-mean.values.lrd)^2)* var(X.lrd)  + sum( (Y.hrd-mean.values.hrd)^2* as.numeric(calc_robust_Var(P))) 
}else{
  v = sum( (Y.lrd-mean.values.lrd)^2)* var(X.lrd)  + sum( (Y.hrd-mean.values.hrd)^2* var(X.hrd)  )
}
score.test.obs=score^2/v  

X1 = X.hrd - mean(X.hrd)
X2 = X.lrd - mean(X.lrd)  

C = NULL
sc=c()
var = c()

for (j in 1:nboot){
  Xhrd = sample(X1[],nhrd,replace=TRUE)
  Xlrd = sample(X2[],nlrd,replace=TRUE)
  
  score=sum((Y.lrd-mean.values.lrd)*Xlrd) + sum((Y.hrd-mean.values.hrd)*Xhrd)
  
  total.var=sum( (Y.lrd-mean.values.lrd)^2)* var(Xlrd)  + sum( (Y.hrd-mean.values.hrd)^2* var(Xhrd)  ) 
  
  score.test=score^2/total.var	
  sc = c(sc, score)
  var = c(var, total.var)
  
  C =c(C,score.test)
}
cc = (sum(C<=score.test.obs)+1)/(nboot+1)
1-cc
return(cc)




#******************************************************************#

#********** REGULAR SCORE TEST FOR COMMON VARIANT ANALYSIS ***********#

#*******************************************************************#


#rd.group is not actively used in regular score calculations but it is an ingredient in "calc_mean" function, than is why I added "rd.group" in "regScore_Rao_cov"

regScore_Rao=function(Y,X,Z,rd.group){
   
  values=calc_mean(Y,X,Z,rd.group)  
  
  mean.values=values$mean.values
  Y=values$Y
  X=values$X
  
  score=sum( (Y-mean.values)*X)
  v = sum( (Y-mean.values)^2)*var(X)    
  score.test=score^2/v
  cc = 1-pchisq(score.test,1)
  return (cc)
}

regScore_Rao_old=function(Y,X){
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  score=sum( (Y-mean(Y))*X)
  v = sum( (Y-mean(Y))^2)*var(X)    #In anova(Rao),if we do, var(X)*(N-1)/N we get the same result.
  score.test=score^2/v
  cc = 1-pchisq(score.test,1)
  return (cc)
}

#******************************************************************#

regScore_perm=function(Y,X,Z,nperm,rd.group){ 
  	
 	values=calc_mean(Y,X,Z,rd.group)  
  
  	mean.values=values$mean.values
  	Y=values$Y
  	X=values$X

  	test = calc_score_test_cov(X,Y,mean.values)
  	C = NULL
 for (j in 1:nperm){
	k = sample(length(Y))
	YY = Y[k]
    a = calc_score_test_cov(X,YY,mean.values)
    C = c(C,a)
  }
  cc = (sum(C>=test)+1)/(nperm+1)
  return (cc)
}

regScore_perm_zeynep_old=function(Y,X,nperm){
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  test = calc_score_test(X,Y)
  C = NULL
 for (j in 1:nperm){
	k = sample(length(Y))
	YY = Y[k]
    a = calc_score_test(X,YY)
    C = c(C,a)
  }
  cc = (sum(C>=test)+1)/(nperm+1)
  return (cc)
}

#******************************************************************#

calc_score_test = function(X,Y,mean.values){
 
  score=sum( (Y-mean.values)*X)
  v = sum( (Y-mean.values)^2)*var(X)    
  score.test=score^2/v
  return (score.test)
}

calc_score_test_old = function(X,Y){
 
  score=sum( (Y-mean(Y))*X)
  v = sum( (Y-mean(Y))^2)*var(X)    
  score.test=score^2/v
  return (score.test)
}

#******************************************************************#

#*********** RVS TEST FOR COMMON VARIANT ANALYSIS -- ASY ***********#

#*******************************************************************#


RVS_asy=function(Y,X,Z,rd.group,P,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_asy, should be True or False!\n')
    return(NULL)
  }
  
  	values=calc_mean(Y,X,Z,rd.group)  
  
  	mean.values=values$mean.values
  	Y=values$Y
  	X=values$X
  	rd.group=values$rd.group

  
    X.hrd = X[rd.group == 1]
    X.lrd = X[rd.group == 0]

    Y.hrd = Y[rd.group == 1]
    Y.lrd = Y[rd.group == 0]

    mean.values.hrd = mean.values[rd.group == 1]
    mean.values.lrd = mean.values[rd.group == 0]

   	nhrd = length(rd.group==1)
   	nlrd = length(rd.group==0)

  
  score=sum( (Y-mean.values)*X)
  
  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = sum( (Y.lrd-mean.values.lrd)^2)* var(X.lrd)  + sum( (Y.hrd-mean.values.hrd)^2* as.numeric(calc_robust_Var(P)))
  }else{
    v = sum( (Y.lrd-mean.values.lrd)^2)* var(X.lrd)  + sum( (Y.hrd-mean.values.hrd)^2* var(X.hrd)  )
  }
   score.test=score^2/v

  cc = 1-pchisq( score.test,1)
  return (cc)
}


RVS_asy_old=function(Y,X,P,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_asy, should be True or False!\n')
    return(NULL)
  }
  ## run anlysis for non-NA X only
    rd.group = rd.group[!is.na(X)]
  	Y = Y[!is.na(X)]
  	X = X[!is.na(X)]
  
    X.hrd1 = X[rd.group == 0]
    X.hrd2 = X[rd.group == 1]
    X.lrd = X[rd.group == 2]

    Y.hrd1 = Y[rd.group == 0]
    Y.hrd2 = Y[rd.group == 1]
    Y.lrd = X[rd.group == 2]

   	nhrd1 = length(rd.group==0)
  	nhrd2 = length(rd.group==1)
  	nlrd = length(rd.group==2)

  
  score=sum( (Y-mean(Y))*X)
  
  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd1-mean(Y))^2* as.numeric(calc_robust_Var(P))) + sum( (Y.hrd2-mean(Y))^2* as.numeric(calc_robust_Var(P))  )
  }else{
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd1-mean(Y))^2* var(X.hrd1)  )  +sum( (Y.hrd2-mean(Y))^2* var(X.hrd2)  ) 
  }
   score.test=score^2/v

  cc = 1-pchisq( score.test,1)
  return (cc)
}

#******************************************************************#

#******* RVS TEST FOR COMMON VARIANT ANALYSIS -- BOOTSTRAP ********#

#*******************************************************************#


RVS_btrap=function(Y,X,Z,rd.values,P,nboot,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_brap, should be True or False!\n')
    return(NULL)
  }
  
   	values=calc_mean(Y,X,Z,rd.group)  
  
  	mean.values=values$mean.values
  	Y=values$Y
  	X=values$X
  	rd.group=values$rd.group

    X.hrd = X[rd.group == 1]
    X.lrd = X[rd.group == 0]

    Y.hrd = Y[rd.group == 1]
    Y.lrd = Y[rd.group == 0]

    mean.values.hrd = mean.values[rd.group == 1]
    mean.values.lrd = mean.values[rd.group == 0]

   	nhrd = length(Y.hrd)
  	nlrd = length( Y.lrd)

    score=sum( (Y-mean.values)*X)

  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = sum( (Y.lrd-mean.values.lrd)^2)* var(X.lrd)  + sum( (Y.hrd-mean.values.hrd)^2* as.numeric(calc_robust_Var(P))) 
  }else{
    v = sum( (Y.lrd-mean.values.lrd)^2)* var(X.lrd)  + sum( (Y.hrd-mean.values.hrd)^2* var(X.hrd)  )
  }
   score.test.obs=score^2/v  
  
 X1 = X.hrd - mean(X.hrd)
 X2 = X.lrd - mean(X.lrd)  

   C = NULL
   sc=c()
   var = c()
   
  for (j in 1:nboot){
    Xhrd = sample(X1[],nhrd,replace=TRUE)
    Xlrd = sample(X2[],nlrd,replace=TRUE)
     
 	score=sum((Y.lrd-mean.values.lrd)*Xlrd) + sum((Y.hrd-mean.values.hrd)*Xhrd)
	
	total.var=sum( (Y.lrd-mean.values.lrd)^2)* var(Xlrd)  + sum( (Y.hrd-mean.values.hrd)^2* var(Xhrd)  ) 

	score.test=score^2/total.var	
  sc = c(sc, score)
  var = c(var, total.var)
  
    C =c(C,score.test)
  }
  cc = (sum(C<=score.test.obs)+1)/(nboot+1)
  1-cc
  return(cc)
}


RVS_btrap_old=function(Y,X,P,nboot,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_brap, should be True or False!\n')
    return(NULL)
  }
  
    rd.group = rd.group[!is.na(X)]
  	Y = Y[!is.na(X)]
  	X = X[!is.na(X)]
  
    X.hrd1 = X[rd.group == 0]
    X.hrd2 = X[rd.group == 1]
    X.lrd = X[rd.group == 2]

    Y.hrd1 = Y[rd.group == 0]
    Y.hrd2 = Y[rd.group == 1]
    Y.lrd = Y[rd.group == 2]

   	nhrd1 = length(Y.hrd1)
  	nhrd2 = length( Y.hrd2)
  	nlrd = length(Y.lrd)
  	
     score=sum( (Y-mean(Y))*X)

  	  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd1-mean(Y))^2* as.numeric(calc_robust_Var(P))) + sum( (Y.hrd2-mean(Y))^2* as.numeric(calc_robust_Var(P))  )
  }else{
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd1-mean(Y))^2* var(X.hrd1)  )  +sum( (Y.hrd2-mean(Y))^2* var(X.hrd2)  ) 
  }
   score.test.obs=score^2/v  
  
 X1 = X.hrd1 - mean(X.hrd1)
 X2 = X.hrd2 - mean(X.hrd2)  
 X3 = X.lrd - mean(X.lrd)  
 
   C = NULL
  for (j in 1:nboot){
    Xhrd1 = sample(X1[],nhrd1,replace=TRUE)
    Xhrd2 = sample(X2[],nhrd2,replace=TRUE)
    Xlrd = sample(X3[],nlrd,replace=TRUE)
     
 	score=sum((Y.lrd-mean(Y))*Xlrd) + sum((Y.hrd1-mean(Y))*Xhrd1) + sum((Y.hrd2-mean(Y))*Xhrd2)
	
	total.var=sum( (Y.lrd-mean(Y))^2)* var(Xlrd)  + sum( (Y.hrd1-mean(Y))^2* var(Xhrd1)  )  +sum( (Y.hrd2-mean(Y))^2* var(Xhrd2)  ) 

	score.test=score^2/total.var	

    C =c(C,score.test)
  }
  cc = (sum(C<=score.test.obs)+1)/(nboot+1)
  return(cc)
}

#*******************************************************************#

