



#******************************************************************#

#********** REGULAR SCORE TEST FOR COMMON VARIANT ANALYSIS ***********#

#*******************************************************************#


regScore_Rao_zeynep=function(Y,X){
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  score=sum( (Y-mean(Y))*X)
  v = sum( (Y-mean(Y))^2)*var(X)    #In anova(Rao),if we do, var(X)*(N-1)/N we get the same result.
  score.test=score^2/v
  cc = 1-pchisq(score.test,1)
  return (cc)
}



regScore_Rao = function(Y,X){
  ## run anlysis for non-NA X only
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  a = (glm(Y~1,family='binomial'))
  b = glm(Y~X,family='binomial')
  x = anova(a,b,test='Rao')
  cc = 1-pchisq(x$Rao[2],1)
  #res<-list(est=coef(summary(b))[2,1],
  #          Rao=x$Rao[2],
  #          p_Rao=x[2,6])
  return (cc)
}



#******************************************************************#

regScore_perm_zeynep=function(Y,X,nperm){
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  test = calc_score_test(X,Y)
  C = NULL
 for (j in 1:nperm){
	k = sample(length(Y))
	YY = Y[k]
    a = calc_score_testp(X,YY)
    C = c(C,a)
  }
  cc = (sum(C>=test)+1)/(nperm+1)
  return (cc)
}

regScore_perm = function(Y,X,nperm){
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  ncase = sum(Y==1)  ## number of case with X!=NA
  ncont = sum(Y==0)  ## number of cont with X!=NA
  X1=X[Y==1]
  X2=X[Y==0]
  test = calc_score_test(X1,X2)
  C = NULL  
    for (j in 1:nperm){
    k = sample(length(X))
    X = X[k]
    X11=X[Y==1]
    X21=X[Y==0]
    a = calc_score_test(X11,X21)
    C = c(C,a)
  }

  cc = (sum(C>=test)+1)/(nperm+1)
  return (cc)
}
#******************************************************************#

calc_score_test_zeynep = function(X,Y){
 
  score=sum( (Y-mean(Y))*X)
  v = sum( (Y-mean(Y))^2)*var(X)    
  score.test=score^2/v
  return (score.test)
}

calc_score_test = function(M1,M2){
  X.new = c(M1,M2)
  Y.new = c(rep(1,length(M1)),rep(0,length(M2)))
  p = length(M1)/length(X.new)
  q = 1 - p
  S = q*sum(M1)-p*sum(M2)
  vs = p*q*(length(X.new))*var(X.new)
  x = S^2/vs
  return (x)
}





#******************************************************************#

#*********** RVS TEST FOR COMMON VARIANT ANALYSIS -- ASY ***********#

#*******************************************************************#

RVS_asy=function(Y,X,P,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_asy, should be True or False!\n')
    return(NULL)
  }
  ## run anlysis for non-NA X only
  
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  a = (glm(Y~1,family='binomial'))
  b = glm(Y~X,family='binomial')
  p = length(X[Y==1])/length(X)
  q = 1 - p
  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = q*as.numeric(calc_robust_Var(P))/var(X) + p*var(X[Y==0],na.rm=TRUE)/var(X)
  }else{
    v = q*var(X[Y==1]) + p*var(X[Y==0])  ## regular estimation for variance on case at var(X[Y==1]) from test_RVS_asy
  }

  x = anova(a,b,test='Rao')
  x_rvs=x$Rao[2]/v  ## var(S_inRao)=#case*#cont/#sample*var(X), var(RVS)=#case*#cont/#sample*(q*var(Xcase)+p*var(Xcont))
  cc = 1-pchisq(x_rvs,1)
  return (cc)
}

RVS_asy_zeynep=function(Y,X,P,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_asy, should be True or False!\n')
    return(NULL)
  }
  ## run anlysis for non-NA X only
    rd.group = rd.group[!is.na(X)]
  	Y = Y[!is.na(X)]
  	X = X[!is.na(X)]
  
    X.hrd = X[Y == 1]
    X.lrd = X[Y == 0]

    Y.hrd = Y[Y == 1]
    Y.lrd = Y[Y == 0]

   	nhrd1 = length(Y==1)
   	nlrd = length(Y==0)

  score=sum( (Y-mean(Y))*X)
  
  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd-mean(Y))^2) * as.numeric(calc_robust_Var(P))
  }else{
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd-mean(Y))^2* var(X.hrd)  )
  }
  
   score.test=score^2/v
  
  cc = 1-pchisq( score.test,1)
  return (cc)
}

X = X[i,]
P = P[i,]

X <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/X.txt", sep="\t", header = F)
Y <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Y.txt", sep="\t", header = F)
P <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/P.txt", sep="\t", header = F)
M <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/M.txt", sep="\t", header = T)
Z <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Z.txt", sep="\t", header = T)

rd.group <- M$hrg


for (i in 1:30){
 # print(RVS_asy(Y,X[i,],P[i,]))
  
  print(RVS_btrap_zeynep(Y,X[,i],P[i,], 3000))
  
}



#******************************************************************#

#******* RVS TEST FOR COMMON VARIANT ANALYSIS -- BOOTSTRAP ********#

#*******************************************************************#
RVS_btrap = function(Y,X,P,nboot,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_brap, should be True or False.\n')
    return(NULL)
  }
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  maf=sum(X)/(length(X)*2); maf0=min(maf,1-maf)
  if(maf0<0.05){cat('Warning: MAF of the SNP is',maf0,', try to group with other rare variants and use RVS_rare1.\n')}
  
  ncase1 = sum(Y==1)
  ncont1 = sum(Y==0)
  p = length(X[Y==1])/length(X)
  q = 1 - p
  #  if(RVS %in% c('TRUE','True','true','T','t')) {
  vcase = as.numeric(calc_robust_Var(P))
  #  } else{
  #    vcase = var(X[Y==1])
  #  }
  vcont = var(X[Y==0])
  if(RVS %in% c('TRUE','True','true','T','t')) {
    Tobs = (mean(X[Y==1]) - mean(X[Y==0]))/sqrt(vcase/ncase1+vcont/ncont1)
  }else{
    X1=X[Y==1]; X2=X[Y==0]
    Tobs=calc_score_test(X1,X2) 
  }  
  X1 = X[Y==1] - mean(X[Y==1])
  X2 = X[Y==0] - mean(X[Y==0])
  C = NULL
  if(RVS %in% c('TRUE','True','true','T','t')) {
    for (j in 1:nboot){
      Xca = sample(X1[],ncase1,replace=TRUE)
      Xco = sample(X2[],ncont1,replace=TRUE)
      vcase = var(Xca)
      vcont = var(Xco)
      C =c(C,(mean(Xca) - mean(Xco))/sqrt(vcase/ncase1+vcont/ncont1))
    }## enf for j
  }else{
    for(j in 1:nboot){
      k=sample(length(X));
      Xnew=X[k]
      X1=Xnew[Y==1];X2=Xnew[Y==0]
      C=c(C,calc_score_test(X1,X2))
    } ## end for j
  }## else if RVS
  cc = (sum(abs(C)>=abs(Tobs))+1)/(nboot+1)
  return(cc)
  
}

RVS_btrap_zeynep=function(Y,X,P,nboot,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_brap, should be True or False!\n')
    return(NULL)
  }
  
    Y = Y[!is.na(X)]
    X = X[!is.na(X)]
  
  
    X.hrd = X[Y == 1]
    X.lrd = X[Y == 0]

    Y.hrd = Y[Y == 1]
    Y.lrd = Y[Y == 0]

   	nhrd = length(Y.hrd)
  	nlrd = length(Y.lrd)
  	
     score=sum( (Y-mean(Y))*X)

  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd-mean(Y))^2* as.numeric(calc_robust_Var(P)))
  }else{
    v = sum( (Y.lrd-mean(Y))^2)* var(X.lrd)  + sum( (Y.hrd-mean(Y))^2* var(X.hrd)  )  
  }
   score.test.obs=score^2/v  
  
 X1 = X.hrd - mean(X.hrd)
 X3 = X.lrd - mean(X.lrd)  
 
   C = NULL
  for (j in 1:nboot){
    Xhrd = sample(X1[],nhrd,replace=TRUE)
    Xlrd = sample(X3[],nlrd,replace=TRUE)
     
 	score=sum((Y.lrd-mean(Y))*Xlrd) + sum((Y.hrd-mean(Y))*Xhrd) 
	
	total.var=sum( (Y.lrd-mean(Y))^2)* var(Xlrd)  + sum( (Y.hrd-mean(Y))^2* var(Xhrd)  )

	score.test=score^2/total.var	

    C =c(C,score.test)
  }
  cc = (sum(C<=score.test.obs)+1)/(nboot+1)
  return(1-cc)
}



#*******************************************************************#


X <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/X.txt", sep="\t", header = F)
Y <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Y.txt", sep="\t", header = F)
P <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/P.txt", sep="\t", header = F)
M <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/M.txt", sep="\t", header = T)


for (i in 1:30){
  #print(RVS_btrap(Y,X[i,],P[i,], 10000))
  
  print(RVS_btrap_zeynep(Y,X[i,],P[i,], 10000))
  
}
