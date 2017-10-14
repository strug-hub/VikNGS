
#******************************************************************#

#********** REGULAR SCORE TEST FOR RARE VARIANT ANALYSIS ***********#

#*******************************************************************#
calc_ScoreV_zeynep=function (X, Y) 
{
    L = ncol(X)
    S = NULL
    for (i in 1:L) {
        Yn = Y[!is.na(X[, i]),]
        Xn = X[!is.na(X[, i]), i]
        s = sum((Yn - mean(Yn)) * Xn)
        S = c(S, s)
    }
    S = t(S)
    return(S)
}
calc_ScoreV = function(X,Y){ 
  L = ncol(X)
  S = NULL
  for (i in 1:L){
    Yn = Y[!is.na(X[,i]),]
    Xn = X[!is.na(X[,i]),i]
    xca = Xn[Yn==1]
    xco = Xn[Yn==0]
    my = mean(Yn)
    s = sum(xca,na.rm=T)*(1-my) - my*sum(xco,na.rm=T)
    S = c(S,s)
  }
  S = t(S)
  return (S)
}


calc_ScoreV_regVar_zeynep=function (X, Y)
{
	Ym=NULL

	for (i in 1:ncol(X)){
    Yn = Y[!is.na(X[,i])]    
    Ym[i] = sqrt( (sum((Yn - mean(Yn))^2)/length(Yn))*length(Y))  #length(Yn))*length(Y) term is like a correction to the different sample sizes in X columns after NAs are removed

  	}
 #  X = check_hom_matrix(X, hom,multiplier)

  	diag_S = diag(Ym) %*% cov(X,use="pairwise.complete.obs")%*%diag(Ym) 

  	return(diag_S)
}

calc_ScoreV_regVar = function(X,Y){

  X1 = X[Y==1,]
  X2 = X[Y==0,]
  l1 = nrow(X1)
  l2 = nrow(X2)
  J = ncol(X1)

  a =colSums(is.na(X1),na.rm=T)
  b =colSums(is.na(X2),na.rm=T)

  ncase = rep(l1,J) - a
  ncont = rep(l2,J) - b
  nn = ncase+ncont
  L= length(Y)
  Xm = cov(X,use="pairwise.complete.obs")
  vs  = sqrt(diag(ncase*ncont/nn^2))
  diag_S  = vs%*%Xm%*%vs*L
  return (diag_S)
}

#******************************************************************#

#************** RVS TEST FOR RARE VARIANT ANALYSIS ***************#

#*******************************************************************#





sample_bootstrap_zeynep=function(X,Y) {
 
 data1=data.frame(X,Y)
 data1.split=split(as.data.frame(data1),data1$V1)
 cluster.no=length(data1.split)

 nX=ncol(X)
 for (i in 1:nX) {
 	 	for (j in 1:cluster.no) 	{	
  		X.cluster=calc_centralize_matrix(as.data.frame(data1.split[[j]][,i]))
  		nX.cluster=length(X.cluster)
  		row.no=sample(1:nX.cluster,nX.cluster,replace=TRUE)
   		data1.split[[j]][,i]=X.cluster[row.no] 	
 			}
 	}
 	data2=unsplit(data1.split,data1$V1)
 	X.sampled=data2[,1:nX]
  return(X.sampled)
}


sample_bootstrap = function(X1,X2){
  case = nrow(X1)
  cont = nrow(X2)
  X1 = calc_centralize_matrix(X1)
  X2 = calc_centralize_matrix(X2)
  ca = sample(1:case,case,replace=TRUE)
  co = sample(1:cont,cont,replace=TRUE)
  Xca = as.matrix(X1[ca,])
  Xco = as.matrix(X2[co,])
  X = rbind(Xca,Xco)
  #return (X)
  return(list('Xca'=Xca,'Xco'=Xco))
}

#*********************************************************************#

calc_ScoreV_likely_zeynep = function(X,Y,hom,multiplier){
  
    X.lrd = X[Y == 0, ]
    X.hrd = X[Y == 1, ]

    Y.lrd = Y[Y == 0]
    Y.hrd = Y[Y == 1]

    X.lrd = check_hom_matrix(X.lrd, hom = hom, multiplier = multiplier)
    X.hrd = check_hom_matrix(X.hrd, hom = hom, multiplier = multiplier)

	Ym.lrd=Ym.hrd=NULL

 for (i in 1:ncol(X)){
    
    Yn = Y[!is.na(X[,i]),]
    Yn.lrd = Y.lrd[!is.na(X.lrd[,i])]    
    Yn.hrd = Y.hrd[!is.na(X.hrd[,i])]    
    
    Ym.lrd[i] = sqrt( (sum( (Yn.lrd - mean(Yn))^2)/length(Yn.lrd))*length(Y.lrd)) 
    Ym.hrd[i] = sqrt( (sum((Yn.hrd - mean(Yn))^2)/length(Yn.hrd))*length(Y.hrd)) 
	
 }
 
  diag_S = diag(Ym.lrd) %*% cov(X.lrd,use="pairwise.complete.obs") %*%diag(Ym.lrd) + diag(Ym.hrd) %*% cov(X.hrd,use="pairwise.complete.obs")%*%diag(Ym.hrd)  
  return (diag_S)
}


calc_ScoreV_likely = function(X,Y,hom,multiplier,snp_loop){
  X1=X[Y==1,]
  X2=X[Y==0,]
  l1 = nrow(X1)
  l2 = nrow(X2)
  J = ncol(X1)

  a =colSums(is.na(X1),na.rm=T)
  b =colSums(is.na(X2),na.rm=T)

  ncase = rep(l1,J) - a
  ncont = rep(l2,J) - b
  nn = ncase+ncont

  Yhat = mean(Y)
  L= length(Y)
  p = l2/l1
  q = l1/l2
  X1=check_hom_matrix(X1,hom=hom,multiplier=multiplier)
  X2=check_hom_matrix(X2,hom=hom,multiplier=multiplier)
  Xm1 = cov(X1,use="pairwise.complete.obs")*(l1-1)
  Xm2 = cov(X2,use="pairwise.complete.obs")*(l2-1)

  vs  = sqrt(diag(ncase*ncont/nn^2))
  diag_S  = vs%*%(p*Xm1+ q*Xm2)%*%vs
  return (diag_S)
}

#*********************************************************************#
calc_ScoreV_RVS_zeynep = function(X,Y,P,hom,multiplier){
  
  X.lrd = X[Y == 0, ]
  X.hrd = X[Y == 1, ]

  Y.lrd = Y[Y == 0]
  Y.hrd = Y[Y == 1]

  L = nrow(P)
	V = rep(NA,L)
	for (i in 1:L){
  	V[i] = calc_robust_Var(P[i,])
	}
	V = diag(sqrt(as.numeric(V)))

  X.lrd = check_hom_matrix(X.lrd, hom = hom, multiplier = multiplier)
  X.hrd = check_hom_matrix(X.hrd, hom = hom, multiplier = multiplier)

	Ym.lrd=Ym.hrd=NULL

 	for (i in 1:ncol(X)){
    
    Yn = Y[!is.na(X[,i]),]
    Yn.lrd = Y.lrd[!is.na(X.lrd[,i])]    
    Yn.hrd = Y.hrd[!is.na(X.hrd[,i])]    

	  Ym.lrd[i] = sqrt( (sum( (Yn.lrd - mean(Yn))^2)/length(Yn.lrd))*length(Y.lrd)) 
    Ym.hrd[i] = sqrt( (sum((Yn.hrd - mean(Yn))^2)/length(Yn.hrd))*length(Y.hrd)) 
  	
  	}
 
  cor.X.hrd = cor.X.function(X.hrd, hom = hom, multiplier = multiplier)
  vargroupHRD = t(V) %*% cor.X.hrd %*% V 

 	diag_S = diag(Ym.lrd) %*% cov(X.lrd,use="pairwise.complete.obs") %*%diag(Ym.lrd) + diag(Ym.hrd) %*% vargroupHRD %*%diag(Ym.hrd)  

  return (diag_S)
}



calc_ScoreV_RVS = function(X,Y,P,hom,multiplier,snp_loop){
  X1=X[Y==1,]
  X2=X[Y==0,]
  L = nrow(P)
  V = rep(NA,L)
  for (i in 1:L){
    V[i] = calc_robust_Var(P[i,])
  }
  V = diag(sqrt(as.numeric(V)))
  l1 = nrow(X1)
  l2 = nrow(X2)
  Yhat = mean(Y)
  L= length(Y)
  p = l2/l1
  q = l1/l2
  
  #zeynep
  X1=check_hom_matrix(X1,hom=hom,multiplier=multiplier)
  X2=check_hom_matrix(X2,hom=hom,multiplier=multiplier)
  cor.X=cor.X.function(X1,hom=hom,multiplier=multiplier,snp_loop)
  
  Sgcase = t(V)%*%cor.X%*%V*l1
  Xm2 = calc_centralize_matrix(X2)
  Xm2[is.na(Xm2)]<-0  ## if not set NA to 0, t(Xm2)%*%Xm2 will have too many NAs
  ## I have checked, t(Xm2)%*%Xm2 is the same as cov(X2) here.
  vs  = var(Y)
  diag_S  = as.numeric(vs)*(p*Sgcase + q*t(Xm2)%*%Xm2)
  return (diag_S)
}




#*********************************************************************#


RVS_rare1= function(Y,X,P,njoint=5,nboot,RVS='TRUE', hom=1,multiplier=1,snp_loop=1){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_rare, should be True or False.\n')
    return(NULL)
  }

  X1 = as.matrix(X[Y==1,])
  X2 = as.matrix(X[Y==0,])
  S = calc_ScoreV(X,Y)
  if(RVS %in% c('TRUE','True','true','T','t')){
    Sigma = calc_ScoreV_RVS(X,Y,P,hom,multiplier,snp_loop)  #(*)
  }else{
        Sigma = calc_ScoreV_likely(X,Y,hom,multiplier,snp_loop) # (*)
  }
  
  
  w = rep(1,ncol(X))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  
  Q = NULL
  L = NULL

  for (i in 1:nboot){
    Xs = sample_bootstrap(X1,X2) # (*)

    Xa = Xs$Xca  # (*)
    Xb = Xs$Xco  # (*) 
    X = rbind(Xa,Xb)
    S = calc_ScoreV(X,Y)
    Sigma = calc_ScoreV_likely(X,Y,hom,multiplier)
    w = rep(1,ncol(X))
    L = c(L,test_CAST(S,Sigma,w))
    Q = c(Q,test_Calpha(S,Sigma))
  }
  pl = (sum(L<=SLobs)+1)/(nboot+1)
  pQ = (sum(Q<=SQobs)+1)/(nboot+1)
  preturn=c(pl,pQ)
  return (preturn)
}

#Changes in "RVS_rare1" : I put (*) beside the lines
# Delete "snp_loop" in functions "calc_ScoreV_RVS" and "calc_ScoreV_likely"

RVS_rare1_zeynep=function(Y,X,P,njoint=5,nboot,RVS='TRUE', hom=1,multiplier=1){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_rare, should be True or False!\n')
    return(NULL)
  }
  
  S = calc_ScoreV(X,Y)
  if(RVS %in% c('TRUE','True','true','T','t')){
    Sigma = calc_ScoreV_RVS_zeynep(X,Y,P,hom,multiplier) #(*)
  }else{
    Sigma = calc_ScoreV_likely_zeynep(X,Y,hom,multiplier)  #(*)
  }
  
  
  w = rep(1,ncol(X))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  Q = NULL
  L = NULL
  
  for (i in 1:nboot){
    Xs = sample_bootstrap_zeynep(X,Y)
    
    S = calc_ScoreV(Xs,Y)
    
    Sigma = calc_ScoreV_likely_zeynep(Xs,Y,hom,multiplier)
    w = rep(1,ncol(Xs))
    L = c(L,test_CAST(S,Sigma,w))
    Q = c(Q,test_Calpha(S,Sigma))
  
    
  }
  pl = (sum(L<=SLobs)+1)/(nboot+1)
  pQ = (sum(Q<=SQobs)+1)/(nboot+1)
  preturn=c(pl,pQ)
  
  #print(mean(L))
  
  return (preturn)
}

for (i in 1:6){
  
  X <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/X.txt", sep="\t", header = F)
  Y <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/Y.txt", sep="\t", header = F)
  P <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/P.txt", sep="\t", header = F)
  M <- read.csv("C:/Users/Scott/Desktop/RVS-master/example/M.txt", sep="\t", header = T)
  
  
  njoint=5
  t1<-(i-1)*njoint+1;
  t2<-i*njoint;
  X <- t(X)
  X <- X[,t1:t2]
  P <- P[t1:t2,]
#  print(RVS_rare1(Y,X,P, nboot=10000))
  print(RVS_rare1_zeynep(Y,X,P, nboot=1000))
  
}

#*********************************************************************#

check_hom_matrix_zeynep=function(X,hom=1,multiplier=1)
{
  if(!hom%in%c(1,2))
  {cat('Wrong input for hom input in check_hom_matrix, hom need be 1 or 2.');
    return;
  }
   var0=apply(X, MARGIN = 2, function(x) {  sum(max(x, na.rm=TRUE)-min(x,na.rm=TRUE)==0)})
 ## column index for those all column are the same, might be more than one column
   var0.index=which(var0==1)
  #check if there are any problematic columns,
  if (sum(var0)>0) { 
    # this statement is used if there are column that the whole column has the same genotype number.
    na.index=apply(X, MARGIN = 1, function(x) { sum(is.na(x))==0})
    qq= which(na.index=='TRUE')
    if(hom==1) {j=qq[1]
    }else{j=sample(qq,1)
    }
    if(any(X[j,var0.index]==0)) {ind.0=which(X[j,var0.index]==0);X[j,var0.index[ind.0]]=10^(-15) 
    }
	ind.non0=which(X[j,var0.index]!=0)
    if(multiplier==1){
    X[j,var0.index[ind.non0]]=X[j,var0.index[ind.non0]]/2  #non-zero columns replace it by half
    }else{
      X[j,var0.index[ind.non0]]=X[j,var0.index[ind.non0]]*2  #non-zero columns replace it by 2 times
     if(any(X[j,var0.index[ind.non0]]>2)) {X[j,var0.index[ind.non0]>2]=2}
    }
  }  
  return(X)
}


check_hom_matrix<-function(X,hom=1,multiplier=1)
{
  if(!hom%in%c(1,2))
  {cat('Wrong input for hom input in check_hom_matrix, hom need be 1 or 2.');
    return;
  }
   var0=apply(X, MARGIN = 2, function(x) {  sum(max(x, na.rm=TRUE)-min(x,na.rm=TRUE)==0)})
 ## column index for those all column are the same, might be more than one column
   var0.index=which(var0==1)
  #check if there are any problematic columns,
  if (sum(var0)>0) { 
    # this statement is used if there are column that the whole column has the same genotype number.
    na.index=apply(X, MARGIN = 1, function(x) { sum(is.na(x))==0})
    qq= which(na.index=='TRUE')
    if(hom==1) {j=qq[1]
    }else{j=sample(qq,1)
    }
    ind.0=which(X[j,var0.index]==0)
    ind.non0=which(X[j,var0.index]!=0)
    X[j,var0.index[ind.0]]=10^(-15)    #if the column has all zeros , replace the first row with 10^-15

    if(multiplier==1){
    X[j,var0.index[ind.non0]]=X[j,var0.index[ind.non0]]/2  #non-zero columns replace it by half
    }else{
      X[j,var0.index[ind.non0]]=X[j,var0.index[ind.non0]]*2  #non-zero columns replace it by 2 times
     if(any(X[j,var0.index[ind.non0]]>2)) {X[j,var0.index[ind.non0]>2]=2}
    }
 #   if(X[j,var0.index[ind.non0]]>2) {X[j,var0.index[ind.non0]]=2}
  }  
  return(X)
}

#*********************************************************************#






