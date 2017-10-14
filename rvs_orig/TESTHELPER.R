calc_robust_Var = function(P){
  Sq = 4*P[3] + P[2]
  Sm = 2*P[3] + P[2]
  S = Sq - Sm^2
  return (S)
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