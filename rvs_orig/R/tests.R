### in this file, include the simplified program for parallel computing, also similar to the setting of lm

#' RVS using asymptotic distribution for score test statistic
#' 
#' This functions actually includes two association tests which differ only at the estimation of the variance of score statistic by option 'method'.
#' It uses regular formula for \eqn{Var_{case}(E(G_{ij}|D_{ij}))} when method = 'Regular' (also called likelihood method in the paper),
#' \eqn{Var_{case}(G_{ij})}  is used for \eqn{Var_{case}(E(G_{ij}|D_{ij}))} when method = 'RVS'. It uses the population frequency P
#' to estimate variance for case if method=RVS, otherwise use the expected genotype probability for case directly.
#'
#' Note, both function scaled the variance in equation (1) in Appendix A by dividing N_case*Ncontrol/N
#' also note the test statistics in anova gives us Rao=s^2/var(s^2), where RVS=s^2/robust var(s^2)
#' so RVS=Rao *var(s^2)/robust var(s^2), s=sum(y_j-mean(y))E(G_ij|D_ij), therefore var(s)//var(E(Gij|Dij)), var(X) in code.
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/22570057}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/24733292}
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return   p-values for the snp
#' @export
RVS_asy = function(Y,X,P,RVS='TRUE'){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_asy, should be True or False.\n')
    return(NULL)
  }
  ## run anlysis for non-NA X only
  
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  a = (glm(Y~1,family='binomial'))
  b = glm(Y~X,family='binomial')
  p = length(X[Y==1])/length(X)
  q = 1 - p
#  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = q*as.numeric(calc_robust_Var(P)) + p*var(X[Y==0],na.rm=TRUE)
 # }else{
 #   v = q*var(X[Y==1]) + p*var(X[Y==0])  ## regular estimation for variance on case at var(X[Y==1]) from test_RVS_asy
 # }
  x = anova(a,b,test='Rao')
  if(RVS %in% c('TRUE','True','true','T','t')) {
    x_rvs=x$Rao[2]*var(X)/v  ## var(S_inRao)=#case*#cont/#sample*var(X), var(RVS)=#case*#cont/#sample*(q*var(Xcase)+p*var(Xcont))
  }else{
    x_rvs=x$Rao[2]
  }
  cc = 1-pchisq(x_rvs,1)
 # res<-list(coefficients=coef(summary(b))[2,1],
#       score=x$Rao[2],
 #      p_Rao=x[2,6],
#       RVS=x_rvs,
#       p_rvs=cc)
  return (cc)
}


#' use RVS to test associaton by bootrtrap, given phenotype (Y), expected values of genotypes for case and controls (X)
#' estimated genotype frequency (P) and number of times of bootstrap (nboot)
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param nboot - number of bootstrap
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return   p-values for the snp
#' @export
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


#' Use regular Score/Trend test for association for given phenotype (Y) and expected values of genotypes for case and controls (X).
#' it is evaluation by asymptotic distribution
#'
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @return p-values the variant (double)
#' @export
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

#' Use regular Score/Trend test for association for given phenotype (Y) and expected values of genotypes for case and controls (X)
#' it is evaluation by permutation distribution
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @param nperm - number of permutations
#' @return p-values the variant (double)
#' @export
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


# RVS analysis with rare variants for one group
# Robust Variance Estimate
# Get p-values from RVS with CAST and C-alpha (Resampling: Bootstrap, Variance Estimate: Robust)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values, it has dimensions n by J, J - number of variants, n - sample size
#' @param nboot - number of bootstrap
#' @param njoint - number of SNPs grouping together for one test, default is 5
#' @param P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return two p-values for CAST and C-alpha
#' @export
RVS_rare1= function(Y,X,P,njoint=5,nboot,RVS='TRUE', hom=1,multiplier=1,snp_loop=1){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_rare, should be True or False.\n')
    return(NULL)
  }

  X1 = as.matrix(X[Y==1,])
  X2 = as.matrix(X[Y==0,])
  S = calc_ScoreV(X,Y)
  if(RVS %in% c('TRUE','True','true','T','t'))
  {     Sigma = calc_ScoreV_RVS(X,Y,P,hom,multiplier,snp_loop)
  }else{
        Sigma = calc_ScoreV_likely(X,Y,hom,multiplier,snp_loop)
  }
  w = rep(1,ncol(X))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  Q = NULL
  L = NULL

  for (i in 1:nboot){
    Xs = sample_bootstrap(X1,X2)
    Xa = Xs$Xca
    Xb = Xs$Xco
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


# RVS analysis with rare variants for one group
# Robust Variance Estimate
# Get p-values from RVS with CAST and C-alpha (Resampling: Bootstrap, Variance Estimate: Robust)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values, it has dimensions n by J, J - number of variants, n - sample size
#' @param nboot - number of bootstrap
#' @param njoint - number of SNPs grouping together for one test, default is 5
#' @param P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return two p-values for CAST and C-alpha
#' @export
RVS_rare2= function(Y,X,P,njoint=5,nboot,RVS='TRUE', hom=1,multiplier=1,snp_loop=1){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_rare, should be True or False.\n')
    return(NULL)
  }
  
  X1 = as.matrix(X[Y==1,])
  X2 = as.matrix(X[Y==0,])
  S = calc_ScoreV(X,Y)
  if(RVS %in% c('TRUE','True','true','T','t'))
  {     Sigma = calc_ScoreV_RVS_test(X,Y,P,hom,multiplier,snp_loop)
  }else{
    Sigma = calc_ScoreV_likely(X,Y,hom,multiplier,snp_loop)
  }
  w = rep(1,ncol(X))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  Q = NULL
  L = NULL
  
  for (i in 1:nboot){
    Xs = sample_bootstrap(X1,X2)
    Xa = Xs$Xca
    Xb = Xs$Xco
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

# RVS analysis with rare variants
# Robust Variance Estimate
# Get p-values from RVS with CAST and C-alpha (Resampling: Bootstrap, Variance Estimate: Robust)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values, it has dimensions n by J, J - number of variants, n - sample size
#' @param nboot - number of bootstrap
#' @param njoint - number of SNPs grouping together for one test, default is 5
#' @param P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return two p-values for CAST and C-alpha
#' @export
RVS_rare= function(Y,X,P,njoint=5,nboot,RVS='TRUE',hom=1,multiplier=1){
  if(!RVS %in% c('TRUE','True','true','T','t','FALSE','False','false','F','f')) {
    cat('Wrong input for option RVS in RVS_rare, should be True or False.\n')
    return(NULL)
  }
  
  nsnp=nrow(P)
  if(nsnp%%njoint!=0) {nloop=(nsnp-nsnp%%njoint)/njoint ## the last couple of snps(<5) combined to the last group
  }else {nloop=nsnp/njoint}
  
  p.RVS.rare=NULL
  for(i in 1:nloop){
    t1<-(i-1)*njoint+1; 
    if(i < nloop) {t2<-i*njoint;
    }else{t2 <- nsnp}
    Xsub=X[,t1:t2]
    
    tt=RVS_rare1(Y,X[,t1:t2],P[t1:t2,],njoint,nboot,RVS,hom=hom,multiplier=multiplier,snp_loop=i)
    p.RVS.rare<-rbind(p.RVS.rare,tt)
  }
  
  row.names(p.RVS.rare)=paste0('loop',1:nloop)
  p.RVS.rare=data.frame(p.RVS.rare)
  colnames(p.RVS.rare)=c('p.CAST','p.Calpha')
  return (p.RVS.rare)
}



# regScore test for J joint rare variants
# Get p-values from regScore with CAST and C-alpha (Resampling: permutation)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values/genotype calls for J variants, it has dimensions n by J, J - number of variants, n - sample size
#' @param njoint number of rare variants to be grouped together
#' @param nperm - number of bootstrap
#' @return two p-values for CAST and C-alpha
#' @export
regScore_rare1 = function(Y,X,njoint=5,nperm){
  S = calc_ScoreV(X,Y)
  Sigma = calc_ScoreV_regVar(X,Y)
  w = rep(1,length(X[1,]))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  Q = NULL
  L = NULL

  for (i in 1:nperm){
    k = sample(length(Y))
    YY = Y[k]
    S = calc_ScoreV(X,YY)
    Sigma = calc_ScoreV_regVar(X,YY)
    w = rep(1,length(X[1,]))
    L = c(L,test_CAST(S,Sigma,w))
    Q = c(Q,test_Calpha(S,Sigma))
  }
  pl = (sum(L<=SLobs)+1)/(nperm+1)
  pQ = (sum(Q<=SQobs)+1)/(nperm+1)
  return (c(pl,pQ))
}


# regScore test for rare variants, each p-value is for J rare varaints
# Get p-values from regScore with CAST and C-alpha (Resampling: permutation)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values/genotype calls for J variants, it has dimensions n by J, J - number of variants, n - sample size
#' @param njoint number of rare variants to be grouped together
#' @param nperm - number of bootstrap
#' @return two p-values for CAST and C-alpha
#' @export
regScore_rare = function(Y,X,njoint,nperm){

nsnp=ncol(X)
if(nsnp%%njoint!=0) {nloop=(nsnp-nsnp%%njoint)/njoint ## the last couple of snps(<5) combined to the last group
}else {nloop=nsnp/njoint}

p.geno.rare=matrix(NA,nrow=nloop,ncol=2)

for(k in 1:nloop){
  t1<-(k-1)*njoint+1;
  if(k < nloop) {t2<-k*njoint;
  }else{t2 <- nsnp}
  Xtmp=X[,t1:t2]
  p.geno.rare[k,]<-regScore_rare1(Y,Xtmp,njoint,nperm)
}
  p.geno.rare<-data.frame(p.geno.rare)
  colnames(p.geno.rare)<-c('p.CAST','p.Calpha')
  rownames(p.geno.rare)<-paste0('loop',1:nloop)
   return (p.geno.rare)
}
