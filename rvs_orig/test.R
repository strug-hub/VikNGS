calc_robust_Var = function(P){
  Sq = 4*P[3] + P[2]
  Sm = 2*P[3] + P[2]
  S = Sq - Sm^2
  return (S)
}


setwd("C:/Users/Scott/Desktop/RVS-master/")
geno = vcf_process()

Y = geno[["Y"]]
X = geno[["Geno"]]
P = geno[["P"]]

X = X[,i]
P = P[i,]

for(i in 1:ncol(X))
  print(RVS_asy(Y,X[,i],P[i,],RVS='TRUE'))


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
    cat('Wrong input for option RVS in RVS_asy, should be True or False!\n')
    return(NULL)
  }
  ## run anlysis for non-NA X only
  
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  maf=sum(X)/(length(X)*2)
  maf0=min(maf,1-maf)
  
  if(maf0<0.05){cat('Warning: MAF of the SNP is',maf0,', try to group with other rare variants and using RVS_rare1.\n')}
  
  a = (glm(Y~1,family='binomial'))
  b = glm(Y~X,family='binomial')
  p = length(X[Y==1])/length(X)
  q = 1 - p
  if(RVS %in% c('TRUE','True','true','T','t')) {
    v = q*as.numeric(calc_robust_Var(P)) + p*var(X[Y==0],na.rm=TRUE)
  }else{
    v = q*var(X[Y==1]) + p*var(X[Y==0])  ## regular estimation for variance on case at var(X[Y==1]) from test_RVS_asy
  }
  x = anova(a,b,test='Rao')
  x_rvs=x$Rao[2]*var(X)/v  ## var(S_inRao)=#case*#cont/#sample*var(X), var(RVS)=#case*#cont/#sample*(q*var(Xcase)+p*var(Xcont))
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

  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  maf=sum(X)/(length(X)*2)
  maf0=min(maf,1-maf)
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

calc_score_test = function(M1,M2){
  X = c(M1,M2)
  Y = c(rep(1,length(M1)),rep(0,length(M2)))
  p = length(M1)/length(X)
  q = 1 - p
  S = q*sum(M1)-p*sum(M2)
  vs = p*q*(length(X))*var(X)
  x = S^2/vs
  return (x)
}

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

library('MASS')
library('CompQuadForm')

setwd("C:/Users/Scott/Desktop/RVS-master/")
geno = vcf_process()

Y = geno[["Y"]]
X = geno[["Geno"]]
P = geno[["P"]]

nsnp = geno[["nsnp"]]
RVS='TRUE'
njoint=5
hom=1
multiplier=1
nboot=100000
i=1
t1<-(i-1)*njoint+1;
t2<-i*njoint;
X <- X[,t1:t2]
P <- P[t1:t2,]

RVS_rare1(Y,X,P,njoint,nboot,RVS,hom,multiplier)
  
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
  
  nsnp=nrow(P)
  if(nsnp%%njoint!=0) {nloop=(nsnp-nsnp%%njoint)/njoint ## the last couple of snps(<5) combined to the last group
  }else {nloop=nsnp/njoint}
  
  p.RVS.rare=NULL
  for(i in 1:nloop){
    t1<-(i-1)*njoint+1; 
    if(i < nloop) {
      t2<-i*njoint;
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
  {     
    Sigma = calc_ScoreV_RVS(X,Y,P,hom,multiplier)
  }else{
    Sigma = calc_ScoreV_likely(X,Y,hom,multiplier,snp_loop)
  }
  w = rep(1,ncol(X))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  

  
  Q = NULL
  L = NULL
  for (i in 1:nboot){
    
    if(i %% 1000 == 0)
      print(i)
    
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

test_CAST= function(S,Sigma,w){
  T = (w)%*%t(S)
  Sg = t(w)%*%Sigma%*%w
  linear = as.numeric(T)/sqrt(as.numeric(Sg))
  p = 2*(1-pnorm(abs(linear)))
  return (p)
}

#' calc_ScoreV_likely: Get Variance of vector S by likelihood method
#'
#' it is called in \code{RVS_rare}
#'@param X  the expected conditional genotype probability for case and controls
#'@param Y   the phenotype
#'@return the variance of score statistic vector S by RVS using regular expected prob
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

#' calc_ScoreV gets score vector S for nsnp number of rare variants
#'
#' it is called in \code{regScore_rare}, \code{RVS_rare}
#'@param X   the expected conditional genotype probability
#'@param Y   the phenotype
#' @return the score statistic vector S for nsnp rare variant
calc_ScoreV = function(X,Y){ # NEED to check!!!
  L = ncol(X)
  S = NULL
  for (i in 1:L){
    Yn = Y[!is.na(X[,i])]
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

#
#' Function to bootstrap sampling
#'
#' It first centralize the input matrix \eqn{X_1, X_2} for case and controls, then take a random sample from each of them with the same case or control size.
#' It calls \code{\link{calc_centralize_matrix}} and called by \code{\link{RVS_rare}}.
#'
#' @param X1   genotype likelihood matrix for case
#' @param X2   genotype likelihood matrix for control
#' @return two matrices which are the  bootstrap sampling of genotype likelihood matrix for case (Xca) and control (Xco).
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

#' Function to centralize a matrix X by column
#'
#'  This function will centeralized each column of the input X by the mean of that column calculated from non-missing data. It is called by \code{\link{sample_bootstrap}} and \code{\link{calc_ScoreV_RVS}}
#' @param X  a matrix. In our package, X is the genotype likelihood matrix, one column for each snp
#' @return A volumn-centeralized matrix for the given matrix X
#' @export
calc_centralize_matrix = function(H){
  l = ncol(H)
  mX = NULL
  for (i in 1:l){
    mm = mean(H[,i],na.rm=T)
    mv = H[,i]-mm
    mX = cbind(mX,mv)
  }
  return (mX)
}

#' Function to obtain p-value by Calpha method (quadratic test) for J rare variant
#'
#' it is called in \code{\link{regScore_rare}} and \code{\link{RVS_rare}}
#' @param S   vector of length J which includes the Scores for each variant
#' @param Sigma   the variance matrix of vector S.
#' @return p-values from Calpha test for J variants (double)
#' @export
#' @importFrom CompQuadForm davies
test_Calpha=function(S,Sigma){
  quad  =  (S)%*%t(S)
  lambda = Re(eigen(Sigma)$values)  ## the real part of the eigenvalue of Sigma will be the weight
  qp = davies(quad,lambda)   ### in CompQuadForm library, to calculate P[Q>quad] using davies's method
  return (abs(qp$Qq))
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






#' Function to remove the homozygosity of a rare SNP in case or control
#' 
#' Given the original case or control genotype matrix, return the edited matrix
#' 
#'@param X the expected condidtional genotype probablity matrix for grouped rare variants
#'@param hom with value 1 or 2. 1 means making changes with the 1st nonNA element, 2 means making changes with a random nonNA element
#'@param multiplier with value 1 or 2, 1 is dividing by 2 and 2 is mulitpling by 2.
#'@return newX the edited expected condidtional genotype probablity matrix
check_hom_matrix<-function(U,hom=1,multiplier=1)
{
  if(!hom%in%c(1,2))
  {cat('Wrong input for hom input in check_hom_matrix, hom need be 1 or 2.');
    return;
  }
  var0=apply(U, MARGIN = 2, function(x) {  sum(max(x, na.rm=TRUE)-min(x,na.rm=TRUE)==0)})
  ## column index for those all column are the same, might be more than one column
  var0.index=which(var0==1)
  #check if there are any problematic columns,
  if (sum(var0)>0) { 
    
    
    # this statement is used if there are column that the whole column has the same genotype number.
    na.index=apply(U, MARGIN = 1, function(x) { sum(is.na(x))==0})
    qq= which(na.index=='TRUE')
    if(hom==1) {
      j=qq[1]
    }else{
      j=sample(qq,1)
    }
    
    
    ind.0=which(U[j,var0.index]==0)
    ind.non0=which(U[j,var0.index]!=0)
    U[j,var0.index[ind.0]]=10^(-15)    #if the column has all zeros , replace the first row with 10^-15
    
    if(multiplier==1){
      U[j,var0.index[ind.non0]]=U[j,var0.index[ind.non0]]/2  #non-zero columns replace it by half
    }else{
      U[j,var0.index[ind.non0]]=U[j,var0.index[ind.non0]]*2  #non-zero columns replace it by 2 times
      if(any(U[j,var0.index[ind.non0]]>2)) {U[j,var0.index[ind.non0]>2]=2}
    }
  
    
  
  } 
  return(U)
}



#' calc_ScoreV_RVS calculates variance of vector S by RVS using robust variance
#'
#'  it is called in \code{RVS_rare}
#'  
#'@param X   the expected conditional genotype probability for case and controls
#'@param Y   the phenotype
#'@param P  genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#' @return the variance of score statistic vector S by  RVS using robust variance
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
  cor.X=cor.X.function(X1,hom=hom,multiplier=multiplier)
  
  Sgcase = t(V)%*%cor.X%*%V*l1
  Xm2 = calc_centralize_matrix(X2)
  Xm2[is.na(Xm2)]<-0  ## if not set NA to 0, t(Xm2)%*%Xm2 will have too many NAs
  ## I have checked, t(Xm2)%*%Xm2 is the same as cov(X2) here.
  vs  = var(Y)
  diag_S  = vs*(p*Sgcase + q*t(Xm2)%*%Xm2)
  return (diag_S)
}


cor.X.function=function(X1,hom=1,multiplier=1){
  hom=hom; multiplier=multiplier
  # rm.index=NA
  cor.X=suppressWarnings(cor(X1,use="pairwise.complete.obs")) 
  
  if (any(is.na(cor.X))) {
    
    t1=data.frame(which(is.na(cor.X),arr.ind=TRUE))
    t1=t1[t1[,1]>t1[,2],]
    n.snp.prob=nrow(t1)
    for(i in 1:n.snp.prob){
      subX=NULL
      k=t1[i,1];w=t1[i,2];
      subX=na.omit(cbind(X1[,k],X1[,w]))
      subX=check_hom_matrix(subX,hom,multiplier) 
      cor.X[k,w]=cor.X[w,k]=cor(subX[,1],subX[,2])
    }
  }
  cor.X
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
