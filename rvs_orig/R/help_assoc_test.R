## In this file, we include all the helper function for association test
## calc_score_test: the regular score test, including the regular variance of score statistic
## calc_robust_Var: the robust variance of score statistic (using prob(G_ij) for case)
## calc_scoreV: the vector of score statistic for J rare variants
## calc_scoreV_regVar(X,Y): the regular variance of score vector.
## calc_scoreV_RVS(X1,X2,Y): the robust variance of score vector (using p(G_ij) for case)
## calc_scoreV_RVS_regVar(X1,X2,Y): the robust variance of score vector using var(X(Y=1)) and var(X(Y=0))
## test_CAST, test_Calpha and test_Hotelling, the three test statistic for rare variants with input score S and variance Sigma
## calc_centralize_matrix, sample_boostrap used in the rare variant associaiton test.

# These packages should be downloaded from R-CRAN
library('MASS')
library('CompQuadForm')

#' Function to calculate the score test for given expected genotype probability for case and controls seperately (M1, M2)
#'
#' This function calcuates the score test \eqn{T=S^2/var(S)}, for a given \eqn{j, S_j=\sum_i (Y_i-bar(Y))E(G_{ij}|D_{ij})=\sum_{case}(1-bar(Y))E(G|D)-\sum_{cont}bar(Y)E(G|D)}. It is called in \code{\link{regScore_perm}}.  
#' @param M1   a vector of the expected  genotype probability for case \eqn{E(G_{ij}|D_{ij})} on one snp 
#' @param M2   a vector of the expected genotype probability for control \eqn{E(G_{ij}|D_{ij})} on one snp
#' @return the score test statistic calculated from M1 and M2
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



#' clac_score_test1 calculates the score test for given expected genotype probability for cases and controls seperately (M1, M2)
#'
#' calc_score_test1 is called in \code{\link{regScore_perm}}
#' @param M1 -  the expected value of genotype for case E(G_ij|D_ij) on one snp (dimension: ncase*1)
#' @param M2 -  the expected value of genotype for control E(G_ij|D_ij) on one snp (dimension: ncont*1)
#' @return  - the score test statistic calculated from M1 and M2 and empirical variance estimation seperated for case and controls.
calc_score_test1 = function(M1,M2){
  X = c(M1,M2)
  ncase <- length(M1)
  ncont <- length(M2)
  Y = c(rep(1,ncase),rep(0,ncont))

  p = ncase/length(X)
  q = 1 - p
  S = q*sum(M1)-p*sum(M2)
  vs = p*ncont*(q*var(X[Y==1])+p*var(X[Y==0]))
  x = S^2/vs
  return (x)
}


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



#' calc_ScoreV_regVar: Get Variance of vector S by regular method
#'
#'  it is called in \code{regScore_rare}
#' @param X -  the expected genotype probability
#' @param Y -  the phenotype
#' @return the variance of score statistic vector S by regular method
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

#' Function to remove the homozygosity of a rare SNP in case or control
#' 
#' Given the original case or control genotype matrix, return the edited matrix
#' 
#'@param X the expected condidtional genotype probablity matrix for grouped rare variants
#'@param hom with value 1 or 2. 1 means making changes with the 1st nonNA element, 2 means making changes with a random nonNA element
#'@param multiplier with value 1 or 2, 1 is dividing by 2 and 2 is mulitpling by 2.
#'@return newX the edited expected condidtional genotype probablity matrix
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


#' Function to edit the correlation matrix of expected genotype probability. It either remove the NA correlation pairs or remove the same pairs (complete LDs).
#' 
#' Given the original case or control genotype matrix, return the edited correlation matrix. In the check_hom_matrix, we are trying to solve the problem that the whole column is homozygous.
#' But this NA correlation pairs happens when one column is not homozygous, but when paired with another column, the non-missing values might be homozygous. And it will return the edited correlation matrix
#' while check_hom_matrix return the edited genotype probability matrix.
#' 
#'@param X the expected condidtional genotype probablity matrix for grouped rare variants
#'@param hom with value 1 or 2. 1 means making changes with the 1st nonNA element, 2 means making changes with a random nonNA element
#'@param multiplier with value 1 or 2, 1 is dividing by 2 and 2 is mulitpling by 2.
#'@param snp_loop will tell the relative location of the probablematic snp.
#'@return cor.X the edited correlation from the expected condidtional genotype probablity matrix

cor.X.function=function(X1,hom=1,multiplier=1,snp_loop){
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
#   index=which(cor.X==1,arr.ind=TRUE)
#   index=data.frame(index[index[,1]>index[,2],])
#   while(nrow(index)>0){
#       k=index[1,1];w=index[1,2];
#       if(sum(is.na(X1[,k]))>sum(is.na(X1[,w]))) {X1=X1[,-k]
#       }else{X1=X1[,-w]
#       }
#       cat('rare variants group',snp_loop,' has complete LD variants, so one of the pair is removed.\n')
#       cor.X=suppressWarnings(cor(X1,use="pairwise.complete.obs")) 
#       index=which(cor.X==1,arr.ind=TRUE)
#       index=data.frame(index[index[,1]>index[,2],])
#   }
#   return(list(cor.X=cor.X,rm.index=rm.index))
  return(cor.X)
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
  cor.X=cor.X.function(X1,hom=hom,multiplier=multiplier,snp_loop)
  
  Sgcase = t(V)%*%cor.X%*%V*l1
  Xm2 = calc_centralize_matrix(X2)
  Xm2[is.na(Xm2)]<-0  ## if not set NA to 0, t(Xm2)%*%Xm2 will have too many NAs
  ## I have checked, t(Xm2)%*%Xm2 is the same as cov(X2) here.
  vs  = var(Y)
  diag_S  = vs*(p*Sgcase + q*t(Xm2)%*%Xm2)
  return (diag_S)
}


#
#' Function to obtain p-value by CAST method (linear test) for J rare variant
#'
#' it is called in \code{regScore_rare} and \code{RVS_rare}
#' @param S   vector of length J which includes the Scores for each variant
#' @param Sigma   the variance matrix of vector S.
#' @param w   the weight for CAST test
#' @return p-values from CAST test for J variants (double)
#' @export
#' @importFrom stats pnorm
test_CAST= function(S,Sigma,w){
  T = (w)%*%t(S)
  Sg = t(w)%*%Sigma%*%w
  linear = as.numeric(T)/sqrt(as.numeric(Sg))
  p = 2*(1-pnorm(abs(linear)))
  return (p)
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

#' Function to obtain pvalue by Hotelling test for J number of rare variants
#' 
#' @param S   vector of length J which includes the Scores for each variant
#' @param Sigma the variance matrix of vector S.
#' @return p-values from Hotelling test for J variants (double)
#' @import MASS
test_Hotelling = function(S,Sigma){
  dd  =  try(ginv(Sigma), silent = TRUE)
  if ( class(dd) == 'try-error'){
    cat('Inverse_error','\n')
    Sigma = diag(diag(Sigma))
  }
  X =ginv(Sigma)
  quad = S%*%X%*%t(S)
  rank = qr(X)$rank
  p = 1-pchisq(quad,rank)
  return (p)
}

#' Function to centralize a matrix X by column
#'
#'  This function will centeralized each column of the input X by the mean of that column calculated from non-missing data. It is called by \code{\link{sample_bootstrap}} and \code{\link{calc_ScoreV_RVS}}
#' @param X  a matrix. In our package, X is the genotype likelihood matrix, one column for each snp
#' @return A volumn-centeralized matrix for the given matrix X
#' @export
calc_centralize_matrix = function(X){
  l = ncol(X)
  mX = NULL
  for (i in 1:l){
    mm = mean(X[,i],na.rm=T)
    mv = X[,i]-mm
    mX = cbind(mX,mv)
  }
  return (mX)
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

