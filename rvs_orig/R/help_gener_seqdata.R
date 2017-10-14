#' Function to produce a sequencing read at one locus with converage=ndepth from given a genotype and the error scores (one error for one read).
#'
#' It uses formula \eqn{p(D=g|G_1=g)=1-e,p(D=other\ call|G_1=g)=e/3}  to generate the read for a given genotype \eqn{G_1} and error rate \eqn{e}.
#' It is called from the \code{\link{generate_seqdata}} functions.
#' @param genotype  vector of two alleles for a given genotype, ie. genotype CT is given as ('C', 'T')
#' @param error  vector of sequencing errors, length of this vector is the same as read depth (ndepth).
#' @param ndepth  the read depth (integer)
#' @return a row vector of sequence reads consisting of 'A','G','T' or 'C' for a single variant with length the same as read depth (ndepth).
#' @export
seq_call <- function(genotype, error, ndepth){
  if(length(error) != ndepth) {
    cat ('dimension error in seq_call.\n')
    return(NULL)
  }
  vect_reads = rep(NA,ndepth)
  all_alleles = c('A','G','T','C')
  for (i in 1:ndepth){
    s = rbinom(1,1,0.5)+1
    g = genotype[s]
    ng = all_alleles[all_alleles!=g]
    e3 = error[i]/3
    e23 = e3*2
    e = e3*3
    r = runif(1)

    if (r <= e3) {
      vect_reads[i] <- ng[1]
    }else if (r > e3 & r <= e23) {
      vect_reads[i] <- ng[2]
    }else if (r > e23 & r <= e) {
      vect_reads[i] <- ng[3]
    }else{
      vect_reads[i] <- g
    }## end if r
  }# end for loop
  return (vect_reads)
}


#' Function to  compute the genotype likelihood of a read given the sequencing error rate \eqn{e} and true genotype \eqn{G_1G_2}: \eqn{P(D=g|G=G_1G_2,e), g, G_1, G_2 } have values of A, C, G or T.
#'
#'This function uses the formula \eqn{P(D=g|G=G_1G_2)=1/2*P(D=g|G_1)+1/2*P(D=g|G_2)} and \eqn{p(D=g|G_1=g)=1-e,p(D=other\ call|G_1=g)=e/3} 
#'to compute the genotype likelihood  \eqn{P(D=g|G=G_1G_2,e)}, more details are in  Appendix B of the paper. It is called from  functions \code{\link{calc_pobs_ndepth}} and \code{\link{calc_Mr}}.
#'
#' @param sing_read  single read, e.g. A.
#' @param A1  allele 1 in the true genotype, e.g. if \eqn{G_1G_2=AC} then A1='A'
#' @param A2   allele 2 in the true genotype, e.g. A2='C'
#' @param error  sequencing error, decimal with value between 0 and 1.
#' @return  likelihood of observing sing_read given genotype (A1A2) and error rate
calc_psingle = function(sing_read,A1,A2,error){
  if (sing_read == A1){
    p1 = 1-error
  }else{
    p1 = error/3
  } ## end if

  if (sing_read == A2){
    p2 = 1-error
  }else{
    p2 = error/3
  } ## end if
  
  p = 1/2*p1 + 1/2*p2
  return (p)
}

#' Function to calculate the likelihood matrix \eqn{P(D|G)} given a vector of error rates and a vector of sequence reads.
#'
#' This function calculates a likelihood matrix \eqn{P(D|G)} given a vector of error rates and a vector of sequence reads, where D is the sequencing data and G is the true genotype. It uses function \code{\link{calc_psingle}} to calculate the likelihood for each read given one possible genotype (AA, Aa or aa).
#'  It is called from function \code{\link{generate_seqdata}}.
#' 
#' @param Error  vector of sequencing error rates. Similarly as error in the function \code{\link{seq_call}}, but concatenate all error values over all samples.
#' @param vect_row_reads  row vector of sequence reads consisting of 'A','G','T' or 'C', similar as output from function \code{\link{seq_call}}, but concatenated over all samples in same fashion as Error.
#' @param ndepth  integer, read depth, should have value  the same as length of error or length of vect_row_reads.
#' @return a matrix with 3 columns representing P(read|AA), P(read|Aa) and P(read|aa), where row numbers equal to ndepth.
calc_pobs_ndepth= function(Error,vect_row_reads,ndepth=1){
  L = length(Error)
  L2=length(vect_row_reads)
  if(L!=ndepth || L2!=ndepth) {
    cat("Lengths mismatch between error, reads and read-depth in cal_pobs_ndepth.\n");
    return(NULL)
  }

  M = NULL
  LG  =c('TT','CT','CC')

  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      G = LG[j]
      A1 = substring(G,1,1)
      A2 = substring(G,2,2)
      m = c(m,calc_psingle(vect_row_reads[i],A1,A2,Error[i]))
    }
    M = rbind(M,m)
  }
  return (M)
}


#' Function to calcualte likelihood \eqn{P(D|G)} for simulating data, more general case of \code{\link{calc_pobs_ndepth}}
#'
#' This function calculates a likelihood matrix \eqn{P(D|G)} given a vector of error rates and a vector of sequence reads, where D is the sequencing data and G is the true genotype. It uses function \code{\link{calc_psingle}} to calculate the likelihood for each read given one possible genotype (AA, Aa or aa).
#' It is a more general case of  \code{\link{calc_pobs_ndepth}} with an additional param rdv, and called from \code{\link{generate_seqdata}}. 
#'
#' @param rdv  integer vector of read depths concatenated for all samples.
#' @param error  vector of base call error rates; Uses the same error input in the \code{seq_call} function but concatenate all error values over all samples.
#' @param vect_row_reads  row vector of sequence reads consisting of 'A','G','T' or 'C', same as output from function \code{seq_call}, but concatenated over all samples in same fashion as error.
#' @param ndepth read depth (integer), it eaquals to the length of error or vect_row_reads, or sum of vector rdv
#' @return a matrix (M) with 3 columns representing  P(read|AA), P(read|Aa) and P(read|aa), where row dimensions equal length(error)=length(vect_row_reads).
calc_Mr = function(error,vect_row_reads,ndepth=1,rdv){
  L = length(error)
  L2=length(vect_row_reads)
  if(L!=ndepth || L2!=ndepth) {cat("Lengths mismatch between error, reads and read-depth in calc_Mr."); return(NULL)}

  L = length(rdv)
  M = NULL
  LG  =c('TT','CT','CC')
  S = 0
  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      LL = 1
      for (kk in 1:rdv[i]){
        G = LG[j]
        A1 = substring(G,1,1)
        A2 = substring(G,2,2)
        LL = LL*calc_psingle(vect_row_reads[S+kk],A1,A2,error[S+kk])
      }
      m =c(m,LL)
    }
    S =S+rdv[i]
    M = rbind(M,m)
  }
  return (M)
}


#' Function to calculate the conditional expected  genotype probability \eqn{E(P(G_{ij}|D_{ij}))} given the genotype likelihoods \eqn{P(D_{ij}|G_{ij}=g)} and frequencies. Simplified version of \code{\link{calc_EG_general}}.
#'
#' Calculate the conditional expected genotype probability given the observed sequencing data \eqn{D_{ij}} from formula \eqn{E(P(G_{ij}|D_{ij}))=\sum_{g=0}^2 gP(G_{ij}=g|D_{ij})}, where \eqn{ P(G_{ij}=g|D_{ij})=P(D_{ij}|G_{ij}=g)*P(G_{ij}=g)/P(D_{ij})}.
#' \eqn{P(D_{ij}|G_{ij}=g)} (input M) is from the VCF file or function \code{\link{calc_pobs_ndepth}}, \eqn{p(G_{ij}=g)} (input p) is output from function \code{\link{calc_EM}}.
#' All values for \eqn{P(G_{ij}=g|D_{ij})} are scaled by \eqn{P(D_{ij})}. 
## It is called by \code{\link{generate_seqdata}}.
## talk to Andriy, he said in the kk loop, he did it for more general case, we only need rdv=1 for our case.
## so I renamed Andriy's get_EG to calc_EG_general, delete the kk loop in calc_EG function.
## all values are rescaled by P(D_ij), (see m below, but it is not probability, so m/sum(m) gives us sth. like prob.
#'
#' @param M   genotype likelihoods with dimension of number of sample times 3 (double), each column is the genotype likelihood \eqn{P(D_{ij}|G=AA, Aa\ or\ aa)} for one locus. It uses output from \code{\link{get_exp_geno}} or \code{\link{get_exp_MAF}} for VCF input.
#' @param p   genotype frequencies \eqn{P(G=AA, Aa\ or\ aa)} or \eqn{p(G=0,1,2\ minor\ allele)} for each SNP, it is a vector of length 3 and uses the output from  function \code{\link{calc_EM}}.
#' @param rdv   read depth for all samples. Dummy variable in \code{\link{calc_EG}}, listed here to be consitent with \code{\link{calc_EG_general}}. It is a vector of integers of one with length equal to the number of samples.
#' @return  a vector with the same length as rdv, containing conditional expectation probability \eqn{ E(P(G_{ij}|D_{ij}))}.
calc_EG <- function(M, p, rdv) {
  LL = length(rdv)

  EG = NULL
  g = c(0,1,2)
  for (i in 1:LL){
    m = NULL

    ### m=(M[i,1]*p[1],M[i,2]*p[2],M[i,3]*p[3])
    for  (j in 1:3){
      L=M[i,j]
      m = c(m,L*p[j])
    }
    pm = sum(m/sum(m)*g)
    EG = c(EG,pm)
  }
  return (EG)
}

#' More general function to calculate the conditional expected  genotype probability \eqn{E(P(G_{ij}|D_{ij}))} given the genotype likelihoods \eqn{P(D_{ij}|G_{ij}=g)} and frequencies, see \code{\link{calc_EG}} for simplified version.
#'
#'
#'  Calculate the conditional expected genotype probability given the observed sequencing data \eqn{D_{ij}} from formula \eqn{E(P(G_{ij}|D_{ij}))=\sum_{g=0}^2 gP(G_{ij}=g|D_{ij})}, where \eqn{ P(G_{ij}=g|D_{ij})=P(D_{ij}|G_{ij}=g)*P(G_{ij}=g)/P(D_{ij})}. It is called from \code{\link{generate_seqdata}}.
#'  
#' @param M   genotype likelihoods with dimension of number of sample by 3 (double), each column is the genotype likelihood \eqn{P(D_{ij}|G=AA, Aa\  or\ aa)} for one locus; uses output from \code{\link{calc_pobs_ndepth}} for simulation data.
#' @param p   genotype frequencies \eqn{P(G=AA, Aa\ or\ aa)} for each variant; output from function \code{\link{calc_EM}}.
#' @param rdv   read depth for all samples, should be a vector of integer with length equal to number of samples.
#' @return  a vector of containing conditional expected genotype probability. Length should be the same as rdv.
#'
calc_EG_general= function(M,p,rdv){
  LL = length(rdv)
  S = 0
  EG = NULL
  g = c(0,1,2)
  for (i in 1:LL){
    m = NULL
    for  (j in 1:3){
      L = 1
      for (kk in 1:rdv[i]){

        L = L*M[S + kk,j]
      }
      m = c(m,L*p[j])

    }
    S = S + rdv[i]
    pm = sum(m/sum(m)*g)
    ## added by Jiafen, simulation found with EG<0.
    if(pm<0) {pm=0}
    if(pm>2) {pm=2}
    EG = c(EG,pm)
  }
  
  return (EG)
}


#'Function to calculate the variance of conditional expected genotype probability \eqn{Var(E(P(G_{ij}|D_{ij})))}.
#'
#' This is the variance for the simplified function \code{\link{calc_EG}}. It used the formula \eqn{var(X)=E(X^2)-E(X)^2}, where \eqn{E(X)} is the same as the result of \code{\link{calc_EG}}.
## talk to Andriy, he said in the kk loop, he did it for more general case, we only need rdv=1 for our case.
#### uses output from \code{calc_pobs_ndepth} function for simulation data and output from \code{\link{get_exp_geno}} or \code{\link{get_exp_MAF}} for VCF input
#' @param M   genotype likelihoods \eqn{P(D_{ij}|G=AA, Aa\ or\ aa)} for one locus, matrix of dimension number of samples by 3 (double);
#' @param p   genotype frequencies \eqn{P(G=AA, Aa\ or\ aa)} (double); output from function \code{\link{calc_EM}}.
#' @param rdv   read depth (vector of integers) for all samples
#' @return the variance of \eqn{E(P(G_{ij}|D_{ij}))}.
calc_EG_Var = function(M,p,rdv){
  LL = length(rdv)
  EG2 = NULL
  g = c(0,1,2)
  for (i in 1:LL){
    m = NULL
    for  (j in 1:3){
      L=M[i,j]
      m = c(m,L*p[j])
    }
    pm = sum(m/sum(m)*g^2) - sum(m/sum(m)*g)^2
    EG2 = c(EG2,pm)
  }
  return (EG2)
}

#' Function to calculate the variance of conditional expected genotype probability \eqn{Var(E(P(G_{ij}|D_{ij})))}, more general case of \code{\link{calc_EG_Var}}
#'
#' This is the variance for function \code{\link{calc_EG_general}}. It used the formula \eqn{var(X)=E(X^2)-E(X)^2}, where \eqn{E(X)} is the same as the result of \code{\link{calc_EG_general}}.
#'
#' @param M   genotype likelihoods \eqn{P(D_{ij}|G=AA, Aa\ or\ aa)} for one locus, matrix of dimension number of samples by 3 (double);
#' @param p   genotype frequencies \eqn{P(G=AA, Aa\ or\ aa)} (double); output from function \code{\link{calc_EM}}.
#' @param rdv read depth (vector of integers) for all samples
#' @return the variance of \eqn{E(P(G_ij|D_ij))}
calc_EG_Var_general = function(M,p,rdv){
  L = length(rdv)
  S = 0
  EG = NULL
  g = c(0,1,2)
  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      L = 1
      for (kk in 1:rdv[i]){

        L = L*M[S + kk,j]
      }## end for kk
      m = c(m,L*p[j])
    }##end for j
    S = S + rdv[i]
    pm = sum(m/sum(m)*g^2) - sum(m/sum(m)*g)^2
    EG = c(EG,pm)
  }
  return (EG)
}

#' Function to use EM algorithm to estimate the genotype frequencies in the sample
#' 
#' Given a matrix (M) with dimension number of sample by 3 consisting of genotype likelihoods for each individual, this function uses the EM algorithm to estimate genotype frequencies.
#'  It returns a vector with three decimals to stands for probability of 0, 1 or 2 minor alleles in the whole sample
#'
#' @param M  matrix of genotype likelihoods derived from Phred-scale likelihood in VCF input for one variant.
#' @return a vector with three decimals to stands for probability of 0, 1 or 2 minor alleles.
#' @export
calc_EM <- function(M){
  p_0 = 0.15
  q_0 = 0.15
  q_n = 1
  p_n = 0
  d_n = 0
  k = 0
  while ((p_n - p_0)^2 + (q_n - q_0)^2>0.000001){
    d_0 = 1-p_0 - q_0
    v = c(p_0,q_0,d_0)
    p_D = M%*%(v)
    E_p = M[,1]*p_0/p_D
    E_q = M[,2]*q_0/p_D
    p_n= p_0
    q_n = q_0
    d_n = 1-q_0-p_0
    p_0 = sum(E_p)/length(E_p)
    q_0 = sum(E_q)/length(E_q)
    k = k+1
    if (k==1000){
 #     cat('hi','\n')
      return ( c(p_0, q_0, 1-p_0-q_0) )
    }
  }
  return ( c(p_0, q_0, 1-p_0-q_0) )
}

#' Function to calculate the probability from given odds
#' 
#' @param odds odds 
#' @return probs probability corresponding to the given odds
odds.2.probs<-function(odds)
{
  probs<-odds/(1+odds)
  return(probs)
}

#' Function to simulate the correlated binary data with given OR.
#'
#' @param N Population size, greater than ncase+ncont
#' @param preval prevalence rate of the disease
#' @param mafco MAF for controls
#' @param OR the odds ratio
#' @param ncase Number of cases in the taken sample, taken from case of the N population 
#' @param ncont Number of controls in the sample, take from the controls of N population
#' @return a dataframe include genotype data (x) and phenotype data (y) in each column
sim.corr.binary.data<-function(N,preval,mafco,OR,ncase,ncont)
{
  Ncont=N-floor(N*preval)
  
  prob.y0<-Ncont/N; prob.y1=1-prob.y0
  
  ## p(x=0,y=0)=p(y=0|x=0)*P(x=0)=p(y=0|x=0)*(P(x=0,y=0)+P(x=0,y=1))=(1/(1+exp(beta0)))*(P(x=0,y=0)+p(x=0,y=1))==>p(x=0,y=1)=exp(beta0)*P(x=0,y=0)
  ## similarly, P(x=1,y=1)=exp(beta0+beta1)*p(x=1,y=0); p(x=2,y=1)=exp(beta0+2beta1)*p(x=2,y=0)
  ### HWE
  ## p(x=0,y=0)=P(x=0|y=0)*p(y=0)=(1-mafco)^2*Ncont/N; p(x=1,y=0)=2*(1-mafco)*mafca*Ncont/N; p(x=2,y=0)=mafco^2*Ncont/N, sum(p(x=i,y=0))=p(y=0)=Ncont/N
  ## p(x=0,y=1)=(1-mafca)^2*Ncase/N; p(x=1,y=1)=2*mafca*(1-mafca)*Ncase/N; p(x=2,y=1)=mafca^2*Ncase/N
  ### so exp(beta0)*P(x=0,y=0)+exp(beta0+beta1)*p(x=1,y=0)+exp(beta0+2beta1)*p(x=2,y=0)=Ncase/N=1-Ncont/N
  ## exp(beta0)=(N-Ncont)/(p(x=0,y=0)+exp(beta1)*p(x=1,y=0)+exp(2beta1)*p(x=2,y=0))=(N-Ncont)/(Ncont*(mafca*(1-exp(beta1))-1)^2)
  
  ebeta0=(N-Ncont)/(Ncont*(mafco*(1-OR)-1)^2)
  ## ebeta0=odds.y1.x0=exp(beta0); OR=ebeta1:=exp(beta1)
  odds.y1.x0<-ebeta0;
  odds.y1.x1<-ebeta0*OR
  odds.y1.x2<- ebeta0*OR^2
  
  #  cat('ebeta0=',ebeta0,'\n')
  Ncase=N-Ncont
  
  ## p(x=0,y=0), p(x=1,y=0), p(x=0,y=1), p(x=1,y=1)
  prob.x0.y0=Ncont/N*(1-mafco)^2; prob.x1.y0=2*mafco*(1-mafco)*Ncont/N;
  prob.x0.y1=odds.y1.x0*prob.x0.y0; prob.x1.y1=odds.y1.x0*OR*prob.x1.y0
  
  Y=c(rep(1,Ncase),rep(0,Ncont))
  tmp=runif(Ncase,min=0,max=prob.y1)  ## prob(x=0,y=1)+prob(x=1,y=1)+prob(x=2,y=1)=prob(y=1)
  x1<-ifelse(tmp>prob.x0.y1,ifelse(tmp<prob.x0.y1+prob.x1.y1,1,2),0)
  tmp=runif(Ncont,min=0,max=prob.y0)  #prob(x=0,y=0)+prob(x=1,y=0)+prob(x=2,y=0)=prob(y=0)
  x0<-ifelse(tmp>prob.x0.y0,ifelse(tmp<prob.x0.y0+prob.x1.y0,1,2),0)
  x=c(x1,x0)
  case.index=sample(1:Ncase,ncase)
  cont.index=sample((Ncase+1):N,ncont)
  x=x[c(case.index,cont.index)]; y=Y[c(case.index,cont.index)]
  dat<-data.frame(x=x,y=y)
  return(dat)
}




#' Function to calculate expected genotype data for simuation data, given simulated reads, error rate and read depth.
#' 
#' @param reads  simulated reads
#' @param error_vector a vector for error rate for each simulated data of each person, one person by one person
#' @param rdepth_vector a vector for read depth for each simulated data of each person, one person by one person
#' @return a list includes expected genotype data (exp_geno) and population frequency (geno_freq)
calc_simu_exp<-function(reads,error_vector,rdepth_vector){
  t_ndepth=sum(rdepth_vector)
  M  = calc_Mr(error_vector,reads,t_ndepth,rdepth_vector)
  Mm = calc_pobs_ndepth(error_vector,reads,t_ndepth)
  p = calc_EM(M)
  #    cat('i=',i,',p=',p,',sum_p=',sum(p),'\n')
  p[p<0]=0
  EG = calc_EG_general(Mm,p,rdepth_vector)
  return(list(exp_geno=EG,geno_freq=p))
}






