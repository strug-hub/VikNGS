#
#' Generate sequence data for one indiviudal with read depth rd 

#' it calls functions \code{seq_call} to generate sequence data with read depth = rd, each read depth has error error[i].
#' @param maf, a numeric to specify the maf for one variant
#' @param rd, a numeric, read depth of the sequence
#' @param error, a vector with length rd, error rate for each squence call
#' @return a list including read depth (rdepth), error rate (error) and a vector of length rdepth showing the sequencing data (seq)
generate_seq_1ind<-function(maf,rd,error)
{
  k=rbinom(1,2,maf)
  gen_vect = c('CC','CT','TT')
  if (k==2){genotype=c('C','C')}
  if (k==1){genotype=c('C','T')}
  if (k==0){genotype=c('T','T')}
  a=seq_call(genotype, error, ndepth=rd)
 return(list('rdepth'=rd,'error'=error,'seq'=a,'geno'=k))
}
  

#' Generate sequence data for ncase and ncontrol indiviudal at one loci

#' it calls functions \code{generate_seq_1ind} to generate sequence data for ncase and ncont at one loci, read depth for two samples (case,cont) follow 
#' normal distribution of N(mdcase,sdcase) and N(mdcont, sdcont) respectively.
#' @param ncase, an integer to specify how many case samples 
#' @param ncont, an integer to specify how many control samples
#' @param mdcase, the mean read depth for case
#' @param sddcase, the stand variance of the read depth for case
#' @param mdcont, the mean read depth for cont
#' @param sddcont, the stand deviation of the read depth for cont
#' @param mecase, the mean of error rate, probability that sequence call is wrong (double) for case>0
#' @param sdecase, standard deviation of error rate (double) for case >0
#' @param mecont,  the mean of error rate, probability that sequence call is wrong (double) for cont >0
#' @param sdecont, standard deviation of error rate (double) for cont >0
#' @param mafcase, a numeric to specify the minor allele freq for the variant in case
#' @param mafcont, a numeric to  specify the minor allele freq for the variant in cont
#' @return a list including four vectors:  concatenated read depth for each individual (vrdepth, leng(vrdepth)=ncase+ncont), error rate of each sequence for each sample (verror, leng(verror)=sum(vrdepth)), the sequencing data (seq, leng=leng(verror)) and true genotype data (vgeno,leng=leng(vrdepth)), 
#' @export
generate_seqdata_1var <- function(ncase, ncont, mdcase, sddcase, mdcont, sddcont, mecase,sdecase,mecont,sdecont, mafcase,mafcont){
    seq=NULL
    vrdepth=NULL
    verror=NULL
    vgeno=NULL
    for (j in 1:ncase){
      rd = round(sddcase*rnorm(1) + mdcase)
      if (rd <= 0){rd=1}
      error = sdecase*rnorm(rd) + mecase
      error[error<0]=abs(error[error<0]) 
      error[error>1]=error[error>1]-floor(error[error>1])
      tmp=generate_seq_1ind(mafcase,rd,error)
      seq=c(seq,tmp$seq)
      vgeno=c(vgeno,tmp$geno)
      vrdepth=c(vrdepth,rd)
      verror=c(verror,error)
    }
    
    for (j in 1:ncont){
      rd = round(sddcont*rnorm(1) + mdcont)
      if (rd <= 0){rd=1}
      error = sdecont*rnorm(rd) + mecont
      error[error<0]=abs(error[error<0]) 
      error[error>1]=error[error>1]-floor(error[error>1])
      tmp=generate_seq_1ind(mafcont,rd,error)
      seq=c(seq,tmp$seq)
      vgeno=c(vgeno,tmp$geno)
      vrdepth=c(vrdepth,rd)
      verror=c(verror,error)
    }
    return(list('seq'=seq,'vgeno'=vgeno,'vrdepth'=vrdepth,'verror'=verror))
}


#' use the simulated data to calculate the expected genotype probability and population frequency
#' 
#' use the output of  \code{generate_seqdata_1var} as input to calculate the expected genotype probability and population frequency
#' It calls subfunction \code{calc_Mr} and \code{calc_EM} for population frequency of P(G=AA,AB,BB),
#' It calls subfunction \code{calc_pobs_ndepth}, \code{calc_EG_general} and the populationfrequency for expected genotype probability
#' @param vgeno, the concatenated true genotype for each sample
#' @param seq, the concatenated sequence data
#' @param verror, the concatenated error rate for each sequence data
#' @param vrdepth, the concatenated error rate for each sample  
#' @return a list including expected genotype probability (exp_prob, leng=nsample), population fequency (pop_frq,leng=3),true genotype(geno),and read depth(vrdepth), two numerics for number of cases (ncase) and controls (ncont)
#' @export
seq_eg_freq=function(vgeno,seq,verror,vrdepth)
  {
    ## the total sequence read depth in all samples for this variant
    t_ndepth=length(verror)
    
    ##general case of calc_pobs_ndepth,get genotype likelihood P(Di|Gi) for simulation data 
    M  = calc_Mr(verror,seq,t_ndepth,vrdepth)
    p = calc_EM(M)
    p[p<0]=0
    #get likelihood matrix (t_ndepth*3) P(Di|G=TT), P(Di|G=CT), P(Di|G=CC), given a vector of error rates and a vector of sequence reads, it is concatenate from genotype likelihood  P(D=g|G=G1G2,e) for a single base read given the sequencing error rate and true genotype
    Mm = calc_pobs_ndepth(verror,seq,t_ndepth)
    EG = calc_EG_general(Mm,p,vrdepth)
  return (list('exp_prob'=EG,'pop_freq'=p,'geno'=vgeno,'vrdepth'=vrdepth,'ncase'=ncase,'ncont'=ncont))
}

#' Generate sequence data for nsnp variants
#' 
#' It calls functions \code{generate_seqdata_1var},\code{seq_eg_freq},\code{rbindlist} to generate sequence data,
#' return a list including expected genotype probability (exp_prob, leng=nsample), population fequency (pop_frq,leng=3),true genotype(geno),concatenated error rate(verror) and read depth(vrdepth)
#'
#' @param nsnp, number of variants (integer) >0
#' @param ncase, number of cases (integer) >0
#' @param ncont, number of controls (integer) >0
#' @param mdcase, average read depth in cases (double) >0
#' @param sddcase, standard deviation of read depth in cases  (double) >0
#' @param mdcont, average read depth in controls (double) >0
#' @param sddcont, standard deviation of read depth in controls (double) >0
#' @param mecase, average error rate, probability that sequence call is wrong (double) for case>0
#' @param sdecase, standard deviation of error rate (double) for case >0
#' @param mecont,  average error rate, probability that sequence call is wrong (double) for cont >0
#' @param sdecont, standard deviation of error rate (double) for cont >0
#' @param mmafcase, a vector of MAF for nsnp variants in cases, HWE is assumed >0
#' @param mmafcont, a vector of MAF for nsnp variants in controls
#' @return a list including expected genotype probability (exp_prob, leng=nsample), population fequency (pop_frq,leng=3),true genotype(geno),and read depth(vrdepth), and ncase and ncont
#' @export
seqdata_nvar<-function(nsnp,ncase, ncont, mdcase, sddcase, mdcont, sddcont, mecase,sdecase,mecont,sdecont, mmafcase,mmafcont){
  for ( i in 1:nsnp)
  {
    tmp=generate_seqdata_1var(ncase, ncont, mdcase, sddcase, mdcont, sddcont, mecase, sdecase,mecont,sdecont, mmafcase[i],mmafcont[i])
    tt=seq_eg_freq(tmp$vgeno,tmp$seq,tmp$verror,tmp$vrdepth)
    if (i==1) {
      seq=tt
    }else{ seq=rbindlist(seq,tt)}
  }
  return(seq)
}

#' rbindlist can rbind combine two lists with the same structure. 
#' 
#' @param list1, the first list,each element is a vector
#' @param list2, the second list,each element is a vector
#' @return a list,  each element is the rbinding of the corresponding element from list1 and list2
#' @export
rbindlist <- function(list1,list2)
{
  list3=mapply(rbind,list1,list2,SIMPLIFY=FALSE)
  return(list3)
}


#' clist can combine two lists with the same structure. 
#' 
#' @param list1, the first list, each element is a vector
#' @param list2, the second list, each element is a vector
#' @return a list with the same structure of list1 or list2, each element is the concatenation of the corresponding element from list1 and list2
#' @export
clist <-function(list1,list2){
   list3=mapply(c,list1,list2,SIMPLIFY=FALSE)
   return(list3)
}


