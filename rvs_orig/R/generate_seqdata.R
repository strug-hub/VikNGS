# Originally include two functions, generate sequence data under null and alternative hypothesis, now only one function, generate_seqdate

#' Function to simulate data  for one variant by given parameters
#' 
#' @param N population size
#' @param preval prevalence rate for the disease
#' @param Ncont number of controls
#' @param mafco MAF for controls
#' @param OR the odds ratio
#' @param ncase the number of case
#' @param ncont the number of controls 
#' @param mdcase, average read depth in cases (double) >0#' @param sdcase, standard deviation of read depth in cases  (double) >0
#' @param mdcont, average read depth in controls (double) >0
#' @param sdcont, standard deviation of read depth in controls (double) >0#' @param me, average error rate, probability that sequence call is wrong (double) >0
#' @param sde, standard deviation of error rate (double) >0
#' @return MM - conditional expected value E(Gij|Dij)
#' @return P -estimated genotype frequencies P(G=0), P(G=1), P(G=2)
#' @return G - true genotype from which sequence data generated
#' @export 
generate_seqdata_OR1<-function(N, preval, ncase, ncont, mafco,OR,mdcase, sdcase, mdcont, sdcont, me, sde){
 Ntotal=ncase+ncont
  pheno.geno<-sim.corr.binary.data(N,preval,mafco, OR, ncase,ncont)
  v = NULL # variant in the reads concatenated person by person
  erv = NULL # error rate for each reads concatenated person by person
  rdv = NULL # vector of read depth for each indvidual (sum(rdv)=length(v)=length(erv))
  
  for(i in 1:Ntotal)
  {
    if(pheno.geno$y[i]==1){
      rd = round(sdcase*rnorm(1) + mdcase)
    }else{
      rd = round(sdcont*rnorm(1) + mdcont)
    }
    if (rd <= 0){rd=1}
    error = sde*rnorm(rd) + rep(me,rd)
    k = pheno.geno$x[i]
    gen_vect = c('TT','CT','CC')
    if (k==0){genotype=c('T','T')}
    if (k==1){genotype=c('C','T')}
    if (k==2){genotype=c('C','C')}
    
    a = seq_call(genotype, error, ndepth=rd) ## read from one individual
    v = c(v,a)
    erv = c(erv,error)
    rdv = c(rdv,rd)
  }
  
  write.table(v,paste0(ncase,'case_',ncont,'control_OR',OR,'_reads_for_simulation.txt',col.names=F,sep=' ',quote=F,append=T))
#  write.table(rdv,'read_depth_for_reads.txt',col.names=F,sep=' ',quote=F,append=T)
  # final=list(pheno=pheno.geno$y,geno<-pheno.geno$x,exp_geno=exp_geno$exp_geno,geno_freq=exp_geno$geno_freq)
  
  exp_geno<-calc_simu_exp(v,erv,rdv)
  ## reads have different length for each simulation, so rbind will have warning.
  final=list(geno=pheno.geno$x,pheno=pheno.geno$y,rdepth_vector=rdv,exp_geno=exp_geno$exp_geno,pop_freq=exp_geno$geno_freq)
  return(final)
}


#' Function to simulate data  for one variant by given parameters
#' 
#' @param N population size
#' @param preval prevalence rate for the disease
#' @param Ncont number of controls
#' @param mafco MAF for controls
#' @param OR the odds ratio
#' @param ncase the number of case
#' @param ncont the number of controls 
#' @param mdcase average read depth in cases (double) >0#' @param sdcase, standard deviation of read depth in cases  (double) >0
#' @param mdcont average read depth in controls (double) >0
#' @param sdcont standard deviation of read depth in controls (double) >0#' @param me, average error rate, probability that sequence call is wrong (double) >0
#' @param sde standard deviation of error rate (double) >0
#' @param nsnp number of variant
#' @return MM - conditional expected value E(Gij|Dij)
#' @return P -estimated genotype frequencies P(G=0), P(G=1), P(G=2)
#' @return G - true genotype from which sequence data generated
#' @export 
#' 
generate_seqdata_OR=function(N, preval, ncase, ncont, mafco,OR,mdcase,sdcase,mdcont,sdcont,me,sde,nsnp){
  seqdata.alt=foreach(i=1:nsnp,.combine=rbindlist)%dopar%{
    tmp<-generate_seqdata_OR1(N, preval, ncase, ncont, mafco,OR,mdcase, sdcase, mdcont, sdcont, me, sde)
    return(tmp)
  }
return(seqdata.alt)
}

