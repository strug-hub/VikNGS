


generate_seqdata_OR1_zeynep<-function(pop.data,ncase,ncont,rd.group,nhrd1,nhrd2,nhrd3,nlrd1,nlrd2,mdhrd1,sdhrd1,mdhrd2,sdhrd2,mdhrd3,sdhrd3,mdlrd1,sdlrd1,mdlrd2,sdlrd2,me,sde){
 Ntotal=ncase+ncont
  pheno.geno<-sim.corr.binary.data(N,preval,mafco, OR, ncase,ncont)
  v = NULL # variant in the reads concatenated person by person
  erv = NULL # error rate for each reads concatenated person by person
  rdv = NULL # vector of read depth for each indvidual (sum(rdv)=length(v)=length(erv))
  
  for(i in 1:Ntotal)
  {
    if(rd.group[i]==0){
      rd = round(sdhrd1*rnorm(1) + mdhrd1)
    }else if(rd.group[i]==1) { 
      rd = round(sdhrd2*rnorm(1) + mdhrd2)
    }else if(rd.group[i]==2) {	
      rd = round(sdhrd3*rnorm(1) + mdhrd3)	
    }else if(rd.group[i]==3) {	
      rd = round(sdlrd1*rnorm(1) + mdlrd1)	
    }else {rd = round(sdlrd2*rnorm(1) + mdlrd2)	}
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
  
#  write.table(v,paste0(ncase,'case_',ncont,'control_OR',OR,'_reads_for_simulation.txt',col.names=F,sep=' ',quote=F,append=T))
#  write.table(rdv,'read_depth_for_reads.txt',col.names=F,sep=' ',quote=F,append=T)
  # final=list(pheno=pheno.geno$y,geno<-pheno.geno$x,exp_geno=exp_geno$exp_geno,geno_freq=exp_geno$geno_freq)
  
  exp_geno<-calc_simu_exp(v,erv,rdv)
  ## reads have different length for each simulation, so rbind will have warning.
  final=list(geno=pheno.geno$x,pheno=pheno.geno$y,rdepth_vector=rdv,exp_geno=exp_geno$exp_geno,pop_freq=exp_geno$geno_freq)
  return(final)
}




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

