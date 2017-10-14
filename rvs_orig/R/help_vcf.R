#' Function to seprates the case and control IDs.
#'
#' This function read the header of the vcf file (rows starting with # in vcf) and split the case and controls according to the caseID_file. It returns two vectors each of which is the position index of case or control in the list of samples in vcf file (the last row starting with #). 
#' @param vcf_file The saving address and filename  of the vcf file in the format of 'address/filename'.
#' @param caseID_file The  file including case IDs in the format of 'address/filename'.
#' @param nhead The number of the lines in the vcf file which started with #. 
#' @param ncol_ID The number of columns before the sample IDs start in the last line of the headers in the vcf file, default=9.
#' @return a list which inclues two vectors that correspond to the position index of the cases or controls in the last head row of the VCF with IDs, but columns are shifted by the ncol_ID columns.
# @examples getIDs("../data/1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf","../data/bam.txt",128,9)
get_IDs <- function(vcf_file,caseID_file,nhead=128,ncol_ID=9){
  filen<-vcf_file
  filecon = file(filen, open='r')
  on.exit(close(filecon))
  
  tt2 = readLines(filecon, n=nhead) ## read the first nheads of the file
  
  S = unlist(strsplit(tt2[nhead],'\t')) #S is now a list of all the ID's including the names of columns before IDs start
  Sl = length(S)
  h1 = ncol_ID+1
  S=S[h1:Sl] #Have to get rid of the first few components which are not the IDs
  
  CaseIDs <- read.table(caseID_file) #CASE ID's
  L = length(S)
  case <- NULL
  cont <- NULL
  #IF theres an identical ID like the Exome ID in the vcf, can call it a case, otherwise call it a control
  case=which(S%in%CaseIDs$V1)
  controls=1:L
  controls=controls[-case]
  output <- list("cases"=case,"controls"=controls)
  return(output)
}


#' Function to filter variants in VCF files based on column 2 (variant location), 3 (reference call), 4 (alternative call), 7 (PASS) and missingness rate. 
#'
#' This function read the whole vcf file, filter out variants and return the genotype likelihood in VCF file for case and control separately.
#' it calls \code{\link{get_likelihood_vcf}} to convert the vcf format into phred score.
#' The keeping criteria are: variants that PASS filters, do not have duplicates, have only one reference allele with missing rate <=missing_th missing rate.
#' @param vcf_file The saving address and filename  of the vcf file in the format of 'address/filename'.
#' @param caseID_file The  file including case IDs in the format of 'address/filename'.
#' @param nvars total number of variants in the vcf file. default=1000
#' @param nread each time the read.table will read nread number of rows in file, default=300
#' @param nhead The number of the lines in the vcf file which started with #, default=128. 
#' @param ncol_ID default =9, the number of columns before each individual's data
#' @param missing_th the missing rate for filter out a variant, default is set=0.5
#' @return The function returns a list including 6 matrices and 2 vectors:
#' @return  6 matrices:case00, case01, case11, cont00, cont01 and cont11 which contain the genotype likelihoods from VCF, 2 vectors for variant chromosome and location.
#  @examples filter_VCF_SNPs('1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf',IDs$cases,IDs$controls)
#' @export
filter_VCF_SNPs<- function(vcf_file,caseID_file,nvars=1000, nread=300, nhead=128, ncol_ID=9,missing_th=0.5){
  
  tt=get_IDs(vcf_file,caseID_file,nhead,ncol_ID)
  case=tt$cases
  cont=tt$controls
  
  
  # Genotype likelihoods for cases
  # A0M -  major homozygous
  # A1M -  minor heterozygous
  # A2M -  minor homozygous
  A0M <- NULL
  A1M <- NULL
  A2M <- NULL
  #
  # Genotype likelihoods for controls
  # B0M -  major homozygous
  # B1M -  minor heterozygous
  # B2M -  minor homozygous
  B0M <- NULL
  B1M <- NULL
  B2M <- NULL
  Loc <- NULL  ## to store the location of returned variants
  Chr <- NULL
  
  
  filen = vcf_file  # ../data/1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf
  r1_wc=pipe(paste("wc -l",filen),open='r')
  true.rows<-scan(r1_wc,what=list(0,NULL),quiet=TRUE)[[1]];
  close(r1_wc)
  true.snps<-true.rows-nhead
  
  true.snps <- nvars
  
  if(true.snps<nvars) {cat('Warning: the file inluding', true.snps,' variants, less than given total variants.\n'); nvars<-true.snps}
  
  mod <- nvars %% nread
  if (mod ==0) { nloop = nvars/nread
  s <- TRUE
  }else {
    nloop <- floor(nvars/nread)+1
    s <- FALSE
  }
  
  filecon = file(filen, open='r')
  on.exit(close(filecon))
  
  loop <- 0
  while ( loop< nloop) {
    loop <- loop+1
    print(loop)
    
    # Start reading vcf file
    # when filecon has tens of thousands rows, it is hard to read them all in one time.
    if(s == TRUE | loop < nloop) {
      FF = read.table(filecon, nrows = nread, sep='\t')
    }
    if (s == FALSE & loop ==nloop) {
      FF = read.table(filecon, nrows = mod, sep='\t')
    }
    
    if ( class(FF) == 'try-error'){break}
    if ( length(FF) == 0){break}
    
    if(loop<nloop) {tt=loop*nread}
    if(loop==nloop) {tt=nvars}
    cat('while loop:',loop,', ',tt,' variants read.\n')
    
    
    chr<-as.character(FF$V1)
    loc<-as.character(FF$V2)
    FF$V4=gsub('TRUE','T',FF$V4)
    FF$V5=gsub('TRUE','T',FF$V5)
    ref<-as.character(FF$V4)
    alt<-as.character(FF$V5)
    filter<-as.character(FF$V7)
    
    
    #
    # select variants have no duplication (uniq_loc), pass the filter (pass_loc)
    # have unique reference allele or alternative allele (uniq_ref_loc and uniq_alt_loc)
    ### see the intersection afterwards
    #
    
    ## find the location of variants have no duplication
    dup=which(duplicated(loc))
    unique_loc=which(!loc%in%loc[dup])   ### location of the unique location
    
    ## find the location of variants that pass the criteria in GATK algorithm
    pass_loc=which(filter%in%'PASS')
    ## find the location of variants that only have one letter for reference allele or alternative allele
    uniq_ref_loc=which(nchar(ref)==1)
    uniq_alt_loc=which(nchar(alt)==1)
    
    
    ## find the location of variants that meet all the above criteria
    t11=intersect(unique_loc,pass_loc)
    t12=intersect(t11,uniq_ref_loc)
    t13=intersect(t12,uniq_alt_loc)
    
    nsnp_keep=length(t13)
    
    if(nsnp_keep==0) {cat('No variants left after filtering bycolumn 2, 3, 4 and 7 in VCF files in loop',loop,'\n'); break;
    }else{
      cat('Num.of.variant_kept_filtby_vcf_col2347=',nsnp_keep,'\n')
      ## variants pass the above criteria
      Fpass1<-FF[t13,]
      nfilt=nrow(FF)-nrow(Fpass1)
      loop_id=rep(loop,nfilt)
      removed=cbind(loop_id,FF[-t13,1:7])
      ## if removed any variant, write it out
      if(nrow(removed)>0) {
        if(loop == 1) {
          write.table(removed,'Variants_removed_by_dup_polymorphism_filter.txt',sep='\t',row.names=F)
        }else{
          write.table(removed,'Variants_removed_by_dup_polymorphism_filter.txt',append=T,sep='\t',row.names=F,col.names=F)
        }
      }   
      
      #### get the likelihood for these first filtered variants
      AA=get_likelihood_vcf(Fpass1,ncol_ID)  ## return AA$SNPs, AA$L_matrix$L00, AA$L_matrix$L01, AA$L_matrix$L11.
      
      L1 = AA$SNPs[,'loc']
      print(L1[1])
      chr= AA$SNPs[,'chr']
      
      A00=data.frame(AA$L_matrix$L00[case,])
      A01=data.frame(AA$L_matrix$L01[case,])
      A02=data.frame(AA$L_matrix$L11[case,])
      
      B00=data.frame(AA$L_matrix$L00[cont,])
      B01=data.frame(AA$L_matrix$L01[cont,])
      B02=data.frame(AA$L_matrix$L11[cont,])
      
      #### select variants that with missing rate smaller than the threshold given by TH for case and controls separately .
      
      TH = missing_th
      
      case1=apply(A00,2,is.na)  ## if A00 is missing for one row, A01, A02 are missing as well
      case2=apply(case1,2,sum)  ## count how many people have missing probability for each location
      mis_case=case2/nrow(A00) ## missing rate for each variant in case
      
      cont1=apply(B00,2,is.na)
      cont2=apply(cont1,2,sum)
      mis_cont=cont2/nrow(B00)
      
      col_keep= (mis_case<TH) & (mis_cont<TH)
      
      if(sum(col_keep)==0) {cat('No variants are kept after missingness filter in loop',loop,'\n');
      }else{
        if(class(A0M)%in%'NULL'){
          A0M =A00[,col_keep];A1M =A01[,col_keep];A2M =A02[,col_keep]
          B0M =B00[,col_keep];B1M =B01[,col_keep];B2M =B02[,col_keep]
        }else{
          A0M = cbind(A0M, A00[,col_keep])
          A1M = cbind(A1M, A01[,col_keep])
          A2M = cbind(A2M, A02[,col_keep])
          
          B0M = cbind(B0M, B00[,col_keep])
          B1M = cbind(B1M, B01[,col_keep])
          B2M = cbind(B2M, B02[,col_keep])
        }## end if class(A0M)%in%NULL
        Loc = c(Loc,L1[col_keep])  ## variants will be analyzed further
        Chr = c(Chr,chr[col_keep])
        if(sum(col_keep)<length(col_keep))
        {
          n_rm=length(col_keep)-sum(col_keep) 
          removed=AA$SNPs[!col_keep,]
          loop_id=rep(loop,n_rm)
          removed=cbind(loop_id,removed)
          cat('n_missing_keep=',sum(col_keep),'\n')
          
          if(loop==1) {
            write.table(removed,'Variants_removed_by_missingness.txt',row.names=F,sep='\t')
          }else{
            write.table(removed,'Variants_removed_by_missingness.txt',row.names=F,col.names=F,sep='\t',append=T)
          }## end if loop==1
        }## end if sum(col_keep)<length(col_keep)
        
        cat('total dimension of case so far: ', dim(data.frame(A0M)),'\n\n')
      }## end if sum col_keep==0
    }# end if nsnp_keep==0
  }## end while
  
  rm('A00','A01','A02','B00','B01','B02','AA','L1')
  if(length(Loc)==0) { cat('No variants are kept in the whole VCF file, check your VCF files.\n'); return(NULL);
  }else{
    snps<-data.frame(cbind(Chr,Loc))
    colnames(snps)=c('chr','loc')
    return(list("SNPs"=snps,"case00"=A0M,"case01"=A1M,"case11"=A2M,"cont00"=B0M,"cont01"=B1M,"cont11"=B2M))
  } ## end if length(Loc)==0
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




#' Function to generates the expected probabilities of the genotypes E(G_ij|D_ij).
#'
#' Using the genotype likelihood for case and controls from VCF file to generates the population frequency by calling function \code{\code{calc_EM}} and then use it to calculate the expected genotype probabilities E(G_ij|D_ij) by calling function \code{\link{calc_EG}}. Variants with homozygous call in the whole sample (the std of their E(G_ij|D_ij) <10^4) will be removed and the expected genotype probability will be reorganized so that the top rows are for case and the rest for controls.
#' @param A0M Genotype likelihood for major homozygous of the cases, output of filter_VCF_SNPs.
#' @param A1M Genotype likelihood for minor heterozygous of the cases, output of filter_VCF_SNPs.
#' @param A2M Genotype likelihood for minor homozygous of the cases, output of filter_VCF_SNPs.
#' @param B0M Genotype likelihood for major homozygous of the controls, output of filter_VCF_SNPs.
#' @param B1M Genotype likelihood for minor heterozygous of the controls, output of filter_VCF_SNPs.
#' @param B2M Genotype likelihood for minor homozygous of the controls, output of filter_VCF_SNPs.
#' @param chr chromosome of the variant, output from the  function filter_VCF_SNPs
#' @param loc location of the variant, output from the function filter_VCF_SNPs.
#' @return a list includes a matrix (exp_cond_prob) contains the expected genotype likelihood for each person (row) at each locus (column); a matrix (pop_freq) for genotype frequency in all the sample at each locus (size: number of variant *3); a data frame (SNPs) contains chr and loc of the returned variants and number of variants (nrm_hom) removed because of homozygous calls in the whole sample.
# @examples
# geneexp <- getgenexp(SNPs$A0M,SNPs$A1M,SNPs$A2M,SNPs$B0M,SNPs$B1M,SNPs$B2M,SNPs$L11)
#' @export
get_exp_geno <- function(A0M,A1M,A2M,B0M,B1M,B2M,chr,loc){
  
  A0 = A0M
  A1 = A1M
  A2 = A2M
  
  B0 = B0M
  B1 = B1M
  B2 = B2M
  loc = loc
  
  MM = NULL
  P = NULL
  
  #number of snps
  Lsnp = ncol(data.frame(A0))
  
  for (i in 1:Lsnp){
    if(Lsnp==1)
    { A=cbind(A0,A1,A2)
    B=cbind(B0,B1,B2)
    }else{
      A = cbind(A0[,i],A1[,i],A2[,i])
      B = cbind(B0[,i],B1[,i],B2[,i])
    }
    M = rbind(A,B)
    
    EG = rep(NA,nrow(M))
    kk0 = !is.na(M[,1])
    M = M[kk0,]
    ## length(M[,1]) is the number of sample, here read depth for each sample =1
    rdv = rep(1,nrow(M))
    p = calc_EM(M)
    #p = c(0.9^2,2*0.9*0.1,0.1^2)
    
    #     if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)>0)){
    #       cat('Unreasonable genotype frequency by EM algorithm.')
    #       return (NULL)
    #     }
    t1<-which(p<(-0.00001)); p[t1]=0;
    t1<-which(p>1.00001); p[t1]=1;
    P = rbind(P,p)
    
    EG[kk0] = calc_EG(M,p,rdv)
    
    #print(calc_EG(M,p,rdv))
    MM = cbind(MM,EG)
  } ## end for loop
  
  
  
  t1<-paste(chr,loc,sep=':')
  
  rownames(P)<-paste('snp',t1,sep='@')
  colnames(P)=c('P00','P01','P11')
  
  colnames(MM)=rownames(P)
  t1<-1:nrow(MM)
  rownames(MM)<-paste('sample',t1,sep='')
  
  SNP<-data.frame(cbind(as.character(chr),as.character(loc)))
  colnames(SNP)=c('chr','loc')
  
  ### to filter out homozyous variants
  nsnp<-nrow(SNP)
  ## check the homozygous variants and remove them
  hom.all_sd=rep(FALSE,nsnp)
  for(i in 1:nsnp)
  { 
    tmp4=sd(MM[,i],na.rm=T)
    if(abs(tmp4-0)<0.0001) {hom.all_sd[i]=TRUE}
  }## end for i in 1:nsnp
  
  
  homs=sum(hom.all_sd)
  if(homs>0) {
    Geno=MM[,!hom.all_sd]
    write.table(SNP[hom.all_sd,],paste0(homs,'Variants_removed_bc_homozygosity_in_whole_sample.txt'),row.names=F,quote=F,sep='\t',col.names=T)
    cat(homs,'variant(s) is(are) removed because of homozygous call in all samples (ie. monomorphic). \n')
    rm_hom=SNP[hom.all_sd,]
    SNP=SNP[!hom.all_sd,]
    P=P[!hom.all_sd,]
  }else{
    Geno=MM
  } ## end if homs>0
  if(nrow(SNP)==0){
    cat('No variants are left after removing homozygous variants (in expected genotype probability) in the whole samples.\n'); return(NULL);
  }else{
    cat('Expected genotype probabilities are calculated for', nrow(SNP),' variants \n');
  }
  return(list("SNPs"=SNP,"pop_frq"=P,"exp_cond_prob"=Geno,'nrm_hom'=sum(hom.all_sd),'ncase'=nrow(A0)))
}


#' Function to calculate the minor allele frequency (MAF) from conditional expected genotype probability.
#'
#' It calculates the MAF from conditional expected genotype probability and returns the matrices of expected genotype probability and genotype frequncy for variants indicated by maf_cut and group which are to be inputed into a test function of the rvs package.
#'
#' @param P genotype frequency for the whole sample which is the output of  function get_exp_geno.
#' @param MM expected genotype probability,  output of  function get_exp_geno.
#' @param chr the chrom info of each variant
#' @param loc  the location info of each variant
#' @param ncase integer, the number of cases
#' @param maf_cut the minor allele frequency cut-off  for common or rare variants
#' @param common Boolean variable to indicate common or rare variants
#' @return matrix of genotype frequency (P), expected genotype probability (Geno) and variant information (SNPs), number of case (ncase), number of control (ncont) and number of variants (nsnp) and a vector for phenotyp (Y). 
#' @export
# @examples
# getMAF(geneexp$P,geneexp$MM)
get_exp_MAF <- function(P,MM,chr,loc,ncase,maf_cut=0.05,common=common){
  if(!common%in%c('TRUE','true','T','t','FALSE','false','F','f')){
    cat('Wrong indicator of common, need be TRUE or FALSE only.\n');
    return(NULL)
  }
  
  
  SNPs<-data.frame(cbind(as.character(chr),as.character(loc)))
  colnames(SNPs)=c('chr','loc')
  SNPs$chr=as.character(SNPs$chr)
  SNPs$loc=as.character(SNPs$loc)
  
  
  maf = P%*%c(0,1,2)/2
  maf[maf>0.5] = 1-maf[maf>0.5]
  SNPs<-cbind(SNPs,maf)
  
  if(common %in% c('TRUE','true','T','t'))
  {  kk=which(maf>=maf_cut)
  }else{
    kk=which(maf<maf_cut)
  }
  n.all=nrow(SNPs)
  n.min=length(kk)
  cat(length(kk),'out of', n.all,'variants satisfies the MAF condition provided.\n')
  if(length(kk)==0){return(NULL)
  }else{
    if(n.all==1){ncont=length(MM)-ncase
    }else{ 
      
      
      n.all=nrow(MM)
      MM = MM[,kk]
      P = P[kk,]
      SNPs=SNPs[kk,]
      maf=maf[kk]
      
      ncont=n.all-ncase
    }## end if n.all==1
  } ## if kk==0
  Y=c(rep(1,ncase),rep(0,ncont))
  nsnp=nrow(SNPs)
  Geno=MM
  return(list("SNPs"=SNPs,"P"=P,"Geno"=MM,"ncase"=ncase,"ncont"=ncont,"nsnp"=nsnp,"Y"=Y))
}

# 
#' Function to get Likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields
#' 
#' This function fetches the genotype likelihoods from the VCF file (wo header) for the analysis.
#' The first step is to read the vcf file as a matrix and outputs an R object with 6 lists including 5 vectors and 1 list including 3 matrix:
#' 5 vectors: chr (1st column of VCF file), position (2nd column of VCF file, POS),ref allele (4th column, REF), alt allele (5th column, ALT), filter indicator (7th column, FILT) 
#' 1 list: likelihood of L(00|ref,alt), L(01|ref,alt), L(11|ref,alt), which are three matrix with row for sample and col for variants.
#' They are calculated by \code{\link{get_lsingle_vcf}} from the PL part of the 10th column and onwards of VCF in format of GT:AD:DP:GQ:PL, example: 0/0:2,0:4:6:0,6,42 
#' GL = genotype likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields in the form of log10-scale,
#' also called as Phred-scaled likelihoods: -log10(L(00|ref,alt)), -log10(L(01|ref,alt)), -log10(L(11|ref,alt)).
#' PL : the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field) (Integers)
#'
#' @param M  The vcf file without the headers in a matrix, each row stands for a variant with number of columns = ncol_ID + number of sampls.
#' @param ncol_ID the number of columns before the data for each individual
#' @return Outputs an R object with
#' 1. vector for chromosome of variants
#' 2. vector for the position of variant location
#' 3. vector for the reference allele
#' 4. ALT allele (vector)\cr
#' 5. FILTER indicator for 'PASS' (vector)\cr
#' 6. likelihoods: L(D|0), L(D|1), L(D|2). List of 3 matrix:
#' L0 - matrix of P(D|0), rows are individuals and columns are SNPS\cr
#' L1 - matrix of P(D|1), rows are individuals and columns are SNPS\cr
#' L2 - matrix of P(D|2), rows are individuals and columns are SNPS\cr
get_likelihood_vcf <- function(M, ncol_ID ){
  ref = as.vector(M[,4])
  alt = as.vector(M[,5])
  loc = as.vector(M[,2])
  filter = as.vector(M[,7])
  chr = as.vector(M[,1])
  
  L00 = NULL
  L01 = NULL
  L11 = NULL
  Lsnp = nrow(M)
  Lind = ncol(M)-ncol_ID
  
  for (i in 1:Lsnp){
    l0=NULL
    l1=NULL
    l2=NULL
    format=unlist(strsplit(as.character(M[i,9]),split=':'))
    t1=which(format%in%'PL')
    for (j in 1:Lind){
      temp<-M[i,j+ncol_ID]
      if(temp =='./.') {
        p_loci=c(NA,NA,NA)
      }else {
        p_loci=get_lsingle_vcf(temp,t1)
      }
      l0=c(l0,p_loci[1])
      l1=c(l1,p_loci[2])
      l2=c(l2,p_loci[3])
    }
    
    
    L00 = cbind(L00,l0)
    L01 = cbind(L01,l1)
    L11 = cbind(L11,l2)
  }
  SNP=cbind(chr,loc,ref,alt,filter)
  return ( list('SNPs'=SNP, "L_matrix"=list( "L00"=L00, "L01"=L01, "L11"=L11) ) )
}


#' Function to parse Phred-scaled likelihoods
#'
#' This function parses the data column (GT:AD:DP:GQ:PL) to get the PL values -10log10(x),
#' and converts them to a vector of likelihoods: L0,L1,L2. It's called within \code{\link{get_likelihood_vcf}}.
#' newer version of workflow of GATK suggested calling the BAM for each individual first then combined them, 
#' so the format column can have different length for each variant like the example Lizhen provided.
#'
#' @param Mij   a single data unit of format'GT:AD:DP:GQ:PL'. 
#' @param n_pl an integer to indicate the PL field. 
#' @return a vector of three values: l0, l1, l2
# @example get_lsingle_vcf('0/0:2,0:4:6:0,6,42')
get_lsingle_vcf <- function( Mij, n_pl ){
  a = unlist(strsplit(as.vector(Mij),':'))
  if(a[n_pl]=='.') {return(c(NA,NA,NA)) }
  l = as.numeric(unlist(strsplit(a[n_pl],',')))
  l0 = l[1]
  l1 = l[2]
  l2 = l[3]
  Pl0 = 10^(-l0/10)
  Pl1 = 10^(-l1/10)
  Pl2 = 10^(-l2/10)
  return ( c(Pl0,Pl1,Pl2) )
}


#' the main function to call the VCF file and return the  expected probabilities of the genotypes E(G_ij|D_ij) from VCF file
#'
#' @param file the VCF file
#' @param caseID_file the caseID files
#' @param nhead the row number of your VCF which contains all sample ID
#' @param ncol_ID the number of columns before the sample ID
#' @param missing_th the missing rate cut off for a SNP to be filtered, if missing>missing_th in case or controls, the SNP will be removed
#' @param nvar_tot should be the number of SNPs in your VCF file, but if VCF is too large, you can choose only read in the first nsnp SNPs
#' @param nread the number of SNPs read each time by R (nread<=nsnp)
#' @param maf_cut the MAF to be used to filter out SNP,
#' @param common Boolean variable to indicate common (TRUE) or rare (FALSE) variants 
#'
#' @return a list  for case or control (depends on group) including expected conditional probabilities, SNP info, Y, ncase/ncont etc 
#' @export
vcf_process<-function(file='example/example_1000snps.vcf',
                      caseID_file='example/caseID_old.txt',
                      nhead=128, ncol_ID=9, missing_th=0.2, nvar_tot=1000, nread=300, maf_cut=0.05, common=TRUE)
{
  ## seperate the case and controls, return their column locations
  IDs <- get_IDs(file,caseID_file,nhead,ncol_ID);
  
  nsample=length(IDs$cases)+length(IDs$controls)
  cat('This VCF includes ',nsample,' samples with ', length(IDs$cases), 'cases.\n')
  cat('cases location index in all samples:\n')
  cat(IDs$cases)
  cat('\n\n')
  
  ## filter out SNPs
  SNPs <- filter_VCF_SNPs(file, caseID_file, nvars=nvar_tot, nread=nread,nhead=nhead, ncol_ID=ncol_ID,missing_th=missing_th);
  
  if(class(SNPs)%in%'NULL') { cat ('No variants left after filtering out your VCF files\n'); return(NULL)
  }
  nsnp.left=nrow(SNPs$SNPs)
  cat(nsnp.left, 'SNPs out of ', nvar_tot, 'SNPs are kept.\n\n')
  
  ## getgenexp1 has second filter on missingness of samples, they have the same result if the missing_th2=missing in filter_SNPs.
  ##  exp_prob <- getgenexp1(SNPs$case00, SNPs$case01, SNPs$case11, SNPs$cont00, SNPs$cont01, SNPs$cont11, SNPs$chr, SNPs$Loc, missing_th2=0.5)
  exp_prob<-get_exp_geno(SNPs$case00, SNPs$case01, SNPs$case11, SNPs$cont00, SNPs$cont01, SNPs$cont11, chr=SNPs$SNPs[,'chr'], SNPs$SNPs[,'loc'])
  if(class(exp_prob)%in%'NULL') { cat ('No variants left after removing homozygous variants. \n'); return(NULL)
  }
  geno<-get_exp_MAF(exp_prob$pop_frq,exp_prob$exp_cond_prob,exp_prob$SNPs[,'chr'],exp_prob$SNPs[,'loc'],ncase=length(IDs$cases),maf_cut=maf_cut,common=common)
  
  if(!class(geno)%in%'NULL') { 
    ngeno=nrow(geno$SNPs)
    if(common==TRUE) {group='common'
    }else{ group='rare'}
    cat('There are ', ngeno, group, 'variants. \n\n')
    
    SNPs=geno$SNPs; P=geno$P; Geno=geno$Geno;ncase=geno$ncase
    ncont=geno$ncont; nsnp=geno$nsnp;Y=geno$Y;
    write.table(SNPs,paste0(nsnp,'snps_left_for',ncase,'ncase',ncont,'ncont.txt'),row.names=F,quote=F,sep='\t',col.names=T)
    return(geno)
  }
}
