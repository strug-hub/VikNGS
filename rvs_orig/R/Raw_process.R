#' Function to add the phenotype in the 6th column of the genotype file recognizable by plink, eg, file with extension 'ped' or 'raw'.
#'
#' @param genofile additive coding genotype files, it is the .raw file in plink
#' @param caseIDfile  the file including the list of case IDs
#' @return a list including- geno: the genotype data with phenotype on 6th column, snps: variant information in bim format (in plink) and ncase: number of case
#' @export
geno_add_2pheno<-function(genofile,caseIDfile)
{
  t1=nchar(genofile)
  if(substr(genofile,t1-2,t1)%in%'raw') {
     tt=TRUE
  }else {
    tt=FALSE}
   geno=read.table(genofile,header=tt)
   caseID=read.table(caseIDfile,header=FALSE)
   geno[geno[,2]%in%caseID$V1,6]=1
   geno[!geno[,2]%in%caseID$V1,6]=0
   snps<-read.table(paste0(substr(genofile,1,t1-3),'bim'),header=FALSE)
   return(list("geno"=geno,"snps"=snps,"ncase"=nrow(caseID)))
}


#' Function to filter out variants in the genotype hard call file 
#' 
#' It will filtler out variants by (1) duplicated position, (2) homozygous genotype for all samples, (3) missing rate greater than given missing-cut and (4) given maf.
#'
#' @param geno  additive coding genotype files, it is the .raw file in plink with 6th columns are phenotype data (output of geno_add_2pheno if no phenotype)
#' @param snps  a matrix of the variant informations including at least chr, position, alternative and reference alleles
#' @param missing_cut  a decimal, the missing rate cut off
#' @param maf_cut  a decimal, the minor allele cut off
#' @param common  a boolean variable, common = T, return the variants with maf>=maf_cut, ow, return variants with maf<maf_cut.
#' @return a list including- geno: a matrix with the genotype data only (size of nsample*nsnp), snps: matrix with variant information, Y: a vector of phenotype corresponding to the order of the sample in geno, rm_call0: a vector including the position of variants with hom reference calls and rm_miss: a vector including the position of variants with missing_rate>missing_cut in snps after removing rm_call0.
#' @export
filter_genocall<-function(geno,snps,missing_cut=0.5,maf_cut=0.05,common=T)
{
  ## remove the duplicated position
  dup1 <- which(duplicated(snps$V4))
  if(length(dup1)>0){
  unique_loc=which(!snps$V4%in%snps$V4[dup1])   ### location of the unique location
  snps=snps[unique_loc,]
  geno1=geno[,c(1:6,6+unique_loc)]
  }else{
  geno1=geno;
  }
  
  cat(length(dup1),'duplciate location.\n')
  ### the homozygous calls, the only genotype is 0 in the whole sample
  
  miscall<-NA
  for (i in 7:ncol(geno1))
  {
    if(length(levels(as.factor(geno1[,i])))==1 && levels(as.factor(geno1[,i]))=="0") {
      miscall<-c(miscall,i)
  }## end if length(levels geno1[,i])
} ## end for
    if(length(miscall)>1)
    {
    miscall<-miscall[!is.na(miscall)]
    geno1=geno1[,-c(1:6,miscall)]
    snp1=snps[-(miscall-6),]
    }else{
    geno1=geno1[,-c(1:6)]
    snp1=snps
  }## enf if length miscall>1
  
#  cat(length(miscall),'variants with homozygous reference calls in all samples.\n')
  
  snp1$V5<-as.character(snp1$V5) #class(snp1$V5)<-'character' will convert the factor into numbers
  snp1$V6<-as.character(snp1$V6)
  indel<-which(nchar(snp1$V5)>1|nchar(snp1$V6)>1)
  if(length(indel)>0){
    snp1.1=snp1[-indel,]
    geno1.1=geno1[,-indel]
  }else{
    snp1.1=snp1
    geno1.1=geno1
  }
  
  cat(length(indel),'indels to be removed.\n')
  
  ## missing_rate
  X1=geno1.1[geno[,6]%in%'1',]
  X2=geno1.1[geno[,6]%in%'0',]
  filt_miss=filter_by_miss(X1,X2,missing_cut=missing_cut)
  misfilter=which(filt_miss$k%in%'FALSE')
  if(length(misfilter)>0){
  geno2=geno1.1[,filt_miss$k]
  snp2=snp1.1[-misfilter,]
  }else{
  geno2=geno1.1
  snp2=snp1.1
  } ## end if length misfilter>0
#  cat(length(misfilter),'variants filtered by high missing rate.\n')
  ## maf
   maf=get_allele_frq(filt_miss$M)
   snp2=cbind(snp2,maf)
   if (common==T)
   { n_snp=which(maf>=maf_cut)
     if(length(n_snp)==0){cat('No SNP with maf greater than given maf_cut',maf_cut,'\n'); return(NULL);
      }else{
        geno.final=geno2[,maf>=maf_cut]
        snp.final=snp2[maf>=maf_cut,]
      }## if length(n_snp)==0)
   }else{
     n_snp1=which(maf<=maf_cut)
     if(length(n_snp1)==0){cat('No SNP with maf smaller than given maf_cut',maf_cut,'\n'); return(NULL);
     }else{
     geno.final=geno2[,maf<maf_cut]
     snp.final=snp2[maf<maf_cut,]
     } ## if length(n_snp1)==0)
   }## end if common==T
   Y=geno[,6]
   geno=geno.final
   snps=snp.final
   rm_call0=miscall-6
   rm_miss=misfilter
 #  save(Y,geno,snps,rm_call0,rm_miss,file='genotype_call_filtered.RData')
  return(list("Y"=Y,"geno"=geno.final,"snps"=snp.final,"rm_call0"=rm_call0,"rm_miss"=misfilter))
}

#' Function to deal with genotype hard call files in additive coding
#' 
#' It will read the additive coding genotype file, filter out variants and then return the phenotype and genotype file necessary for association test.
#' It calls functions \code{\link{geno_add_2pheno}} and \code{\link{filter_genocall}}
#'
#' @param genofile  additive coding genotype files, it is the .raw file in plink
#' @param caseIDfile  the file including the list of case IDs
#' @param missing_cut  a decimal, the missing rate cut off
#' @param maf_cut  a decimal, the minor allele cut off
#' @param common  a boolean variable, common = T, return the variants with maf>=maf_cut, ow, return variants with maf<maf_cut.
#' @return a list including- geno: a matrix with the genotype data only (size of nsample*nsnp), snps: matrix with variant information, Y: a vector of phenotype corresponding to the order of the sample in geno, rm_call0: a vector including the position of variants with hom reference calls and rm_miss: a vector including the position of variants with missing_rate>missing_cut in snps after removing rm_call0.
#' @export
geno_process<-function(genofile,caseIDfile,missing_cut,maf_cut,common){
  t1=geno_add_2pheno(genofile,caseIDfile)
  cat('There are',nrow(t1$snps),'variants in the inputed genotype files.\n')
  t2=filter_genocall(t1$geno,t1$snps,missing_cut=missing_cut,maf_cut=maf_cut,common=common)
  cat(length(t2$rm_call0),'monomorphic variants removed.\n')
  
  cat(length(t2$rm_miss),'remaining variants are removed because missing rate higher than',missing_cut,'\n')
  if(common==T) {
      cat(nrow(t2$snps),'remaining variants have MAF greater than',maf_cut,'\n')
  }else{
    cat(nrow(t2$snps),'remaining variants have MAF less than',maf_cut,'\n')
  }
  return(list("Y"=t2$Y,"geno"=t2$geno,"snps"=t2$snp,'rm_call0'=t2$rm_call0,"rm_miss"=t2$rm_miss))
}

 

#' Function to take a subset the samples from a larger genotype data
#' 
#' This function will take a subset of samples' genotype data from a larger dataset, only return the genotype for the wanted individuals.
#' Similar as the option --keep in plink.
#' 
#' @param path1  path and file name of the larger genotype files where subset is taken from, ped file is prefered
#' @param path2  the file including the list of IDs for individuals to be kept
#' @return a list including the genotype data matrix (geno1) and a interger for the number of samples (n1).
subset_ind<-function(path1,path2)
{
  raw<-read.table(path1,header=F,stringsAsFactors=F)
  keep<-read.table(path2,header=F,stringsAsFactors=F)

  geno1=raw[raw[,2]%in%keep$V1,]
  n1<-nrow(geno1)

  return( list('geno1'=geno1,'n1'=n1))
}


#' Function to find the variants presenting in two datasets.
#' 
#' This function try to find the variants presenting in two datasets, it will return their location index. Variants presenting in both dataset but with duplicated record will not be returned.
#' @param cord_set1  the coordinate of the variants in file1
#' @param cord_set2  the coordinate of the variants in file2
#' @return a list include two vectors for the location of variants presenting in both files in file1 (cord1) and  file2 (cord2).
#' @export
common_coord = function(cord_set1,cord_set2){
  set1 = NULL
  set2 = NULL
  L = length(cord_set1)
  for (i in 1:L){
    k = which(cord_set2 == cord_set1[i])
    if (length(k)==1){
      set1 = c(set1,i)
      set2 = c(set2,k)
    }
  }
  com_cord = list(cord1 = set1,cord2 = set2)
  return (com_cord)
}

#' Function to extract the genotype columns for chosen variants from a ped format file.
#'
#' Given the full set of genotype file in plink format and variants index would like to keep, this function will return the genotype columns for the wanted variants (the first 6 columns are removed). Similar to the option --extract in plink software, but the first 6 columns are removed in this function.
#' @param M the genotype file in ped format,Note, each snp have two columns in the ped file and the first 6 column need to be included.
#' @param set the return of common_coord, the vector of integers representing the location of the common snps in the two sets.
#' @return the genotype file for the snps common in two sets only, each snp has two columns.
#' @export
subset_snp = function(M, set){ ## Note, each snp have two columns in the ped file

  set_k = 2*set + 5 ## 2*(set-1)+7, the start column Num of the genotype for the set-th snp
  set_k1 = set_k + 1  ## 2*set+6, the end column Num of the genotype for the set-th snp

  L = length(set_k)
  all_colns = NULL #1:6
  for (i in 1:L){
    all_colns  = c(all_colns,set_k[i],set_k1[i])
  }
  Msb = M[,all_colns]
  return (Msb)
}


#' Function to change all TRUE genotype into T 
#' 
#' In R, if all the individual has genotype 'T' in one strand, R automatically change them into 'TRUE'. This function will change all 'TRUE' genotye back into T.
#'
#' @param M the genotype matrix
#' @return M the genotype matrix after change TRUE into T
#' @export
snp_true_2T <-function(M)
{
  kk=which(M[1,] == 'TRUE')
  for(i in kk){ M[,i]=factor('T')}
  return(M)
}

#'  Function to combine two ped files for case and controls and transform them into additive coding
#'
#'   This function remove variants that are missing in one sample, monomorphic in whole sample, monomorphic in each sample with diff alleles, non biallelic
#'  Alternatively, for those know plink, all these can be done by plink to merge the two files and generate the additive coding with option --recodeA
#'
#' @param M1 the genotype data for the first sample (normally case)
#' @param M2 the genotype data for the secondsample (normally controls)
#' @param ncols  the number of column before the genotype data start
#' @return A list including - Ad_all: the additive coding for all snps left, ncase: the first ncase rows are for cases;
#'  minor: the minor allele for each variants; rmed, the index of the removed snps from the input
#' @export
ped_2_additive = function(M1,M2,ncols){
  L = (length(M1[1,])-ncols)/2
  A1 = NULL
  A2 = NULL
  Minor=NULL
  rmed=NULL
  for (i in 1:L){
    c1=2*i-1+ncols
    c2=2*i+ncols
#    c1=1
#    c2=2
    samp1=M1[,c(c1,c2)]
    samp2=M2[,c(c1,c2)]
    l1 = levels(unlist(c(list(M1[,c1],M1[,c2]))))
    l2 = levels(unlist(c(list(M2[,c1],M2[,c2]))))

    set1=l1[l1 != '0']
    set2=l2[l2 != '0']
    d1 = length(set1)
    d2 = length(set2)

    if(d1 ==0 && d2 > 0)
    { cat('snp',i,' is missing in the first sample; removed. \n')
      rmed=c(rmed,i)
    }else if(d1>0 && d2 == 0){
      cat('snp',i,' is missing in the second sample; removed. \n')
      rmed=c(rmed,i)
    }else if(d1 == 1 && d2 == 1 && set1 == set2){
      cat('snp',i,' is monomorphic in the whole samples; removed. \n')
      rmed=c(rmed,i)
    }else if(d1 == 1 && d2 == 1 && set1 != set2){
      cat('snp',i,' is monomorphic in each sample, and different samples with different calls; removed. \n')
      rmed=c(rmed,i)
    }else if (d1 + d2>=3 && d1 + d2 <=4 && length(unique(c(set1,set2))) > 2)
    {cat('snp',i,' is biallelic in each sample, overall there are more than 2 alleles; removed. \n')
      rmed=c(rmed,i)
    }else if ( d1>2 || d2>2){cat('snp ',i,' has more than 2 alleles; removed. \n')
        rmed=c(rmed,i)
    }else{
     ma=minor_allele(M1[,c(c1,c2)],M2[,c(c1,c2)])
     A1=cbind(A1,convert(samp1,ma))
     A2=cbind(A2,convert(samp2,ma))
     Minor=c(as.character(Minor),as.character(ma))
    }
  }
    Mall=rbind(A1,A2)
    nsamp1=nrow(A1)
     addi = list("Ad_all"=Mall,"nsamp1"=nsamp1,Minor=Minor,rmed=rmed)
     return (addi)
}

#' Function to find the minor allele in the whole sample
#' 
#' @param M1  the 2 column genotype for a variant from sample1's ped file
#' @param M2  the 2 column genotype for a variant from sample2's ped file
#' @return  the minor allele from the two inputed samples (character).
minor_allele <- function(M1, M2){
  c1 =unlist(c(list(M1[,1],M1[,2],M2[,1],M2[,2])))
  c1=c1[!c1%in%'0']
  t1=data.frame(table(c1))
  tt=which(t1[,1]%in%'0')
  if(length(tt>0)) {t1=t1[-tt,] }
  t1=t1[order(t1$Freq),]
  minor=(t1[1,1])
  return(minor)
}


#' convert two columns of a variant from ped file into additive coding
#'
#' @param M  2 columns of a  variant from ped file
#' @param ref the first allele for the  variant in map file 
#' @return a vector representing additive coding with length the same as sample size for that variant.
convert = function(M,ref){
  L = length(M[,1])
  Ad = rep(0,L)
  for (i in 1:L){
    a  = M[i,]
    k = which(a=='0')
    if (length(k)>0){
      Ad[i] = NA
    }else{
      k = which(a==as.character(ref))
      Ad[i] = length(k)
    }
  }
  return (Ad)
}


#' Function to filter out variants by missing rate p
#'
#' It removes all variants that missing rate are smaller than p in either sample.
#' @param M1 the additive coding of genotype for sample1
#' @param M2 the additive coding of genotype for sample2
#' @param missing_cut the missing rate cut-off, only keep variants with missing rate <p in both samples
#' @return a list including row combined genotype from both sample (M)  and a vector of index whether a variant pass the missingness filter (k).
#' @export
filter_by_miss = function(M1,M2,missing_cut){
  L = length(M1[1,])
  J1 = length(M1[,1])
  J2 = length(M2[,1])
  p=missing_cut
  k = rep(FALSE,L)
  for (i in 1:L){
    ll1 = sum(is.na(M1[,i]))
    ll2 = sum(is.na(M2[,i]))
    if ((ll1<p*J1) & (ll2<p*J2)){
        k[i] = TRUE
    }
  }
  tt=which(k>0)
#  cat('filter_keep=',tt)
  M = rbind(M1[,k],M2[,k])
  return (list(M=M, k =k))
}

#' Function to calculate the allele frequency from additive coding genotype
#'
#'Function to calculate the allele frequency from additive coding genotype, the number of allele is not necessarily for minor allele. For people know plink, there is another option to calcualte the MAF.
#'
#' @param M the additive coding genotype, a matrix with column for each variant
#' @return a vector of allele frequency calculated based on the non-missing data
#' @export
get_allele_frq = function(M){
  L = length(M[1,])
  maf = rep(0,L)
  for (i in 1:L){
    LL = sum(!is.na(M[,i]))
    X = sum(M[!is.na(M[,i]),i])
    maf[i] = X/(2*LL)
  }
  return (maf)
}

#' Function to reset the additive coding genotype into nubmer of minor alleles.
#'
#' This is not needed if your genotype coding is based on the number of minor allele. It calls \code{\link{get_allele_frq}} for the allele frequency then reset the genotype matrix if frequency >0.5.
#' @param M the additive coding genotype matrix for the whole sample
#' @return a matrix for additive genotype in terms of number of minor allels. 
reset_minor_additive= function(M){
  L = length(M[1,])
  maf = get_allele_frq(M)
  for (i in 1:L){
    if (maf[i]>0.5){
      a = abs(M[,i] - 2)
      M[,i] = a
    }
  }
  return (M)
}


#' Function to combine two sets of ped pairs for two samples
#'
#' @param file1 the file name for the first sample before the extension. Eg. file1=sample1 if the file names are sample1.ped/map
#' @param file2  the file name for the second sample before the extension. Eg. file1=sample2 if the file names are sample2.ped/map
#' @param keep1 if only a subset of sample1 is needed in the final analyis, keep1=1, ow, keep1=NA
#' @param keep2 if only a subset of sample2 is needed in the final analyis, keep2=1, ow, keep2=NA
#' @param missing_cut the missing rate to remove the high-missing variants
#' @return a list include data sets 'common' and 'rare' for common  or rare variants and 'snp.Nomiss' for the SNPs kept in common and rare.
#' @export
combine_twogeno<-function(file1,file2,keep1=NA,keep2=NA,missing_cut=0.5)
{
  ## find the common SNPs in two datasets
  map1=read.table(paste0(file1,'.map'))
  map2=read.table(paste0(file2,'.map'))
  snpsub=common_coord(map1$V4,map2$V4)
  
  #### ped file to take the common variants in both datasets
  ped1=read.table(paste0(file1,'.ped'))
  ped2=read.table(paste0(file2,'.ped'))
  
  ped1.com=subset_snp(ped1,snpsub$cord1) ## only the genotype data are kept.
  ped2.com=subset_snp(ped2,snpsub$cord2)
  ## convert the subsetted ped file into additive coding
  ped1.com=snp_true_2T(ped1.com) ## if one column contains allele 'T', R think it is TRUE, need to change it into character T.
  ped2.com=snp_true_2T(ped2.com)
  
  pedall<-ped_2_additive(ped1.com,ped2.com,0) ## two columns geno into 1 column additive coding.
  snp.all<-map1[snpsub$cord1,][-pedall$rmed,]
  
#  save(snp.all,file='genotype_allsample_ped2additive')
#  cat('keep1=',keep1,'keep2=',keep2,'\n')
  if(!is.na(keep1))
  {   tmp1<-read.table(paste0(file1,'_keep.txt'))
  index1<-which(ped1$V2%in%tmp1$V2)
  nkeep1 <- length(index1)
  }
  if(!is.na(keep2))
  {   tmp2<-read.table(paste0(file2,'_keep.txt'))
      index2<-which(ped2$V2%in%tmp2$V2)
  }
  if(exists("index1") && exists("index2")) {
    ped.final <- pedall$Ad_all[c(index1,index2+pedall$nsamp1),]
  }else if (exists("index1") && (!exists("index2"))){
    rm=(1:pedall$nsamp1)[-index1]
    ped.final <- pedall$Ad_all[-rm,]
  }else if (!exists("index1") && (exists("index2")) ){
    ped.final <- pedall$Ad_all[c(1:pedall$nsamp1,index2+pedall$nsamp1),]
  }else{
    ped.final <- pedall$Ad_all
  }
  if(!exists("nkeep1")) nkeep1<-pedall$nsamp1

  nall<-nrow(ped.final)
 # cat('ncase=',ncase,'nall=',nall,'\n')
  geno2test=filter_by_miss(ped.final[1:nkeep1,],ped.final[(nkeep1+1):nall,],missing_cut=missing_cut)
  snp.Nomiss=snp.all[geno2test$k,]
  Y=c(rep(1,nkeep1),rep(0,nall-nkeep1))
  
  maf=get_allele_frq(geno2test$M)
  common=geno2test$M[,maf>=0.05]
  rare=geno2test$M[,maf<0.05]
  snp.Nomiss=cbind(snp.Nomiss,maf)
  final=list('common'=common,'rare'=rare,'snp.Nomiss'=snp.Nomiss,"Y"=Y)
  cat('\nThere are', ncol(common), 'common SNPs and ',ncol(rare),' rare SNPs.\n')
  cat('All data are saved in list final, written to  data genotype_data.RData.\n') 
  #save(final,file='genotype_data.RData')
  return(final)
}  



