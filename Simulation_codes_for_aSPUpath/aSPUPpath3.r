###################### Adaptive Pathway analysis based on the
###################### sum of powered score (SPU) tests and aSPU test
###################### based on permutations (and simulation may NOT work
###################### for "smalln, large p",
###################### hence unable to adjust for possible covariates);
###################### ADAPTIVELY combine the gene-specific test
###################### statistics (again!) by powered gene-specific stats!!!
###################### Version 3: modified from version 2, aSPUPpath2.r
###################### difference: un-standardizing gene-spec stat as
######################   (sum_i SNP_i^r)^(1/r),
###################### and standardizing gene-spec stat as
######################   (sum_i SNP_i^r/nSNPs)^(1/r),
######################  while in V2 the former is (sum_i SNP_i^r);
###################### as V2, V3 should avoid over-penalizing large genes 
###################### since as r->infty, Lr-norm-->max, which is NOT
######################  inversely proportional to 1/nSNPs.
######################Wei Pan, weip@biostat.umn.edu, 11/25/11 (Yes, Black Friday at home after going to gym!)

#library("mvtnorm")

# Introduction to input parameters:

# Y: disease lables; =0 for controls, =1 for cases;
# X: genotype data; each row for a subject, and each column for an SNP;
#    the value of each element is the # of the copies for an allele;

#Input:
#      Y: a vector of 0-1 response, nx1;
#      X: design matrix for genotypes;
#      B=200: number of permutations to obtain p-values;
#Goal:
#      fit model: logit(Pr(Y=1)) = a0+ a1*X ,
#      test H0: a1=0.

############################################
############################################
############################################

###extract the first fw PCs such that they explain at least cutoff*100%
### of total variation from X;
### Input: X: nObs by nSNPs;
###        cutoff: threshold to determine how many first few PCs to use;
### Output: a matrix consisting of the first few PCs;
   
extractPCs<-function(X, cutoff=0.95){

if (is.null(ncol(X)) || ncol(X)<2) Xpcs<-X
else{
X.pca<-prcomp(X, center=FALSE)
##find the min k s.t. the first k PCs explain >= cutoff variations:
Xevs<-X.pca$sdev^2
for(k in 1:ncol(X))
  if (sum(Xevs[1:k])/sum(Xevs) >= cutoff) break;
Xpcs<-X %*% X.pca$rot[,1:k]
} 

Xpcs
}

#calculate a permuttaion p-value matrix based on a permuted stat matrix:
PermPvs<-function(T0s){
B=nrow(T0s); n=ncol(T0s)
P0s<-matrix(1, nrow=B, ncol=n)
for(j in 1:n)
  for(b in 1:B)
    P0s[b,j] = sum( abs(T0s[b,j]) < abs(T0s[-b,j]) )/(B-1)
return(P0s)
}

# Main function for the pathway analysis based on the SPU tests and the adaptive SPU (aSPU) test:
##############################################
# Input: Y: a vector of 0 or 1 for controls or cases respectively;
#        X: SNP codings; each row for one subject and each column for one SNP;
#           the SNPs (with the number stored in nSNPs) from one gene are 
#           stored CONSECUTIVELY, starting from gene 1;
#        nSNPs: a vector; each element is the number of the SNPs in the 
#               corresponding gene;
#        pow: gamma values for the SPU tests; i.e. jth SPU test stat is
#             the sum of the powered score stat's with power=pow[j];
#        pow2: gamma values for gene-specific test stats;
#        B: # of permuttaions used to calculate p-values;
#        usePCs: indicating whether to extract PCs and then use PCs of X;
#        varprop: the proportion of the variations explained (cutoff) that 
#                 determines how many first PCs to use. 
# Output: for gene-unstandardized, then for gene-standardized (by #SNPs),
#          for each value of pow2, give p-values for the SPU tests with the 
#           gamma values in pow; then for each value of pow2,
#           give the p-value for the aSPU test, followed by aSPUpapth. 
#         specifically,
#####0) component unstdPsUnnorm as a vector of
#####p's for un-normed gene-spec test stats and gene-unstandardized data: 
#####   SPU(pow[1], pow2[1]), SPU(pow[2], pow2[1]),..., SPU(pow[Npow], pow2[1]),
#####   SPU(pow[1], pow2[2]), SPU(pow[2], pow2[2]),..., SPU(pow[Npow], pow2[2]),
#####   ...
#####   SPU(pow[1], pow2[Npow2]), SPU(pow[2], pow2[Npow2]),..., SPU(pow[Npow], pow2[Npow2]),
#####   aSPU(, pow2[1]), aSPU(, pow2[2]), ..., aSPU(, pow2[Npow2]), aSPUpath(,);
#####1) component unstdPs as a vector of
#####p's for gene-unstandardized data: 
#####   SPU(pow[1], pow2[1]), SPU(pow[2], pow2[1]),..., SPU(pow[Npow], pow2[1]),
#####   SPU(pow[1], pow2[2]), SPU(pow[2], pow2[2]),..., SPU(pow[Npow], pow2[2]),
#####   ...
#####   SPU(pow[1], pow2[Npow2]), SPU(pow[2], pow2[Npow2]),..., SPU(pow[Npow], pow2[Npow2]),
#####   aSPU(, pow2[1]), aSPU(, pow2[2]), ..., aSPU(, pow2[Npow2]), aSPUpath(,);
#####2) component stdPs as a vector of
#####p's for gene-standardized data: 
#####   SPU(pow[1], pow2[1]), SPU(pow[2], pow2[1]),..., SPU(pow[Npow], pow2[1]),
#####   SPU(pow[1], pow2[2]), SPU(pow[2], pow2[2]),..., SPU(pow[Npow], pow2[2]),
#####   ...
#####   SPU(pow[1], pow2[Npow2]), SPU(pow[2], pow2[Npow2]),..., SPU(pow[Npow], pow2[Npow2]),
#####   aSPU(, pow2[1]), aSPU(, pow2[2]), ..., aSPU(, pow2[Npow2]), aSPUpath(,);


aPathSPUPtest<-function(Y, X, nSNPs, pow=1:8, pow2=1, B=200, 
                        usePCs=F, varprop=0.95 ){

##get rid of the genes with 0 SNP:
nSNPs = nSNPs[nSNPs>=1]

nSNPs0<-rep(0, length(nSNPs))
if (usePCs){
  Xg<-NULL
  for(iGene in 1:length(nSNPs)){
    if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs[1:(iGene-1)])+1
    indx=(SNPstart:(SNPstart+nSNPs[iGene]-1))
    Xpcs<-extractPCs(X[, indx], cutoff=varprop)
    Xg<-cbind(Xg, Xpcs)
    if (is.null(ncol(Xpcs))) nSNPs0[iGene]=1
    else nSNPs0[iGene]=ncol(Xpcs)
    }
  } else { Xg=X; nSNPs0=nSNPs}

n<-length(Y)
k<-ncol(Xg)

#######construction of the score vector:
U<-t(Xg) %*% (Y-mean(Y))

#??????????????????????
nGenes=length(nSNPs0)

   # test stat's:
   Ts<-StdTs<-rep(0, length(pow)*nGenes)
   TsUnnorm<-Ts<-StdTs<-rep(0, length(pow)*nGenes)
   for(j in 1:length(pow))
    for(iGene in 1:nGenes){
      if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
      indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
      if (pow[j] < Inf){
        a= (sum(U[indx]^pow[j]))
        TsUnnorm[(j-1)*nGenes+iGene] = a
        Ts[(j-1)*nGenes+iGene] = sign(a)*((abs(a)) ^(1/pow[j]))
        StdTs[(j-1)*nGenes+iGene] = sign(a)*((abs(a)/nSNPs0[iGene]) ^(1/pow[j]))
        # (-1)^(1/3)=NaN!
        #Ts[(j-1)*nGenes+iGene] = (sum(U[indx]^pow[j]))^(1/pow[j])
        }
      else {
           TsUnnorm[(j-1)*nGenes+iGene] = Ts[(j-1)*nGenes+iGene] = StdTs[(j-1)*nGenes+iGene] =max(abs(U[indx]))
           }
      }

   # Permutations: 
   T0s = StdT0s = matrix(0, nrow=B, ncol=length(pow)*nGenes)
   T0sUnnorm=T0s = StdT0s = matrix(0, nrow=B, ncol=length(pow)*nGenes)
   for(b in 1:B){
     Y0 <- sample(Y, length(Y))
     #########Null score vector:
     U0<-t(Xg) %*% (Y0-mean(Y0))

     # test stat's:
     for(j in 1:length(pow))
      for(iGene in 1:nGenes){
        if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
        indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
        if (pow[j] < Inf){
          a = (sum(U0[indx]^pow[j]))
          T0sUnnorm[b, (j-1)*nGenes+iGene] = a
          T0s[b, (j-1)*nGenes+iGene] = sign(a)*((abs(a)) ^(1/pow[j]))
          StdT0s[b, (j-1)*nGenes+iGene] = sign(a)*((abs(a)/nSNPs0[iGene]) ^(1/pow[j]))
          }
        else T0sUnnorm[b, (j-1)*nGenes+iGene] = T0s[b, (j-1)*nGenes+iGene] = StdT0s[b, (j-1)*nGenes+iGene] = max(abs(U0[indx]))
        }
     }

   # obtian means/sds for each gene-level stat based on permutations:
#   mu0s<-apply(T0s, 2, mean)
#   sd0s<-apply(T0s, 2, sd)
#   sd0s<-ifelse(sd0s<1e-20, 1e-20, sd0s)
#   StdTs<-(Ts - mu0s)/sd0s
#   StdTs<-Ts/nSNPs0
#   StdT0s<-T0s
#   for(b in 1:B)
#     StdT0s[b,] = (T0s[b,])/nSNPs0
#     StdT0s[b,] = (T0s[b,] - mu0s)/sd0s

   #combine gene-level stats to obtain pathway-lelev stats:
   Ts1<-Ts2<-TsU<-rep(0, length(pow)*length(pow2))
   T0s1<-T0s2<-T0sU<-matrix(0, nrow=B, ncol=length(pow)*length(pow2))
   for(j2 in 1:length(pow2)){
     for(j in 1:length(pow)){
       TsU[(j2-1)*length(pow) +j] = sum(TsUnnorm[((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
       Ts1[(j2-1)*length(pow) +j] = sum(Ts[((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
       Ts2[(j2-1)*length(pow) +j] = sum(StdTs[((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
       for(b in 1:B){
         T0sU[b, (j2-1)*length(pow) +j] = sum(T0sUnnorm[b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
         T0s1[b, (j2-1)*length(pow) +j] = sum(T0s[b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
         T0s2[b, (j2-1)*length(pow) +j] = sum(StdT0s[b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
         }
       }  
      } 

   # permutation-based p-values:
   pPerm1 = pPerm2 = pPermU = rep(NA, length(pow)*length(pow2));
   pvs = NULL;

      for(j in 1:(length(pow)*length(pow2))) {
        pPermU[j] = sum( abs(TsU[j]) < abs(T0sU[,j]))/B
        pPerm1[j] = sum( abs(Ts1[j]) < abs(T0s1[,j]))/B
        pPerm2[j] = sum( abs(Ts2[j]) < abs(T0s2[,j]))/B
        }
      
      P0s1 = PermPvs(T0s1)
      P0s2 = PermPvs(T0s2)
      P0sU = PermPvs(T0sU)
      minP0s1 = apply(P0s1, 1, min)
      minP0s2 = apply(P0s2, 1, min)
      minP0sU = apply(P0sU, 1, min)
      minP1 =  sum( min(pPerm1) > minP0s1 )/B
      minP2 =  sum( min(pPerm2) > minP0s2 )/B
      minPU =  sum( min(pPermU) > minP0sU )/B
      minP1s<-minP2s<-minPUs<-rep(NA, length(pow2))
      for(j2 in 1:length(pow2)){
        minP0s1 = apply(P0s1[, ((j2-1)*length(pow)+1):(j2*length(pow))] , 1, min)
        minP0s2 = apply(P0s2[, ((j2-1)*length(pow)+1):(j2*length(pow))], 1, min)
        minP0sU = apply(P0sU[, ((j2-1)*length(pow)+1):(j2*length(pow))], 1, min)
        minP1s[j2] =  sum( min(pPerm1[((j2-1)*length(pow)+1):(j2*length(pow))]) > minP0s1 )/B
        minP2s[j2] =  sum( min(pPerm2[((j2-1)*length(pow)+1):(j2*length(pow))]) > minP0s2 )/B
        minPUs[j2] =  sum( min(pPermU[((j2-1)*length(pow)+1):(j2*length(pow))]) > minP0sU )/B
        }
      pvs<-list(unstdPs=c(pPerm1, minP1s, minP1), stdPs=c(pPerm2, minP2s, minP2), unstdPsUnnorm=c(pPermU, minPUs, minPU) )
 return(pvs) 
}  

 


 

