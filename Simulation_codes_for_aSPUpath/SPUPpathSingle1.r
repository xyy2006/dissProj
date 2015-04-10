###################### Pathway analysis based on the
###################### sum of powered score (SPU) tests and aSPU test
###################### based on permuttaions (and simulation may NOT work
###################### for "smalln, large p",
###################### hence unable to adjust for possible covariates);
###################### BUR test each gene SEPARATELY (and combine their test
###################### statistics by Max!!!
###################### modified from aSPU/prog/wSPUP1.r
######################Wei Pan, weip@biostat.umn.edu, 10/2/11

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
#Output: p-values of the wSPU tests in the order of supplied pow values;
#        finally, give the p-value of the awSPU test (that combines the
#           SPUs tests with pow by taking their min P-value and adjust for
#           multiple testing.

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
#        B: # of permuttaions used to calculate p-values;
#        usePCs: indicating whether to extract PCs and then use PCs of X;
#        varprop: the proportion of the variations explained (cutoff) that 
#                 determines how many first PCs to use. 
# Output: p-values for the SPU tests with the gamma values in pow, then
#         the p-value for the aSPU test.

PathSPUPtestSingle<-function(Y, X, nSNPs, pow=1:8, B=200, usePCs=F, varprop=0.95 ){

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
   Ts<-rep(0, length(pow)*nGenes)
   for(j in 1:length(pow))
    for(iGene in 1:nGenes){
      if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
      indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
      if (pow[j] < Inf)
        Ts[(j-1)*nGenes+iGene] = (sum(U[indx]^pow[j]))
        # (-1)^(1/3)=NaN!
        #Ts[(j-1)*nGenes+iGene] = (sum(U[indx]^pow[j]))^(1/pow[j])
      else Ts[(j-1)*nGenes+iGene] = max(abs(U[indx]))
      }

   # Permutations: 
   T0s = matrix(0, nrow=B, ncol=length(pow)*nGenes)
   for(b in 1:B){
     Y0 <- sample(Y, length(Y))
     #########Null score vector:
     U0<-t(Xg) %*% (Y0-mean(Y0))

     # test stat's:
     for(j in 1:length(pow))
      for(iGene in 1:nGenes){
        if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
        indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
        if (pow[j] < Inf)
          T0s[b, (j-1)*nGenes+iGene] = (sum(U0[indx]^pow[j]))
        else T0s[b, (j-1)*nGenes+iGene] = max(abs(U0[indx]))
        }
     }

   # obtian means/sds for each gene-level stat based on permutations:
   mu0s<-apply(T0s, 2, mean)
   sd0s<-apply(T0s, 2, sd)
   sd0s<-ifelse(sd0s<1e-20, 1e-20, sd0s)
   StdTs<-(Ts - mu0s)/sd0s
   StdT0s<-T0s
   for(b in 1:B)
     StdT0s[b,] = (T0s[b,] - mu0s)/sd0s

   #combine gene-level stats to obtain pathway-lelev stats:
   Ts1<-Ts2<-rep(0, length(pow))
   T0s1<-T0s2<-matrix(0, nrow=B, ncol=length(pow))
   for(j in 1:length(pow)){
       Ts1[j] = max(abs(Ts[((j-1)*nGenes+1):(j*nGenes)]))
       Ts2[j] = max(abs(StdTs[((j-1)*nGenes+1):(j*nGenes)]))
       for(b in 1:B){
         T0s1[b, j] = max(abs(T0s[b, ((j-1)*nGenes+1):(j*nGenes)]))
         T0s2[b, j] = max(abs(StdT0s[b, ((j-1)*nGenes+1):(j*nGenes)]))
         }
       } 

   # permutation-based p-values:
   pPerm1 = pPerm2 = rep(NA, length(pow));
   pvs = NULL;

      for(j in 1:length(pow)) {
        pPerm1[j] = sum( abs(Ts1[j]) < abs(T0s1[,j]))/B
        pPerm2[j] = sum( abs(Ts2[j]) < abs(T0s2[,j]))/B
        }
      P0s1 = PermPvs(T0s1)
      P0s2 = PermPvs(T0s2)
      minP0s1 = apply(P0s1, 1, min)
      minP0s2 = apply(P0s2, 1, min)
      minP1 =  sum( min(pPerm1) > minP0s1 )/B
      minP2 =  sum( min(pPerm2) > minP0s2 )/B
      pvs<-list(unstdPs=c(pPerm1, minP1), stdPs=c(pPerm2, minP2) )
 return(pvs) 
}  

 

# Main function for the SPU tests and the adaptive SPU (aSPU) test:

SPUPtest<-function(Y, X, pow=1:8, B=200){

n<-length(Y)
k<-ncol(X)

Xg<-X
#calculate weights based on ONLY controls:
#Xg0<-Xg[Y==0,]
#freq<-apply(Xg0, 2, sum)
#p0s<-(freq+1)/(2*nrow(Xg0)+2)
#wts<-1/sqrt(n*p0s*(1-p0s))
#for(i in 1:n)
#   Xg[i,]<-Xg[i,]*wts


#######construction of the score vector:
U<-t(Xg) %*% (Y-mean(Y))

   # test stat's:
   Ts<-rep(0, length(pow))
   for(j in 1:length(pow)){
     if (pow[j] < Inf)
        Ts[j] = sum(U^pow[j]) else Ts[j] = max(abs(U))
     }

   # Permutations: 
   T0s = matrix(0, nrow=B, ncol=length(pow))
   for(b in 1:B){
     Y0 <- sample(Y, length(Y))
#     Xg <- X
#     #calculate weights based on ONLY controls:
#     Xg0<-Xg[Y0==0,]
#     freq<-apply(Xg0, 2, sum)
#     p0s<-(freq+1)/(2*nrow(Xg0)+2)
#     wts<-1/sqrt(n*p0s*(1-p0s))
#     for(i in 1:n)
#       Xg[i,]<-Xg[i,]*wts
     #########Null score vector:
     U0<-t(Xg) %*% (Y0-mean(Y0))

     # test stat's:
     for(j in 1:length(pow))
       if (pow[j] < Inf)
         T0s[b, j] = sum(U0^pow[j]) else T0s[b, j] = max(abs(U0))

     }

   # permutation-based p-values:
   pPerm = pPerm0 = rep(NA, length(pow));
   pvs = NULL;

      for(j in 1:length(pow)) {
        pPerm0[j] = sum( abs(Ts[j]) < abs(T0s[,j]))/B
        }
      P0s = PermPvs(T0s)
      minP0s = apply(P0s, 1, min)
      minP =  sum( min(pPerm0) > minP0s )/B
      pvs<-c(pPerm0, minP)
 return(pvs) 
}  

 

