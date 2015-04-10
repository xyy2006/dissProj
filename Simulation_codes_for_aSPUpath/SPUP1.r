###################### Sum of powered score (SPU) tests and aSPU test
###################### based on permuttaions (and simulation may NOT work
###################### for "smalln, large p",
###################### hence unable to adjust for possible covariates);
###################### modified from aSPU/prog/wSPUP1.r
######################Wei Pan, weip@biostat.umn.edu, 9/23/11

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

#calculate a permuttaion p-value matrix based on a permuted stat matrix:
PermPvs<-function(T0s){
B=nrow(T0s); n=ncol(T0s)
P0s<-matrix(1, nrow=B, ncol=n)
for(j in 1:n)
  for(b in 1:B)
    P0s[b,j] = sum( abs(T0s[b,j]) < abs(T0s[-b,j]) )/(B-1)
return(P0s)
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

 

