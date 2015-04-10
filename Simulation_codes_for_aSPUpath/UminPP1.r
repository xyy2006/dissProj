#############permutation-based UminP test.
############
######################Wei Pan, weip@biostat.umn.edu, 9/30/11

#Input:
#      Y: a vector of 0-1 response, nx1;
#      X: design matrix for covariates (e.g. gene-gene interactions) of
#           interest (i.e. to be tested), nxk;
#      Z: design matrix for nuisance covariates, nxk2;
#Goal:
#      fit model: logit(Pr(Y=1)) = a0+ a1*X+a2*Z ,
#      test H0: a1=0.

#Examples: suppose that one has two groups of variables in design matrices
#          X1 and X2 respectively, and Y is a binary response vector,
#          1) if one would like to test the marginal effects of X1, then
#             Tests5V1(Y, X=X1);
#          2) if one would like to test the effects of X1 after ajusting for X2,
#             Tests5V1(Y, X=X1, Z=X2);
#          3) if one would like to test the effects of X1 in the presence of
#             possible interactions with X2, then
#             Tests5V1(Y, X=cbind(X1, X12), Z=X2),
#             where X12 is the design matrix for interactions;
#          4) if one would like to test on the interaction between X1 and X2
#             after adjusting for the main effects of X1 and X2, then
#             Tests5V1(Y, X=X12, Z=cbind(X1,X2)),
#             where X12 is the design matrix for interactions;

############################################
############################################
############################################
# Main function for the SPU tests and the adaptive SPU (aSPU) test:
#Input:
#      Y: a vector of 0-1 response, nx1;
#      X: design matrix for genotypes, nxk1;
#      Z: design matrix for (nuisance) covariates, nxk2;
#      pow: the power's used to construct SPU tests and obatining the aSPU test;
#      B=200: number of permutations to obtain p-values;
#Goal:
#      fit model: logit(Pr(Y=1)) = a0+ a1*X ,
#      test H0: a1=0.
#Output: p-value of the UminP test based on permutations

UminPP<-function(Y, X, Z=NULL, B=200 ){

if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
n<-length(Y)
k<-ncol(X)
#k2<-ncol(Z)

#######construction of the score vector and its cov matrix:
if (is.null(Z)){
## NO nuisance parameters:
   Xg <- X
   Xbar<-apply(Xg, 2, mean)
   Xgb<-Xg
   for(i in 1:nrow(Xg))
      Xgb[i,]<-Xg[i,]-Xbar
   ##score vector:
   U<-t(Xg) %*% (Y-mean(Y))
   ##cov of the score stats:
   CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)

} else {
## with nuisance parameters:
tdat1<-data.frame(trait=Y, Z)
fit1<-glm(trait~.,family="binomial",data=tdat1)
pis<-fitted.values(fit1)
Us<-XUs<-matrix(0, nrow=n, ncol=k)
Xmus = X
for(i in 1:k){
    tdat2<-data.frame(X1=X[,i], Z)
    fit2<-glm(X1~.,data=tdat2)
    Xmus[,i]<-fitted.values(fit2)
    XUs[, i]<-(X[,i] - Xmus[,i])
    }
U<-t(XUs) %*% (Y - pis)

CovS<-matrix(0, nrow=k, ncol=k)
for(i in 1:n)
  CovS<-CovS + pis[i]*(1-pis[i])*XUs[i,] %*% t(XUs[i,])
}

#test stat:
if (k==1) T = U^2/(CovS[1,1]+1e-20)
else T = max(U^2/(diag(CovS)+1e-20))

T0s<-rep(0, B)
   # bootstrap:
   Y0=Y
   for(b in 1:B){
     if (is.null(Z)) {
       Y0 <- sample(Y, length(Y))
       #########Null score vector:
       U0<-t(Xg) %*% (Y0-mean(Y0))
       }
     else{
       ## with nuisance parameters:
       for(j in 1:n)
         Y0[j] <- sample(c(1,0), 1, prob=c(pis[j], 1-pis[j]) )
       tdat0<-data.frame(trait=Y0, Z)
       fit0<-glm(trait~.,family="binomial",data=tdat0)
       pis0<-fitted.values(fit0)
       U0<-t(XUs) %*% (Y0 - pis0)
       }
   if (k==1) T0s[b] = U0^2/(CovS[1,1]+ 1e-20)
   else T0s[b] = max(U0^2/(diag(CovS) +1e-20))
   }
pv =sum(T0s > T)/B
return(pv)
}
