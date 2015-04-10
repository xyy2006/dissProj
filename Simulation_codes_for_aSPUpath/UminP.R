#' permutation-based UminP test.
#'
#' fit model: logit(Pr(Y=1)) = a0+ a1*X+a2*Z, test H0: a1=0.
#'
#' @param Y Response or phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele. Matrix with dimension n by g (n : number of observation, p : number of genotype data)
#'
#' @param Z covariates. Matrix with dimension n by k (n :number of observation, k : number of covariates)
#'
#' @param n.perm number of permutations.
#'
#' @export
#' @return p-value of the UminP test based on permutations
#'
#' @author Wei Pan
#'
#' @examples
#'
#' data(exdat)
#' out <- UminP(Y = exdat$Y, X = exdat$X, Z = NULL, n.perm = 200)
#' out
#'


UminP<-function(Y, X, Z=NULL, n.perm=200 ){

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

T0s<-rep(0, n.perm)
   # bootstrap:
   Y0=Y
   for(b in 1:n.perm){
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
pv =sum(T0s > T)/n.perm
return(pv)
}
