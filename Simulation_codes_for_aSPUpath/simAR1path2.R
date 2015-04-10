######## Similar to simAR1path.R except that teh causal SNPs are removed;
######## simulated a pathway with multiple genes, each one with
######## SNPs from a latent multivariate Gaussian variable with
######## an AR1 correlation structure. The causal SNP0 is the
######## firsdt SNP in each causal gene.
######## Wei Pan, weip@biostat.umn.edu
######## 9/26/11

########Input:
#               nGenes: # of the genes in a pathway;
#               nGenes1: # of the genes with causal SNPs;
#               nSNPlow, nSNPup: # of SNPs inside each gene is drawn from
#                                Unif(nSNPlow, nSNPup);
#               nSNP0: # of causal SNPs in each gene;
#               LOR: association in log OR between a causal SNP and outcome;
#               MAF0: MAF for a causal SNP0;
#               n: # of cases (= # of controls);
#               MAFslow, MAFsup: MAF's of the SNPs are drawn from
#                                 Unif(MAFslow, MAFsup);
#               rholow, rhoup: the SNPs in eahc gene are from a latent Normal
#                    variable with a AR(rho) corr structure, rho's are drawn
#                    from Unif(rholow, rhoup); the SNPs in diff genes are indep;
#               p0: background disease prevalence;i.e. intercept=log(p0/(1-p0))
########Output: a list of the binary outcome Y (=0 or 1) and SNPs (=0, 1 or 2);
#               Y is a vector of length 2n; X is a matrix of 2n by nSNP.

simPathAR1SNP2<-function(nGenes=20, nGenes1=1, nSNPs=NULL,
                 nSNPlow=1, nSNPup=20, nSNP0=1,
                 LOR=0,
                 n=500, MAFslow=0.05, MAFsup=0.4, rholow=0, rhoup=0.8, p0=0.05){

if (is.null(nSNPs)) nSNPs=sample(nSNPlow:nSNPup, nGenes, replace=T)
## total # of SNPs:
q<-sum(nSNPs)

rhos=runif(nGenes, rholow, rhoup)

###get the overall corr matrix for the laten var's:
R<-matrix(0, nrow=q, ncol=q)
for(iGene in 1:nGenes){
  if(iGene==1) Rstart=0 else Rstart=sum(nSNPs[1:(iGene-1)])
  for(i in 1:nSNPs[iGene])
    for(j in 1:nSNPs[iGene])
    R[Rstart+i, Rstart+j]<-rhos[iGene]^(abs(i-j))
  }
svd.R<-svd(R)
R1<-svd.R$u %*% diag(sqrt(svd.R$d))
#R2<- diag(sqrt(svd.R$d)) %*% t(svd.R$v)
#Note above: R1 %*% t(R1)= R

##background disease prev:
b0<-log(p0/(1-p0))
##LOR:
b1<-rep(0, q)
for(i in 1:nGenes1)
  if (i==1) b1[1:nSNP0]=LOR
  else b1[sum(nSNPs[1:(i-1)])+(1:nSNP0)]=LOR

MAFs<-runif(q, MAFslow, MAFsup)
cutoff<-qnorm(MAFs)

X<-matrix(0, nrow=n+n, ncol=q)
Y<-rep(0, n+n); Y[(n+1):(2*n)]<-1
i<-1
#sampling controls:
while ( i <= n){
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X2<-ifelse(X1<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X3<-ifelse(X1<cutoff, 1, 0)
  X4<-X2+ X3
  pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==0){
    X[i, ]<-X4
    i<-i+1
    }
  }
#sampling cases:
while ( i <= 2*n){
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X2<-ifelse(X1<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X3<-ifelse(X1<cutoff, 1, 0)
  X4<-X2+ X3
  pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==1){
    X[i, ]<-X4
    i<-i+1
    }
  }

######exclude causal SNPs:
X = X[, abs(b1)<1e-10]

#######
snp.chrom=snp.loc=NULL
for(i in 1:nGenes1){
  snp.chrom=c(snp.chrom, rep(i, nSNPs[i]-nSNP0))
  if (nSNPs[i]-nSNP0>=1)
    snp.loc=c(snp.loc, 1:(nSNPs[i]-nSNP0))
  }
if (nGenes1 < nGenes){
for(i in (nGenes1+1):nGenes){
  snp.chrom=c(snp.chrom, rep(i, nSNPs[i]))
  snp.loc=c(snp.loc, 1:nSNPs[i])
  }
}
snp.info = cbind(1:(q-nGenes1*nSNP0), snp.chrom, snp.loc)

#gene.start<-gene.end<-rep(0, nGenes)
#for(i in 1:nGenes)
#  if (i==1) {gene.start[i]=1; gene.end[i]=nSNPs[i];}
#  else  {gene.start[i]=sum(nSNPs[1:(i-1)])+1; gene.end[i]=sum(nSNPs[1:i]) ;}
#gene.info=cbind(1:nGenes, 1:nGenes, gene.start, gene.end)
nSNPsNone0=nSNPs
for(i in 1:nGenes1)
  nSNPsNone0[i] = nSNPs[i]-nSNP0
gene.info=cbind(1:nGenes, 1:nGenes, rep(0, nGenes), nSNPsNone0)

pathway = list(pathway1=as.character(1:nGenes))

list(Y=Y, X=X, snp.info=snp.info, gene.info=gene.info, pathway=pathway,
     nSNPs=nSNPsNone0)

}


######## Similar to simAR1path.R except that teh causal SNPs are removed;
######## simulated a pathway with multiple genes, each one with
######## SNPs from a latent multivariate Gaussian variable with
######## an AR1 correlation structure. The causal SNP0 is the
######## firsdt SNP in each causal gene.
######## Wei Pan, weip@biostat.umn.edu
######## 9/26/11

########Input:
#               nGenes: # of the genes in a pathway;
#               nGenes1: # of the genes with causal SNPs;
#               nSNPlow, nSNPup: # of SNPs inside each gene is drawn from
#                                Unif(nSNPlow, nSNPup);
#               nSNP0: # of causal SNPs in each gene;
#               LOR: association in log OR between a causal SNP and outcome;
#               MAF0: MAF for a causal SNP0;
#               n: # of cases (= # of controls);
#               MAFslow, MAFsup: MAF's of the SNPs are drawn from
#                                 Unif(MAFslow, MAFsup);
#               rholow, rhoup: the SNPs in eahc gene are from a latent Normal
#                    variable with a AR(rho) corr structure, rho's are drawn
#                    from Unif(rholow, rhoup); the SNPs in diff genes are indep;
#               p0: background disease prevalence;i.e. intercept=log(p0/(1-p0))
########Output: a list of the binary outcome Y (=0 or 1) and SNPs (=0, 1 or 2);
#               Y is a vector of length 2n; X is a matrix of 2n by nSNP.

simPathAR1SNP2v2<-function(nGenes=10, nGenes1=5, nSNPs=NULL, ncSNPs = NULL,
                 nSNPlow=1, nSNPup=20, nSNP0=1:3,
                 LOR=0.3,
                 n=100, MAFslow=0.05, MAFsup=0.4, rholow=0, rhoup=0.8, p0=0.05){

if (is.null(nSNPs)) nSNPs=sample(nSNPlow:nSNPup, nGenes, replace=T)
if (is.null(ncSNPs)) ncSNPs=pmin(nSNPs[1:nGenes1], sample(nSNP0, nGenes1, replace=T) )
## total # of SNPs:
q<-sum(nSNPs)

rhos=runif(nGenes, rholow, rhoup)

###get the overall corr matrix for the laten var's:
R<-matrix(0, nrow=q, ncol=q)
for(iGene in 1:nGenes){
  if(iGene==1) Rstart=0 else Rstart=sum(nSNPs[1:(iGene-1)])
  for(i in 1:nSNPs[iGene])
    for(j in 1:nSNPs[iGene])
    R[Rstart+i, Rstart+j]<-rhos[iGene]^(abs(i-j))
  }
svd.R<-svd(R)
R1<-svd.R$u %*% diag(sqrt(svd.R$d))
#R2<- diag(sqrt(svd.R$d)) %*% t(svd.R$v)
#Note above: R1 %*% t(R1)= R

##background disease prev:
b0<-log(p0/(1-p0))
##LOR:
b1<-rep(0, q)
for(i in 1:nGenes1)
  if (i==1) b1[1:ncSNPs[1]]=LOR  else b1[sum(nSNPs[1:(i-1)])+(1:ncSNPs[i])]=LOR

MAFs<-runif(q, MAFslow, MAFsup)
cutoff<-qnorm(MAFs)

X<-matrix(0, nrow=n+n, ncol=q)
Y<-rep(0, n+n); Y[(n+1):(2*n)]<-1
i<-1
#sampling controls:
while ( i <= n){
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X2<-ifelse(X1<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X3<-ifelse(X1<cutoff, 1, 0)
  X4<-X2+ X3
  pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==0){
    X[i, ]<-X4
    i<-i+1
    }
  }
#sampling cases:
while ( i <= 2*n){
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X2<-ifelse(X1<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
  X3<-ifelse(X1<cutoff, 1, 0)
  X4<-X2+ X3
  pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==1){
    X[i, ]<-X4
    i<-i+1
    }
  }

######exclude causal SNPs:
X = X[, abs(b1)<1e-10]

#######
snp.chrom=snp.loc=NULL
for(i in 1:nGenes1){
  snp.chrom=c(snp.chrom, rep(i, nSNPs[i]-ncSNPs[i]))
  if (nSNPs[i]-ncSNPs[i]>=1)
    snp.loc=c(snp.loc, 1:(nSNPs[i]-ncSNPs[i]))
  }
if (nGenes1 < nGenes){
for(i in (nGenes1+1):nGenes){
  snp.chrom=c(snp.chrom, rep(i, nSNPs[i]))
  snp.loc=c(snp.loc, 1:nSNPs[i])
  }
}
snp.info = cbind(1:(q-sum(ncSNPs)), snp.chrom, snp.loc)

#gene.start<-gene.end<-rep(0, nGenes)
#for(i in 1:nGenes)
#  if (i==1) {gene.start[i]=1; gene.end[i]=nSNPs[i];}
#  else  {gene.start[i]=sum(nSNPs[1:(i-1)])+1; gene.end[i]=sum(nSNPs[1:i]) ;}
#gene.info=cbind(1:nGenes, 1:nGenes, gene.start, gene.end)
nSNPsNone0=nSNPs
for(i in 1:nGenes1)
  nSNPsNone0[i] = nSNPs[i]-ncSNPs[i]
gene.info=cbind(1:nGenes, 1:nGenes, rep(0, nGenes), nSNPsNone0)

pathway = list(pathway1=as.character(1:nGenes))

list(Y=Y, X=X, snp.info=snp.info, gene.info=gene.info, pathway=pathway,
     nSNPs=nSNPsNone0)

}
