######## simulated a pathway with multiple genes, each one with
######## SNPs from a latent multivariate Gaussian variable with
######## an Ind correlation structure. The causal SNP0 is the
######## first SNP in each causal gene.
######## Wei Pan, weip@biostat.umn.edu
######## 9/23/11

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
#               p0: background disease prevalence;i.e. intercept=log(p0/(1-p0))
########Output: a list of the binary outcome Y (=0 or 1) and SNPs (=0, 1 or 2);
#               Y is a vector of length 2n; X is a matrix of 2n by nSNP.

simPathIndSNP<-function(nGenes=20, nGenes1=5, nSNPs=NULL,
                 nSNPlow=1, nSNPup=20, nSNP0=1,
                 LOR=0,
                 n=500, MAFslow=0.05, MAFsup=0.4, p0=0.05){

if (is.null(nSNPs)) nSNPs=sample(nSNPlow:nSNPup, nGenes, replace=T)

## total # of SNPs:
q<-sum(nSNPs)

##background disease prev:
b0<-log(p0/(1-p0))
##LOR:
b1<-rep(0, q)
for(i in 1:nGenes1)
  if (i==1) b1[1:nSNP0]=LOR   else b1[sum(nSNPs[1:(i-1)])+(1:nSNP0)]=LOR

MAFs<-runif(q, MAFslow, MAFsup)
cutoff<-qnorm(MAFs)

X<-matrix(0, nrow=n+n, ncol=q)
Y<-rep(0, n+n); Y[(n+1):(2*n)]<-1
i<-1
#sampling controls:
while ( i <= n){
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X2<-ifelse(X0<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X3<-ifelse(X0<cutoff, 1, 0)
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
  X2<-ifelse(X0<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X3<-ifelse(X0<cutoff, 1, 0)
  X4<-X2+ X3
  pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==1){
    X[i, ]<-X4
    i<-i+1
    }
  }

#######
snp.chrom=snp.loc=NULL
for(i in 1:nGenes){
  snp.chrom=c(snp.chrom, rep(i, nSNPs[i]))
  snp.loc=c(snp.loc, 1:nSNPs[i])
  }
snp.info = cbind(1:q, snp.chrom, snp.loc)

#gene.start<-gene.end<-rep(0, nGenes)
#for(i in 1:nGenes)
#  if (i==1) {gene.start[i]=1; gene.end[i]=nSNPs[i];}
#  else  {gene.start[i]=sum(nSNPs[1:(i-1)])+1; gene.end[i]=sum(nSNPs[1:i]) ;}
#gene.info=cbind(1:nGenes, 1:nGenes, gene.start, gene.end)
gene.info=cbind(1:nGenes, 1:nGenes, rep(0, nGenes), nSNPs)

pathway = list(pathway1=as.character(1:nGenes))

list(Y=Y, X=X, snp.info=snp.info, gene.info=gene.info, pathway=pathway,
     nSNPs=nSNPs)

}


simPathIndSNPv2<-function(nGenes=10, nGenes1=5, nSNPs=NULL, ncSNPs = NULL,
                 nSNPlow=1, nSNPup=20, nSNP0=1:3,
                 LOR=.5,
                 n=50, MAFslow=0.05, MAFsup=0.4, p0=0.05){

if (is.null(nSNPs)) nSNPs=sample(nSNPlow:nSNPup, nGenes, replace=T)
if (is.null(ncSNPs)) ncSNPs=pmin(nSNPs, sample(nSNP0, nGenes1, replace=T) )

## total # of SNPs:
q<-sum(nSNPs)

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
  X2<-ifelse(X0<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X3<-ifelse(X0<cutoff, 1, 0)
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
  X2<-ifelse(X0<cutoff, 1, 0)
  X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
  X3<-ifelse(X0<cutoff, 1, 0)
  X4<-X2+ X3
  pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==1){
    X[i, ]<-X4
    i<-i+1
    }
  }

#######
snp.chrom=snp.loc=NULL
for(i in 1:nGenes){
  snp.chrom=c(snp.chrom, rep(i, nSNPs[i]))
  snp.loc=c(snp.loc, 1:nSNPs[i])
  }
snp.info = cbind(1:q, snp.chrom, snp.loc)

#gene.start<-gene.end<-rep(0, nGenes)
#for(i in 1:nGenes)
#  if (i==1) {gene.start[i]=1; gene.end[i]=nSNPs[i];}
#  else  {gene.start[i]=sum(nSNPs[1:(i-1)])+1; gene.end[i]=sum(nSNPs[1:i]) ;}
#gene.info=cbind(1:nGenes, 1:nGenes, gene.start, gene.end)
gene.info=cbind(1:nGenes, 1:nGenes, rep(0, nGenes), nSNPs)

pathway = list(pathway1=as.character(1:nGenes))

list(Y=Y, X=X, snp.info=snp.info, gene.info=gene.info, pathway=pathway,
     nSNPs=nSNPs)

}
