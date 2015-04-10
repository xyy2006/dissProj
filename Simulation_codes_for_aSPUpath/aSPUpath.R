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


#' Pathway based Adaptive Sum of powered score tests (SPUpath and aSPUpath)
#'
#' It gives p-values of the SPUpath tests and aSPUpath test.
#'
#' @param Y Response or phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#' or It can be any quantitative traits. Vector with length n (number of observations)
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele. Matrix with dimension n by g (n : number of observation, p : number of genotype data)
#'
#' @param cov covariates. Matrix with dimension n by k (n :number of observation, k : number of covariates)
#'
#' @param model Use "gaussian" for quantitative trait, and use "binomial" for binary trait.
#'
#' @param nSNPs a vector; each element is the number of the SNPs in the
#'              corresponding gene
#'
#' @param pow SNP specific power(gamma values) used in SPUpath test.
#'
#' @param pow2 GENE specific power(gamma values) used in SPUpath test.
#'
#' @param n.perm number of permutations.
#'
#' @param usePCs indicating whether to extract PCs and then use PCs of X
#'
#' @param varprop the proportion of the variations explained (cutoff) that
#'                 determines how many first PCs to use.
#'
#' @export
#' @return Test Statistics and p-values for SPUpath tests and aSPUpath test. There are three versions. "std" version is the exactly same version with the paper (Pan, Kwak and Wei 2015). "unnorm" is unnormalized one which erased the power of 1/gamma in outside of bracket and didn't divide the number of SNP for each Gene in gene level SPU statistics. equation (3) on the paper, outside power of 1/gamma is set to 1 and it didn't divide the number of gene k_g. "unstd" version it didn't divide the number of SNP for each Gene. In paper equation (3) we didn't divide the sum of weighted scores with k_g(the number of SNP).
#'
#' @examples
#'
#' dat1<-simPathAR1SNP(nGenes=20, nGenes1=5, nSNPlim=c(1, 20), nSNP0=1,
#'                     LOR=.2, n=100, MAFlim=c(0.05, 0.4), p0=0.05 )
#' p.pathaspu<- aSPUpath(dat1$Y, dat1$X, nSNPs = dat1$nSNPs,
#'          model = "binomial", pow=1:8, pow2=c(1, 2, 4, 8), n.perm=100)
#'
#' @seealso \code{\link{simPathAR1SNP}} \code{\link{aSPUpathSingle}}


aSPUpath <- function (Y, X, cov = NULL, model = c("binomial", "gaussian"),
    nSNPs, pow = 1:8, pow2 = 1, n.perm = 200, usePCs = F, varprop = 0.95)
{
    model = match.arg(model)
    nSNPs = nSNPs[nSNPs >= 1]
    nSNPs0 <- rep(0, length(nSNPs))
    if (usePCs) {
        Xg <- NULL
        for (iGene in 1:length(nSNPs)) {
            if (iGene == 1)
                SNPstart = 1
            else SNPstart = sum(nSNPs[1:(iGene - 1)]) + 1
            indx = (SNPstart:(SNPstart + nSNPs[iGene] - 1))
            Xpcs <- extractPCs(X[, indx], cutoff = varprop)
            Xg <- cbind(Xg, Xpcs)
            if (is.null(ncol(Xpcs)))
                nSNPs0[iGene] = 1
            else nSNPs0[iGene] = ncol(Xpcs)
        }
    }
    else {
        Xg = X
        nSNPs0 = nSNPs
    }
    n <- length(Y)
    k <- ncol(Xg)
    if (is.null(cov)) {
        XUs <- Xg <- X
        r <- Y - mean(Y)
        U <- as.vector(t(Xg) %*% r)
    }
    else {
        tdat1 <- data.frame(trait = Y, cov)
        fit1 <- glm(trait ~ ., family = model, data = tdat1)
        pis <- fitted.values(fit1)
        XUs <- matrix(0, nrow = n, ncol = k)
        Xmus = X
        for (i in 1:k) {
            tdat2 <- data.frame(X1 = X[, i], cov)
            fit2 <- glm(X1 ~ ., data = tdat2)
            Xmus[, i] <- fitted.values(fit2)
            XUs[, i] <- (X[, i] - Xmus[, i])
        }
        r <- Y - pis
        U <- t(XUs) %*% r
    }
    nGenes = length(nSNPs0)
    Ts <- StdTs <- rep(0, length(pow) * nGenes)
    TsUnnorm <- Ts <- StdTs <- rep(0, length(pow) * nGenes)
    for (j in 1:length(pow)) for (iGene in 1:nGenes) {
        if (iGene == 1)
            SNPstart = 1
        else SNPstart = sum(nSNPs0[1:(iGene - 1)]) + 1
        indx = (SNPstart:(SNPstart + nSNPs0[iGene] - 1))
        if (pow[j] < Inf) {
            a = (sum(U[indx]^pow[j]))
            TsUnnorm[(j - 1) * nGenes + iGene] = a
            Ts[(j - 1) * nGenes + iGene] = sign(a) * ((abs(a))^(1/pow[j]))
            StdTs[(j - 1) * nGenes + iGene] = sign(a) * ((abs(a)/nSNPs0[iGene])^(1/pow[j]))
        }
        else {
            TsUnnorm[(j - 1) * nGenes + iGene] = Ts[(j - 1) *
                nGenes + iGene] = StdTs[(j - 1) * nGenes + iGene] = max(abs(U[indx]))
        }
    }
    T0s = StdT0s = matrix(0, nrow = n.perm, ncol = length(pow) *
        nGenes)
    T0sUnnorm = T0s = StdT0s = matrix(0, nrow = n.perm, ncol = length(pow) *
        nGenes)
    for (b in 1:n.perm) {
        Y0 <- sample(Y, length(Y))
        U0 <- t(XUs) %*% (Y0 - mean(Y0))
        for (j in 1:length(pow)) for (iGene in 1:nGenes) {
            if (iGene == 1)
                SNPstart = 1
            else SNPstart = sum(nSNPs0[1:(iGene - 1)]) + 1
            indx = (SNPstart:(SNPstart + nSNPs0[iGene] - 1))
            if (pow[j] < Inf) {
                a = (sum(U0[indx]^pow[j]))
                T0sUnnorm[b, (j - 1) * nGenes + iGene] = a
                T0s[b, (j - 1) * nGenes + iGene] = sign(a) *
                  ((abs(a))^(1/pow[j]))
                StdT0s[b, (j - 1) * nGenes + iGene] = sign(a) *
                  ((abs(a)/nSNPs0[iGene])^(1/pow[j]))
            }
            else T0sUnnorm[b, (j - 1) * nGenes + iGene] = T0s[b,
                (j - 1) * nGenes + iGene] = StdT0s[b, (j - 1) *
                nGenes + iGene] = max(abs(U0[indx]))
        }
    }
    Ts1 <- Ts2 <- TsU <- rep(0, length(pow) * length(pow2))
    T0s1 <- T0s2 <- T0sU <- matrix(0, nrow = n.perm, ncol = length(pow) *
        length(pow2))
    for (j2 in 1:length(pow2)) {
        for (j in 1:length(pow)) {
            TsU[(j2 - 1) * length(pow) + j] = sum(TsUnnorm[((j -
                1) * nGenes + 1):(j * nGenes)]^pow2[j2])
            Ts1[(j2 - 1) * length(pow) + j] = sum(Ts[((j - 1) *
                nGenes + 1):(j * nGenes)]^pow2[j2])
            Ts2[(j2 - 1) * length(pow) + j] = sum(StdTs[((j -
                1) * nGenes + 1):(j * nGenes)]^pow2[j2])
            for (b in 1:n.perm) {
                T0sU[b, (j2 - 1) * length(pow) + j] = sum(T0sUnnorm[b,
                  ((j - 1) * nGenes + 1):(j * nGenes)]^pow2[j2])
                T0s1[b, (j2 - 1) * length(pow) + j] = sum(T0s[b,
                  ((j - 1) * nGenes + 1):(j * nGenes)]^pow2[j2])
                T0s2[b, (j2 - 1) * length(pow) + j] = sum(StdT0s[b,
                  ((j - 1) * nGenes + 1):(j * nGenes)]^pow2[j2])
            }
        }
    }
    pPerm1 = pPerm2 = pPermU = rep(NA, length(pow) * length(pow2))
    pvs = NULL
    for (j in 1:(length(pow) * length(pow2))) {
        pPermU[j] = sum(abs(TsU[j]) < abs(T0sU[, j]))/n.perm
        pPerm1[j] = sum(abs(Ts1[j]) < abs(T0s1[, j]))/n.perm
        pPerm2[j] = sum(abs(Ts2[j]) < abs(T0s2[, j]))/n.perm
    }
    P0s1 = PermPvs(T0s1)
    P0s2 = PermPvs(T0s2)
    P0sU = PermPvs(T0sU)
    minP0s1 = apply(P0s1, 1, min)
    minP0s2 = apply(P0s2, 1, min)
    minP0sU = apply(P0sU, 1, min)
    minP1 = sum(min(pPerm1) > minP0s1)/n.perm
    minP2 = sum(min(pPerm2) > minP0s2)/n.perm
    minPU = sum(min(pPermU) > minP0sU)/n.perm
    minP1s <- minP2s <- minPUs <- rep(NA, length(pow2))
    for (j2 in 1:length(pow2)) {
        minP0s1 = apply(P0s1[, ((j2 - 1) * length(pow) + 1):(j2 *
            length(pow))], 1, min)
        minP0s2 = apply(P0s2[, ((j2 - 1) * length(pow) + 1):(j2 *
            length(pow))], 1, min)
        minP0sU = apply(P0sU[, ((j2 - 1) * length(pow) + 1):(j2 *
            length(pow))], 1, min)
        minP1s[j2] = sum(min(pPerm1[((j2 - 1) * length(pow) +
            1):(j2 * length(pow))]) > minP0s1)/n.perm
        minP2s[j2] = sum(min(pPerm2[((j2 - 1) * length(pow) +
            1):(j2 * length(pow))]) > minP0s2)/n.perm
        minPUs[j2] = sum(min(pPermU[((j2 - 1) * length(pow) +
            1):(j2 * length(pow))]) > minP0sU)/n.perm
    }
    pvs <- list(unstdPs = c(pPerm1, minP1s, minP1), stdPs = c(pPerm2,
        minP2s, minP2), unstdPsUnnorm = c(pPermU, minPUs, minPU))
    return(pvs)
}
