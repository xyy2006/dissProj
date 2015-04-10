## calculate logit p-values for each SNPs
##

#library(postgwas)

getlogitp <- function(dat) {
    n.Y <- length(dat$Y)
    nc.X <- ncol(dat$X)

    logitPs <- NULL
    for(i in 1:nc.X) { # i  =1
        tdat1<-data.frame(trait=dat$Y, X = dat$X[,i])
        fit1<-glm(trait~.,family="binomial",data=tdat1)
        logitPs <- c(logitPs, summary(fit1)$coeff[2,4])
    }
    logitPs
}


pvalcorr <- function(rmat) {
    d <- ncol(rmat)
    PcorMat <- matrix(0, nrow = d, ncol = d)
    for(i in 1:d)
        for(j in 1:d) {
            x <- abs(rmat[i,j])
            PcorMat[i,j] = 0.2982 * x^6 - 0.0127 * x^5 + 0.0588 * x^4 + 0.0099 * x^3 +
                0.6281 * x^2 - 0.0009 * x
        }
    PcorMat
}

GatesSimesnHyst <- function(pvec, X, snp.info, gene.info) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) {
        if(sum(snp.info[,2]==gene.info[g,1]) == 1 ) {
            pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
            M <- length(pvecG)
            Me <- 1

        } else {
            pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
            XG <- X[,snp.info[,2]==gene.info[g,1]]
            M <- length(pvecG)

            rmat <- cor(XG)
            pcmat <- pvalcorr(rmat)
            pcmat.svd <- svd(pcmat)

            egval <- pcmat.svd$d
            Me <- NULL
            for(i in 1:M) # i = 1
                Me <- c(Me, i - sum( as.numeric( egval[1:i] > 1 ) * (egval[1:i] - 1)  ) )

            if( sum(Me<0) > 0 )
                Me[Me < 0] <- which(Me < 0)
        }

        PG <- min(sort(pvecG) * Me[M] / Me)
        keyG <- which(sort(pvecG) * Me[M] / Me == PG)

        PGs <- c(PGs, PG)
        keyGs <- c(keyGs, gpos + keyG)
        gpos = gpos + M
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(cov(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)

    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)

    pval <- c(PgatesSime, Physt)
    names(pval) <- c("PgatesSimes", "Physt")
    pval
}


GatesSimesnHyst2 <- function(pvec, X, snp.info, gene.info) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) {

        pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
        M <- length(pvecG)
        Me <- 1

        if(M > 1) {
            order.pvecG <- order(pvecG)
            XG <- X[,snp.info[,2]==gene.info[g,1]] ## SNP data for given gene
            for(i in 2:M) { #i = 2
                XGn <- XG[,order.pvecG[1:i]]
                rmat <- cor(XGn)
                ## estimate correlation of p-values
                pcmat <- pvalcorr(rmat)

                ## save eigen values of corr matrix
                pcmat.svd <- svd(pcmat)
                egval <- pcmat.svd$d

                Me <- c(Me, i - sum( as.numeric( egval[1:i] > 1 ) * (egval[1:i] - 1) ))
            }
        }
        PG <- min(sort(pvecG) * Me[M] / Me)
        keyG <- which(sort(pvecG) * Me[M] / Me == PG)

        PGs <- c(PGs, PG)
        keyGs <- c(keyGs, gpos + keyG)
        gpos = gpos + M
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(cov(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)
    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)
    Physt2 = 1 - pchisq( Hyst , df = 2*n.gene)


    pvec
    M2 <- length(pvec)
    Me2 <- 1

    order.pvec <- order(pvec)
    for(i in 2:M2) { #i = 2
        XGn <- X[,order.pvec[1:i]]
        rmat <- cor(XGn)
        ## estimate correlation of p-values
        pcmat <- pvalcorr(rmat)

        ## save eigen values of corr matrix
        pcmat.svd <- svd(pcmat)
        egval <- pcmat.svd$d

        Me2 <- c(Me2, i - sum( as.numeric( egval[1:i] > 1 ) * (egval[1:i] - 1) ))
    }

    Pgates <- min(sort(pvec) * Me2[M2] / Me2)

    pval <- c(PgatesSime, Physt, Physt2, Pgates)
    names(pval) <- c("PgatesSimes", "Physt", "Physt2", "Pgates")
    pval
}



GatesSimesnHyst3 <- function(pvec, X, snp.info, gene.info, controls ) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) {

        pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
        M <- length(pvecG)
        Me <- 1

        if(M > 1) {
            order.pvecG <- order(pvecG)
            XG <- X[controls, snp.info[,2]==gene.info[g,1]] ## SNP data for given gene
            for(i in 2:M) { #i = 2
                XGn <- XG[, order.pvecG[1:i]]
                rmat <- cor(XGn)
                ## estimate correlation of p-values
                pcmat <- pvalcorr(rmat)

                ## save eigen values of corr matrix
#                pcmat.svd <- svd(pcmat)
                egval <- eigen(pcmat, only.values=TRUE)$values

                Me <- c(Me, i - sum( as.numeric( egval[1:i] > 1 ) * (egval[1:i] - 1) ))
            }
        }
        PG <- min(sort(pvecG) * Me[M] / Me)
        keyG <- which(sort(pvecG) * Me[M] / Me == PG)

        PGs <- c(PGs, PG)
        keyGs <- c(keyGs, gpos + keyG)
        gpos = gpos + M
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(cov(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)
    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)
    Physt2 = 1 - pchisq( Hyst , df = 2*n.gene)


    pvec
    M2 <- length(pvec)
    Me2 <- 1

    order.pvec <- order(pvec)
    for(i in 2:M2) { #i = 2
        XGn <- X[controls, order.pvec[1:i]]
        rmat <- cor(XGn)
        ## estimate correlation of p-values
        pcmat <- pvalcorr(rmat)

        ## save eigen values of corr matrix
        pcmat.svd <- svd(pcmat)
        egval <- pcmat.svd$d

        Me2 <- c(Me2, i - sum( as.numeric( egval[1:i] > 1 ) * (egval[1:i] - 1) ))
    }

    Pgates <- min(sort(pvec) * Me2[M2] / Me2)

    pval <- c(PgatesSime, Physt, Physt2, Pgates)
    names(pval) <- c("PgatesSimes", "Physt", "Physt2", "Pgates")
    pval
}



GatesSimesnHyst4 <- function(pvec, X, snp.info, gene.info, controls ) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) {

        pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
        M <- length(pvecG)
        Me <- 1

        if(M > 1) {
            order.pvecG <- order(pvecG)
            XG <- X[controls, snp.info[,2]==gene.info[g,1]] ## SNP data for given gene
            for(i in 2:M) { #i = 2
                XGn <- XG[, order.pvecG[1:i]]
                rmat <- cor(XGn)
                ## estimate correlation of p-values
                pcmat <- 0.7723 * rmat^6 - 1.5659 * rmat^5 + 1.201 *
                    rmat^4 - 0.2355 * rmat^3 + 0.2184 * rmat^2 +
                        0.6086 * rmat

                ## save eigen values of corr matrix
#                pcmat.svd <- svd(pcmat)
                egval <- eigen(pcmat, only.values=TRUE)$values

                Me <- c(Me, i - sum( as.numeric( egval[1:i] > 1 ) * (egval[1:i] - 1) ))
            }
        }
        PG <- min(sort(pvecG) * Me[M] / Me)
        keyG <- which(sort(pvecG) * Me[M] / Me == PG)

        PGs <- c(PGs, PG)
        keyGs <- c(keyGs, gpos + keyG)
        gpos = gpos + M
    }

#    cat("keyGs :", keyGs);
#    cat("\n");
#    cat("PGs :", PGs);

    PgatesSime <- min( PGs * n.gene / rank(PGs) )

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(cov(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)
    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)
    Physt2 = 1 - pchisq( Hyst , df = 2*n.gene)

    pval <- c(PgatesSime, Physt, Physt2)
    names(pval) <- c("PgatesSimes", "Physt", "Physt2")
    pval
}




GATES2 <- function (ldmatrix, p)
{
    snpcount <- length(p)
    if (!all(dim(ldmatrix) == snpcount))
        stop("function GATES: Argument 'ldmatrix' is not rectangular or does not match the length of argument vector 'p'.\n")
    if (any(is.na(ldmatrix)))
        stop("function GATES: Argument 'ldmatrix' may not contain NA values.\n")
    if (snpcount < 1)
        stop("function GATES: No SNP provided.\n")
    if (snpcount == 1)
        return(p)
    ldmatrix <- 0.7723 * ldmatrix^6 - 1.5659 * ldmatrix^5 + 1.201 *
        ldmatrix^4 - 0.2355 * ldmatrix^3 + 0.2184 * ldmatrix^2 +
        0.6086 * ldmatrix
    eff.snpcount.fun <- function(ldmat) {
        ldmat <- as.matrix(ldmat)
        snpcount.local <- dim(ldmat)[1]
        if (snpcount.local <= 1)
            return(1)
        ev <- eigen(ldmat, only.values = TRUE)$values
        ev <- ev[ev > 1]
        snpcount.local - sum(ev - 1)
    }
    eff.snpcount.global <- eff.snpcount.fun(ldmatrix)

    candid <- sapply(1:snpcount, function(i) {
        (eff.snpcount.global * p[i]) / eff.snpcount.fun(ldmatrix[1:i,1:i])
    }
                     )

    Pg <- min(candid)
    keyGloc <- which(min(candid) == candid)[1]

    return(c(Pg, keyGloc))
}




GatesSimes <- function(pvec, X, snp.info, gene.info, controls ) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) { # g = 1

        if( sum(snp.info[,2]==gene.info[g,1]) == 1 ) {
            Pg <- pvec[snp.info[,2]==gene.info[g,1]]
        } else {
            pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
            Xn <- X[controls, snp.info[,2]==gene.info[g,1]]
            o.pv <- order(pvecG)
            ldmat <- cor(Xn[,o.pv])
            Pg <- GATES2(ldmatrix = ldmat, p = sort(pvecG) )[1]
        }

        PGs <- c(PGs, Pg)
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )
    PgatesSime
}



#    out.H <- Hyst(pvec = pvec1, X = dat1$X, snp.info = dat1$snp.info,
#                  gene.info = dat1$gene.info, controls = dat1$Y==0 )


Hyst <- function(pvec, X, snp.info, gene.info, controls ) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) { #g = 6
        if( sum(snp.info[,2]==gene.info[g,1]) == 1 ) {
            Pg <- pvec[snp.info[,2]==gene.info[g,1]]
            keyG <- 1
            M <- 1
        } else {
            pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
            M <- length(pvecG)
            Xn <- X[controls, snp.info[,2]==gene.info[g,1]]
            o.pv <- order(pvecG)
            ldmat <- cor(Xn[,o.pv])
            out <- GATES2(ldmatrix = ldmat, p = sort(pvecG) )
            keyG <- out[2]
            Pg <- out[1]
        }

        PGs <- c(PGs, Pg)
        keyGs <- c(keyGs, gpos + keyG)
        gpos = gpos + M

    }
#    cat("keyGs :", keyGs);
#    cat("\n");
#    cat("PGs :", PGs);

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(cor(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)
    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)
#    Physt2 = 1 - pchisq( Hyst , df = 2*n.gene)
}
