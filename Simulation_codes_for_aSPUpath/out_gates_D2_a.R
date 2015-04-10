source("../programs/simINDpath.R")
source("../programs/simAR1path.R")
source("../programs/simAR1path2.R")
source("../programs/gates3.R")
source("../programs/UminPP1.r")
source("../programs/aSPUpath.R")

LORs=c(0.1, 0.15, 0.2, 0.25, 0.3, .35, 0.4)
nLOR <- length(LORs)
sims <- 1:1000

Res <- NULL;

for(iLOR in 1:nLOR){

    pvs <- NULL;
    for(isim in sims){
        set.seed(isim)

        dat1<-simPathIndSNPv2(nGenes=20, nGenes1=10, nSNPlow=1, nSNPup=20, nSNP0=1:3,
                            LOR=LORs[iLOR],
                            n=500, MAFslow=0.05, MAFsup=0.4, p0=0.05)

        p.aspup <- aSPUpath(Y = dat1$Y, X = dat1$X, nSNPs=dat1$nSNPs,
                            pow = 1:8, pow2= c(1,2,4,8),
                            model="binomial", n.perm = 500 )$stdPs[37]

        p.uminp <- UminPP(Y = dat1$Y, X = dat1$X, B=500)

        pvec1 <- getlogitp(dat1)
        o.pv <- order(pvec1)
        ldmat <- cor(dat1$X[dat1$Y==0, o.pv])

        out.Gs <- GatesSimes(pvec = pvec1, X = dat1$X, snp.info = dat1$snp.info,
                             gene.info = dat1$gene.info, controls = dat1$Y==0 )

        out.H <- Hyst(pvec = pvec1, X = dat1$X, snp.info = dat1$snp.info,
                      gene.info = dat1$gene.info, controls = dat1$Y==0 )

        out1 <- c(p.aspup, p.uminp, out.Gs, out.H)

        pvs <- rbind(pvs, out1 )

    }
    Res[[iLOR]] = pvs
}
names(Res) <- LORs

save(Res, file = "out_gates_D2_a.RData")

