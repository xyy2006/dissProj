source("../programs/simINDpath.R")
source("../programs/simAR1path.R")
source("../programs/simAR1path2.R")
source("../programs/gates3.R")
source("../programs/UminPP1.r")
library(SNPath)

LORs=c(0.1, 0.15, 0.2, 0.25, 0.3, .35, 0.4, .45, .5)
nLOR <- length(LORs)
sims <- 1:1000

Res <- NULL;

for(iLOR in 1:nLOR){

    pvs <- NULL;
    for(isim in sims){
        set.seed(isim)

        dat1<-simPathAR1SNPv2(nGenes=20, nGenes1=10, nSNPlow=1, nSNPup=20, nSNP0=1:3,
                            LOR=LORs[iLOR],
                            n=500, MAFslow=0.05, MAFsup=0.4, p0=0.05)

        p.grass <- grass(snp.dat=t(dat1$X), snp.info=dat1$snp.info,
                         gene.info=dat1$gene.info, gene.set=dat1$pathway,
                         y=dat1$Y, gene.def="rel", B = 200)

        pvs <- c(pvs, p.grass[[1]] )

    }
    Res[[iLOR]] = pvs
}
names(Res) <- LORs

save(Res, file = "out_grass_D2_c.RData")

