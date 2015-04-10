library(SNPath)
#library(aSPU)
source("../programs/simINDpath.R")
source("../programs/simAR1path.R")
source("../programs/simAR1path2.R")
source("../programs/UminPP1.r")
source("../programs/SPUP1.r")
source("../programs/SPUPpathSingle1.r")
source("../programs/aSPUPpath3.r")

LORs=c(0.1, 0.15, 0.2, 0.25, 0.3, .35, 0.4, .45, .5)
nLOR <- length(LORs)
sims <- 1:1000

Res <- NULL;

for(iLOR in 1:nLOR){

    pvs <- NULL;
    for(isim in sims){
        set.seed(isim)

        dat1<-simPathIndSNPv2(nGenes=20, nGenes1=10, nSNPlow=1, nSNPup=100, nSNP0=1:3,
                            LOR=LORs[iLOR],
                            n=500, MAFslow=0.05, MAFsup=0.4, p0=0.05)

#        p.aspu <- aSPU(dat1$Y, dat1$X, pow=c(1:8, Inf), n.perm=500)

        p.aspu2 <- SPUPtest(dat1$Y, dat1$X, pow=c(1:8, Inf), B=500 )

        p.pathaspu<- aPathSPUPtest(dat1$Y, dat1$X, dat1$nSNPs,
                                   pow=1:8, pow2=c(1, 2, 4, 8),B=500)

        p.pathaspu.pc<- aPathSPUPtest(dat1$Y, dat1$X, dat1$nSNPs,
                                      pow=1:8, pow2=c(1, 2, 4, 8),B=500, usePCs = T)

        p.pathaspu.sg <- PathSPUPtestSingle(dat1$Y, dat1$X, dat1$nSNPs,
                                    pow=1:8, B=500, usePCs = F )

        p.pathaspu.sg.pc <- PathSPUPtestSingle(dat1$Y, dat1$X, dat1$nSNPs,
                                    pow=1:8, B=500, usePCs = T)

        pminP<-UminPP(dat1$Y, dat1$X, B=200)

#        p.grass<-grass(snp.dat=t(dat1$X), snp.info=dat1$snp.info,
#                       gene.info=dat1$gene.info, gene.set=dat1$pathway,
#                       y=dat1$Y, gene.def="rel", B=200)

#        p.plink<-plinkSet(snp.dat=t(dat1$X), snp.info=dat1$snp.info,
#                          gene.info=dat1$gene.info, gene.set=dat1$pathway,
#                          y=dat1$Y, snp.method="logiReg", gene.def="rel",
#                          snp.pcut=0.05, snp.r2cut=0.5)

        pvs <- rbind(pvs, c(
      p.pathaspu$unstdPs[37], p.pathaspu$stdPs[37], p.pathaspu$unstdPsUnnorm[37],
      p.pathaspu.pc$unstdPs[37], p.pathaspu.pc$stdPs[37], p.pathaspu.pc$unstdPsUnnorm[37],
      p.pathaspu.sg$unstdPs[9], p.pathaspu.sg$stdPs[9],
      p.pathaspu.sg.pc$unstdPs[9], p.pathaspu.sg.pc$stdPs[9],
#      p.aspu$pvs[10], p.aspu$pvs[2],
      p.aspu2[10], p.aspu2[2],
      pminP) )

    }
    colnames(pvs) <-
        c("p.pathaspu.unstd","p.pathaspu.std", "p.pathaspu.unstdunnorm",
          "p.pathaspu.unstd.pc","p.pathaspu.std.pc", "p.pathaspu.unstdunnorm.pc",
          "p.pathaspu.unstd.sg","p.pathaspu.std.sg",
          "p.pathaspu.unstd.sg.pc","p.pathaspu.std.sg.pc",
 #         "p.aspu", "p.ssu",
          "p.aspu2", "p.ssu2",
          "p.uminp" )
    Res[[iLOR]] = pvs
}
names(Res) <- LORs



save(Res, file = "out_paspu_D2_b.RData")

