bim.to.map <- function(bim, outfile) {
    if (is.character(bim)) {
        bim <- read.table(bim, sep="\t", header=F, stringsAsFactors=F)
        for (chrm in unique(bim[,1])) {
            b <- bim[bim[,1]==chrm,]
            ppos <- b[,4]
            gpos <- b[,3]
            
            pdiff <- diff(ppos) / 1000000
            gdiff <- diff(gpos)
            rrate <- gdiff / pdiff
            # the rrate vector is one shoreter than p/gpos vectors, so we just
            # duplicated the rrate of the first marker
            rrate <- c(rrate[1], rrate)
            
            write.table(data.frame(pposition=ppos, rrate=rrate, gposition=gpos),
                paste0(outfile, chrm, '.txt'), sep=" ", row.names=F, col.names=T, quote=F)
        })
    }
}