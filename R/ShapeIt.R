#' Convert a bim file to the genetic map format used by ShapeIt.
#' Recombiantion rate is estimated as (cM1 - cM0) / ((pos1 - pos0) / 1000000)
bim.to.map <- function(bim, outfile, chrms=NULL) {
    if (is.character(bim)) {
        bim <- read.table(bim, sep="\t", header=F, stringsAsFactors=F)
    }
    if (is.null(chrms)) {
        chrms <- unique(bim[,1])
    }
    scipen <- getOption('scipen')
    options(scipen=10)
    for (chrm in chrms) {
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
    }
    options(scipen=scipen)
}

merge.shapeit <- function(sample.files, hap.files, sample.out, hap.out) {
    stopifnot(length(sample.files)==length(hap.files))
    sample.header <- read.table(sample.files[1], nrow=2, header=F, stringsAsFactors = F)
    samples <- do.call(rbind, lapply(sample.files, read.table, skip=2, header=F, stringsAsFactors=F))
    haps <- read.table(hap.files[1], header=F, stringsAsFactors=F)
    for (i in 2:length(hap.files)) {
        h <- read.table(hap.files[i], header=F, stringsAsFactors=F)
        stopifnot(all(haps[,2] == h[,2]))
        haps <- cbind(haps, h[,-c(1:5)])
    }
    write.table(sample.header, sample.out, row.names=F, col.names=F, quote=F)
    write.table(samples, sample.out, row.names=F, col.names=F, quote=F, append = T)
    write.table(haps, hap.out, col.names=F, row.names=F, quote=F)
}