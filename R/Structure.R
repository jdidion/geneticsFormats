# Split a structure file into chromosomes. Optionally, split into smaller files
# by number of SNPs or basepairs. If replace.map is true, the second line of the
# input file will be replaced by distance between markers in cM from the mapfile.
split.structure.file <- function(infile, outdir, ext=".str", 
        mapfile=NULL, replace.map=F, snp.window=NULL, bp.window=NULL, strict=F) {
    anno<-read.table(infile, nrow=2, header=F)
    geno<-read.table(infile, skip=2, header=F)
    samp <- geno[,1:2]
    # TODO: Strict STRUCTURE format requires the first two columns of the header lines to be tabs
    if (ncol(anno) < ncol(geno)) {
        samp<-rbind(data.frame(V1=c("SNPID","SNPID"),V2=c("POP","POP")), samp)
    }
    geno<-geno[,-c(1:2)]
    colnames(geno) <- colnames(anno)
    chrm.ixs<-which(anno[2,] == -1)
    if (replace.map || (!is.null(snp.window) || !is.null(bp.window))) {
        map <- read.table(mapfile)
    }
    for (i in 1:length(chrm.ixs)) {
        beg <- chrm.ixs[i]
        end <- ifelse(i == length(chrm.ixs), ncol(anno), chrm.ixs[i+1]-1)
        print(paste(i,beg,end))
        a <- anno[,beg:end]
        g <- geno[,beg:end]
        m <- map[match(a[1,], map[,2]),]
        chrm <- m[1,1]
        if (chrm == 23) {
            chrm <- "X"
        }
        if (replace.map) {
            d <- diff(m[,3])
            a[2,] <- c(-1, d)
            anno[2,beg:end] <- c(-1, d)
        }
        if (!is.null(snp.window)) {
            for (win in seq(1, nrow(m), snp.window)) {
                aa <- a[,win:(win+snp.window-1)]
                aa[2,1] <- -1
                gg <- g[,win:(win+snp.window-1)]
                fname <- paste("chr", chrm, "_", win, "-", min(win+snp.window-1, nrow(m)), ext, sep="")
                write.table(cbind(samp, rbind(aa, gg)), file.path(outdir, fname), row.names=F, col.names=F, quote=F)
            }
        }
        else if (!is.null(bp.window)) {
            for (win in seq(3000000, max(m[,4]), bp.window)) {
                w <- (m[,4] >= win) & (m[,4] < (win + bp.window))
                aa <- a[,w]
                gg <- g[,w]
                if (sum(w) == 0) {
                    next
                }
                if (sum(w) == 1) {
                    aa <- matrix(aa, nrow=2, ncol=1)
                    gg <- matrix(gg, nrow=nrow(g), ncol=1)
                }
                aa[2,1] <- -1
                fname <- paste("chr", chrm, "_", format(win, scientific=F), "-", 
                    format(min(win+bp.window, max(m[,4])), scientific=F), ext, sep="")
                write.table(cbind(samp, rbind(aa, gg)), file.path(outdir, fname), row.names=F, col.names=F, quote=F)
            }
        }
        else {
            write.table(cbind(samp, rbind(a, g)), file.path(outdir, paste("chr", chrm, ext, sep="")), row.names=F, col.names=F, quote=F)
        }
    }
    if (replace.map) {
        write.table(cbind(samp, rbind(anno, geno)), infile, row.names=F, col.names=F, quote=F)
    }
}

summarize.fs.windows <- function(K, indir, outdir, chromosomes=c(1:19,"X"), bp.size=1000000) {
    results <- NULL
    for (chrm in chromosomes) {
        pos <- 3000000
        while (pos < 200000000) {
            test <- file.path(indir, paste("K", K[1], sep=""),
                paste("chr", chrm, "_", format(pos, scientific=F), "-", format(pos+bp.size, scientific=F), ".", K[1], ".meanP", sep=""))
            if (!file.exists(test)) {
                pos <- pos + bp.size
                next
            }
            fname <- file.path(indir, "K{0}",
                paste("chr", chrm, "_", format(pos, scientific=F), "-", format(pos+bp.size, scientific=F), ".{0}", sep=""))
            out <- file.path(outdir, 
                paste("chr", chrm, "_", format(pos, scientific=F), "-", format(pos+bp.size, scientific=F), sep=""))
            cmd <- paste("/Users/johndidion/Library/Enthought/Canopy_64bit/User/bin/python ~/software/fs/fastStructure/chooseK.py --input=", fname, " --K=", paste(K, collapse="-"), " --out=", out, sep="")
            print(cmd)
            system(cmd)
            temp <- summarize.fs(out, F)
            results <- rbind(results, data.frame(chrm=chrm, start=pos, 
                end=pos+bp.size, complexity=temp[1], components=temp[2]))
            pos <- pos + bp.size
        }
    }
    results
}

plot.fs.windows <- function(s, ...) {
    plot.genome.data(data.frame(chrm=paste("mm",s$chrm,sep=""), x=s$start+500000, y=s$components,pch=20,size=0.3),
        type=c("meanline","line","point"), gps=list(line=gpar(col='darkgray'),point=gpar(col='black',alpha=0.6), meanline=gpar(col="red",lty=2)),
        ylim=c(0,50), yscale=T, ...)
}

summarize.fs <- function(prefix, plot=T, main=NULL) {
    compl<-read.table(paste(prefix,".complexity.csv",sep=""), sep=",")
    compl<-compl[order(compl[,1]),]
    maxml <- which(compl[,2]==max(compl[,2]))
    compo<-read.table(paste(prefix,".components.csv",sep=""), sep=",")
    compo<-compo[order(compo[,1]),]
    tab <- table(compo[,2])
    compo.freq <- min(as.integer(names(tab)[tab==max(tab)]))
    
    if (plot) {
        plot(compl,type="l",xlab="K", ylab="Marginal Likelihood", 
            main=paste("FastStructure Model Complexity,", main))
        points(compl[maxml,], col="red", pch=20)
        segments(compo.freq,0,compo.freq,min(compl[,2]),col="blue")
    }
    
    list(max.complexity=compl[maxml,1], max.components=compo.freq, complexity=compl, components=compo)
}

summarize.fs2 <- function(outdir, fams, K=1:25, plot=T) {
    data <- lapply(fams, function(fam) {
        sum_file <- file.path(outdir, fam, 'k_summary.csv')
        summary <- read.table(sum_file, sep=",", header=T, stringsAsFactors=F)
        summary <- summary[summary$K %in% K,]
        summary <- summary[order(summary$K),]
        summary <- data.frame(fam=fam, summary)
        
        summary$size1 <- 0
        summary$size2 <- 0
    
        maxML <- summary[min(which(summary$ML == max(summary$ML))),"K"]
        summary$size1[summary$K == maxML] <- 2.5
    
        tab <- table(summary$bestK)
        bestK <- as.integer(names(tab)[min(which(tab==max(tab)))])
        summary$size2[summary$bestK==bestK] <- 2.5
        
        list(summary=summary, maxML=maxML, bestK=bestK)   
    })
    
    summary <- do.call(rbind, lapply(data, "[[", 'summary'))
    summary$fam <- factor(summary$fam, levels=fams)
    params <- do.call(rbind, lapply(data, function(x) data.frame(maxML=x$maxML, bestK=x$bestK)))
    rownames(params) <- fams
    
    if (plot)
        my.grid.arrange(list(
            ggplot(summary, aes(x=K, y=ML, colour=fam)) + geom_line(size=1) + 
                geom_point(aes(size=size1), shape=21, fill='black') + 
                scale_size_identity() + 
                theme_bw(),
            ggplot(summary, aes(x=K, y=bestK, colour=fam)) + geom_line(size=1) + 
                geom_point(aes(size=size2), shape=21, fill='black') + 
                scale_size_identity() + 
                theme_bw()
            ), ncol=1
        )
    
    list(summary=summary, params=params)
}

my.grid.arrange <- function(plots, ...) {
    grobs <- lapply(plots, ggplotGrob)
    max.width <- as.list(do.call(grid::unit.pmax, lapply(grobs, function(x) x$widths[2:5])))
    for (i in 1:length(grobs)) {
        grobs[[i]]$widths[2:5] <- max.width
    }
    do.call(grid.arrange, c(grobs, list(...)))
}

summarize.structure <- function(plink.prefix, result.dir, name, k, chromosomes=c(1:19,"X"), sbs=F, locus.af=F) {
    if (length(k) > 1) {
        Ks <- lapply(k, function(kk) summarize.structure(plink.prefix, result.dir, name, kk, chromosomes, sbs, locus.af))
        names(Ks) <- k
        return(Ks)
    }
    
    indivs.file <- paste(plink.prefix, ".fam", sep="")
    indivs <- read.table(indivs.file)[,1:2]
    colnames(indivs) <- c("Pop","Name")
    
    if (sbs) {
        map.file <- paste(plink.prefix, ".bim", sep="")
        map <- read.table(map.file)
        colnames(map) <- c("Chr","Name","cM","Pos")
    }
    
    chrms <- lapply(chromosomes, function(chrm) {
        print(paste("Chromosome", chrm))
        snps <- NULL
        if (sbs) {
            snps <- map[map$Chr == chrm,]
        }
        result.files <- paste(file.path(result.dir, paste("K", k, sep=""), paste(name, k, paste("chr", chrm, sep=""), "out_", sep=".")), c("f", "ss"), sep="")
        summarize.structure.chrm(indivs, snps, result.files[1], result.files[2], k, chrm, locus.af)
    })
    names(chrms) <- chromosomes
    chrms
}

# Summarize site-by-site STRUCTURE results for a single K and chromosome.
# indivs: two-column matrix (FAMID, INDIVID) in same order as in STRUCTURE input files, or path to PLINK .fam file
# snps: four-column matrix (chrm, ID, cM, pos) in same order as in STRUCTURE input files, or path to PLINK .map file
# name: prefix for STRUCTURE output files
# dir: parent directory; sub-directories/files should be of the form K{K}/{name}.{K}.chr{chrm}.out_[f/ss]
# k: K value (number of populations)
# chrm: chromosome to summarize
summarize.structure.chrm <- function(indivs, snps, summary.file, sbs.file, k, chrm, locus.af=F) {
    K.labels <- paste("K", 1:k, sep="")
    
    # open whole-chromosome output file
    con <- file(summary.file, "r")
    # population allele frequency
    skip(con, 35)
    afdiv <- read.table(con, nrows=k, strip.white=T, header=T, row.names=1, na.strings="-", check.names=F)
    colnames(afdiv) <- K.labels
    rownames(afdiv) <- K.labels
    # expected within-population heterozygosity
    skip(con, 2)
    het <- read.table(con, sep=":", nrows=k, strip.white=T, header=F, row.names=1)
    # summary stats
    skip(con, 2)
    stats <- read.table(con, sep="=", nrows=k+6, strip.white=T, header=F, row.names=1)
    # predicted ancestry of each sample
    skip(con, 4)
    anc <- read.table(con, nrows=nrow(indivs), strip.white=T, header=F, row.names=1)[,-4]
    colnames(anc) <- c("Name", "PctMissing", "Pop", K.labels)
    
    af <- NULL
    
    if (!is.null(snps) && locus.af) {
        skip(con, 6)
        af <- do.call(rbind, lapply(1:nrow(snps), function(i) {
            lines <- readLines(con, 3)
            id <- substr(lines[1],regexpr(":", lines[1])+2,nchar(lines[1])-1)
            N <- as.integer(substr(lines[2], 1, 1))
            tab <- read.table(con, nrows=N, strip.white=T, header=F)
            tab <- cbind(data.frame(snp=id, allele=tab[,1], freq=as.numeric(gsub("[\\(\\)]", "", tab[,2], perl=T))), tab[,3:ncol(tab)])
            colnames(tab)[4:ncol(tab)] <- K.labels
            skip(con, 1)
            tab
        }))
    }
    
    close(con)
    
    result <- list(afdiv=afdiv, het=het, stats=stats, anc=anc)
    
    if (!is.null(snps) && !is.null(sbs.file)) {
        # table of joint probabilities
        sbs <- read.table(sbs.file, skip=2, header=F)
        # compress to table of population probabilities
        sbs <- cbind(sbs[,1:2], do.call(cbind, 
            lapply(seq(3, ncol(sbs), k), function(ki) apply(sbs[,ki:(ki+k-1)], 1, sum))))
        colnames(sbs) <- c("Indiv", "Locus", K.labels)
    
        # one site-by-site file per individual
        result$indiv <- by(sbs, sbs[,1], function(x) cbind(snps[x[,2],], x[,3:ncol(x)]))
        # one covariate file per locus
        result$locus <- by(sbs, sbs[,2], function(x) cbind(indivs[x[,1],], x[,3:ncol(x)]))
    }
    
    result
}

skip <- function(con, n) return(invisible(readLines(con, n, warn=F)))

# Use the Delta K method to choose the best K.
# LnProb is an Nx2 matrix, with one row per K, each row being the
# LnProbMean,LnProbStdev, rownames = K values.
deltaK <- function(LnProb) {
  LnProb <- LnProb[order(as.integer(rownames(LnProb))),]
  x <- 2:(nrow(LnProb)-1)
  dK <- sapply(x, function(i)
      abs(LnProb[i+1,1] - 2.0 * LnProb[i,1] + LnProb[i-1,1]) / LnProb[i,2])
  names(dK) <- rownames(LnProb)[x]
  dK
}

plot.probs <- function(Ks, outfile=NULL) {
    mat1 <- do.call(rbind, lapply(Ks, function(k) unlist(lapply(k, function(chrm) chrm$stats[1,1]))))
    mat2 <- mat1
    for (i in 1:nrow(mat2)) {
        for (j in 1:ncol(mat2)) {
            mat2[i,j] <- exp(mat1[i,j]) / sum(exp(mat1[i,-j]))
        }
    }
    sums1 <- data.frame(K=as.integer(rownames(mat1)), LnProb=apply(mat1, 1, sum))
    sums2 <- data.frame(K=as.integer(rownames(mat2)), LnProb=apply(mat2, 1, sum))
    if (!is.null(outfile)) jpeg(outfile, width=7, height=7, units="in", res=300)
    plot(sums2, type="l") 
    if (!is.null(outfile)) dev.off()
    mat2
}

plot.chrm.structure <- function(chrm.result, sample.names, outfile, k=10) {
    n <- length(chrm.result)
    jpeg(outfile, width=11, height=n+2, units="in", res=300)
    par(mfrow=c(n,1), oma=c(3,0,3,0), mar=c(0,0,0.8,0))
    cols <- brewer.pal(k,"Set3")
    for (i in 1:n) {
        data <- t(chrm.result[[i]][,5:(4+k)])
        barplot(data, space=0, col=cols, border=NA, xlab="", ylab="", xaxt="n", yaxt="n", main=sample.names[i], cex.main=0.8)
        if (i == 1) {
            legend(ncol(data)/2, 1.5, paste("K",1:k,sep=""), fill=brewer.pal(10,"Set3"), xpd=NA, horiz=T, xjust=0.5)
        }
        else if (i == n) {
            axis(1)
        }
    }
    dev.off()
}

write.cov.files <- function(result, outdir, exclude.indiv=NULL) {
    j <- 1
    for (chrm in 1:20) {
        if (chrm==20) {
            ch <- "X"
        }
        else {
            ch <- as.character(chrm)
        }
        print(ch)
        r <- result[[chrm]]$locus
        s <- result[[chrm]]$indiv[[1]]$Name
        for (i in 1:length(s)) {
            outfile <- file.path(outdir, paste("covariates_", sprintf("%06d", j), ".txt", sep=""))
            rr <- r[[i]]
            rr <- cbind(rr[,1:2], rep(1, nrow(rr)), rr[,3:ncol(rr)])
            if (!is.null(exclude.indiv)) {
                rr <- rr[-exclude.indiv,]
            }
            write.table(rr, outfile, row.names=F, col.names=F, quote=F)
            j <- j + 1
        }
    }
    print(j-1)
}
