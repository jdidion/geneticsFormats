trim <- function (x) gsub("^\\s+|\\s+$", "", x)

read.assignments <- function(path, pos) {
    lapply(readLines(path), function(line) {
        x <- c("1", unlist(strsplit(trim(line), "[: ]", perl=TRUE)))
        starts <- pos[ as.integer(x[seq(1, length(x) - 2, 2)]) ]
        ends <- pos[ as.integer(x[seq(3, length(x), 2)]) ]
        anc <- x[seq(2, length(x) - 1, 2)]
        anc1 <- names(pops)[ as.integer(substr(anc, 1, 1)) + 1 ]
        anc2 <- names(pops)[ as.integer(substr(anc, 2, 2)) + 1 ]
        cbind(starts, ends, anc1, anc2)
    })
}

load.assignments <- function(lamp.dir, samples, pops, chromosomes=c(1:19, "X")) {
    assignments <- list()
    
    for (chrm in chromosomes) {
        chrm.dir <- paste(lamp.dir, "/chr", chrm, "/LAMP/", sep="")
        assn.file <- paste(chrm.dir, "chr", chrm, ".50.15.out", sep="")
        pos.file <- paste(chrm.dir, "chr", chrm, ".pos", sep="")
        
        pos <- as.integer(readLines(pos.file, warn=FALSE))
        assn <- read.assignments(assn.file, pos)

        n <- length(assn)
        if (n == 0) {
            stop(paste("No assignments for chromosome", chrm))
        }
        
        for (i in 1:n) {
            name <- samples[i]
            assignments[[name]] <- rbind(assignments[[name]], data.frame(
                chrm=rep(paste("mm", chrm, sep=""), nrow(assn[[i]])), 
                x1=as.integer(assn[[i]][,1]), x2=as.integer(assn[[i]][,2]),
                y=rep(1, nrow(assn[[i]])), anc1=assn[[i]][,3], anc2=assn[[i]][,4],
                stringsAsFactors=FALSE))
        }
    }
    
    summary <- do.call(cbind, lapply(assignments, function(x) {
        size <- x$x2 - x$x1 + 1
        tab <- tapply(c(size, size), c(x$anc1, x$anc2), sum)
        tot <- sum(tab)
        (tab / tot)[names(pops)]
    }))
    summary[is.na(summary)] <- 0
    
    list(assignments=assignments, summary=summary)
}

plot.lamp.genome <- function(assignments, pops, legend=TRUE, ...) {
    data <- list()
    for (name in names(assignments)) {
        assn <- assignments[[name]]
        data[[paste(name, "1")]] <- list(data=cbind(assn[,1:4], fill=unlist(pops[assn[,5]]), stringsAsFactors=FALSE))
        data[[paste(name, "2")]] <- list(data=cbind(assn[,1:4], fill=unlist(pops[assn[,6]]), stringsAsFactors=FALSE))
    }
    if (legend) {
        legend <- data.frame(col=unlist(pops), row.names=names(pops), stringsAsFactors=FALSE)
    }
    else {
        legend <- NULL
    }
    plot.genome.data(data, type="blocks", ylim=c(0, 1), legend=legend, ...)
}

plot.lamp.indivs <- function(assignments, output.dir, ...) {
    for (name in names(assignments)) {
        assn <- assignments[[name]]
        data <- list(
            list(data=cbind(assn[,1:4], fill=unlist(pops[assn[,5]]), stringsAsFactors=FALSE)),
            list(data=cbind(assn[,1:4], fill=unlist(pops[assn[,6]]), stringsAsFactors=FALSE)))
        fname <- gsub("[/ ]", "_", name)
        plot.genome.data(data, type="blocks", ylim=c(0, 1), title=name,
            pdf.file=paste(output.dir, "/", fname, ".pdf", sep=""), ...)
    }
}

plot.lamp.summary <- function(summary, pops, output.dir, ...) {
    pdf(paste(output.dir, "/pop_summary.pdf", sep=""))
    par(mar=c(5, 5, 2, 2))
    mp <- barplot(summary, col=unlist(pops[rownames(summary)]), border=NA, yaxt="n", horiz=TRUE)
    axis(2, mp, colnames(summary), las=1, tick=FALSE, cex.axis=0.7)
    #legend(0, ncol(summary)+6, legend=rownames(summary), fill=unlist(pops[rownames(summary)]), horiz=TRUE, cex=0.7, xpd=NA)
    dev.off()
}

assn.to.tract.length <- function(assn, col) {
    tab <- do.call(rbind, by(assn, assn$chrm, function(y) {
        w <- which(diff(as.integer(factor(y$anc1)))!=0)
        if (length(w) == 0) {
            c(y[1, col], y[nrow(y), 3] - y[1, 2] + 1)
        }
        else {
            w <- c(0, w, nrow(y))
            t(sapply(2:length(w), function(i) c(y[w[i], col], y[w[i], 3] - y[w[i-1] + 1, 2] + 1)))
        }
    }))
    data.frame(ss=tab[,1], size=as.integer(tab[,2]))
}

plot.tract.lengths <- function(assn, pops, window.size=5000000, pdf.file=NULL) {
    tab <- do.call(rbind, lapply(assn, function(x) {
        rbind(assn.to.tract.length(x, 5), assn.to.tract.length(x, 6))
    }))
    print(tab)
    pop.names <- unique(tab[,1])
    breaks<-seq(0,195000000, window.size)
    counts <- do.call(rbind, lapply(pop.names, function(p) {
        hist(tab[tab[,1] == p, "size"], breaks=breaks, plot=FALSE)$counts
    }))
    colnames(counts) <- breaks[-1] / 1000000
    if (!is.null(pdf.file)) {
        pdf(pdf.file, pdf.width, pdf.height)
        on.exit(dev.off())
    }
    barplot(counts, col=unlist(pops[as.character(pop.names)]), xlab="Tract size (Mb)")
}

# Plot the results of analyzing a single chromosome (or sub-chromosomal region) using LAMP.
# * assn_file = LAMP-LD output file
# * pos_file = pos file supplied to LAMP
# * pops = list of pop_name=color, in the order that the population haplotype files were
#   supplied to LAMP
# * samples = list of sample names in the order that they appear in assn_file
plot.lamp.chrm <- function(assn.file, pos.file, pops, samples, title=NULL, pdf.file=NULL, pdf.width=10, pdf.height=10) {
    pos <- as.integer(readLines(pos_file))
    assns <- read.assignments(assn.file, pos)
    
    n <- length(assns)
    if (n == 0) {
        return(NULL)
    }
    
    names(assns) <- samples
    
    if (!is.null(pdf.file)) {
        pdf(pdf.file, pdf.width, pdf.height)
        on.exit(dev.off())
    }
    
    hsize <- c(5, 1, 1)
    vsize <- c(0, 1, 5)
    if (!is.null(title))
        vsize[1] <- 2

    pushViewport(viewport(layout=grid.layout(3, 3, 
        heights=unit(vsize, c('lines', 'null', 'lines')),
        widths=unit(hsize, c('lines', 'null', 'lines')))))
    if (!is.null(title))
        grid.text(title, y=unit(convertY(unit(1, 'npc'), 'lines', TRUE) - 1, 'lines'))
    grid.text(paste("Chromosome Position (Mb)"), y=unit(3, 'lines'))

    # draw sample labels
    xmin <- pos[1]
    xmax <- pos[length(pos)]
    pushViewport(viewport(layout=grid.layout(n+1, 1), xscale=c(0, 1), layout.pos.row=2, layout.pos.col=1, name='data'))
    for (i in 1:n) {
        pushViewport(viewport(layout=grid.layout(1, 1, heights=1), 
            layout.pos.row=i, layout.pos.col=1, xscale=c(0, 1), yscale=c(0, 1)))
        name <- samples[i]
        grid.text(name, y=0.25, gp=gpar(cex=0.7))
        upViewport()
    }
    upViewport()

    # draw sample data
    pushViewport(viewport(layout=grid.layout(n+1, 1), xscale=c(xmin, xmax), layout.pos.row=2, layout.pos.col=2, name='data'))
    for (i in 1:n) {
        pushViewport(viewport(layout=grid.layout(1, 1, heights=1), 
            layout.pos.row=i, layout.pos.col=1, xscale=c(xmin, xmax), yscale=c(0, 1)))
        data <- assns[[i]]
        for (j in 1:nrow(data)) {
            start <- as.integer(data[j, 1])
            end <- as.integer(data[j, 2])
            col1 <- pops[[data[j, 3]]]
            col2 <- pops[[data[j, 4]]]
            grid.rect(
                x=unit(start, "native"), width=unit(end - start, "native"), 
                y=unit(0.05, "native"), height=unit(0.45, "native"),
                just='left', gp=gpar(col=NA, fill=col1))
            grid.rect(
                x=unit(start, "native"), width=unit(end - start, "native"), 
                y=unit(0.5, "native"), height=unit(0.45, "native"),
                just='left', gp=gpar(col=NA, fill=col2))
        }
        upViewport()
    }
    
    # draw axis
    pushViewport(viewport(layout.pos.row=n+1, layout.pos.col=1, xscale=c(xmin, xmax)))
    xticks <- seq(xmin, xmax, length.out=8)
    grid.xaxis(xticks, as.character(round(xticks / 1000000)), gp=gpar(fontsize=9))
    upViewport()
    
    upViewport()
    
    # draw legend
    midpts <- seq(0.15, 0.8, length.out=length(pops))
    for (i in 1:length(pops)) {
        grid.rect(x=unit(midpts[i] - 0.01, "npc"), width=unit(0.02, "npc"),
                  y=unit(1, "lines"), height=unit(0.02, "npc"),
                  gp=gpar(col=NA, fill=pops[[i]]))
        grid.text(names(pops)[i], x=unit(midpts[i] + 0.05, "npc"), y=unit(1, "lines"))
    }
    
    assns
}

summarize.chrm <- function(indir, outdir, chrm, pops, samples, leg, windows=c(25,50,100,300), founders=c(15,150)) {
    indir <- paste(indir, "/chr", chrm, "/LAMP", sep="")
    pos <- paste(indir, "/chr", chrm, ".pos", sep="")
    result.mixed <- list()
    result.single <- list()
    ancs <- NULL
    
    for (w in windows) {
        for (f in founders) {
            infile <- paste(indir, "/chr", chrm, ".", w, ".", f, ".out", sep="")
            print(infile)
            if (file.exists(infile)) {
                title <- paste("LAMP Ancestry Assignment, Chrm. ", chrm, ", ", w, " SNP Windows, ", f, " Founders", sep="")
                outfile <- paste(outdir, "/lamp_", w, "_", f, "_chr", chrm, ".pdf", sep="")
                assns <- plot.lamp.chrm(infile, pos, pops, samples, title, outfile)
                if (!is.null(assns)) {
                    result.mixed[[paste(chrm, w, f, sep="_")]] <- lapply(assns, function(x) {
                        anc <- apply(force.matrix(x[,3:4], 1), 1, function(y) paste(sort(y), collapse="_"))
                        ancs <<- c(ancs, anc)
                        tot <- as.integer(x[nrow(x),2]) - as.integer(x[1,1]) + 1
                        sapply(unique(anc), function(a) {
                            w <- anc == a
                            sum(as.integer(x[w,2]) - as.integer(x[w,1]) + 1) / tot
                        })
                    })
                    result.single[[paste(chrm, w, f, sep="_")]] <- lapply(assns, function(x) {
                        tot <- (as.integer(x[nrow(x),2]) - as.integer(x[1,1]) + 1) * 2
                        sapply(names(pops), function(p) {
                            w1 <- x[,3] == p
                            s1 <- sum(as.integer(x[w1,2]) - as.integer(x[w1,1]) + 1)
                            w2 <- x[,4] == p
                            s2 <- sum(as.integer(x[w2,2]) - as.integer(x[w2,1]) + 1)
                            (s1 + s2) / tot
                        })
                    })
                }
            }
        }
    }

    # plot summary of all ancestry combintations
    n <- unique(ancs)
    s <- sapply(samples, function(s) {
        x <- rep(0, length(n))
        names(x) <- n
        for (r in result.mixed) {
            x[ names(r[[s]]) ] <- x[ names(r[[s]]) ] + r[[s]]
        }
        x <- x / length(result.mixed)
        x
    })
    m <- match(names(leg), rownames(s))
    s <- s[m[!is.na(m)],]
    pdf(paste(outdir, "/lamp_chr", chrm, "_anc_summary.pdf", sep=""))
    mp <- barplot(s, col=leg[rownames(s)], border=NA, xaxt="n")
    axis(1, mp, colnames(s), las=3, tick=FALSE, cex.axis=0.7)
    legend(1, 1.15, legend=rownames(s), fill=leg[rownames(s)], ncol=4, cex=0.7, xpd=NA)
    dev.off()
    
    # plot summary of each population
    s <- sapply(samples, function(s) {
        x <- rep(0, length(pops))
        names(x) <- names(pops)
        for (r in result.single) {
            x[ names(r[[s]]) ] <- x[ names(r[[s]]) ] + r[[s]]
        }
        x <- x / length(result.single)
        x
    })
    m <- match(names(pops), rownames(s))
    s <- s[m[!is.na(m)],]
    pdf(paste(outdir, "/lamp_chr", chrm, "_pop_summary.pdf", sep=""))
    par(mar=c(5, 5, 2, 2))
    mp <- barplot(s, col=unlist(pops[rownames(s)]), border=NA, yaxt="n", horiz=TRUE)
    axis(2, mp, colnames(s), las=1, tick=FALSE, cex.axis=0.7)
    legend(0, ncol(s)+6, legend=rownames(s), fill=unlist(pops[rownames(s)]), horiz=TRUE, cex=0.7, xpd=NA)
    dev.off()
}
