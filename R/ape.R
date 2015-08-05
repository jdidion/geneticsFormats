# This function is equivalent to running SEQBOOT + DNADIST. It creates B different samplings of
# a data set and calls dist.dna on each. seq must be a matrix with rows = samples and cols = SNPs,
# and cells must be A/C/G/T or an ambiguous character, and must be lower case. The output is a list
# of distance matrices (of class dist). This can be written to a PHYLIP-compatible file using
# write.dist. If use.unsampled=TRUE, one of the distance matrices will be computed using seq
# directly. NOTE: dist.dna can use upwards of 6 GB for a large matrix (150 x 650k).
multidist <- function(seq, B=100, use.unsampled=FALSE, cores=3, model="K80") {
    if (use.unsampled) {
        my.dist(seq, FALSE)
    }
    mclapply(1:100, my.dist.wrapper, seq, model, mc.cores=cores, mc.cleanup=TRUE)
}

my.dist.wrapper <- function(i, seq, model) {
    print(paste("Bootstrap", i, "\n"))
    my.dist(seq, model)
}

my.dist <- function(seq, model="K80", sample=TRUE) {
    if (sample) {
        snps <- ncol(seq)
        cols <- sample(snps, replace=TRUE)
        seq <- as.DNAbin(seq[,cols])
    }
    else {
        seq <- as.DNAbin(seq)
    }
    dist.dna(seq, model)
}

# Convert a matrix of snps x samples to a DNAbin. snps is a two-column matrix of A,B alleles.
geno.as.DNAbin <- function(snps, geno) {
    as.DNAbin(t(tolower(g.n.het.as.ambig(snps, geno))))
}

# Writes a list of distance matricies to a file in PHYLIP format.
write.dists <- function(d, path, rename.file=NULL) {
    f <- file(path, 'w')
    new.names <- NULL

    for (i in 1:length(d)) {
        dd <- d[[i]]
        if (!is.matrix(dd)) {
            dd <- as.matrix(dd)
        }
        if (!is.null(rename.file)) {
            if (is.null(new.names)) {
                full.names <- rownames(dd)
                new.names <- as.character(1:nrow(dd))
                write(paste(new.names, full.names, sep='='), rename.file, 1)
            }
            rownames(dd) <- new.names
        }
        write.dist(dd, f)
    }
}

write.dist <- function(d, f) {
    names <- sprintf('%10-s', substr(rownames(d), 1, 10))
    writeln(f, nrow(d))
    for (i in 1:nrow(d)) {
        writeln(f, paste(names[i], paste(formatC(round(d[i,], 6), digits=6, format='f'), collapse=' '), sep=' '))
    }
}

relabel.tree <- function(infile, outfile, names.file, sd=0) {
    p <- read.tree(infile)
    p$tip.label <- read.table(names.file, sep='=', header=FALSE, colClasses='character', row.names=1)[p$tip.label, 1]
    if (!is.na(sd)) {
        p$edge.length <- round(p$edge.length, sd)
    }
    write.tree(p, outfile)
}

reroot.tree <- function(infile, new.root, outfile=infile) {
    write.tree(root(read.tree(infile), new.root), outfile, resolve.root=T)
}
