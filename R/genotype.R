# Sequence conversion functions for SNPs in generic table format. Functions are named according
# to the input and output formats. So for example, g.n.as.i converts genotypes in nucleotide 
# format to genotypes in integer format. The formats are:
# 1. Genotype (g): A vector of genotypes. Genotypes may be encoded as nucleotides (n),
# alleles (a) or integers (i). Numeric genotypes are coded A=1, H=2, B=3, V=4, N=5.
# 2. Haplotype (h): A pair of vectors (typically a list, but may also be rows or columns in a
# matrix). Again, characters may be nucleotides or integers, but H and V are not allowed.
# Additionally, these formats may be converted to VINO (v) format, which is binary (non-VINO calls
# are converted to one number/character and VINO calls are converted to another). Methods that
# convert to v format *cannot* be passed to the general-purpose convert() function.

# Loads a tab-delimited file with N columns, where N > 2 (sample, snp, ...),
# and returns a list of N-2 data frames. Each data frame is MxN (M=samples, N=snps).
load.long.format <- function(path, snpIDs, sep="\t", col.names=c("sample","snp","call","intensity"), fill=list("N",0)) {
    header <- (col.names[1] == TRUE)
    data <- read.table(path, sep=sep, header=header, stringsAsFactors=FALSE, 
        check.names=FALSE, comment.char="")
    n <- ncol(data)
    if (!header && length(col.names) == n) {
        colnames(data) <- col.names
    }
    result <- list()
    for (i in 3:n) {
        name <- colnames(data)[i]
        result[[name]] <- first.col.as.names(
            dcast(data, snp~sample, fill=fill[[i-2]], value.var=name), snpIDs)
    }
    result
}

first.col.as.names <- function(df, snpIDs) {
    data.frame(df[,-1], row.names=df[,1], stringsAsFactors=FALSE, check.names=FALSE)[snpIDs,]
}

# Given a three-column table of snp,genotype,count (where each snp may appear more than onece),
# create a consensus table of snp,genotype,freq where each snp appears exactly once and freq is
# the fraction of calls that were the consensus genotype. If two or more genotypes have the same
# frequency, the consensus is N and the value of freq is NA.
make.consensus <- function(counts) {
    cons <- do.call(rbind, by(counts, counts[,1], function(x) {
        if (nrow(x) == 1) {
            geno <- x[1,2]
            freq <- 1.0
        }
        else {
            x <- x[order(x[,3], decreasing=TRUE),]
            if (x[1,3] == x[2,3]) {
                geno <- "N"
                freq <- NA
            }
            else {
                geno <- x[1,2]
                freq <- x[1,3] / sum(x[,3])
            }
        }
        list(geno=geno, freq=freq, snp=x[1,1])
    }))
    data.frame(geno=unlist(cons[,1]), freq=unlist(cons[,2]), row.names=unlist(cons[,3]), stringsAsFactors=FALSE)
}

# Get H and N values for a set of SNPs (margin=1) or samples (margin=2).
as.HN <- function(geno, margin=2) {
    t(apply(geno, margin, function(x) {
        N <- sum(x == "N")
        l <- length(x)
        H <- ifelse(N == l, 0, sum(x == "H") / (l - N))
        c(H, N / l)
    }))
}

### Conversion of het calls to ambiguous characters ###

ambig <- c(
    "A", # 1 = AA
    "C", # 2 = CC
    "M", # 3 = AC
    "G", # 4 = GG
    "R", # 5 = AG
    "S", # 6 = CG
    "V", # 7 = ACG
    "T", # 8 = TT
    "W", # 9 = AT
    "Y", # 10 = CT
    "H", # 11 = ACT
    "K", # 12 = GT
    "D", # 13 = AGT
    "B", # 14 = CGT
    "N"  # 15 = ACGT
    )

# Given an Nx2 matrix (rows=SNPs, cols=ref and var alleles), return a vector
# of ambiguous characters.
n.as.ambig <- function(a) {
    if (!is.matrix(a)) a <- as.matrix(a)
    m <- matrix(15, nrow(a), ncol(a))
    m[a=="A"] <- 1
    m[a=="C"] <- 2
    m[a=="G"] <- 4
    m[a=="T"] <- 8
    ambig[bitOr(m[,1], m[,2])]
}

# Converts het calls in genotypes to ambiguous characters.
g.n.het.as.ambig <- function(a, g, H="H") {
    h <- n.as.ambig(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    for (i in 1:nrow(g)) {
        g[i,g[i,]==H] <- h[i]
    }
    g
}

### Converters to VINO format ###

# these SNP frequencies are based on Sanger data
snp.freq <- matrix(c(
    "A","C",0.09,
    "A","G",0.33,
    "A","T",0.09,
    "C","G",0.07,
    "C","T",0.33,
    "G","T",0.09),
    ncol=3, byrow=TRUE)

# create a two column matrix (Allele A, Allele B) based on observed SNP frequencies.
simulate.snps <- function(n) {
    cs <- c(0, cumsum(snp.freq[,3]))
    ns <- runif(n)
    ns[ns == 1.0] <- 0.99
    m <- matrix(NA, nrow=n, ncol=2)
    for (i in 2:length(cs)) {
        w <- ns >= cs[i-1] & ns < cs[i]
        if (any(w)) {
            m[w,1] <- snp.freq[i-1,1]
            m[w,2] <- snp.freq[i-1,2]
        }
    }
    m
}

# Returns a matrix of VxN, where V is the set of rows with any V calls. The resulting genotypes are
# trinary: A,H,B are converted to one genotype, V to another, and N's remain the same.
g.i.as.v.n <- function(g, N="N", id.suffix="V")  .g.i.as.v(g, N, id.suffix, FALSE)
g.i.as.v.a <- function(g, N="NN", id.suffix="V") .g.i.as.v(g, N, id.suffix, TRUE)

.g.i.as.v <- function(g, N, id.suffix, diploid) {
    if (!is.matrix(g)) g <- as.matrix(g)
    w <- apply(g, 1, function(x) any(x==4))
    if (!any(w)) {
        return(NULL)
    }
    g <- g[w,,drop=FALSE]
    rownames(g) <- paste(rownames(g)[w], id.suffix, sep="")
    n <- nrow(g)
    a <- simulate.snps(n)
    if (diploid) {
        a <- cbind(paste(a[,1], a[,1], sep=""), paste(a[,2], a[,2], sep=""))
    }
    dimnames(a)<-list(rownames(g),c("A","B"))
    m <- matrix(N, nrow=n, ncol=ncol(g), dimnames=dimnames(g))
    for (i in 1:n) {
        m[i, g[i,] <= 3] <- a[i,1]
        m[i, g[i,] == 4] <- a[i,2]
    }
    list(anno=a, vino=m)
}

g.n.as.v.n <- function(g, V.in="V", N.in="N", N.out="N", id.suffix="V") {
    if (!is.matrix(g)) g <- as.matrix(g)
    w <- apply(g, 1, function(x) any(x==V.in))
    if (!any(w)) {
        return(NULL)
    }
    else {
        ids <- paste(rownames(g)[w], id.suffix, sep="")
        g <- g[w,,drop=FALSE]
    }
    d <- dim(g)
    a <- simulate.snps(d[1])
    colnames(a) <- c("A","B")
    m <- matrix(rep(a[,1], d[2]), nrow=d[1], ncol=d[2], dimnames=list(ids, colnames(g)))
    m[g == N.in] <- N.out
    for (i in 1:d[1]) {
        m[i, g[i,] == V.in] <- a[i,2]
    }
    list(anno=a, vino=m)
}

# Creates synthetic markers for all markers with any V calls. The resulting genotypes are
# trinary: A,B are converted to one genotype, V to another, and N's remain the same.
# Returns a Mx2N matrix, where each pair of columns is a set of haploid genotypes.
g.n.as.v.h.matrix <- function(g, ...) {
    v <- g.n.as.v.n(g, ...)
    m <- v$vino
    d <- dim(m)
    m2 <- matrix(NA, nrow=d[1], ncol=d[2] * 2, dimnames=list(rownames(m),
        as.vector(rbind(paste(colnames(m), "1", sep="."), paste(colnames(m), "2", sep=".")))))
    m2[,seq(1,ncol(m2),2)] <- m
    m2[,seq(2,ncol(m2),2)] <- m
    list(anno=v$a, vino=m2)
}

### Other conversion functions ###

# General-purpose conversion function that takes a matrix with both annotation and sequence
# columns, breaks it into separate annotation and sequence matricies, calls the specified converter,
# and returns one of the following:
# A) a matrix of identical dimensions as the input matrix
# B) for converters that produce haplotype matricies, a matrix with twice as many sample columns
# C) for converters that produce list of sample haplotypes, the list is returned directly.
convert <- function(m, fn, A.col=1, B.col=2, first.geno.col=3, ...) {
    a <- m[,c(A.col, B.col)]
    g <- m[,first.geno.col:ncol(m)]
    newg <- fn(a, g, ...)
    if (is.list(newg)) {
        newg
    }
    else {
        data.frame(m[1:(first.geno.col-1)], newg)
    }
}

# Convert genotypes in nucleotide format to numeric format.
g.n.as.i <- function(a, g, H="H", V="V", V.out=4, N.out=5) {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    m <- matrix(N.out, nrow(g), ncol(g), dimnames=dimnames(g))
    for (i in 1:nrow(m)) {
        m[i, g[i,]==a[i,1]] <- 1
        m[i, g[i,]==H] <- 2
        m[i, g[i,]==a[i,2]] <- 3
        m[i, g[i,]==V] <- V.out
    }
    mode(m) <- "integer"
    m
}

# Convert genotypes in nucleotide format to allele format. If sort=TRUE, alleles for heterozygous
# calls will always appear in (alphabetically) sorted order, otherwise the reference allele will
# be first followed by the variant alle.
g.n.as.a <- function(a, g, H.in="H", V.in="V", N.in="N", V.out="VV", N.out="NN", sort=FALSE) {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    H <- paste(H.in,H.in,sep="")
    V <- paste(V.in,V.in,sep="")
    N <- paste(N.in,N.in,sep="")
    m <- t(apply(g, 1, function(x) paste(x, x, sep="")))
    dimnames(m) <- dimnames(g)
    if (sort) {
        h <- apply(a, 1, function(x) paste(sort(x), collapse=""))
    }
    else {
        h <- apply(a, 1, paste, collapse="")
    }
    for (i in 1:nrow(m)) {
        m[i, m[i,]==H] <- h[i]
        m[i, m[i,]==V] <- V.out
        m[i, m[i,]==N] <- N.out
    }
    m
}

# Convert a MxN matrix of genotypes to a Mx2N matrix of alleles (every pair of
# columns is a set of haploid (unphased) genotypes).
g.n.as.h.matrix <- function(a, g, H="H", V="V", N.in="N", N.out="N") {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    g[g==V] <- N.out
    if (N.in != N.out) {
        g[g==N.in] <- N.out
    }
    g1 <- g
    g2 <- g
    for (i in 1:nrow(g)) {
        w <- g[i,]==H
        g1[i, w] <- a[i,1]
        g2[i, w] <- a[i,2]
    }

    m <- matrix(N.out, nrow(g), ncol(g) * 2)
    rownames(m) <- rownames(g)
    colnames(m) <- as.vector(rbind(paste(colnames(g), "1", sep="."), paste(colnames(g), "2", sep=".")))
    m[,seq(1, ncol(m), 2)] <- g1
    m[,seq(2, ncol(m), 2)] <- g2
    m
}

# Convert a MxN matrix of numeric genotypes to a Mx2N matrix of alleles (every pair of
# columns is a set of haploid (unphased) genotypes). VINOs are converted to N's.
g.i.as.h.matrix <- function(a, g, N.out="N") {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    
    m1 <- matrix(N.out, nrow(g), ncol(g))
    m2 <- matrix(N.out, nrow(g), ncol(g))
    
    for (i in 1:nrow(g)) {
        m1[i, g[i,]<=2] <- a[i,1]
        m1[i, g[i,]==3] <- a[i,2]
        m2[i, g[i,]==1] <- a[i,1]
        m2[i, g[i,]==2 | g[i,]==3] <- a[i, 2]
    }

    m <- matrix(N.out, nrow(g), ncol(g) * 2)
    rownames(m) <- rownames(g)
    colnames(m) <- as.vector(rbind(paste(colnames(g), "1", sep="."), paste(colnames(g), "2", sep=".")))
    m[,seq(1, ncol(m), 2)] <- m1
    m[,seq(2, ncol(m), 2)] <- m2
    m
}

# Returns a list of two-element vectors (one for each sample), where each element is a haplotype.
# For heterozygous sites, the first haplotype has the reference allele and the second has the
# variant allele.
g.n.as.h.list <- function(a, g, H.in="H", V.in="V", N.in="N", V.out=NA, N.out=NA) {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    if (!is.na(V.out)) g[g==V.in] <- V.out
    if (!is.na(N.out)) g[g==N.in] <- N.out
    gg1 <- g
    gg2 <- g
    for (i in 1:nrow(g)) {
        w <- g[i,]==H.in
        gg1[i, w] <- a[i,1]
        gg2[i, w] <- a[i,2]
    }
    sapply(colnames(gg1), function(n) {
        c(paste(gg1[,n], collapse=""), paste(gg2[,n], collapse=""))
    }, simplify=FALSE, USE.NAMES=TRUE)
}

# Convert genotypes in integer format to nucleotide format. If H=Null, ambiguous characters
# will be used for heterozygous calls.
g.i.as.n <- function(a, g, H="H", V="V", N="N") {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (is.null(H)) {
        H <- n.as.ambig(a)
    }
    else if (!is.na(H)) {
        H <- rep(H, nrow(g))
    }
    m <- matrix(N, nrow(g), ncol(g), dimnames=dimnames(g))
    for (i in 1:nrow(g)) {
        m[i,g[i,]==1] <- a[i,1]
        m[i,g[i,]==2] <- H[i]
        m[i,g[i,]==3] <- a[i,2]
        m[i,g[i,]==4] <- V
    }
    m
}

# Convert genotypes in integer format to allelic format. If sort=TRUE, alleles for heterozygous
# calls will always appear in (alphabetically) sorted order, otherwise the reference allele will
# be first followed by the variant alle.
g.i.as.a <- function(a, g, V="VV", N="NN", sort=FALSE) {
    if (!is.matrix(a)) a <- as.matrix(a)
    if (!is.matrix(g)) g <- as.matrix(g)
    m <- matrix(N, nrow(g), ncol(g), dimnames=dimnames(g))
    A <- paste(a[,1], a[,1], sep="")
    B <- paste(a[,2], a[,2], sep="")
    if (sort) {
        H <- apply(a, 1, function(x) paste(sort(x), collapse=""))
    }
    else {
        H <- apply(a, 1, paste, collapse="")
    }
    for (i in 1:nrow(g)) {
        m[i, g[i,]==1] <- A[i]
        m[i, g[i,]==2] <- H[i]
        m[i, g[i,]==3] <- B[i]
        m[i, g[i,]==4] <- V
    }
    list(anno=cbind(A,B), geno=m)
}

# Convert a MxN matrix of numeric genotypes (rows = SNPs, cols = samples) to a Mx2N matrix of
# alleles (every pair of columns is a set of haploid (unphased) genotypes).
g.i.as.h.n <- function(a, g, V="V") {
    if (!is.matrix(a)) a <- as.matrix(a)
    m <- matrix(N, nrow(g), ncol(g) * 2)
    rownames(m) <- rownames(g)
    colnames(m) <- as.vector(rbind(paste(colnames(g), "1", sep="."), paste(colnames(g), "2", sep=".")))
    for (i in 1:nrow(g)) {
        n <- c(a[i,], V)
        m[i,] <- n[as.vector(sapply(g[i,], function(x) {
            if (x == 2) {
                c(1, 2)
            }
            else if (x > 1) {
                c(x, x) - 1
            }
            else {
                c(x, x)
            }
        }))]
    }
    m
}

# Returns a list (one per element) of vectors, where each vector has two strings, one
# for each haplotype. For heterozygous sites, the first haplotype has the reference
# allele and the second has the variant allele.
g.i.as.h.n.list <- function(a, g, V.out="V", N.out="N") {
    g <- g.i.as.n(a, g)
    g.n.as.h.list(a, g, V.out=V.out, N.out=N.out)
}

# Convert from "copy" format to nucleotide format. In copy format, genotypes are encoded
# as 0, 1, or 2 copies of the variant allele or -1 for an N.
g.c.as.n <- function(a, g, N="N") {
    if (!is.matrix(g)) g <- as.matrix(g)
    g <- g + 1
    g[g == 0] <- 5
    g.i.as.n(a, g, N=N)
}

# Convert haplotypes in integer format to genotypes. Haplotypes have one line
# per SNP and two columns per sample (0 = A, 1 = B).
h.i.as.n <- function(a, h, col.names=NULL) {
    g.i <- sapply(seq(1, ncol(h), 2), function(i) {
        apply(h[,i:(i+1)] + 1, 1, function(x) bitOr(x[1], x[2]))
    })
    m <- matrix(-1, nrow(g.i), ncol(g.i), dimnames=dimnames(g.i))
    if (!is.null(col.names)) {
        colnames(m) <- col.names
    }
    for (i in 1:nrow(g.i)) {
        m[i,g.i[i,]==1] <- a[i,1]
        m[i,g.i[i,]==2] <- a[i,2]
        m[i,g.i[i,]==3] <- "H"
    }
    m
}

## Down-sampling methods
# Reduce a set of markers by removing those with the
# highest missingness and lowest MAF. m is a MxN matrix
# of numeric genotypes.
downsample.pp <- function(m) {
    print ("H")
    H <- m
    H[H!=2] <- 0
    H[H==2] <- 1

    print ("A")
    A <- m
    A[A!=1] <- 0
    A[A==1] <- 2
    A <- A + H
    
    print ("B")
    B <- m
    B[B!=3] <- 0
    B[B==3] <- 2
    B <- B + H
    
    print ("N")
    N <- m
    N[N!=5] <- 0
    N[N==5] <- 1
    
    list(A=A, B=B, N=N)
}

# Downsample by randomly iteratively removing the markers
# with the highest missingness and lowest MAF unitl the 
# standard deviation of individual missingness falls below 
# a specified value.
downsample.sd <- function(d, max.N.sd=0.01, per.iter=100) {
    with(d, {
        ix <- 1:nrow(A)
    
        rem <- NULL
        i <- 1
        while (nrow(A) > 0) {
            sample.miss <- apply(N,2,sum)/nrow(N)
            sample.miss.sd <- round(sd(sample.miss),3)
            if (sample.miss.sd < max.N.sd) {
                break
            }
            print(paste("Iteration",i,"SD =", sample.miss.sd))
            a<-apply(A,1,sum)
            b<-apply(B,1,sum)
            freq <- (pmin(a,b)/(a+b)) * 2
            miss <- apply(N,1,sum)/ncol(N)
            worst <- order(freq*(1-miss))[1:per.iter]
            rem <- c(rem, ix[worst])
            A <- A[-worst,]
            B <- B[-worst,]
            N <- N[-worst,]
            ix <- ix[-worst]
            i <- i + 1
        }
    
        sort(rem)
    })
}

# Downsample by randomly iteratively removing a specified 
# fraction of N markers from the sample with the highest 
# missingness unitl the standard deviation of individual
# missingness falls below a specified value.
downsample.worst <- function(d, max.N.sd=0.01, per.iter=0.1) {
    with(d, {
        ix <- 1:nrow(A)
    
        rem <- NULL
        iters <- NULL
        i <- 1
        while (nrow(A) > 0) {
            sample.miss <- apply(N,2,sum)/nrow(N)
            sample.miss.sd <- round(sd(sample.miss),3)
            if (sample.miss.sd <= max.N.sd) {
                break
            }
            print(paste("Iteration",i,"N =",nrow(N),"SD =", sample.miss.sd, ))
            iters <- rbind(iters, data.frame(I=i, N=nrow(N), sd=sample.miss.sd))
            w <- which(sample.miss == max(sample.miss))
            if (length(w) > 1) {
                N.snps <- which(apply(N[,w],1,function(x) any(x>0)))
            }
            else {
                N.snps <- which(N[,w] > 0)
            }
            
            snp.miss <- apply(N[N.snps,],1,sum)/length(N.snps)
            
            to.rem <- round(length(N.snps) * per.iter)
            if (to.rem < 1) {
                break
            }
            
            #sel <- sample(which(ww), to.rem)
            sel <- N.snps[order(snp.miss, decreasing=T)[1:to.rem]]
            
            rem <- c(rem, ix[sel])
            N <- N[-sel,]
            ix <- ix[-sel]
            i <- i + 1
        }
    
        list(rem=sort(rem), iters=iters)
    })
}