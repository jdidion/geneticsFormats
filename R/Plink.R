# Read PLINK PED and MAP files.

read.plink <- function(file) {
    ped <- read.table(paste(file, ".ped", sep=""), header=F, stringsAsFactors=F, colClasses="character", comment.char="")
    map <- read.table(paste(file, ".map", sep=""), header=F, stringsAsFactors=F, colClasses=c("character","character","numeric","integer"), comment.char="")

    colnames(map) <- c("Chromosome", "SNPID", "Morgans", "BP")
    map$Chromosome <- factor(map$Chromosome, levels=unique(map$Chromosome))

    fam <- ped[,1:6]
    colnames(fam) <- c("FamilyID","IndividualID","PaternalID","MaternalID","Sex","Phenotype")
    # factor Family ID and Phenotype
    fam$FamilyID <- factor(fam$FamilyID)
    fam$Phenotype <- factor(fam$Phenotype)
    # convert Sex to an integer, with 0 = other
    other <- !(fam$Sex %in% c("1","2"))
    if (any(other)) {
        fam$Sex[other] <- "0"
    }
    fam$Sex <- as.integer(fam$Sex)
    
    gen <- ped[,7:ncol(ped)]
    temp <- lapply(seq(1, ncol(gen), 2), function(i) {
        a1 <- gen[,i]
        a2 <- gen[,i+1]
        
        g1 <- rep(0, length(a1))
        g2 <- rep(0, length(a2))
        
        u <- setdiff(unique(c(a1, a2)), "0")
        for (j in 1:length(u)) {
            g1[a1 == u[j]] <- j
            g2[a2 == u[j]] <- j
        }
        if (length(u) < 2) {
            u <- c(u, rep(NA, 2-length(u)))
        }
        names(u) <- c("A","B")
        
        list(alleles=u, geno=bitOr(g1, g2))
    })
    
    map <- cbind(map, do.call(rbind, lapply(temp, function(x) x$alleles)))
    gen <- do.call(cbind, lapply(temp, function(x) x$geno))
    
    list(fam=fam, map=map, gen=gen)
}

# Given a "full" PLINK .genome file (made by using the --genome --genome-full flags),
# produce an NxN pairwise matrix of IBS2* values.
genome.to.matrix <- function(g) {
    if (is.character(g)) {
        g <- read.table(g, header=T)
    }
    u <- unique(c(g$IID1,g$IID2))
    m <- matrix(0, nrow=length(u), ncol=length(u), dimnames=list(u,u))
    w <- cbind(ix1=match(g$IID1, u), ix2=match(g$IID2, u))
    IBS2 <- g$HETHET/(g$HOMHOM+g$HETHET)
    for (i in 1:nrow(w)) {
        m[w[i,1],w[i,2]] <- m[w[i,2],w[i,1]] <- IBS2[i]
    }
    m
}

# Given a matrix of IBS2* values and a threshold, remove the
# fewest samples that eliminate all related pairs in the matrix.
remove.related <- function(m, max.ibs2=0.8) {
    while (max(m[upper.tri(m)]) > max.ibs2) {
        n <- apply(m, 1, function(x) sum(x > max.ibs2))
        w <- which(n == max(n))
        if (length(w) > 1) {
            w <- sample(w,1)
        }
        m <- m[-w,-w]
    }
    m
}

subsample.indivs.hwe <- function(prefix, N, extract.file, outdir=".", niters=10, 
        threshold=0.05, num.signif=0) {
    # load the fam file
    fam <- read.table(paste(prefix, "fam", sep="."))
    fam <- fam[fam[,6] == 1,]
    # perform niters subsamplings of individuals, perform HWE tests, and extract
    # holm-corrected p-values
    snps <- NULL
    tab <- sapply(1:niters, function(i) {
        outfile <- file.path(outdir, paste("subsample", i, sep="."))
        write.table(fam[sort(sample(1:nrow(fam),N)),1:2],outfile,row.names=F,col.names=F,quote=F)
        system(paste("plink --bfile", prefix, "--keep", outfile, "--hardy --out", outfile))
        hwe.file <- paste(outfile, "hwe", sep=".")
        hwe <- read.table(hwe.file, header=T, stringsAsFactors=F)
        hwe <- hwe[hwe$TEST=="UNAFF",]
        if (i == 1) {
            snps <<- hwe[,2]
        }
        p.adjust(hwe$P, "holm")
    })
    # for each marker, test whether all subsamplings yielded
    # a p-value greater than threshold
    w <- apply(tab, 1, function(x) sum(x <= threshold) <= num.signif)
    # write a SNP list containing the markers that pass
    write.table(snps[w], extract.file, row.names=F, col.names=F, quote=F)
    return(sum(w))
}
