# Write genotypes as haplotypes in fastPHASE format.
# anno: four-column matrix: chrm, pos, ref, var
# geno: matrix of genotypes (snps x samples) in either integer (i) or nucleotide (n) format
write.fastphase <- function(anno, geno, outdir, format="i") {
    dir.create(outdir, showWarnings=FALSE)
    for (chrm in unique(anno[,1])) {
        print(chrm)
        w <- anno[,1] == chrm
        nsnps <- sum(w)
        print("Converting to haplotypes")
        if (format == "n") {
            g <- g.n.as.h.list(anno[w,3:4], geno[w,], V.out="?", N.out="?")
        }
        else if (format == "i") {
            g <- g.i.as.h.n.list(anno[w,3:4], geno[w,], V.out="?", N.out="?")
        }
        else {
            stop(paste("Unsupported format", format))
        }

        path <- paste(outdir, paste("chr", chrm, ".txt", sep=""), sep="/")
        f <- file(path, "w")
        writeln(f, ncol(geno))
        writeln(f, nsnps)
        writeln(f, paste("P", paste(anno[w,2], collapse=" ")))
        writeln(f, rep("S", nsnps))

        for (n in names(g)) {
            writeln(f, paste("#", n, sep=""))
            writeln(f, g[[n]][1])
            writeln(f, g[[n]][2])
        }

        close(f)
    }
}