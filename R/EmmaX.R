read.emmax.projects <- function(dir, proj.names, phenos) {
    projects <- lapply(proj.names, function(n) {
        print(n)
        proj.results <- lapply(phenos, function(p) {
            print(p)
            d <- file.path(dir, n, "result", p)
            read.emmax.results(d)
        })
        names(proj.results) <- phenos
        proj.results
    })
    names(projects) <- proj.names
    projects
}

# Collate the results of N emmax runs, where N is the number of SNPs.
read.emmax.results <- function(dir, n=NULL, pad=6) {
    if (is.null(n)) {
        results <- do.call(rbind, lapply(Sys.glob(file.path(dir, "out_*.ps")), function(f) {
            s <- nchar(f) - (pad+2)
            i <- as.integer(substr(f, s, s+pad-1))
            print(i)
            tab <- read.table(f, sep="\t", header=F, stringsAsFactors=F)
            rownames(tab) <- i
            tab
        }))
        results <- results[order(as.integer(rownames(results))),]
    }
    else {
        results <- do.call(rbind, lapply(1:n, function(i) {
            read.table(file.path(dir, paste("out_", sprintf(paste("%0",pad,"d",sep=""), i), ".ps", sep="")), 
                sep="\t", header=F, stringsAsFactors=F)
        }))
    }
    colnames(results) <- c("SNP", "SE", "p")
    results
}

plot.emmax.projects <- function(proj.results, maps, outdir, adjust="holm", sig=c(0.1, 0.05)) {
    for (proj.name in names(proj.results)) {
        phenos <- proj.results[[proj.name]]
        for (pheno in names(phenos)) {
            print(paste(proj.name, pheno))
            proj <- phenos[[pheno]]
            p <- proj$p
            if (!is.null(adjust)) {
                p <- p.adjust(p, adjust)
            }
            path <- file.path(outdir, paste(proj.name, pheno, "jpg", sep="."))
            map <- maps[[proj.name]]
            m <- match(proj$SNP, map[,2])
            d <- data.frame(CHR=map[m,1], BP=map[m,4], P=p)
            manhattan.freq(d, threshold=-log10(sig), dev.file=path, 
                main=paste("Project:", proj.name, "Pheno:", pheno, "Correction:", adjust))
        }
    }
}
