# Generate a PennCNV signal intensity file, which has 6 columns: 
# marker name, chromosome, position, genotype, LLR, BAF.
make.penncnv.file <- function(outfile, ...) {
    m <- get.copy.number.metrics.from.database(...)
    d <- data.frame(Name=rownames(m$snps), Chr=m$snps$chrm, Position=m$snps$pos, m$metrics[,c(1,3,2)],
        stringsAsFactors=FALSE)
    colnames(d)[4:6] <- paste(m$name, c("GType", "Log R Ratio", "B Allele Freq"), sep=".")
    write.table(d, outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}