read.fastlmm <- function(f) {
    m <- read.table(f,sep="\t",header=T,stringsAsFactors=F)
    m <- m[order(m$Chr,m$ChrPos),]
    m$P <- p.adjust(m$PValue,"holm")
    data.frame(CHR=m$Chr,BP=m$ChrPos,P=m$P,BETA=m$SnpWeight)
}