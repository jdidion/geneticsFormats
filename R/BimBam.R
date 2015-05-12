# Reads a single.txt file output by BIMBAM and prepares a data frame for plotting.
read.bimbam <- function(path) {
    bb <- read.delim(path, row.names=1, header=FALSE, skip=2)
    data.frame(rsnum=rownames(bb), chr=bb[,2], pos=bb[,1], bf=bb[,3], pv=bb[,7])
}

# Calculate the prior probability of association
ppa <- function(bf, prior=10^-6) {
    bf <- 10^bf # assumption that bf is log10(BF)
    po <- bf * prior / (1 - prior)
    po / (1 + po)
}

bfdp <- function(bf, prior=10^-6) {
    bf <- 10^bf
    1 - (bf * prior) / ((1 + bf)* prior)
}

# Calculate the log10(BF) that has the given PPA using the given prior
bf <- function(ppa=0.99, prior=10^-6) {
    p <- prior / (1 - prior)
    bf <- ppa / ((1 - ppa) * p)
    log10(bf)
}
