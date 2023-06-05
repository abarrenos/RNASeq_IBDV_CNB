## Plot the Coverage of multiple Coverage tables.
## Coverage is normalized between 0 and 1.

prefixes=c("DF1_", "DF1-I_", "DF1-P_", "DF1-P-I_", "DF1-PC_", "DF1-PC-I_")
for (prefix in prefixes) {

    infiles <- list.files(pattern=paste(prefix, ".*.depth$", sep=''))

    ch <- c(rep('A', 3261), rep('B', 2827))
    nt <- c(seq(3261), seq(2827))	# 3261+2827 (ChrA + ChrB)
    a <- data.frame(ch, nt)
    i=0
    for (f in infiles) {
        i=i+1
        print(f)
        t <- read.table(f)
        M <- max(t$V3)
        a <- merge(a, t, by=c(1,2), all=TRUE)
    }
    a[is.na(a)] <- 0
    print(head(a))

    chromosomes <- levels(as.factor(a$ch))
    for (c in chromosomes) {
        png(paste(paste(prefix, ".all.depth", sep=''), c, "png", sep='.'), 
            width=1000, height=1000)
        d <- a[a$ch==c,]
        head(d)
        typ='l'	# l for line, p for points
        plot(d[ ,3], type=typ)
        for (i in 3:5) {
          if (typ == 'p')
            points(d[,i], col=i)
          else
            lines(d[, i], col=i)
        }
        legend(100,300, legend=infiles, col=3:5, lty=2:5)
        dev.off()
    }
    
}

# All together now

infiles <- list.files(pattern=paste("[Bl].*.depth$", sep=''))

ch <- c(rep('A', 3261), rep('B', 2827))
nt <- c(seq(3261), seq(2827))	# 3261+2827 (ChrA + ChrB)
a <- data.frame(ch, nt)
i=0
for (f in infiles) {
    i=i+1
    print(f)
    t <- read.table(f)
    # we'll normalilze per experiment between 0:max -> 0:1
    # We should probably normalize by chromosome below.
    M <- max(t$V3)
    t$V3 = t$V3/M
    a <- merge(a, t, by=c(1,2), all=TRUE)
}
a[is.na(a)] <- 0
print(head(a))

chromosomes <- levels(as.factor(a$ch))
ncols <- dim(a)[2]
for (c in chromosomes) {
    png(paste(paste("all.depth", sep=''), c, "png", sep='.'), 
        width=1000, height=1000)
    d <- a[a$ch==c,]
    head(d)
    typ='l'	# l for line, p for points
    plot(d[ ,3], type=typ)
    for (i in 3:ncols) {
      if (typ == 'p')
        points(d[,i], col=i)
      else
        lines(d[, i], col=i)
    }
    legend(100,300, legend=infiles, col=3:ncols, lty=2:ncols)
    dev.off()
}
