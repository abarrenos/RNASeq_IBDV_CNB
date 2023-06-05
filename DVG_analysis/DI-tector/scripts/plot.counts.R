#!/usr/bin/Rscript

dat <- read.table('counts/rec.bad.cnt', sep="\t", header=T)
sum <- dat[seq(4, 24, by=4), 2:7]
rownames(sum) <- dat$ID[seq(4, 24, by=4)]
sum$tot <- sum$A.AB + sum$B.AB
prop <- sum[,c(1:6)]/sum$tot

tsum <- dat[seq(4, 24, by=4), ]
rownames(tsum) <- c(1:6)
tsum$tot <- tsum$A.AB + tsum$B.AB
tprop <- tsum[ ,c(2:7)]/tsum$tot

stsum <- cbind(tsum[1], stack(tsum[2:8]))

toplot <- cbind(ID=tsum$ID, tprop[ , c(1,2,4,5)])
toplot <- toplot[2:5,]
rownames(toplot) <- 1:4
stoplot <- cbind(toplot[1], stack(toplot[2:5]))
colors <- c("red", "green", "blue", "yellow")
png("counts/barplot-proportion-from-to.png", width=1024, height=1024)
par(mar=c(6.5,4,4,6))	# bottom, left, top, right
barplot(values ~ ind + ID, data=stoplot, col=colors,
        legend.text=c("A.A", "A.B", "B.A", "B.B"),
        args.legend=list(x="right", inset=c(-0.05,0)),
        ylab="Proportion of reads From/To", beside=F)
dev.off()



typ <- read.table('counts/rec.typed.cnt', sep="\t", header=T)

png("counts/barplot-rec-by-type.png", width=1024, height=1024)
par(mar=c(6.5,4,4,6))	# bottom, left, top, right
barplot(count ~ type + ID, data=typ, col=colors,
        legend.text=T,
        args.legend=list(x="right", inset=c(-0.06,0)),
        ylab="Frequency of recombinations by type", beside=F)
dev.off()

tots <- rowsum(typ$count, typ$ID)
for (i in 1:length(typ$count)) {
    typ$prop[i] <- typ$count[i] / tots[ floor((i - 1) / 4) + 1 ]
    cat(typ$ID[i], typ$type[i], typ$prop[i], typ$count[i], tots[ floor((i - 1) / 4) + 1 ], floor((i - 1) / 4) + 1, '\n' )

}

png("counts/barplot-proportion-rec-by-type.png", width=1024, height=1024)
par(mar=c(6.5,4,4,6))	# bottom, left, top, right
barplot(prop ~ type + ID, data=typ, col=colors,
        legend.text=T,
        args.legend=list(x="right", inset=c(-0.06,0)),
        ylab="Proportion of recombinations by type", beside=F)
dev.off()
