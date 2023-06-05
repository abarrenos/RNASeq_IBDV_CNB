#!/bin/bash
# this is to be run from the base directory after 'analyze_sample.sh' has
# been run from within each sample subdirectory.

here=`realpath .`
base=$1		# experiment name (followed by -# (number of sample)

# we will use the counts file because that is what the web page uses

out=summary/counts
plt=summary/plots
prd=summary/period

mkdir -p $out
mkdir -p $plt
mkdir -p $prd

r1_counts=$out/${base}_R1_counts.txt
r2_counts=$out/${base}_R2_counts.txt


# join the three experiments per PFU level
cat ${base}-[123]/out/R1counts.txt \
| grep -v "^None" \
| sort -k1 -n -k3 -t'	' \
| tr -d "'" \
> $r1_counts

cat ${base}-[123]/out/R2counts.txt \
| grep -v "^None" \
| sort -k1 -n -k3 -t'	' \
| tr -d "'" \
> $r2_counts

if [ ! -s $r1_counts -a ! -s $r2_counts ] ; then 
    echo "NO DVGS FOUND FOR $base"
    exit 
fi

# separate the three experiments collated by segment(s) affected by the DVG
grep 'A|A' $r1_counts > $out/${base}_R1counts_AA.txt
grep 'A|B' $r1_counts > $out/${base}_R1counts_AB.txt
grep 'B|A' $r1_counts > $out/${base}_R1counts_BA.txt
grep 'B|B' $r1_counts > $out/${base}_R1counts_BB.txt

grep 'A|A' $r2_counts > $out/${base}_R2counts_AA.txt
grep 'A|B' $r2_counts > $out/${base}_R2counts_AB.txt
grep 'B|A' $r2_counts > $out/${base}_R2counts_BA.txt
grep 'B|B' $r2_counts > $out/${base}_R2counts_BB.txt

# join both read directions
cat $out/${base}_R1counts_AA.txt $out/${base}_R2counts_AA.txt > $out/${base}_R12counts_AA.txt
cat $out/${base}_R1counts_AB.txt $out/${base}_R2counts_AB.txt > $out/${base}_R12counts_AB.txt
cat $out/${base}_R1counts_BA.txt $out/${base}_R2counts_BA.txt > $out/${base}_R12counts_BA.txt
cat $out/${base}_R1counts_BB.txt $out/${base}_R2counts_BB.txt > $out/${base}_R12counts_BB.txt

# separate by DVG type
grep "^3" $r1_counts > $out/${base}_R1counts_3cbsb.txt
grep "^5" $r1_counts > $out/${base}_R1counts_5cbsb.txt
grep "^D" $r1_counts > $out/${base}_R1counts_deletions.txt
grep "^I" $r1_counts > $out/${base}_R1counts_insertions.txt

grep "^3" $r2_counts > $out/${base}_R2counts_3cbsb.txt
grep "^5" $r2_counts > $out/${base}_R2counts_5cbsb.txt
grep "^D" $r2_counts > $out/${base}_R2counts_deletions.txt
grep "^I" $r2_counts > $out/${base}_R2counts_insertions.txt


cat $out/${base}_R1counts_deletions.txt $out/${base}_R2counts_deletions.txt > $out/${base}_R12counts_deletions.txt
cat $out/${base}_R1counts_insertions.txt $out/${base}_R2counts_insertions.txt > $out/${base}_R12counts_insertions.txt

# prepare data with headers for convenience
#	data will be stored as TSV with values enclosed in double quotes
echo '"DVG.s.type"	"Length"	"BP_Pos"	"RI_Pos"	"Delta_Positions"	"Ref"	"Counts"	"Pct_to_Virus"' \
> $out/${base}_R1counts+header.txt
cat $r1_counts \
| sed -e 's/^/"/g' -e 's/$/"/g' -e 's/	/"	"/g' \
>> $out/${base}_R1counts+header.txt

echo '"DVG.s.type"	"Length"	"BP_Pos"	"RI_Pos"	"Delta_Positions"	"Ref"	"Counts"	"Pct_to_Virus"' \
> $out/${base}_R2counts+header.txt
cat $r2_counts \
| sed -e 's/^/"/g' -e 's/$/"/g' -e 's/	/"	"/g' \
>> $out/${base}_R2counts+header.txt

cat $out/${base}_R1counts+header.txt > $out/${base}_R12counts+header.txt
cat $r2_counts \
| sed -e 's/^/"/g' -e 's/$/"/g' -e 's/	/"	"/g' \
>> $out/${base}_R12counts+header.txt



# Make plots
#R --no-save --no-restore <<END
cat > summary/analyze-$base.R <<END
    out <- "$out"
    base <- "$base"
    plt <- "$plt"
    prd <- "$prd"

    hist.col.by <- function(pngfile, dfile, col, by, max) {
        # make a histogram of column 'col' classified by factor column 'by'
	# the x-axis will go from zero to 'max'
	# the data will be read from 'dfile' and the plot, stored in 'pngfile'
	df <- read.table(dfile, sep="\t", stringsAsFactors=TRUE,strip.white=TRUE)
    	colnames(df) <- c("DVG_Type", "Length", "BP_Pos", "RI_Pos", 
                          "Delta_pos", "Ref", "Counts", "Pct_to_Virus")
        #df <- df[order(df[ , by]), ]
        binsize <- 20
        bins <- seq(0, max+binsize, by=binsize)
        png(pngfile, width=1024, height=1024)
        print(paste("plotting to", pngfile))
        groups <- levels(df[ , by])
        #print(groups)
        ngroups <- length(groups)
        plt.dim <- ceiling(sqrt(ngroups))
        #print(paste("Plot will be", plt.dim, "x", plt.dim))
        par(mfrow=c(plt.dim, plt.dim))
        color <- c("cyan", "red", "green", "orange", "blue", "yellow", "white")
        i <- 0
        for (group in groups) {
            i <- i + 1
            grp <- df[df[ ,by] == group, ]
            print(paste("    plotting", group, "[", dim(grp)[1], dim(grp)[2], "]"))
	    hist( grp[ , col], 
            	  breaks=bins, col=color[i], xlab=col, main=group)
        }
        #mtext(paste(pngfile, col, "by", by), side=3, outer=T, line=-3)
        dev.off()
	# repeat separately
        i <- 0
        for (group in groups) {
            i <- i + 1
            grp <- df[df[ ,by] == group, ]
            print(paste("    plotting", group, "[", dim(grp)[1], dim(grp)[2], "]"))
	    par(mfrow=c(1,1))
	    png1 <- paste( file_path_sans_ext(pngfile), i, 'png', sep='.')
	    png(png1, width=1024, height=1024)
	    hist( grp[ , col], 
            	  breaks=bins, col=color[i], xlab=col, main=group)
            dev.off()
	}
        #mtext(paste(pngfile, col, "by", by), side=3, outer=T, line=-3)
        dev.off()
    }

    # plot BP-RI connected by lines
    plot.cuts.1 <- function(pngfile, dfile, max) {
        # this function plots the break-points+reinitiation points in the 
	# case where both lie in the same segment (*_AA or *_BB)
	# Data is stored in 'dfile'
	# x-axis goes from 0 to 'max' (which should be the segment size)
	# The plot is arranged by DVGtype
	#
        df <- read.table(dfile, sep="\t", stringsAsFactors=TRUE,strip.white=TRUE)
    	colnames(df) <- c("DVG_Type", "Length", "BP_Pos", "RI_Pos", 
                          "Delta_pos", "Ref", "Counts", "Pct_to_Virus")
        df\$size <- df\$RI_Pos - df\$BP_Pos
        df\$color <- ifelse(df\$size > 0, "darkgreen", "red")
    	png(pngfile, width=1024, height=1024)
          print(paste("Plotting to", pngfile))
          #print(head(df))
          groups <- levels(df\$DVG_Type)
          #print(groups)
          ngroups <- length(groups)
          plt.dim <- ceiling(sqrt(ngroups))
          print(paste("Plot will be", plt.dim, "x", plt.dim))
          par(mfrow=c(plt.dim, plt.dim))
          color <- c("cyan", "red", "green", "orange", "blue", "yellow", "white")
          i <- 0
          for (group in groups) {
              i <- i + 1
              grp <- df[df\$DVG_Type == group, ]
              #print(head(grp))
              grp <- grp[order(grp\$BP_Pos, grp\$RI_Pos), ]
              #print("After sorting")
              #print(head(grp))
              print(paste("    plotting", group, "[", dim(grp)[1], dim(grp)[2], "]"))
              npoints <- dim(grp)[1]
    	      vspacing <- floor(1000 / npoints)	# make sure it fits
    	      grp\$y <- seq(vspacing, length=npoints, by=vspacing)
              #print(head(grp\$y))
              #print(npoints)
              #print(vspacing)
              for (k in 1:npoints) {
                  grp\$y[k] <- k * 1000 / npoints
              }
    	      xlim=c(0, max)
    	      ylim=c(0, 1000)
    	      plot(0, 0, xlim=xlim, type="n", ylim=ylim, 
        		 xlab='', ylab='',  yaxt='n',
                	 main=paste(group, "jumps (", npoints, "/", sum(grp\$Counts), ")"))
              segments(x0=grp\$BP_Pos, x1=grp\$RI_Pos, y0=grp\$y, y1=grp\$y, col=grp\$color)
          }
        dev.off()
	#
	# now repeat individually
        i <- 0
        for (group in groups) {
            i <- i + 1
            grp <- df[df\$DVG_Type == group, ]
            #print(head(grp))
            grp <- grp[order(grp\$Delta_pos, grp\$BP_Pos, grp\$RI_Pos), ]
            #print("After sorting")
            #print(head(grp))
            print(paste("    plotting", group, "[", dim(grp)[1], dim(grp)[2], "]"))
	    par(mfrow=c(1,1))
	    png1 <- paste( file_path_sans_ext(pngfile), i, 'png', sep='.')
	    png(png1, width=1024, height=1024)
	      npoints <- dim(grp)[1]
    	      vspacing <- floor(1000 / npoints)	# make sure it fits
    	      grp\$y <- seq(vspacing, length=npoints, by=vspacing)
              #print(head(grp\$y))
              #print(npoints)
              #print(vspacing)
              for (k in 1:npoints) {
                  grp\$y[k] <- k * 1000 / npoints
              }
    	      xlim=c(0, max)
    	      ylim=c(0, 1000)
    	      plot(0, 0, xlim=xlim, type="n", ylim=ylim, 
        		 xlab='', ylab='',  yaxt='n',
                	 main=paste(group, "jumps (", npoints, "/", sum(grp\$Counts), ")"))
              segments(x0=grp\$BP_Pos, x1=grp\$RI_Pos, y0=grp\$y, y1=grp\$y, col=grp\$color)
	    dev.off()
        }
    }

    plot.cuts.2 <- function(pngfile, dfile, max1, max2) {
        # this function plots the break-points+reinitiation points in the 
	# case where the DVGs span two segments: *_AB or *_BA
	# Data is stored in 'dfile'
	# x-axis goes from 0 to 'max' (which should be the segment size)
	# The plot is arranged by DVGtype
	#
    	df <- read.table(dfile, sep="\t", stringsAsFactors=TRUE,strip.white=TRUE)
    	colnames(df) <- c("DVG_Type", "Length", "BP_Pos", "RI_Pos", 
                          "Delta_pos", "Ref", "Counts", "Pct_to_Virus")

        df\$end <- df\$RI_Pos + max1	# add length of first segment to get coords in second segment
        df\$color <- rep("darkgreen", dim(df)[1])
    	png(pngfile, width=1024, height=1024)
        print(paste("Plotting to", pngfile))
        #print(head(df))
        groups <- levels(df\$DVG_Type)
        #print(groups)
        ngroups <- length(groups)
        plt.dim <- ceiling(sqrt(ngroups))
        print(paste("Plot will be", plt.dim, "x", plt.dim))
        par(mfrow=c(plt.dim, plt.dim))
        color <- c("cyan", "red", "green", "orange", "blue", "yellow", "white")
        i <- 0
        for (group in groups) {
            i <- i + 1
            grp <- df[df\$DVG_Type == group, ]
            #print(head(grp))
            grp <- grp[order(grp\$BP_Pos, grp\$RI_Pos), ]
            #print("After sorting")
            #print(head(grp))
            print(paste("    plotting", group, "[", dim(grp)[1], dim(grp)[2], "]"))
            npoints <- dim(grp)[1]
    	    vspacing <- floor(1000 / npoints)	# make sure it fits
    	    grp\$y <- seq(vspacing, length=npoints, by=vspacing)
            for (k in 1:npoints) {
                grp\$y[k] <- k * 1000 / npoints
            }
            xlim=c(0, max1 + max2)
    	    ylim=c(0, 1000)
    	    plot(0, 0, xlim=xlim, type="n", ylim=ylim, 
        	       xlab='', ylab='', 
                       xaxt='n', yaxt='n',
                       main=paste(group, "jumps (", npoints, "/", sum(grp\$Counts), ")"))
            at <- c(seq(0, max1, by=1000), max1, seq(max1, max1+max2, by=1000), max1+max2)
            att <- c(seq(0, max1, by=1000), max1, seq(0, max2, by=1000), max2)
            bottom <- 1
            axis(bottom, at=at, labels=as.character(att), las=2)
            segments(x0=grp\$BP_Pos, x1=grp\$end, 
                     y0=grp\$y, y1=grp\$y, 
                     col=grp\$color)
            abline(v=max1, col="black")
        }
        dev.off()

        # Now repeat individually
        i <- 0
        for (group in groups) {
            i <- i + 1
            grp <- df[df\$DVG_Type == group, ]
            #print(head(grp))
            grp <- grp[order(grp\$Delta_pos, grp\$BP_Pos, grp\$RI_Pos), ]
            #print("After sorting")
            #print(head(grp))
            print(paste("    plotting", group, "[", dim(grp)[1], dim(grp)[2], "]"))
	    par(mfrow=c(1,1))
	    png1 <- paste( file_path_sans_ext(pngfile), i, 'png', sep='.')
	    png(png1, width=1024, height=1024)
	      npoints <- dim(grp)[1]
    	      vspacing <- floor(1000 / npoints)	# make sure it fits
    	      grp\$y <- seq(vspacing, length=npoints, by=vspacing)
              for (k in 1:npoints) {
                  grp\$y[k] <- k * 1000 / npoints
              }
              xlim=c(0, max1 + max2)
    	      ylim=c(0, 1000)
    	      plot(0, 0, xlim=xlim, type="n", ylim=ylim, 
        		 xlab='', ylab='', 
                	 xaxt='n', yaxt='n',
                	 main=paste(group, "jumps (", npoints, "/", sum(grp\$Counts), ")"))
              at <- c(seq(0, max1, by=1000), max1, seq(max1, max1+max2, by=1000), max1+max2)
              att <- c(seq(0, max1, by=1000), max1, seq(0, max2, by=1000), max2)
              bottom <- 1
              axis(bottom, at=at, labels=as.character(att), las=2)
              segments(x0=grp\$BP_Pos, x1=grp\$end, 
                       y0=grp\$y, y1=grp\$y, 
                       col=grp\$color)
              abline(v=max1, col="black")
	    dev.off()
        }
    }


    # we'll start by stating the length of each genome segment
    lenA <- 3261
    lenB <- 2827
    r1 <- read.table("$out/${base}_R1counts+header.txt", header=TRUE, 
      stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    r2 <- read.table("$out/${base}_R2counts+header.txt", header=TRUE, 
      stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    r12 <- read.table("$out/${base}_R12counts+header.txt", header=TRUE, 
      stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)


    par(mfrow=c(1,1))
    library(lattice)
    library(tools)
    # BP by DVG
    png("$plt/${base}_R1_BPbyDVG.png", width=1024, height=1024)
    xyplot(Counts ~ BP_Pos | DVG.s.type, groups=Ref, 
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r1)
    dev.off()
    png("$plt/${base}_R2_BPbyDVG.png", width=1024, height=1024)
    xyplot(Counts ~ BP_Pos | DVG.s.type, groups=Ref, 
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r2)
    dev.off()
    png("$plt/${base}_R12_BPbyDVG.png", width=1024, height=1024)
    xyplot(Counts ~ BP_Pos | DVG.s.type, groups=Ref, 
      type="p", pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r12)
    dev.off()

    # BP by Segment
    png("$plt/${base}_R1_BPbySegment.png", width=1024, height=1024)
    xyplot(Counts ~ BP_Pos | Ref, groups=DVG.s.type, 
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r1)
    dev.off()
    png("$plt/${base}_R2_BPbySegment.png", width=1024, height=1024)
    xyplot(Counts ~ BP_Pos | Ref, groups=DVG.s.type,
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r2)
    dev.off()
    png("$plt/${base}_R12_BPbySegment.png", width=1024, height=1024)
    xyplot(Counts ~ BP_Pos | Ref, groups=DVG.s.type,
      type="p", pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r12)
    dev.off()

    # RI by DVG
    png("$plt/${base}_R1_RIbyDVG.png", width=1024, height=1024)
    xyplot(Counts ~ RI_Pos | DVG.s.type, groups=Ref, 
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r1)
    dev.off()
    png("$plt/${base}_R2_RIbyDVG.png", width=1024, height=1024)
    xyplot(Counts ~ RI_Pos | DVG.s.type, groups=Ref, 
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r2)
    dev.off()
    png("$plt/${base}_R12_RIbyDVG.png", width=1024, height=1024)
    xyplot(Counts ~ RI_Pos | DVG.s.type, groups=Ref, 
      type="p", pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r12)
    dev.off()
    
    # RI by Segment
    png("$plt/${base}_R1_RIbySegment.png", width=1024, height=1024)
    xyplot(Counts ~ RI_Pos | Ref, groups=DVG.s.type, 
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r1)
    dev.off()
    png("$plt/${base}_R2_RIbySegment.png", width=1024, height=1024)
    xyplot(Counts ~ RI_Pos | Ref, groups=DVG.s.type,
      type="p", 
      pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r2)
    dev.off()
    png("$plt/${base}_R12_RIbySegment.png", width=1024, height=1024)
    xyplot(Counts ~ RI_Pos | Ref, groups=DVG.s.type,
      type="p", pch=16, 
      auto.key=list(border=TRUE), 
      par.settings=simpleTheme(pch=16), 
      scales=list(x=list(relation='free'), y=list(relation='free')), 
      col=c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(0,0,0,1/4), rgb(1,0,0,1/4)),
      data=r12)
    dev.off()


    # make histograms with a binsize of 20 nt.
    binsize <- 20
    binsA <- seq(0, lenA+binsize, by=binsize)
    binsB <- seq(0, lenB+binsize, by=binsize)

    # now plot breakpoint and reinitiation site positions by segment
    # First R1
    
    dfile <- "$out/${base}_R1counts_AA.txt"
    hist.col.by("$plt/${base}_R1_AA_bp.png", dfile, "BP_Pos", "DVG_Type", lenA) 
    hist.col.by("$plt/${base}_R1_AA_ri.png", dfile, "RI_Pos", "DVG_Type", lenA)
    plot.cuts.1("$plt/${base}_R1_AA_frags.png", dfile, lenA)

    dfile <- "$out/${base}_R1counts_BB.txt"
    hist.col.by("$plt/${base}_R1_BB_bp.png", dfile, "BP_Pos", "DVG_Type", lenB) 
    hist.col.by("$plt/${base}_R1_BB_ri.png", dfile, "RI_Pos", "DVG_Type", lenB)
    plot.cuts.1("$plt/${base}_R1_BB_frags.png", dfile, lenB)
      
    # and R2
    dfile <- "$out/${base}_R2counts_AA.txt"
    hist.col.by("$plt/${base}_R2_AA_bp.png", dfile, "BP_Pos", "DVG_Type", lenA) 
    hist.col.by("$plt/${base}_R2_AA_ri.png", dfile, "RI_Pos", "DVG_Type", lenA)
    plot.cuts.1("$plt/${base}_R2_AA_frags.png", dfile, lenA)

    dfile <- "$out/${base}_R2counts_BB.txt"
    hist.col.by("$plt/${base}_R2_BB_bp.png", dfile, "BP_Pos", "DVG_Type", lenB) 
    hist.col.by("$plt/${base}_R2_BB_ri.png", dfile, "RI_Pos", "DVG_Type", lenB)
    plot.cuts.1("$plt/${base}_R2_BB_frags.png", dfile, lenB)

    plot.cuts.2("$plt/${base}_R1_AB_frags.png", "$out/${base}_R1counts_AB.txt", lenA, lenB)
    plot.cuts.2("$plt/${base}_R1_BA_frags.png", "$out/${base}_R1counts_BA.txt", lenB, lenA)
    
    plot.cuts.2("$plt/${base}_R2_AB_frags.png", "$out/${base}_R2counts_AB.txt", lenA, lenB)
    plot.cuts.2("$plt/${base}_R2_BA_frags.png", "$out/${base}_R2counts_BA.txt", lenB, lenA)

    # plot together (R1+R2)
    hist.col.by("$plt/${base}_R12_AA_bp.png", "$out/${base}_R12counts_AA.txt", "BP_Pos", "DVG_Type", lenA) 
    hist.col.by("$plt/${base}_R12_AA_ri.png", "$out/${base}_R12counts_AA.txt", "RI_Pos", "DVG_Type", lenA)
    plot.cuts.1("$plt/${base}_R12_AA_frags.png", "$out/${base}_R12counts_AA.txt", lenA)

    hist.col.by("$plt/${base}_R12_BB_bp.png", "$out/${base}_R12counts_BB.txt", "BP_Pos", "DVG_Type", lenB) 
    hist.col.by("$plt/${base}_R12_BB_ri.png", "$out/${base}_R12counts_BB.txt", "RI_Pos", "DVG_Type", lenB)
    plot.cuts.1("$plt/${base}_R12_BB_frags.png", "$out/${base}_R12counts_BB.txt", lenB)

    hist.col.by("$plt/${base}_R12_AB_bp.png", "$out/${base}_R12counts_AB.txt", "BP_Pos", "DVG_Type", lenA) 
    hist.col.by("$plt/${base}_R12_AB_ri.png", "$out/${base}_R12counts_AB.txt", "RI_Pos", "DVG_Type", lenB)
    plot.cuts.2("$plt/${base}_R12_AB_frags.png", "$out/${base}_R12counts_AB.txt", lenA, lenB)

    hist.col.by("$plt/${base}_R12_BA_bp.png", "$out/${base}_R12counts_BA.txt", "BP_Pos", "DVG_Type", lenB) 
    hist.col.by("$plt/${base}_R12_BA_ri.png", "$out/${base}_R12counts_BA.txt", "RI_Pos", "DVG_Type", lenA)
    plot.cuts.2("$plt/${base}_R12_BA_frags.png", "$out/${base}_R12counts_BA.txt", lenB, lenA)

    # compute periodicity
    #install.packages("EMD")
    library(EMD)
    
    plot.emd <- function(pdffile, dfile, max) {
    	# same segment (*_AA or *_BB)
        df <- read.table(dfile, sep="\t", stringsAsFactors=TRUE,strip.white=TRUE)
    	colnames(df) <- c("DVG_Type", "Length", "BP_Pos", "RI_Pos", 
                          "Delta_pos", "Ref", "Counts", "Pct_to_Virus")
        df\$size <- df\$RI_Pos - df\$BP_Pos
        df\$color <- ifelse(df\$size > 0, "darkgreen", "red")
    	pdf(pdffile)
        print(paste("Plotting to", pdffile))
        #print(head(df))
        groups <- levels(df\$DVG_Type)
        #print(groups)
        ngroups <- length(groups)
        par(mar=c(1,1,1,1))	# to avoid error margins too large
        for (group in groups) {
            print(paste("    Plotting", dfile, group))
            grp <- df[df\$DVG_Type == group, ]
            
            par(mfrow=c(1,1))
            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 8, dfile)
            text(5, 7, group)
            text(5, 6, "BP_Pos EMD")
            xbp <- rep(0, max)
            xbp[ grp\$BP_Pos ] <- grp\$Counts
            xbp.emd <- emd(xbp, max.imf=15) # plot.imf=F
            nplots <- as.integer(ceiling(xbp.emd\$nimf / 2)) * 2
            par(mfrow=c(nplots/2, 2))            
            for (ni in 1:xbp.emd\$nimf) {
                plot(xbp.emd\$imf[, ni], type='l', xlab=group)
		#print(ni)
            }
            #print(paste("ni =", ni, "xbp.emd\$nimf =", xbp.emd\$nimf, "nplots =", nplots))
            if (xbp.emd\$nimf < nplots) { 
              for (ni in (xbp.emd\$nimf+1):nplots) { 
                plot.new(); #print(paste("fill",ni)) 
              }
            }
            
            par(mfrow=c(1,1))
            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 8, dfile)
            text(5, 7, group)
            text(5, 6, "RI_Pos EMD")
            xri <- rep(0, max)
            xri[ grp\$RI_Pos ] <- grp\$Counts
            xri.emd <- emd(xri, max.imf=15) # plot.imf=F
            xri.emd <- emd(xri, max.imf=15) # plot.imf=F
            nplots <- ceiling(xri.emd\$nimf / 2) * 2
            par(mfrow=c(ceiling(xri.emd\$nimf / 2), 2))
            for (ni in 1:xri.emd\$nimf) {
                plot(xri.emd\$imf[, ni], type='l', xlab=group)
                #print(ni)
            }
            #print(paste("ni =", ni, "xri.emd\$nimf =", xri.emd\$nimf, "nplots =", nplots))
            if (xri.emd\$nimf < nplots) { 
              for (ni in (xri.emd\$nimf+1):nplots) { 
                plot.new(); #print(paste("fill", ni)) 
              } 
            }
        }
        dev.off()
    }

    plot.emd.simple <- function(pdffile, dfile, max) {
    	# same segment (*_AA or *_BB)
        df <- read.table(dfile, sep="\t", stringsAsFactors=TRUE,strip.white=TRUE)
    	colnames(df) <- c("DVG_Type", "Length", "BP_Pos", "RI_Pos", 
                          "Delta_pos", "Ref", "Counts", "Pct_to_Virus")
        df\$size <- df\$RI_Pos - df\$BP_Pos
        df\$color <- ifelse(df\$size > 0, "darkgreen", "red")
    	pdf(pdffile)
        print(paste("Plotting to", pdffile))
        #print(head(df))
        groups <- levels(df\$DVG_Type)
        #print(groups)
        ngroups <- length(groups)
        ngroups <- length(groups)
        plt.dim <- ceiling(sqrt(ngroups))
        print(paste("Plot will be", plt.dim, "x", plt.dim))
        #par(mfrow=c(plt.dim, plt.dim))
        #par(mar=c(1,1,1,1))	# to avoid error margins too large
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        par(mfrow=c(ngroups, 1))
        for (group in groups) {
            print(paste("    Plotting", dfile, group))
            grp <- df[df\$DVG_Type == group, ]
            
            #par(mfrow=c(1,1))
            #plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            #text(5, 8, dfile)
            #text(5, 7, group)
            #text(5, 6, "BP_Pos EMD")
            x <- rep(0, max)
            x[ grp\$BP_Pos ] <- grp\$Counts
            x.emd <- emd(x, max.imf=15) # plot.imf=F
            # Keep 5 components -- we may need more, or less.
            m <- x.emd\$nimf/2 ; M <- x.emd\$nimf
            print(paste('No. IMF', x.emd\$nimf, 'low IMF', m, 'high IMF', M))
            y <- apply( x.emd\$imf[,m:M], 1, sum ) + mean(x.emd\$residue)
            #y <- apply( x.emd\$imf[,5:10], 1, sum ) + mean(x.emd\$residue)
            # plot the original data
            plot(x, type="l", col="grey", xlab=paste(group, "BP_Pos"), main=dfile)
            # add residuals plus mean as wide line
            lines( y, type="l", lwd=2) 
            # add point at the maxima
            n <- length(y)
            i <- y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] > y[3:n]
            tryCatch(
	      { points( which(i), y[i], pch=15, col="blue" ) },
              error=function(cond) { return() }
            )
        }
        for (group in groups) {
            print(paste("    Plotting", dfile, group))
            grp <- df[df\$DVG_Type == group, ]
            
            #par(mfrow=c(1,1))
            #plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            #text(5, 8, dfile)
            #text(5, 7, group)
            #text(5, 6, "BP_Pos EMD")
            x <- rep(0, max)
            x[ grp\$RI_Pos ] <- grp\$Counts
            x.emd <- emd(x, max.imf=15) # plot.imf=F
            # Keep 5 components -- we may need more, or less.
            m <- x.emd\$nimf/2 ; M <- x.emd\$nimf
            #print(paste('No. IMF', x.emd\$nimf, 'low IMF', m, 'high IMF', M))
            y <- apply( x.emd\$imf[,m:M], 1, sum ) + mean(x.emd\$residue)
            #y <- apply( x.emd\$imf[,5:10], 1, sum ) + mean(x.emd\$residue)
            # plot the original data
            plot(x, type="l", col="grey", xlab=paste(group, "RI_Pos"), main=dfile)
            # add residuals plus mean as wide line
            lines( y, type="l", lwd=2) 
            # add point at the maxima
            n <- length(y)
            i <- y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] > y[3:n]
            tryCatch(
	      { points( which(i), y[i], pch=15, col="blue" ) },
              error=function(cond) { return() }
            )
        }
        dev.off()
    }

    plot.emd("$prd/${base}_R12_EMD_AA.pdf", "$out/${base}_R12counts_AA.txt", lenA)
    plot.emd("$prd/${base}_R12_EMD_BB.pdf", "$out/${base}_R12counts_BB.txt", lenB)

    plot.emd.simple("$prd/${base}_R12_EMDsum_AA.pdf", "$out/${base}_R12counts_AA.txt", lenA)
    plot.emd.simple("$prd/${base}_R12_EMDsum_BB.pdf", "$out/${base}_R12counts_BB.txt", lenB)

END



