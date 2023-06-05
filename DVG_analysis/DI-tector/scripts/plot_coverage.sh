#!bin/bash

## PLOT COVERAGE DEPTH FROM SAMTOOLS OR BEDTOOLS ##
#
#  Arguments:
#	CoverageTable	- Input file containing coverage depth table
#	OutputName	- Output name for the png file (Default Input file name)
#


COV=$1
OUT=$2

: ${OUT:~`basename $COV .txt`}

R --vanilla <<END
    cov.table <- read.table("$COV")
    png(paste("$OUT",".A.png", sep=""), width=1000, height=1000)
    fragment.A <- cov.table[cov.table\$V1=='A',3]
    plot(fragment.A)
    dev.off()
    png(paste("$OUT",".B.png", sep=""), width=1000, height=1000)
    fragment.B <- cov.table[cov.table\$V1=='B',3]
    plot(fragment.B)
    dev.off()
END
R --vanilla <<END
    cov.table <- read.table("$COV")
    png(paste("$OUT",".A.l.png", sep=""), width=1000, height=1000)
    fragment.A <- cov.table[cov.table\$V1=='A',3]
    plot(fragment.A, type='l')
    dev.off()
    png(paste("$OUT",".B.l.png", sep=""), width=1000, height=1000)
    fragment.B <- cov.table[cov.table\$V1=='B',3]
    plot(fragment.B, type='l')
    dev.off()
END
  
