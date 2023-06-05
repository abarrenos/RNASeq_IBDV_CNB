#!/usr/bin/bash

SAMTOOLS=NO
BEDTOOLS=NO
MOSDEPTH=NO
SAMBAMBA=YES

np=`nproc`
np=$((np / 2))

bam=$1
if [[ "$bam" == *"srt"* || "$bam" == *"sorted"* ]] ; then
    sorted_bam=$bam
else
    # we assume the bam file is unsorted
   echo "sorting $bam"
   samtools sort $bam -o ${bam%bam}srt.bam
   sorted_bam=${bam%bam}srt.bam
fi

if [ ! -e ${sorted_bam}.bai ] ; then
    echo "indexing $bam"
    samtools index $sorted_bam
fi

reads_prefix=${sorted_bam%.bam}

echo "Computing coverage for $reads_prefix"

mosdepth=~/work/jrrodriguez/birnavirus-ngs/src/mosdepth/mosdepth-master
sambamba=~/work/jrrodriguez/birnavirus-ngs/src/sambamba

# we no trust no stinkin' nobody !

if [ "$SAMTOOLS" = "YES" ] ; then
    if [ ! -e ${reads_prefix}_sam_cov.depth ] ; then
	# without -d0 this cuts depth over 8000, probably adjusting
	# those positions (maybe as gatk DepthOfCoverage?
	# with -d0 it counts all matches but mate-pair overlaps
	# are counted twice
	echo "Computing depth with samtools"
	samtools depth -d 0 ${sorted_bam} \
	> ${reads_prefix}_sam_cov.depth
    fi
fi


if [ "$BEDTOOLS" = "YES" ] ; then
    if [ ! -s ${reads_prefix}_bed_cov.depth ] ; then
	# gives results similar to samtools -d0 (likely counting
	# twice positions in mate-pair overlaps)
	echo "Computing depth with bedtools"
	bedtools genomecov -d -ibam ${sorted_bam} \
	> ${reads_prefix}_bed_cov.depth
    fi
fi


if [ "$MOSDEPTH" = "YES" ] ; then
    if [ ! -s ${reads_prefix}_mosdepth_cov.depth ] ; then
	# to get the same results as samtools and bedtools, we should
	# remove mate-pair overlap detection (effectively double-counting
	# mate-pair overlap regions).
	echo "Computing depth with mosdepth"
	$mosdepth/mosdepth ${reads_prefix} ${sorted_bam}
	zcat ${reads_prefix}.per-base.bed.gz \
	| cut -f1,3,4 > ${reads_prefix}_mosdepth_cov.depth
	python3 $mosdepth/scripts/plot-dist.py \
            ${reads_prefix}.mosdepth.global.dist.txt
	mv dist.html ${reads_prefix}.mosdepth.cum.cov.plot.html
    fi
fi


if [ "$SAMBAMBA" = "YES" ] ; then
    if [ ! -s ${reads_prefix}_sambamba_cov.depth ] ; then
	# without --fix-mate-overlaps, this gives results similar
	# to samtools and bedtools
	# with --fix-mate-overlaps, this gives results similar to
	# mosdepth
	echo "Computing depth with sambamba"
	$sambamba/sambamba depth base \
            --min-coverage=0 ${comment# same as -c 0 } \
            --fix-mate-overlaps ${comment# same as -m } \
            -o ${reads_prefix}_sambamba_cov.stats \
            ${sorted_bam} \
            >& /dev/null	# otherwise there will be a zillion warnings
	tail -n +2 ${reads_prefix}_sambamba_cov.stats\
	| while read r p0 d o ; do 
            echo "$r	$((p0+1))	$d"
	done > ${reads_prefix}_sambamba_cov.depth
    fi
fi


for meth in sam bed mosdepth sambamba ; do
    
    if [ ! -s ${reads_prefix}_${meth}_cov.depth ] ; then continue ; fi
    
    for seg in A B ; do
        # Plot using feedgnuplot
        if [ ! -s ${reads_prefix}_${meth}_cov.$seg.png ] ; then
            cat ${reads_prefix}_${meth}_cov.depth \
            | grep "^$seg" \
            | tr ' ' '	' \
            | cut -f2,3 \
            | feedgnuplot --lines --title "$meth $reads_prefix ($seg)" \
                           --unset grid \
                           --domain \
                           --ylabel "coverage depth" \
                           --xlabel "genome position (segment $seg)" \
                           --hardcopy ${reads_prefix}_${meth}_cov.$seg.png \
                           --exit
        fi
    done

    if [ ! -s ${reads_prefix}_${meth}_cov.AB.png ] ; then
        # Make a joint plot of both segments.
        # For this, we will generate an on the fly new data
        # file using 'join(1)' to merge each data set. We
        # use the field 2 (genome coordinate) as common [j]oin
        # field and include unpairable lines  ([a]ll values)
        # from files 1 and 2 using a value of 0 for [e]mpty 
        # values and produce as [o]utput the join field (0) 
        # and the 3rd column of each input data set (hence 
        # stripping the segment ID).
        # Since both datasets /for segments A and B) are in the
        # same file, we separate them using two subsidiary
        # commands with $(), each grep(1)ping for one segment.
        # Then we simply plot the data with feedgnuplot(1).
        # NOTE: for mosdepth we should probably remove "-a1 -a2"
        join -j2 -a1 -a2 -e0 -o0,1.3,2.3  \
        <( grep '^A' ${reads_prefix}_${meth}_cov.depth )  \
        <( grep '^B' ${reads_prefix}_${meth}_cov.depth )  \
        | feedgnuplot --lines --title "$meth $reads_prefix (A & B)" \
                      --unset grid \
                      --domain \
                      --ylabel "coverage depth" \
                      --xlabel "genome position" \
                      --hardcopy ${reads_prefix}_${meth}_cov.A.B.png \
                      --exit
   fi

    if [ ! -e ${reads_prefix}_${meth}_cov.depth.A.png ] ; then
        OUT_PFX=$reads_prefix
        echo "Plotting $meth coverage depth"
        echo "    using ${OUT_PFX}_${meth}_cov.depth"
        #ls -l ${OUT_PFX}_${meth}_cov.depth
        R --vanilla -q -s <<END
        t <- read.table("${OUT_PFX}_${meth}_cov.depth")
        # t$V1 contains the chromosome ID
        # t$V2 contains the position
        # t$V3 contains the depth at that position
        # get the list of chromosomes
        chromosomes <- levels(as.factor(t\$V1))
        # for each chromosome, plot its sequencing depth
        for (c in chromosomes) {
            png(paste("${OUT_PFX}_${meth}_cov.depth", c, "png", sep='.'), 
                width=1000, height=1000)
            d <- t[t\$V1==c,3]
            plot(d)
            dev.off()
            png(paste("${OUT_PFX}_${meth}_cov.depth", c, "l.png", sep='.'), 
                width=1000, height=1000)
            d <- t[t\$V1==c,3]
            plot(d, type='l')
            dev.off()
        }
END
    fi

done

exit
if [ ! -e ${reads_prefix}_bed_cov.depth.A.png ] ; then
    echo "Plotting bedtools coverage depth"
    OUT_PFX=$reads_prefix
    R --vanilla -q -s <<END
    t <- read.table("${OUT_PFX}_bed_cov.depth")
    # t$V1 contains the chromosome ID
    # t$V2 contains the position
    # t$V3 contains the depth at that position
    # get the list of chromosomes
    chromosomes <- levels(as.factor(t\$V1))
    # for each chromosome, plot its sequencing depth
    for (c in chromosomes) {
        png(paste("${OUT_PFX}_bed_cov.depth", c, "png", sep='.'), 
            width=1000, height=1000)
        d <- t[t\$V1==c,3]
        plot(d)
        dev.off()
        png(paste("${OUT_PFX}_bed_cov.depth", c, "l.png", sep='.'), 
            width=1000, height=1000)
        d <- t[t\$V1==c,3]
        plot(d, type='l')
        dev.off()
    }
END
fi


if [ ! -e ${reads_prefix}_unmapped_R1.fq ] ; then
    echo -n "Generating table of unmapped reads"
    # -f 4 means select the read if it is not mapped, even if its mate is
    # -f 12 means select if both the read and its mate are not mapped
    # -b means output as BAM file
    if [ ! -s $OUT_PFX.unmapped.bam ] ; then
    echo -n ":"
    samtools view -b -f 4 ${reads_prefix}_mapped_sorted.bam \
        > ${reads_prefix}_unmapped.bam
    fi
    if [ ! -s $$OUT_PFX.unmapped.sorted.bam ] ; then
        echo -n "."
        samtools sort -n ${reads_prefix}_unmapped.bam \
            -o ${reads_prefix}_unmapped_sorted.bam
    fi
    echo -n "_"
    bamToFastq -i ${reads_prefix}_unmapped_sorted.bam \
        -fq ${reads_prefix}_unmapped_R1.fq \
        -fq2 ${reads_prefix}_unmapped_R2.fq
    echo ""
fi
