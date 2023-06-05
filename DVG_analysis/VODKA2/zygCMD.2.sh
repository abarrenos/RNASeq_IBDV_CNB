#!/bin/bash

export PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin
VODKA=~/work/jrrodriguez/birnavirus-ngs/src/VODKA/VODKA-master/scripts/VODKA

type=filtered		# use this to employ ALL filtered reads			# containing both Host and Virus sequences
host=Ggallus
type=unmapped
host=gg			# to refer to the unaligned fastq directory


for i in ../fastq-qc/*[1-3] ; do mkdir -p `basename $i` ; done

for i in DF*[1-3] ; do 
    cd $i 
    if [ ! -e ibdv.fna ] ; then 
        ln -s ../ref/*fna . 
    fi
    cd .. 
done

if [ "$type" == "filtered" ] ; then
    for i in *-[1-3] ; do 
	cd $i 
	if [ ! -e ${i}_R1_001_${type}.fastq ] ; then
            #echo \-
            ln -s \
        	`realpath ../../fastq-qc/$i/IlluQC_${type}_files/R1_001.fastq.gz_${type}` \
        	${i}_R1_001_${type}.fastq 
	fi
	if [ ! -e ${i}_R2_001_${type}.fastq ] ; then
            #echo \-
            ln -s \
        	`realpath ../../fastq-qc/$i/IlluQC_${type}_files/R2_001.fastq.gz_${type}` \
        	${i}_R2_001_${type}.fastq 
	fi
	cd ..
    done

    cd ref
    if [ ! -e ${host}.fna.bwt ] ; then
	bwa index ${host}.fna
    fi
    cd ..
    hostopt='-g ../ref/${host}.fna'

elif [ "$type" == "unmapped" ] ; then
      for i in *-[1-3] ; do 
	cd $i 
	echo $i
	if [ ! -e ${i}_R1_001_${type}.fastq ] ; then
            #echo \-
            ln -s \
        	`realpath ../../${host}-unaligned/${i}_R1_${type}.fastq` \
        	${i}_R1_001_${type}.fastq 
	fi
	if [ ! -e ${i}_R2_001_${type}.fastq ] ; then
            #echo \-
            ln -s \
        	`realpath ../../${host}-unaligned/${i}_R2_${type}.fastq` \
        	${i}_R2_001_${type}.fastq 
	fi
	cd .. 
    done
  
    # no sense in using a host reference since host sequences are
    # already filtered out
    hostopt=''
else
    echo "unknown option"
fi


#for i in *[1-3] ; do
for i in DF1-PC_sample-2 ; do
    echo "Processing $i"
    (
        cd $i
        # create input file list
        rm infiles
        touch infiles
        for f in ${i}_R[12]_001_${type}.fastq ; do realpath -e -P $f >> infiles ; done

        perl $VODKA \
            --genomeFA ../ref/ibdv.fna \
            --queryInput infiles \
            --outDir `pwd` \
            --bt2Dir /usr/bin/ \
            --webLogo ~/.local/bin/weblogo
        
        cd ..
    )
done
