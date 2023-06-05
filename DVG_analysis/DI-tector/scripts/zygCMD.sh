#!/bin/bash

export PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin

di_tector=~/work/jrrodriguez/birnavirus-ngs/src/DI-tector/DI-tector_06.py

type=filtered		# use this to employ ALL filtered reads
			# containing both Host and Virus sequences
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


cd ref	# we always need a viral reference
if [ ! -e ibdv.fna.bwt ] ; then
    bwa index ibdv.fna
fi
cd ..


#for i in *-[1-3] ; do
for i in DF*-[1-3] ; do
#for i in DF1-P-MOI-1.0_sample-1 ; do
    if [ ! -d $i ] ; then continue ; fi
    echo "Processing $i"
    (
        cd $i

        echo "    ---> forward reads"
        if [ ! -e ${i}_R1_DI-tector_counts.txt ] ; then
	    echo "Running DI-tector on $i R1"
            python3 $di_tector -d -p 0 -k \
                ${hostopt} \
                -x 4  \
                -t ${i}_R1_DI-tector \
                ../ref/ibdv.fna \
                ${i}_R1_001_${type}.fastq \
                > ../log/$i.R1.log 2>&1
        fi
	
	# sanity check
	tmpfiles_kept=`ls ${i}_R1_DI-tector_temp_* | wc -l`
	if [ "$tmpfiles_kept" -gt 0 ] ; then
	    if [ -s ${i}_R1_DI-tector_temp_file_Virus_Coverage.txt ] ; then
	        # so we can safety check
		cp ${i}_R1_DI-tector_temp_file_Virus_Coverage.txt \
		   ${i}_R1_DI-tector_Virus_Coverage.depth
	        join -j2 -a1 -a2 -e0 -o0,1.3,2.3 \
		<( grep '^A' ${i}_R1_DI-tector_Virus_Coverage.depth ) \
		<( grep '^B' ${i}_R1_DI-tector_Virus_Coverage.depth ) \
		| feedgnuplot --lines \
		    --title "R1 (A & B)" \
		    --unset grid \
		    --domain \
		    --ylabel "coverage depth" \
		    --xlabel "genome position" \
		    --hardcopy coverage_R1.A.B.png
	    fi
	    # clean up to make room (some files may take too much space)
	    rm ${i}_R1_DI-tector_temp_*
	fi
        echo "    <--- reverse reads"
        if [ ! -e ${i}_R2_DI-tector_counts.txt ] ; then
	    echo "Running DI-tector on $i R2"
            python3 $di_tector -d -p 0 -k \
                ${hostopt} \
                -x 4  \
                -t ${i}_R2_DI-tector \
                ../ref/ibdv.fna \
                ${i}_R2_001_${type}.fastq \
                > ../log/$i.R2.log 2>&1
        fi
	
	# sanity check
	tmpfiles_kept=`ls ${i}_R2_DI-tector_temp_* | wc -l`
	if [ "$tmpfiles_kept" -gt 0 ] ; then
	    if [ -s ${i}_R2_DI-tector_temp_file_Virus_Coverage.txt ] ; then
	        # so we can safety check
		cp ${i}_R2_DI-tector_temp_file_Virus_Coverage.txt \
		   ${i}_R2_DI-tector_Virus_Coverage.depth
	        join -j2 -a1 -a2 -e0 -o0,1.3,2.3 \
		<( grep '^A' ${i}_R2_DI-tector_Virus_Coverage.depth ) \
		<( grep '^B' ${i}_R2_DI-tector_Virus_Coverage.depth ) \
		| feedgnuplot --lines \
		    --title "R1 (A & B)" \
		    --unset grid \
		    --domain \
		    --ylabel "coverage depth" \
		    --xlabel "genome position" \
		    --hardcopy coverage_R2.A.B.png
	    fi
	    # clean up to make room (some files may take too much space)
	    rm ${i}_R2_DI-tector_temp_*
	fi

        cd ..
    ) #> log/$i.log 2>&1 #& # XXX JR XXX
    #exit	# for debugging
done
echo -n "working... "
wait  # XXX JR XXX use if & backgrounding is eplicited above
echo "...finished"
