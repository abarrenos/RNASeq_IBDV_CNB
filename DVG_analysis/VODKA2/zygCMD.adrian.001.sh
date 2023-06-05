#!/bin/bash

export PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin

maxproc=1	# maximum number of simultaneous processes
nproc=8 	# number of cores for each the alignments
proc_0=0 	# Set current number of running processes to 0

# Input files
#host_ref="ref/Gg_GRGc7b.fna"
viral_ref="ref/ibdv.fna"
fastq_dir="./unaligned_reads"
paired=true
type=unmapped
host=gg			# to refer to the unaligned fastq directory

# Scripts
#VODKA=~/work/jrrodriguez/birnavirus-ngs/src/VODKA/VODKA-master/scripts/VODKA
VODKA_dir=./VODKA2-main


##--------------------------------------------------##
## ALING READS AGAINS HOST (OPTIONAL) AND RUN VODKA ##
##--------------------------------------------------##

banner "vodka-2"

## Index Reference genomes if required

if [ ! -e $viral_ref.1.bt2 ] ; then
    echo -e "\nIndexing viral reference genome ..."
    bowtie2-build $viral_ref $viral_ref
fi

if [ ! -e $host_ref.1.bt2 ] && [ ! "$host_ref" == "" ]; then
    echo -e "\nIndexing host reference genome ..."
    bowtie2-build $host_ref $host_ref
fi

for i in $fastq_dir/*sample-*R1* ; do
#for i in ./unaligned_reads/DD_DF1-P_sample-1_R1_unmapped.fastq ; do

(	R1=$i
	R2=${i/R1/R2}
	
	sample=`basename $i | sed -e 's/_R1.*//g'`
	
	echo -e "\n\tPROCESSING SAMPLE $sample"

	# Create output directories and symbolic links 
	
	mkdir -p $sample $sample/vodka_scripts $sample/${sample}_vodka_results

	ln -sf $R1 $sample/R1_001.fastq
	ln -sf $R2 $sample/R2_001.fastq
	ln -sf $viral_ref $sample
    if [ ! "$host_ref" == "" ] ; then ln -s $host_ref $sample; fi
	
	
	###- RUN VODKA -###

    # Define output file basename
	outbase=${sample}_3000.150

	## 1) Generate VODKA2 multifata copy-back viral genome database
	if [ ! -e $sample/$outbase.fasta ] ; then 
        echo -e "Generating copy-back DVG database\n"
		echo -e "=================================\n"
		perl $VODKA_dir/genomefa_to_newfasta_cb_v2.pl $viral_ref 3000 150 $sample/$outbase.fasta 
	fi
	# ARGUMENTS
    # <genomefa> Strandard virus reference genome in FASTA format (1 line sequence)
    # <nb_from_right> Number of nt from right end of the sequence (use genome size for a full length analysis)
    # <read_length> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
    # <newfasta_name> Name of the multifata file to be generated. Should be as .<nb_from_right>.<read_length>.fasta
	
	
	## 2) Build bowtie2 index
 	mkdir -p $sample/bwt2-idx
	if [ ! -e $sample/bwt2-idx/$outbase.1.bt2l ] ; then
        echo -e "\nIndexing DVG database\n"
		echo -e   "=====================\n"
		bowtie2-build --large-index $sample/$outbase.fasta $sample/bwt2-idx/$outbase
		#ln -s $sample/bwt2-idx/$outbase* $VODKA_dir
	fi
	# ARGUMENTS
    # <newfasta_name> Name of the multifata file to be generated. Should be as .<nb_from_right>.<read_length>.fasta
    # <bt2_idx> Name of VODKA2 DB index for Bowtie2

	
	# Create input file list
	rm $sample/infiles
	#for f in ${sample}/R[12]_001.fastq ; do realpath -e -P $f >> $sample/infiles ; done
	for f in ${sample}/R[12]_001.fastq ; do echo $f >> $sample/infiles ; done
	
	## 3) Run VODAK2 analysis setup script
	
	echo -e "\nSetting up VODKA" ; echo "y" > yes # Automate interaction
	echo -e   "================\n"
	bash $VODKA_dir/VODKA2_analysis_setup.sh -f $sample/infiles -d $sample/bwt2-idx/$outbase \
	-p $sample -n 5 -v $viral_ref < yes

	mv *$sample*sh $sample/vodka_scripts	# Move vodka scripts to their respective folder
	
	
	## 3) Run VODKA2 analysis
	echo -e "\nRunning VODKA"
	echo -e   "=============\n"
	bash $sample/vodka_scripts/cbDVG_analysis_$sample.sh

) 2>&1 | tee log/${sample}_run_vodka.log &

done 2>&1 | tee log/do_all.log

wait ; rm yes

exit
	# Create input file list
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
	
	
	
	


done



for i in *[1-3] ; do
#for i in DF1-PC_sample-1 ; do
    echo "Processing $i"
    if [ -e $i/FINISHED ] ; then continue ; fi
    date
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
	
	
    touch $i/FINISHED
    date
do   echo "Processing $i"
    (
        cd $i
        # create input file list
        rm infiles
        touch infiles
        for f in ${i}_R[12]_001_${type}.fastq ; do realpath -e -P $f >> infiles ; done

      ###-- VODKA v.1 --### Only run 1 script
	  #
	  #  perl $VODKA \
      #      --genomeFA ../ref/ibdv.fna \
      #      --queryInput infiles \
      #      --outDir `pwd` \
      #      --bt2Dir /usr/bin/ \
      #      --webLogo ~/.local/bin/weblogo
      #  
        cd ..
    )
done
