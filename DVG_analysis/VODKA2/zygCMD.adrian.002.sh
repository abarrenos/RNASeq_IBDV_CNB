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

#for i in $fastq_dir/*sample-*R1* ; do
for i in ./unaligned_reads/DD_DF1-P_sample-1_R1_unmapped.fastq ; do
(	
	R1=`realpath $i`
	R2=`realpath ${i/R1/R2}`	
	sample=`basename $i | sed -e 's/_R1.*//g'`
	
	echo -e "\n\tPROCESSING SAMPLE $sample"

	# Create output directories
	mkdir -p $sample $sample/vodka_scripts $sample/bwt2-idx #$sample/${sample}_vodka_results
	
	# Create symbolic links in output directory
	ln -sf $R1 $sample/${sample}_R1_001.fastq
	ln -sf $R2 $sample/${sample}_R2_001.fastq
	ln -sf `realpath $viral_ref` $sample/`basename $viral_ref`
    if [ ! "$host_ref" == "" ] ; then
	    ln -sf `realpath $host_ref` $sample/`basename $host_ref`
	fi
	ln -sf `realpath $VODKA_dir/*` $sample
	
	# Set the working directory	
	cd $sample
	
	###- RUN VODKA -###

    # Define output file basename
	outbaseCB="${sample}_CB.3000.150"
	outbaseDEL="${sample}_DEL.3000.150"

	## 1) Generate VODKA2 multifata copy-back viral genome database
	if [ ! -e $outbaseCB.fasta ] ; then 
        echo -e "Generating copy-back DVG database\n"
		echo -e "=================================\n"
		perl ../$VODKA_dir/genomefa_to_newfasta_cb_v2.pl `basename $viral_ref` 3000 150 $outbaseCB.fasta 
	fi
	# ARGUMENTS
    # <genomefa> Strandard virus reference genome in FASTA format (1 line sequence)
    # <nb_from_right> Number of nt from right end of the sequence (use genome size for a full length analysis)
    # <read_length> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
    # <newfasta_name> Name of the multifata file to be generated. Should be as .<nb_from_right>.<read_length>.fasta
	
	
	## 2) Build bowtie2 index
	if [ ! -e bwt2-idx/$outbaseCB.1.bt2l ] ; then
        echo -e "\nIndexing DVG database\n"
		echo -e   "=====================\n"
		bowtie2-build --large-index $outbaseCB.fasta bwt2-idx/$outbaseCB
	fi
	# ARGUMENTS
    # <newfasta_name> Name of the multifata file to be generated. Should be as .<nb_from_right>.<read_length>.fasta
    # <bt2_idx> Name of VODKA2 DB index for Bowtie2

	
	# Create input file list
	rm infiles
	#for f in ${sample}/R[12]_001.fastq ; do realpath -e -P $f >> $sample/infiles ; done
	for f in ${sample}_R[12]_001.fastq ; do echo $f >> infiles ; done
	
	## 3) Run VODAK2 analysis setup script
	
	#if [ ! -e cbDVG_analysis_$sample.sh ] ; then
		echo -e "\nSetting up VODKA" ; echo "y" > yes # Automate interaction
		echo -e   "================\n"
		
		bash ../$VODKA_dir/VODKA2_analysis_setup.sh -f infiles -d bwt2-idx/$outbaseCB \
		-p $sample -n 5 -v `basename $viral_ref` < yes ; rm yes

		#mv *$sample*sh vodka_scripts	# Move vodka scripts to their respective folder
	#fi

	## 3) Run VODKA2 analysis
	echo -e "\nRunning VODKA"
	echo -e   "=============\n"
	bash cbDVG_analysis_$sample.sh

	cd ..

) 2>&1 &#| tee log/${sample}_run_vodka.log &

done 2>&1 #| tee log/do_all.log

wait

date 

banner "done"

	

