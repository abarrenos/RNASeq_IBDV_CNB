#!/bin/bash

maxproc=5	# maximum number of simultaneous processes
nproc=8 	# number of cores for each the alignments
proc_0=0 	# Set current number of running processes to 0

# Input files
#host_ref="ref/Gg_GRGc7b.fna"
viral_ref="ref/ibdv.fna"
fastq_dir="./unaligned_reads"

# Scripts
ditector=scripts/DI-tector_06.py
plot_coverage=scripts/plot_coverage.sh
analyze_sample=scripts/analyze_sample.sh
analyze_group=scripts/analyze_group.sh



##------------------------------------------------------##
## ALING READS AGAINS HOST (OPTIONAL) AND RUN DI-TECTOR ##
##------------------------------------------------------##

## Index Reference genomes if required

if [ ! -s "$host_ref.amb" ] && [ ! "$host_ref" == "" ] ; then
    echo -e "\nIndexing host reference genome ..."
    bwa index -a bwtsw "$host_ref"
fi
if [ ! -s "$viral_ref.amb" ] ; then
    echo -e "\nIndexing viral reference genome ..."
    bwa index -a is "$viral_ref"
fi       
    
## Configure alignment against host genome

if [ "$host_ref" == "" ] ; then
   filter_host=""
else
   filter_host="-g $host_ref"
fi

## Create output directories

mkdir -p summary summary/dvg_counts summary/dvg_sequences summary/dvg_sequences summary/stats summary/vir_coverage_plots summary/wo_virus log

## Save sample names in a list
sample_list=()

banner "di-tector v0.6"

## Iterate through the fastq-filtered samples
for i in $fastq_dir/DD_*sample*_R1*; do

    R1=`realpath $i`
    R2=`realpath ${i/R1/R2}`

    if [ ! -s $fastq_in ] ; then
    	echo "File $i does not exist or its size is null"
    fi
    
    sample=`basename $i | sed -e 's/_R1.*//g'`
    group=`echo $sample | sed -s 's/_DF1.*//g'`    
    sample_list+=("$sample")
    mkdir -p $sample	# DI-tector output dir

     ## Create a Queuing system to run a maximum of maxproc simultaneous processes
    
    if (( proc_0++ > maxproc )); then
        
	# If we already have more processes than desired, 
        # wait for one to finish before running a new one
        
        wait -n   # wait for one job to complete. Bash 4.3
	
	(( proc_0-- ))	# Reduce by 1 the number of processes
    fi
    
     ## Run the multiple processes in parallel
    (  
       if [ ! -s ${sample}/${sample}_R1_counts.txt ]; then
           
           # Run DI-tector
	   echo -e "\nProcessing $group sample $sample R1 ..."
	   python3.8 $ditector $filter_host -o $sample\
	   -t "${sample}_R1" -d -k -p 0 -x $nproc \
	   $viral_ref $R1 2>&1 | tee log/${sample}_R1_di-tector.log

	   # Remove heavy temp files
	   rm ${sample}/${sample}_R1_temp_seqment.fq
	   rm ${sample}/${sample}_R1_temp_aln.sam
       fi
       
       if [ ! -s ${sample}/${sample}_R2_counts.txt ]; then
           
           # Run DI-tector
	   echo -e "\nProcessing $group sample $sample R2 ..."
	   python3.8 $ditector $filter_host -o $sample\
	   -t "${sample}_R2" -d -k -p 0 -x $nproc \
	   $viral_ref $R2 2>&1 | tee log/${sample}_R2_di-tector.log
	
	   # Remove heavy temp files
	   rm ${sample}/${sample}_R2_temp_seqment.fq
  	   rm ${sample}/${sample}_R2_temp_aln.sam
       fi
       
  ## Convert SAM files to BAM files to reduce size
       for sam in $sample/*.sam ; do
          samtools view -hb -o `echo $sam | sed -s 's/sam$/bam/g'` $sam && rm $sam
       done

  ## Plot Coverage Depth of the Viral Genome and save png
        
       echo -e "\nGenerating Coverage Plots ..."	
          
       bash $plot_coverage $sample/${sample}_R1_temp_file_Virus_Coverage.txt \
       $sample/${sample}_R1_coverage_depth > log/${sample}_R1_plot_coverage.log 2>&1
       
       bash $plot_coverage $sample/${sample}_R2_temp_file_Virus_Coverage.txt \
       $sample/${sample}_R2_coverage_depth > log/${sample}_R2_plot_coverage.log 2>&1
    )\
    2>&1 | tee $sample/$sample.log &
       
done

wait	# Wait for all the processes to finish

## Copy output files to the summary directory
       
ln -sf `realpath *sample*/*fasta.fa` summary/dvg_sequences/
ln -sf `realpath *sample*/*summary.txt` summary/stats/
ln -sf `realpath *sample*/*counts.txt` summary/dvg_counts/
ln -sf `realpath *sample*/*output_sorted.txt` summary/dvg_output
ln -sf `realpath *sample*/*coverage_depth*png` summary/vir_coverage_plots
ln -sf `realpath *sample*/*woVirus*.bam` summary/wo_virus

echo "### --- DI-TECTOR FINISHED! --- ###"


##--------------------------##
## ANALYZE DI-TECTOR OUTPUT ##
##--------------------------##     

banner "analisis"
		
# Analysing output files for each sample			     
for i in ${sample_list[@]} ; do

    $analyze_sample $i ;

done |& tee log/analyze_sample.log

# Analyze output files for each condition
for i in ${sample_list[@]} ; do

    if [[ $i =~ .*sample-1.* ]] ; then	
	condition=`echo $i | sed -e 's/_sample.*//g'`
    	echo $condition
	$analyze_group $condition
    fi
    
done |& tee log/analyze_group.log

echo "### --- ANALISIS FINISHED! --- ###"






