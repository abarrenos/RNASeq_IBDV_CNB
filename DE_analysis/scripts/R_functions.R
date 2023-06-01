## Functions required for RNA-Seq Data Processing and Differential Expression Analysis
## 
## Reads are aligned and features are quantified using Rsubread package.
## For DE Analysis, DESeq2 and edgeR-limma packages are used following the documentation.
## Some other tools and libraries are used for data clustering, gene annotation and
## and data visualizatuion.
##

##############################
#   Load Required Packages   #
##############################

library(ensembldb)		# needs to be first to avoid S4 inconsistencies
library(Rsubread)		# for read mapping
library(AnnotationHub)		# to seek annotation packages
library(AnnotationForge)	# to build our annotation package
library(GO.db)
library(PFAM.db)

library("biomaRt")		# an alternate approach to retrieve annotation

library(tibble)			# general tools
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(keypress)

library(edgeR)			# RNAseq with edgeR
library(limma)
library(RColorBrewer)
library(gplots)

library(DESeq2)			# RNAseq with DESeq2
library(IHW)			# for p-value adjustment with IHW
library(ggplot2)
library(ashr)

library(GSEABase)
library(fgsea)
library(reactome.db)
library(pathview)

library(gprofiler2)
library(clusterProfiler)
require(DOSE)
library(enrichplot)
library(europepmc)

library(cluster) 
library(factoextra)
library(fpc)
library("NbClust")

library(tcltk)
library(gWidgets2)


#################################################
#   Functions to Generate Directory Hierarchy   #
#################################################

create.out.dirs <- function(rnaseq.out, paired.end) {

	# Define the output folder of the analysis
	if (paired.end == T) {
		out.dir <<- paste(rnaseq.out, "paired_end", sep='/')	# Whether both ends should match or just any end
		logfile <<- paste(out.dir, "log", "RNAseq.log", sep='/')
	} else {
		paired.end <- F
		out.dir <<- paste(rnaseq.out, "single_end", sep='/')
		logfile <<- paste(out.dir, "log", "RNAseq.log", sep='/')
	}

	# Create needed directory hierarchy
	dir.create(rnaseq.out, showWarnings=FALSE)
	dir.create(file.path(rnaseq.out, "counts"), showWarning=FALSE)
        dir.create(out.dir, showWarnings=FALSE)
        dir.create(file.path(out.dir, "log"), showWarning=FALSE)
	dir.create(file.path(out.dir, "img"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/img"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/go"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/pfam"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/cluster"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/cluster/kmeans"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/cluster/pam"), showWarning=FALSE)
	dir.create(file.path(out.dir, "edgeR/cluster/dbscan"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/raw"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/shrunk"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/signif"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/go"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/pfam"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/img"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/cluster"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/cluster/kmeans"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/cluster/pam"), showWarning=FALSE)
	dir.create(file.path(out.dir, "DESeq2/cluster/dbscan"), showWarning=FALSE)
	return(list(out.dir = out.dir, logfile = logfile))
}


#########################################
#   Functions to Generate Otput Files   #
#########################################

# Create new log files
openlog <- function(name=logfile) {
    i <- 0
    name <- sprintf("%s.%03d", logfile, i)
    # check the existence of increasing numbers starting at 0 of log file versions
    while (file.exists(name)) {
        i <- i + 1
        name <- sprintf("%s.%03d", logfile, i)
    }
    # when we find a free number, use it
    # split=T means that the output will go to the screen and to a file
    sink(name, split=T, append=T)
}

closelog <- function() {
    # close ALL sinks. CAUTION CAUTION CAUTION
    while (sink.number()) { sink() }
}

as.png <- function(PLOT=NULL, file='out.png', width=1024, height=1024, overwrite=FALSE) {
    
    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite | ! file.exists(file)  ) {
        if (VERBOSE)
            cat("as.png(): creating", file, '\n')
        tryCatch( {
                png(file, width=width, height=height)
                print(PLOT)
            },
            finally=dev.off()
        )
    }
    return ()
}


###########################################
#   Functions to Interact with the user   #
###########################################

continue.on.key <- function() {
    if (! INTERACTIVE)
        return('')

    cat("\nPress any key to continue: ")
    key <- keypress()
    return(key)
}

continue.on.enter <- function(prompt='Press ENTER to continue: ') {
    if (! INTERACTIVE)
        return('')
    
    return(readline(prompt))
}

#library(keypress)
more.columns <- function (data, columns=c(1:dim(data)[2]), lines=20, header='', prompt="More? ") {
    maxln <- dim(data)[1]
    #print(dim(data))
    if (maxln < lines) lines <- maxln
    cont <- TRUE
    st <- 1; en <- lines
    #print(paste("start",st, "end", en, sep=' '))
    ans <- ''
    while (cont) {
	#print(paste(">>> start",st, "end", en, sep=' '))
        #cat(" ", data[ st:en, columns], '\n', sep='')
        #print(data[ st:en, columns])
        if (header != '') cat(header, '\n', sep='')
        for (i in st:en) {
            cat("<", i, "> ", sep="")
            for (c in columns) {
                cat('\t"', data[i, c], '"', sep='')
            }
            cat('\n')  
        }
        cat(prompt)
        ans <- keypress()
	#print(paste("|", ans, "|", sep=''))
        if (ans == ' ') {
            st <- st + lines
            en <- en + lines
	#print(paste("start",st, "end", en, sep=' '))
        } else if (ans == 'b' | ans == 'B') {
            st <- st - lines
            en <- en - lines
        } else if (ans == '>') {
            en <- maxln
            st <- en - (lines - 1)
        } else if (ans == '<') {
            st <- 1
            en <- lines
        } else if (ans == '' | ans == 'down') {
            st <- st + 1
            en <- en + 1
        } else if (ans == 'up') {
            st <- st - 1
            en <- en - 1
        } else if (ans == 'q' | ans == 'Q') {
            cont <- FALSE
        } else cont <- FALSE
        # check limits
        if (en > maxln) {
            en <- maxln
            st <- en - (lines - 1)
        }
        if (st < 1) {
            st <- 1
            en <- lines
        }
        cat('\n')
	#print(paste("start",st, "end", en, sep=' '))
    }
    return (ans)
}

show.data.frame <- function(df, window.name='data.frame', visible=TRUE) {
    #library(tcltk)
    #library(gWidgets2)
    window.name <- window.name
    window <- gwindow(title=window.name, visible=FALSE)
    tab <- gtable(df,
           container=window)
    if (INTERACTIVE)
        visible(window) <- TRUE	# setting it to FALSE removes window
    return(window)
}


##########################################
#   Functions to Load an Process Reads   #
##########################################

# Align fastQ reads into the reference genome
align.fastq <- function(fastq.dir, reference, aln.dir, save.dir=aln.dir) {
   
	# Get the fastq file names
    R1.fastq.files <- list.files(path=fastq.dir, pattern='R1', full.names=TRUE)
    R2.fastq.files <- list.files(path=fastq.dir, pattern='R2', full.names=TRUE)
    print(R1.fastq.files)
    print(R2.fastq.files)
    
    # we'll check if the output files exist to avoid repeating
    # work already done
	
	ref.fasta <- reference
	ref.name <- sub("\\.[[:alnum:]]+$", "", basename(reference))

	# Go to the alignment directory   
	wd <- getwd()	
	setwd(aln.dir)
	dir.create('ref', showWarnings=FALSE)
	system(paste('ln -sf `realpath ', ref.fasta, '` ./ref', sep=''))
	
	# Check wether sorted bam files already exist.
	if (length(list.files(pattern='.sorted.bam$', ignore.case=T, recursive=T)) > 0){
		bam.files <<- list.files(pattern='.sorted.bam$', ignore.case=T, recursive=T)
		cat('\nUSING EXISTING SORTED BAM FILES:')
		cat('\t', bam.files, sep='\n\t')
	
	# Check wether unsorted bam files exist.
	} else if (length(list.files(pattern='.bam$', ignore.case=T, recursive=T)) > 0){	
		bam.files <<- list.files(pattern='.bam$', ignore.case=T, recursive=T)
		cat('\nUSING EXISTING BAM FILES:\n')
		cat(paste(bam.files), sep='\n\t')
	} else {

	# If Bam files do not exist, perform alignment with Rsubread
		if (! file.exists(paste('ref/', ref.name, '.00.b.tab', sep=''))) {

			# build the reference index inside the 'aln.dir' directory
		    cat('\nBUILDING INDEX\n')
		    buildindex(basename=ref.name,reference=ref.fasta)
		    #dir()
		} else {
			cat('\nUSING INDEXED REFERENCE\n')
		}    
			
		if (! file.exists(paste(basename(R1.fastq.files[1]), '.subread.BAM', sep=''))) {
				
			# Align the reads
			# IMPORTANT NOTE: WE NEED TWO FILES LISTING ALL THE FASTQ FILES TO ALIGN
			#	R1.fastq.files and R2.fastq.files
				
			cat('\nALIGNING USING R_SUBREAD\n')
			align(index=paste('ref', ref.name, sep='/'), readfile1=R1.fastq.files,
				readfile2=R2.fastq.files)        
				
			# Align will generate the output in the fastq directory, we
			# so we move the alignment results to the output directory
			system(paste("mv ", fastq.dir, "/*.BAM .", sep=""))
			system(paste("mv ", fastq.dir, "/*.vcf .", sep=""))
			system(paste("mv ", fastq.dir, "/*.summary .", sep=""))
		}
	
		bam.files <<- list.files(pattern='.subread.bam$', ignore.case=T)
	}
    
	# Alignment Statistics:
	# Get the bam file names and inspect them to count the number
	# of reads that map to each genome position
    
	if ( ! file.exists('bam_files_stats.txt') & ! length(list.files(pattern = '.idxstats$', recursive = T)) == length(bam.files)) {
        
		cat('\nGENERATING BAM FILES STATISTICS\n')
        props <- propmapped(files=bam.files)
        props
        write.table(props, file='bam_files_stats.txt', 
    	    row.names=T, col.names=T, sep='\t')
    }
	cat('READ ALIGNMENT - DONE')
	setwd(wd)
}

# Compute the feature counts from alignment files using a reference GFF/GTF file
compute.feature.counts <- function(aln.dir, annotation, paired.end, count.dir=paste(rnaseq.out, "counts", sep="/")) {
	
    bam.files <- list.files(path=aln.dir, full.names=T, pattern='.bam$', ignore.case=T)
    ref.name  <- sub("\\.[[:alnum:]]+$", "", basename(annotation))
    ref.gtf   <- paste(ref.name, ".gtf", sep='')
    ref.gff   <- paste(ref.name, ".gff", sep='')
    
    count.files <- list.files(path = count.dir, pattern = '.counts$', full.names = T)
    
    if ( ! file.exists(paste(count.dir, 'featureCounts.csv', sep='/')) & ! length(count.files) == length(bam.files)) { 
    	
        cat('\nCOMPUTING FEATURE COUNTS\n')
        fc <<- featureCounts(files=bam.files, 
			annot.ext=ref.gtf, 
		        isGTFAnnotationFile=T,
			isPairedEnd=paired.end, 
		        requireBothEndsMapped=paired.end, 
		        primaryOnly=T, 
		        ignoreDup=T, 
		        useMetaFeatures=T)
		
	#   This means: process all BAM files
	#	Use as annotation the external file reference.gtf which is GTF
	#	BAM files contain paired end reads and we will only consider those
	#	where both ends match against the genome
	#    -------------------------------------------------------------
	#          ----->        <-----
	#	We will filter matches so that if a read may match more than one
	#	place in the genome, we will only consider the primary match and
	#	ignore other, duplicate matches
	#	We use meta-features
	#	We match against any feature in the GTF file, not only genes
	#	(this implies we will need to check the annotation carefully later)
	#	We do not remove chimeric fragments, but do actually count them too

	# SAVE FEATURE COUNTS
	# -------------------
	# now the variable fc contains 4 columns: annotation, target, counts and stats
	# but they exist only in the RAM memory, they are not stored somewhere safe,
	# so, we save them and create new variables so as to make it easier for us to 
	# manipulate the data
		write.csv(fc$counts, file=paste(count.dir, 'featureCounts.csv', sep='/'),
			  row.names=T, col.names=T, quote=F)
		write_delim(fc$stat, file=paste(count.dir, 'featureCounts_stat.txt', sep='/'), 
				delim='      ')

	# save all the contents of 'fc' in an RDS file and in an Rdata file

		saveRDS(fc, file=paste(count.dir, '/featureCounts.rds', sep=''))
		save(fc, file=paste(count.dir, '/featureCounts.RData', sep=''))
		
	# 'fc' can later be recovered with:
	# 		fc <- readRDS(file='featureCounts.rds')  
	# 		fc <- load(file='featureCounts.RData')
		
		cat('\nFEATURE COUNT FINISHED\n')
        
    } else {
    	
        cat(paste('\nUSING ALREADY EXISTING COUNTS (dir: ', count.dir, ')\n', sep=""))
    }
}

# Merge the different feature count files generated with HTSeq and read them into a single
# data frame.
merge.HTSeq.counts <- function(count.dir = paste(rnaseq.out, "counts", sep="/")) {

	count.files <- list.files('analysis.out/counts', pattern = '*.counts$', full.names=T)
	labels <- sub('\\..*$', '', basename(count.files))
	
	#cat('\nMerging count files:\n')
	#cat(count.files, sep='\n')

	# Each file is saved as an element of the list. Each element is a table with two coulmns:
	#	Gene_ID and "Sample name" (counts).

	count.list <- list()
	for (i in labels) {
		
		filename <- list.files('analysis.out/counts', pattern = i, full.names=T)
		count.list[[i]] <- read.table(filename, row.names='V1', header=F, stringsAsFactors=F)  #row.names='V1'
		colnames(count.list[[i]]) <- i
	}

	# Bind all htseq count files into a unique data frame by GeneID.
	all.counts <- bind_cols(count.list)
	
	# Save the merged dataframe into a file
	if (! file.exists(paste(count.dir, 'htseq_counts_merged.tab', sep='/'))){
		write.table(all.counts, paste(count.dir, 'htseq_counts_merged.txt', sep='/'))
	}
	return(all.counts)
}

# Plot counts-counts per million (cpm) correlation in order to define a count threshold.
counts.cpm.plot <- function(counts, cpm, save.png = TRUE, out.png = paste(out.dir, 'edgeR/img', 'edgeR_CPMcounts.png', sep = '/')) {

	plot(x = as.matrix(counts), y = as.matrix(my.cpm), xlim = c(0,20), ylim = c(0,5),
			pch = 1, cex = 0.8, xlab = 'Counts', ylab = 'Counts per Million (CPM)')
	abline(v = 10, col = 2, lwd = 2)
	abline(v = 15, col = 2, lwd = 2)
	#abline(h = 0.3, col = 2)
	arrows(10, 3, 15, 3, angle = 20, code = 3, col = 1, length = 0.2, lwd = 3, lty = 1)
	text(x=12.5, y=3.2, 'Expression\nthreshold', cex= 1)

	if (save.png) {
		as.png({	plot(x = as.matrix(counts), y = as.matrix(my.cpm), xlim = c(0,20), ylim = c(0,5),
							pch = 1, cex = 0.8, xlab = 'Counts', ylab = 'Counts per Million (CPM)')
					abline(v = 10, col = 2, lwd = 2)
					abline(v = 15, col = 2, lwd = 2)
					#abline(h = 0.3, col = 2)
					arrows(10, 3, 15, 3, angle = 20, code = 3, col = 1, length = 0.2, lwd = 3, lty = 1)
					text(x=12.5, y=3.2, 'Expression\nthreshold', cex= 1)
		}, out.png, overwrite = T)
	}
}

