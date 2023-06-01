### MOST RECENT

#
# We start from data that has already been cleaned by performing a quality
# check with FastQC and subsequent edge trimming.
#
# The next step is to align the reads and calculate the counts of reads
# that map to each gene. This has been done as well previously, obtaining
# alignments in BAM format: the .bam files (binary files containing the
# reads aligned to the reference genome sequence) and the .bam.bai files
# containing the indexes for each one). This step is mostly a matter of
# time and has already been done:
# 
#   Paired-end Illumina short-reads were aligned against Coturnix japonica
#   genome (v2.0 primary assembly) using RNA-STAR (1) (--outReadsUnmapped
#   Fastx; --alignIntronMax 10000; -- alignMatesGapMax 10000). PCR and optical
#   duplicates were marked using the MarkDuplicates function of Picard-Tools 
#   (GATK)(2) (TAGGING_POLICY=ALL). Alignment results, saved as BAM files, 
#   were sorted and indexed using samtools (3).
# 
# We take the alignments produced by the Bioinformatics for Proteomics and
# Genomics Service of CNB directly, which are stored in the folder
# 'Alignments_Coturnix', and the reference genome used by them which is in
# the folder 'refGenomes/Coturnix_Japonica'. The reference files correspond
# to the 2.0 primary assembly. This is important for we will need their
# indexes for the next step.
#
# Now we need to calculate the counts of the reads that map to gene regions,
# summing up the reads that match each gene. For this, we use each .bam file
# and process it with the function featureCounts from R package Rsubread, and
# the reference genome data, which is in file Cjaponica.gtf. This last file
# contains the information about the features annotated in the genome of
# Coturnix japonica, including the coordinates of each feature in the reference
# genome. The function featureCounts will use these coordinates to know to
# which feature each read maps.
#

## Load required dependencies and functions for RNA-Seq analysis.

#setwd("/home/jr/work/adrian/rnaseq")
source("scripts/R_functions.R", echo = TRUE)
source("scripts/DESeq_functions.R", echo = TRUE)
source("scripts/Annotation_functions.R", echo = TRUE)

ALIGN <- FALSE					# Do we need to align reads or do we start from the aligned data?
HTSEQ <- TRUE					# Features counted using HTSEQ?
paired.end <- TRUE				# Paired-end reads?
use.online.annotation <- TRUE	# Use online ENSEMBL Db for annotation or create a Db from a GFF file?
VERBOSE <- TRUE
INTERACTIVE <- FALSE


# We will move inside './data/.' to work, so all paths are relative to it.
#reference <- 'Cjaponica'               # Coturnix japonica
#release <- 'Coturnix_japonica_2.0'
#target.organism <- 'Coturnix_japonica'
#ncbi.taxid <- '93934'
#ncbi.genus <- 'Coturnix'
#ncbi.species <- 'japonica'
#ncbi.version <- '2.0'
#aln.dir <- 'Alignments_Coturnix'	# Directory containing BAM files
#ens.db.pkg <- 'EnsDb.Cjaponica'
#ensembl.version <- '105'
#mart.name <- 'cjaponica_gene_ensembl'
#rnaseq.out <- 'analysis'				# Output directory of the analysis
#org.package <- 'org.Cjaponica.eg.db'
#KEGG_org <- 'cjo'

#reference <- "GRCg6a"
#release <- "GRCg6a"
#version <- "6a"	
#org.package <- "org.Ggallus.eg.db"

#reference <- 'Gg_GRGc7b'	# Gallus gallus
#ens.db.pkg <- 'EnsDb.Ggallus'

target.organism <- 'Gallus gallus'
genus <- 'Gallus'
species <- 'gallus'

release <- 'GRCg7b'
version <- '7b'
ensembl.version <- '105' #"106"

org.package <- 'org.Gg.eg.db'
ncbi.taxid <- '9031'

mart.name <- 'ggallus_gene_ensembl'
KEGG_org <- 'gga'
db.passwd <- '1_4_All+All_4_1'


#fastq.dir <- '/mnt/jr/work/jrrodriguez/birnavirus-ngs/gallus_gallus/gg-aln-Rsubread-7b/fastq-qc-paired'
#data <- fastq.dir
#aln.dir <- '~/work/adrian/rnaseq/test-aln'
#aln.dir <- '/mnt/jr/work/jrrodriguez/birnavirus-ngs/gallus_gallus/gg-aln-Rsubread-7b'
#reference <- '~/work/adrian/rnaseq/test-aln/Gg_GRGc7b.fna'



rnaseq.out <- 'analysis.out'
n.genes <- 1000		# Maximum number of top genes to revise

fastq.dir <- './fastq-qc-paired/ribo_depl'	# Directory containing FastQ files
reference.dir <- './reference_genome'
aln.dir <- './alignments'	# Directory containing BAM files
reference <- './reference_genome/Gg_GRGc7b.fna'
#annotation <- '/home/jr/work/adrian/reference_genome/Gg_GRGc7b.gff' # version106
annotation <- './annotation/ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf' # Version109



###############################################################################
###########################     ADRIAN     ####################################
###############################################################################


## GENERATE OUTPUT FILES AND DIRECTORIES
out <- create.out.dirs(rnaseq.out, paired.end)
out.dir <- out$out.dir
logfile <- out$logfile

align.fastq(fastq.dir, reference, aln.dir)

compute.feature.counts(aln.dir, annotation, paired.end)




# Hierarchical Clustering using hclust()
h.cluster.changes <- function(data.table, annotated.data, 
                              distance="euclidean",	# see dist()
                              clusters=1,
                              method="complete",	# see hclust()
                              estimate=T
                              ) {
                              
    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    # Calculate distance matrix  
    distan = dist(nor, method=distance)

    # We will try several clustering strategies

    # Hierarchical agglomerative clustering  
    cat("
    H I E R A R C H I C A L   C L U S T E R I N G
    =============================================
    \n\n")
    
    # We should consider using hcut() instead which cuts and allows fviz use
    
    hc = hclust(distan)
#    plot(mydata.hclust)
#    plot(mydata.hclust,hang=-1)
    print(plot(hc,labels=rownames(data.table),main='Default from hclust'))
    continue.on.enter("Press [ENTER] to continue ")
    
    # Cluster membership
    if ((clusters <= 1) & (estimate == T)) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    member <- cutree(hc, nclust)	# cut to 4 groups
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), " elements)\n")
        
        #View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
        #     title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
        
        # or, using tcltk and gWidget2
        #library(tcltk)
        #library(gWidgets2)
        data.to.show <- data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")]
        # clean up for showing
        data.to.show[ is.na(data.to.show) ] <- 'NA'
        window.name <- paste("Cluster no.", i, "(", length(clus.i), " elements)")
        #window <- gwindow(title=window.name, visible=TRUE)
        #tab <- gtable(data.to.show,
        #       container=window)
        window <- show.data.frame(data.to.show, window.name)
        visible(window) <- F
        #visible(window) <- TRUE	# setting it to FALSE removes window
#       keypress()
    }

    cat("
    
    Silhouette Plot for hierarchical clustering with normalized data
    ----------------------------------------------------------------
    Measure similarity of each object to its own cluster (cohesion) compared to
    other clusters (dispersion). Values range from -1 to +1. Large values
    indicate objects well matched to their own cluster and badly to neighboring
    clusters. If many points have low or negative value, the number of clusters
    is too low or too high.
    \n")
    print(plot(silhouette(cutree(hc, nclust), distan)))

    return(hc)
}


# Hierarchical Clustering using hcut()
hcut.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="hclust",		# see hcut()
                                    distance="euclidean", 	# see hcut()
                                    method="ward.D2",		# see hcut()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# ignored
                                    gap_bootstrap=500,
                                    estimate=FALSE
                                    ) {

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat("
    H c u t
    -------
    \n")

	# Estimate the optimal number of clusters
    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee\n")
        print(fviz_nbclust(nor, hcut, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, hcut, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=hcut, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }

    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    # Clustering (N groups)
    set.seed(123)
    hc <- hcut(nor, k=nclust, hc_func=algorithm, hc_method=method,
				hc_metric=distance, is_diss=FALSE)
    print(head(hc))
    member <- hc$cluster
    cat("Cluster membership counts\n")
    print(table(member))

    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        
		# sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(hc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(hc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



k.means.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="Hartigan-Wong",	# see kmeans()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# n. of random start sets to choose
                                    estimate=F,
                                    gap_bootstrap=500
                                    ) {

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat("
    K - m e a n s
    -------------
    \n")
    
    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        #
        # Scree Plot
        cat("

        Scree plot using normalized data

        This allows us to evaluate how much variation we account for as we
        consider more clusters and decide what a reasonable number of clusters
        might be.
        We draw here the variances accounted for using up to 20 K-means clusters
        \n") 
        # compute variances by row
        wss <- (nrow(nor)-1)*sum(apply(nor, by.row, var))
        for (i in 2:20) wss[i] <- sum(kmeans(nor, centers=i)$withinss)
        print(plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") )
        continue.on.enter("Press [ENTER] to continue ")

        # The scree plot will allow us to see the variabilities in clusters, 
        # we expect that if we increase the number of clusters, then the 
        # within-group sum of squares would come down. 
        #

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee\n")
        print(fviz_nbclust(nor, kmeans, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, kmeans, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=kmeans, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }

    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    # K-means clustering (N groups)
    kc <- kmeans(nor, centers=nclust, nstart=nstart, algorithm=algorithm)
    print(head(kc))
    member <- kc$cluster
    
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(kc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(kc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



pam.cluster.changes <- function(data.table, annotated.data,
                                    distance="euclidean", 	# see dist()
                                    metric="euclidean",		# see pam()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# n. of random start sets to choose
                                    estimate=T,
                                    gap_bootstrap=500
                                    ) {
                                    
    if (nstart != 1) medoids="random" else medoids=NULL
        
    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat("
    Partition around medoids
    ------------------------
    \n")

    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        cat("    Plot by within-cluster sums of squares\n\n")
        cat("    Elbow method: look at the knee\n")
        print(fviz_nbclust(nor, pam, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, pam, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=pam, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }


    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters

    # cluster with partition around medioids for k=N clusters 
    # (data is not dissimilarity but distance)
    # compute distance (we have dissimilarity to a common reference but not
    # between the considered groups).
    cat("Computing clusters (may take some time)...\n")
    eucldist <- dist(nor, method=distance) 
    #cluster
    ### JR ### NOTE: may be worth trying to cluster separately with diss=TRUE)
    pam.clus <- pam(eucldist, k=nclust, diss=TRUE)
    print(head(pam.clus))
    member <- pam.clus$clustering

    cat("Cluster membership information\n")
    print(pam.clus$clusinfo)
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))
    
    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot the partitioning
    print(clusplot(pam.clus, shade = FALSE,labels=F,
	    col.clus="blue",col.p="red",
            span=FALSE,
            main="PAM Cluster Mapping",cex=1.2))
    continue.on.enter("Press [ENTER] to continue ")

    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(pam.clus, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    # this seemingly ignores the data argument!!!
    print(fviz_cluster(pam.clus, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



cluster.changes <- function(data.table, annotated.data,
                            FUN=hcut,
                            clusters=1,			# n. of clusters
                            algorithm="default",	# see below	
                            distance="default", 	# see below
                            method="default",		# see see below
                            nstart=1,		# n. of ran dom start sets to choose
                            eps=1.0,		# epsilon for DBScan
                            gap_bootstrap=100,
                            normalize=TRUE,
                            estimate=TRUE,	# only if clusters > 1
                            output.folder=NULL	# NULL => no output desired
                                    ) {
# FUN -- one of c(hcut, kmeans, pam)
# algorithm -- for 'hcut' one of c(_"hclust"_, "agnes", "diana")
#              for 'kmeans' one of c(_"Hartigan-Wong"_, "Lloyd", "Forgy", "MacQueen")
#              for 'pam' one of c(_"original"_, "o_1", "o_2", "f_3", "f_4", "f_5", 
#                       "faster")
# distance -- for 'hcut' one of c(_"euclidean"_, "manhattan", "maximum", 
#                       "canberra", "binary", "minkowski", "pearson", "spearman", 
#                       "kendall")
#             for 'kmeans' it is ignored
#             for 'pam' one of c(_"euclidean"_, manhattan")
#             for 'dbscan' it is ignored
# method -- for 'hcut' one of c("ward.D"', "ward.D2", "single", "complete", 
#                       "average" (= UPGMA), "mcquitty" (= WPGMA), 
#                       "median" (= WPGMC), "centroid" (= UPGMC), 
#                       "weighted" (=WPGMA), "flexible", "gaverage")
#             for 'kmeans' ignored
#             for 'pam' ignored
#             for 'dbscan' ignored

    # Normalize the data
    by.row <- 1
    by.col <- 2
    if (normalize == TRUE) {
        means <- apply(data.table, by.col, mean)
        sds <- apply(data.table, by.col, sd)
        nor <- scale(data.table,center=means,scale=sds)
    } else {
        nor <- data.table
    }

    fun.name <- deparse(substitute(FUN))
    cat("
    C L U S T E R I N G    W I T H :   ", fun.name, "
    ---------------------------------------------
    \n")

    # re-activate next line to save computation time
    #if (clusters > 1) estimate <- FALSE		# we already know the number
    if (clusters < 1) estimate <- TRUE
    # we may have clusters == 1 and estimate == FALSE, e.g. to find outliers
    
    if (estimate == TRUE) {
        cat("Suggestions for the best number of clusters using", fun.name, "\n")

        if (normalize == TRUE) {
            cat("Plot by within-cluster sums of squares\n")
            cat("    Elbow method: look for a knee (normalized data)\n")
            print(fviz_nbclust(nor, FUN, method="wss"))
            out.png <- sprintf("%s/DESeq2_%s_wss.png",
                output.folder, fun.name)
            as.png(fviz_nbclust(nor, FUN, method="wss"), out.png, overwrite = TRUE)
            continue.on.enter("Press [ENTER] to continue ")
        }

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee (raw data)\n")
        print(fviz_nbclust(data.table, FUN, method="wss"))
        out.png <- sprintf("%s/DESeq2_%s_wss.png",
            output.folder, fun.name)
        as.png(fviz_nbclust(data.table, FUN, method="wss"), out.png, overwrite = TRUE)
        continue.on.enter("Press [ENTER] to continue ")

        if (fun.name != "dbscan") {
            cat("
            Average Silhouette Method

            The average silhouette approach measures the quality of a clustering. It
            determines how well each observation lies within its cluster.

            A high average silhouette width indicates a good clustering. The average
            silhouette method computes the average silhouette of observations for
            different values of k.
            \n")
            print(fviz_nbclust(nor, FUN, method="silhouette"))
            out.png <- sprintf("%s/DESeq2_%s_silhouette.png",
                output.folder, fun.name)
            as.png(fviz_nbclust(nor, FUN, method="silhouette"), out.png, overwrite = TRUE)

            continue.on.enter("Press [ENTER] to continue ")
        }

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (e.g.
        K-means, hierarchical clustering, partition around medoids).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        We are looking for the highest peak identified.

        \n")
        gap_stat <- clusGap(nor, FUN=FUN, K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        out.png <- sprintf("%s/DESeq2_%s_gap_stat.png",
            output.folder, fun.name)
        as.png(fviz_gap_stat(gap_stat), out.png, overwrite = TRUE)
        continue.on.enter("Press [ENTER] to continue ")
    }
    
    # Select the number of clusters
    if ((clusters < 1) & (fun.name != 'dbscan')) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else {
        nclust <- clusters
    }
    
    # Do the clustering (N groups). This is necessarily method-specific
    if (fun.name == "hcut") {
        if (algorithm == "default") algorithm <- 'hclust'
        if (distance == "default") distance <- "euclidean"
        if (method == "default") method <- "complete"
        
        cl <- hcut(nor, k=nclust, 
               hc_func=algorithm, hc_method=method, hc_metric=distance, is_diss=FALSE)
        print(head(cl))
        member <- cl$cluster
    } else if (fun.name == "kmeans") {
        if (algorithm == "default") algorithm <- "Hartigan-Wong"
        if (distance == "default") distance <- "euclidean"
        
        cl <- kmeans(nor, centers=nclust, nstart=nstart, algorithm=algorithm)
        print(head(cl))
        member <- cl$cluster
    } else if (fun.name == "pam") {
        if (algorithm == "default") algorithm <- "faster"
        if (distance == "default") distance <- "euclidean"
        
        # cluster with partition around medioids for k=N clusters 
        # (data is not dissimilarity but distance)
        # compute distance (we have dissimilarity to a common reference but not
        # between the considered groups).
        cat("Computing clusters (may take some time)...\n")
        use_dist_matrix <- FALSE	### JR ### fviz_cluster fails, why?
        if (use_dist_matrix == TRUE) {
            ### JR ### NOTE: may be worth trying to cluster with diss=TRUE)
            # calculate dissimilarity matrix and cluster it
            #distm <- get_dist(nor, method=distance, stand=TRUE) 
            cl <- pam(get_dist(nor, method=distance, stand=TRUE),
                      k=nclust, diss=TRUE, 
                      nstart=nstart, metric=distance, variant=algorithm)
        } else {
            cl <- pam(nor, k=nclust, diss=FALSE, 
                      nstart=nstart, metric=distance, variant=algorithm)
        }

        print(head(cl))
        member <- cl$cluster
        cat("Summary\n")
        print(cl$clusinfo)
    }  else if (fun.name == "dbscan") {
        if (algorithm == "defaut") algorithm <- "hybrid"
        if (distance == "default") distance <- "manhattan"
        cl <- dbscan(nor, eps=eps, MinPts=4, showplot=1)
        
	# showplot=1 makes it produce a movie plot
        continue.on.enter("Press [ENTER] to continue ")
        print(head(cl))
        member <- cl$cluster
        names(member) <- annotated.data$ensembl.gene.id
        nclust <- max(member)
        plot(cl, nor, main="DBScan")
        out.png <- sprintf("%s/DESeq2_%s_plot.png",
            output.folder, fun.name)
        as.png(plot(cl, nor, main="DBScan"), out.png, overwrite = TRUE)
        continue.on.enter("Press [ENTER] to continue ")
    }
        
    # this may take very long and is of little use for now
    if (FALSE) {
        cat("Distance plot\n")
        # plot distances
        distm <- get_dist(nor, distance, stand=TRUE)
        print(fviz_dist(distm,
                  gradient=list(low="blue", mid="white", high="red")))
        out.png <- sprintf("%s/DESeq2_%s_distances.png",
            output.folder, fun.name)
         as.png(fviz_dist(distm,
                    gradient=list(low="blue", mid="white", high="red")),
                    out.png, overwrite = TRUE)
        continue.on.enter("Press [ENTER] to continue ")
    }
    
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(data.table, list(member), mean))

    for (i in (min(member):max(member))) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id 
                                 %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        #View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
        #     title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
        d <- data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")]
        d[ is.na(d) ] <- 'NA'
        t <- paste("Cluster no.", i, "(", length(clus.i), ") elements")
        w <- show.data.frame(d, t, TRUE)
        visible(w) <- FALSE
        #visible(w) <- TRUE
#        keypress()
         clus.file <- sprintf("%s/DESeq2_%s_nc=%03d_c=%03d.tab",
             output.folder, fun.name, max(member), i)
         write.table(data.i, file=clus.file,
             row.names=T, col.names=T, sep='\t')
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(cl, data=nor, geom="point", ellipse.type="convex"))
    out.png <- sprintf("%s/DESeq2_%s_nc=%03d_PCA.png",
            output.folder, fun.name, max(member))
    as.png(fviz_cluster(cl, data=nor, geom="point", ellipse.type="convex"),
           out.png, overwrite = TRUE)
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    rcl <- cl
    if (fun.name == 'hcut') {
        cat("Setting names to gene names\n")
        dimnames(rcl$data)[[1]] <- rownames(gor)
    } else if (fun.name == 'pam') {
        names(cl$cluster <- rownames(gor))
    }
    print(fviz_cluster(rcl, data=gor, 
              ellipse.type="convex", show.clust.cen=FALSE, labelsize=6))
    out.png <- sprintf("%s/DESeq2_%s_nc=%03d_PCA_genes.png",
            output.folder, fun.name, max(member))
    as.png(fviz_cluster(rcl, data=gor, 
              ellipse.type="convex", show.clust.cen=FALSE, labelsize=6),
           out.png, overwrite = TRUE)
    continue.on.enter("Press [ENTER] to continue ")
    
    return(cl)
}


openlog(logfile)
# when a script ends or when stop() is called, there are a number
# of housekeeping tasks to do. on.exit() allows us to add additional
# tasks so that they, too, are executed at the end of the script
# this way we do not need to keep track of all the sink()s.
on.exit(closelog, add=T)




##############################################################################
#
#  PREPARE ANNOTATION SO IT IS READY WHEN WE GET THE RESULTS
#
##############################################################################
#  THIS NEEDS GENERALIZING ### JR ###
#---------------------------------------------------------------------------

# Now is time to prepare to connect all the results we get with the existing
# information  we know from the literature, so we will retrieve infos from
# ENSEMBL and connect with the genes we have kept. We are interested only in
# genes and transcriptomes (non-characterized genes) from the organism
# Gallus galllus. there is no database available for this organism so we
# have to create it by ourselves

# we create two different databases one from the gtf file (small archive)
# and one from the data we extracted from ENSEMBL and stored it in a sqlite 
# file

# we need to load the package that enables us to connect a file as an
# external database ensembldb. After that we create a variable for the DB
# storing it in the memory simple process: we extract the files from ENSEMBL
# and store the data to MySQL account and then with the script
# 'generate-EnsDBs.R' we create a sqlite file that has all the information
# needed to continue

# Define directory to store annotation files

annotation.dir <- './annotation'
target.organism <- 'Gallus gallus'
dir.create(annotation.dir, showWarnings=FALSE)

cat("\n================================================\n")
cat(  "     B U I L D I N G   A N N O T A T I O N\n")
cat(  "================================================\n")

prepareEnsemblDb <- function(	target.organism,
								annotation.dir,
								use.online.annotation = TRUE,
								annotation,
								ensembl.version){

	    ############################################
	    ## USE ONLINE ENSEMBL ANNOTATION DATABASE ##
	    ############################################
	    
	# SELECT ENSEMBL ANNOTATION DATABASE AVAILABLE ON-LINE
    if (use.online.annotation == TRUE) {
		
		cat("\n							\n")
		cat("\tFETCHING ENSEMBL DATABASE\n")
	    cat("\t=========================\n")

		# Create output directory
		net.ens.out <- paste(annotation.dir, 'net.ens.Db', sep = '/')
	    dir.create(net.ens.out, showWarnings=FALSE)

	    #if ( ! file.exists(paste(net.ens.out, 'net.ens.db.rds', sep='/')) ) {

		# Get ENSEMBL annotations for our target
            ah <- AnnotationHub()
            qr <- query(ah, c("EnsDb", target.organism))

		# Choose the most recent (last) one: AH98040 Ensembl Db version 105
		last.ref <- tail(names(qr), n = 1)
		net.ens.db <- qr[[last.ref]]

		# Save online ensemble annotation
		saveRDS(net.ens.db, file=paste(net.ens.out, 'net.ens.db.rds', sep='/'))
		save(net.ens.db, file=paste(net.ens.out, 'net.ens.db.RData', sep='/'))


        ## ERROR LOADING:
        ## Extarnal pointer is not valid. ¿Can online DBs be saved as objects? 
        #} else {
	    #
        #    net.ens.db <- readRDS(file=paste(net.ens.out, 'net.ens.db.rds', sep='/'))
        #}

	    # Check it
	    print(net.ens.db)
	    columns(net.ens.db)
	    #head(keys(net.ens.db, 'GENEID'))

	    return(net.ens.db)
    }


	    #######################################
	    ## BUILD ENSEMBL ANNOTATION DATABASE ##
	    #######################################

    ### BUILD DATABASE FROM GTF/GFF,

    if (use.online.annotation == FALSE) {
		
		cat("\n							\n")
		cat("\tBUILDING ENSEMBL DATABASE\n")
	    cat("\t=========================\n")

	    # Load required packages
        source(system.file("scripts/generate-EnsDBs.R", package = "ensembldb"))

        # Use our local GTF file to generate ENSEMBL Database
	    gtf <- paste(sub('\\.g(f|t)f$', '', annotation), 'gtf', sep='.')
	    gff <- paste(sub('\\.g(f|t)f$', '', annotation), 'gff', sep='.')
	    ref.base <- sub('\\.g(f|t)f$', '', basename(annotation))

	    # Create output directory	
        ens.pkg.out <- paste(annotation.dir, 'ens.pkg.Db', sep='/')
	    dir.create(ens.pkg.out, showWarnings = FALSE)

	    ##############
	    ##	ADRIAN  ##
	    ##############

	    # WARNING MESSAGES:
	    # I'm missing column(s): 'gene_name','entrezid'. The corresponding database column(s) will be empty!
	    # No column 'exon_id' present, created artificial exon IDs by concatenating the transcript ID and the exon number.
	    # Could not determine length for all seqnames.  Unable to retrieve sequence lengths from Ensembl.

	    sqlite <- paste(basename(gtf), 'sqlite', sep = '.')
	    #sqlite <- paste(target.organism, basename(gtf), ensembl.version, 'sqlite', sep='.')
	    ens.pkg.db <- paste(ens.pkg.out, sqlite, sep = '/') 
		
        if ( ! file.exists(ens.pkg.db) ) {

		    # Generate SQLite database in place from GTF
		    #gtf.edb <- ensDbFromGtf(gtf=gtf, 
	        #	    organism=target.organism,
            #       genomeVersion=basename(gtf),
            #       version=ensembl.version,
            #       destDir=ens.pkg.out)			# lacks entrezid

            ensDbFromGtf(	gtf = gtf, 
        				    organism= target.organism, 
        				    genomeVersion = basename(gtf),
        				    version = ensembl.version,
        				    outfile = ens.pkg.db)

		    # Move the database to the output directory
            # system(paste("mv", sqlite, ens.pkg.db))	

		    # Generate Ensembl DB Package
            makeEnsembldbPackage(	ensdb = ens.pkg.db,
        						    version = ensembl.version, 
    		             		    maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             		    author="J. R. Valverde",
                         		    destDir=ens.pkg.out,
                         		    license="Artistic-2.0")
	    }

	    # If building DB from GTF file fails, we may use the GFF3 file (which should contain
	    # the same information

	    ##############
	    ##  ADRIAN  ##
	    ##############

	    # DB FROM GFF DOES NOT WORK. ERROR MESSAGE:
	    # Required columns/fields gene_id;exon_id;biotype not present in the GFF file!

	    sqlite <- paste(basename(gff), 'sqlite', sep = '.')
	    #sqlite <- paste(target.organism, basename(gff), ensembl.version, 'sqlite', sep='.')
	    ens.pkg.db <- paste(ens.pkg.out, sqlite, sep = '/') 
		

        if ( ! file.exists(ens.pkg.db) ) {

		    #gff.edb <- ensDbFromGff(gff=gff, 
		    #	organism=target.organism,
		    #	genomeVersion=basename(gff),
		    #	version=ensembl.version,
		    #	destDir=ens.pkg.out)			# lacks entrezid

		    ensDbFromGff(	gff = gff, 
        				    organism= target.organism, 
        				    genomeVersion = basename(gff),
        				    version = ensembl.version,
        				    outfile = ens.pkg.db)

	    # Move the database to the output directory
            # system(paste("mv", sqlite, ens.pkg.db))	

	    # Generate Ensembl DB Package
            makeEnsembldbPackage(	ensdb = ens.pkg.db,
        		            	version = ensembl.version, 
    		             	    maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             	    author="J. R. Valverde",
                         	    destDir=ens.pkg.out,
                         	    license="Artistic-2.0")


            # Move the database to the output directory
            system(paste("mv", sqlite, ens.pkg.db))				

            makeEnsembldbPackage(ens.pkg.db, version=ensembl.version, 
							    maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
							    author="J. R. Valverde",
							    destDir=ens.pkg.out, license="Artistic-2.0")
        }
		return(ens.pkg.db)
    }
}
    ### BUILD DATABASE FROM LOCAL MYSQL DATABASE
    #
    #    # If building DB from GFF fails, we can still build it from a MySQL file
    #
    #	sqlite <- paste(target.organism, release, ensembl.version, 'sqlite', sep='.')
    #	ens.pkg.db <- paste(folder, '/', sqlite)
    #
    #   if ( ! file.exists(ens.pkg.db) ) {
    #
    #        # This shouln't be needed because we have already done it above
    #        # using the GTF/GFF3 file.
    #        
    #        # This is an alternate way to do generate the EnsDb package from
    #        # MySQL data downloaded from ENSEMBL
    #
    #        local.mysql.db <- paste("mysql/", tolower(target.organism), 
    #								"_core_", ensembl.version, "_6", sep='')
    #        #local.mysql.db <- "mysql/coturnix_japonica_core_104_2"
    #        #local.mysql.db <- "mysql/gallus_gallus_core_105_6"
    #
    #        createEnsDbForSpecies(ens_version = ensembl.version,
    #                user = "sci", pass=db.passwd,
    #                host = "localhost", 
    #                local_tmp=local.mysql.db, 
    #                #species="gallus_gallus", 
    #                species=target.organism,
    #                dropDb=FALSE)
    #
    #        # This should create a .sqlite file like ensDbFromGtf above,
    #        # named
    #        #	Gallus_gallus.GRCg6a.106.sqlite
    #        #       ^
    #        # with similar contents to the one produced from ensDBFromGtf/Gff;
    #        # if none of these does work, then tweak the process by hand in
    #        # source("scripts/generate-EnsDBs.R")
    #
    #        system(paste("mv", sqlite, ens.pkg.db))
    #
    #        makeEnsembldbPackage(ens.pkg.db, version=ensembl.version, 
    #						maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    #						author="J. R. Valverde",
    #						destDir=folder, license="Artistic-2.0")
    #    }


	###################################################
	## EXTRACT GENE ANNOTATION FROM ENSEMBL DATABASE ##
	###################################################
	
# Prepare the annotation of genes we will use.

annotateGenesFromENSEMBL <- function(ensembl.db, gene.names = '', gene.ids = '', save = TRUE, out.dir) {
	
	if (table(gene.names != '') & table(gene.ids == '')) {
		# Retrieve Information by GeneName
		ens.ann <- ensembldb::select(	ensembl.db, 
				      		keytype= 'GENENAME', keys = gene.names, 
				      		columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
						'GENENAME', 'GENEID', 'ENTREZID', 
						'TXNAME', 'TXBIOTYPE', 
						'PROTEINID', 'UNIPROTID'))
		
	} else if (table(gene.names == '') & table(gene.ids != '')) {
		# Retrieve Information by GeneID
		ens.ann <- ensembldb::select(	ensembl.db, 
						keytype= 'GENEID', keys = gene.ids, 
						columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
						'GENENAME', 'GENEID', 'ENTREZID', 
						'TXNAME', 'TXBIOTYPE', 
						'PROTEINID', 'UNIPROTID'))
				                  		
	} else { error("Provide either GeneIds OR GeneNames as input") }
	
	if (save == TRUE){

	# SAVE ANNOTATION
	# ---------------
		write.table(ens.ann, file = paste(out.dir, 'ensembl.annotation.txt', sep = '/'), 
					sep = '\t', row.names = T, col.names = T)

		# check if the amount of genes we have is the same as the number of 
		# the annotations that we have extracted
		
		if ( ! table(ens.ann$GENEID == gene.ids) ) {
			cat("Annotation does not match genes\n")
			cat("Using only one entry per gene\n")
			ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]
		
		} else {
			ens.ann.1 <- ens.ann
		}
		
		write.table(ens.ann.1, file = paste(out.dir, 'ensembl.annotation.uniq.txt', sep = '/'), 
			sep = '\t', row.names = T, col.names = T)
	}
	
	return(ens.ann)
}


ens.db <- prepareEnsemblDb(	target.organism = target.organism,
							annotation.dir = annotation.dir,
							use.online.annotation = TRUE,
							annotation = annotation,
							ensembl.version = ensembl.version)
ens.db

	
###################################################
###   A L T E R N A T E   A N N O T A T I O N   ###
###################################################


	######################################################
	## CREATE ORG PACKAGE FROM NCBI ANNOTATION DATABASE ##
	######################################################

org.package <- 'org.Gg.eg.db'

prepareOrgDb <- function(	org.package = org.package,
							annotation.dir = annotation.dir,
							use.online.annotation = TRUE,
							ncbi.taxid  = ncbi.taxid, # from NCBI Taxonomy browser (ncbi:txid9031)
							genus   = genus, 
							species = species, 
							version = version	){

	## RETRIEVE ONLINE ORG ANNOTATION VIA ANNOTATION HUB

	if (use.online.annotation == TRUE) {

		cat("\n						\n")
		cat("\tFETCHING ORG DATABASE\n")
		cat("\t=====================\n")

		# Use AnnotationHub to seek a suitable package
		ah <- AnnotationHub()
		qo <- query(ah, c('Orgdb', org.package))
		
		# Select the most recent annotation
		last.ref <- tail(names(qo), n = 1)
		org.db <- qo[[last.ref]]

		if ( ! exists('org.db')) {

	## RETRIEVE ONLINE ORG ANNOTATION VIA ORG PACKAGE
			
			if ( ! require(org.package, character.only = TRUE, quietly = TRUE)) {
				BiocManager::install(org.package) }
			
			# Character.only allows to require package form variable name
			require(org.package, character.only = TRUE)
			org.db <- get(org.package)
		}		
	}

	if ( use.online.annotation == FALSE | ! exists('org.db') ) {

		# Then try to build off-line annotation
		    
		#-----------------------------------------------
		# Create Org Package
		#-----------------------------------------------
		#
		# The datasets were first downloaded by hand from NCBI
		# and SwissProt into org.Gg.eg.db.
		# Then the following command had to be used:
		#
		
		# Create output directory
		org.out <- paste(annotation.dir, 'org.Db', sep = '/')
		ncbi.out <- paste(annotation.dir, 'org.Db', 'ncbi', sep = '/')
		dir.create(org.out, showWarnings=FALSE)	
		dir.create(ncbi.out, showWarnings=FALSE)		
		
		if ( ! file.exists(paste(ncbi.out, 'NCBI.sqlite', sep='/'))){
		
			makeOrgPackageFromNCBI(
									author  = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
									maintainer = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
									ncbi.taxid  = ncbi.taxid, # from NCBI Taxonomy browser (ncbi:txid9031)
									genus   = genus, 
									species = species, 
									version = version, 
									outputDir = org.out, 
									#NCBIFilesDir=".", rebuildCache=FALSE)
									NCBIFilesDir = ncbi.out #, 
									#rebuildCache=FALSE
			)
			
		    # We specify a directory to save locally the files used (and 
		    # retrieved from NCBI) to create the Org database.
		    # Normally, if the files are older than one day, they will be
		    # downloaded again. Since they are very large, the download may
		    # take too long and be interrupted frequently. This implies that
		    # if we are unable to download everything in a single day, we are
		    # doomed. The 'rebuildCache' option should avoid the re-downloading,
		    # but it is described as an option for "internal use only and for
		    # testing", and it may imply that no files are downloaded at all
		    # and only local files in the cache directory are used, so we will
		    # not use it unless strictly necessary.
		    # 
		    # GET
		    # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/*
		    # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ARCHIVE/gene2unigene
		    # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
		    # if that fails try using https://... instead  <<< PREFERRED!!!
		    # if that fails try curl

		    # We can download the files by hand first and then run this command
		    # so it uses the already downloaded files.
		    # install the package
		    
		    install.packages(org.package, repos = NULL)
		    
		    if (require(org.package, character.only = T)) {
		    	
		    	# org.package has already been loaded
				# Convert the packahe name to variable name and save the variable
		        
		        #org.db <- eval(parse(text = org.package))
		        org.db <- get(org.package)
		        columns(org.db)
		        
		    } else {
		    
		        org.db <- NULL
		    	cat("UNABLE TO GENERATE ORG ANNOTATION FROM NCBI")
		    
		    }
		    
		    # This fails in Gallus gallus due to download failures from NCBI.
		    # Should use the latest GitHub version installed with
		    # library(devtools)
		    # install_github("Bioconductor/AnnotationHub")
		}
	}
	return(org.db)
}

org.db <- prepareOrgDb(	org.package = org.package,
						annotation.dir = annotation.dir,
						use.online.annotation = TRUE,
						ncbi.taxid  = ncbi.taxid, # from NCBI Taxonomy browser (ncbi:txid9031)
						genus   = genus, 
						species = species, 
						version = version	)
org.db

## At this point ens.db should contain the ENSEMBL data and 
## org.db the Org type data.


# ---------------------------------------------------------------
# O B T A I N   B I O M A R T   A N N O T A T I O N
# ---------------------------------------------------------------

# Create output directory
bm.out <- paste(annotation.dir, 'biomaRT.Db', sep = '/')
dir.create(bm.out, showWarnings = FALSE)

### biomart.ds.name <- 'cjaponica_gene_ensembl'
biomart.ds.name <- 'ggallus_gene_ensembl'

getBiomaRtAnnotation <- function(biomart.ds.name, annotation.dir, save=TRUE) {

		cat("\n									\n")
		cat("\tFETCHING BIOMART ENSEMBL DATABASE\n")
		cat("\t=================================\n")

	bm.out <- paste(annotation.dir, 'biomaRT.Db', sep = '/')
	
	if ( ! file.exists(paste(bm.out, 'biomaRt.annotation.txt', sep='/')) ) {

		## Set up connection to ensembl database
		marts <- listMarts()
		head(marts)
		
		#bm.ensembl <- useMart("ensembl")
		bm.ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
		bm.ens.datasets <- listDatasets(bm.ensembl)
		head(bm.ens.datasets)
		
		## List the available datasets (species)
		listDatasets(bm.ensembl)
		
		## Get BiomaRt dataset name using RegExp
		desired.ds <- subset(bm.ens.datasets, grepl(biomart.ds.name, dataset))
		biomart.ds.name <- desired.ds$dataset		
		
		cat("\n Loading dataset", biomart.ds.name, "\n")
		if ( ! nrow(desired.ds) > 0) { stop("BiomaRt Dataset Not Found!") }

		## Use already known dataset name	
		mart.db <- useDataset(biomart.ds.name, mart = bm.ensembl)
		#mart.db <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = biomart.ds.name)
		
		## Check the available "attributes" - things you can retreive
		attributes <- listAttributes(mart.db)
		head(attributes)
		
		## Check the available "attributes" - things you can filter by
		filters <- listFilters(mart.db)
		head(filters, 15)

	# we cannot get all the annotation at once because it times out
	# full.annot <- getBM(attributes=
	#                       c("ensembl_gene_id", "ensembl_transcript_id", 
	#			   "start_position", "end_position", 
	#                          "chromosome_name", "gene_biotype", 
	#                          "description", 
	#                          "entrezgene_id", "entrezgene_accession", "entrezgene_description", 
	#                          "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003", 
	#                          "goslim_goa_accession", "goslim_goa_description", 
	#                          "pdb", 
	#                          "reactome", "uniprotswissprot"), 
	#                       mart=mart.db)

		# We will retrieve the data in pieces, including ensembl_gene_id in
		# each piece so we can use it as key for merging the annotation
		
		if ( ! file.exists( paste(bm.out, '/biomart.ensembl.tab', sep='')) ) {
		
			bm.ensembl.annot <- getBM(mart = mart.db, attributes = c(	
									"ensembl_gene_id", 
									"ensembl_transcript_id", 
									"start_position",
									"end_position", 
									"chromosome_name",
									"gene_biotype", 
									"description"))
			if (save == TRUE){
			write.table(bm.ensembl.annot, paste(bm.out, '/biomart.ensembl.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.ensembl.annot <- read.table(paste(bm.out, '/biomart.ensembl.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.entrez.tab', sep='')) ) {
			bm.entrez.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id", 
									"entrezgene_id",
									"entrezgene_accession",
									"entrezgene_description"))
			if (save == TRUE){
			write.table(bm.entrez.annot, paste(bm.out, '/biomart.entrez.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.entrez.annot <- read.table(paste(bm.out, '/biomart.entrez.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.go.tab', sep='')) ) {
			bm.go.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id", 
									"go_id",
									"name_1006",
									"definition_1006",
									"go_linkage_type",
									"namespace_1003"))
			if (save == TRUE){
			write.table(bm.go.annot, paste(bm.out, '/biomart.go.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.go.annot <- read.table(paste(bm.out, '/biomart.go.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.goslim.tab', sep='')) ) {
			bm.goslim.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id", 
									"goslim_goa_accession", 								"goslim_goa_description"))
			if (save == TRUE){
			write.table(bm.goslim.annot, paste(bm.out, '/biomart.goslim.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.goslim.annot <- read.table(paste(bm.out, '/biomart.goslim.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.fam.tab', sep='')) ) {
			bm.fam.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id",
									"pfam",
									"pirsf",
									"prints",
									"tigrfam"
									))
			if (save == TRUE){
			write.table(bm.fam.annot, paste(bm.out, '/biomart.fam.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.fam.annot <- read.table(paste(bm.out, '/biomart.fam.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.prosite.tab', sep='')) ) {
			bm.prosite.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id",
									"scanprosite",
									"pfscan" 
									))
			if (save == TRUE){
			write.table(bm.prosite.annot, paste(bm.out, '/biomart.prosite.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.prosite.annot <- read.table(paste(bm.out, '/biomart.prosite.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.sfam.tab', sep='')) ) {
			bm.sfam.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id",
									"superfamily"
									))
			if (save == TRUE){
			write.table(bm.sfam.annot, paste(bm.out, '/biomart.sfam.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.sfam.annot <- read.table(paste(bm.out, '/biomart.sfam.tab', sep=''), 
			header=T, sep='\t')
		}

		if ( ! file.exists(paste(bm.out, '/biomart.extra.tab', sep='')) ) {
			bm.extra.annot <- getBM(mart=mart.db, attributes=c(
									"ensembl_gene_id", 
									"pdb",
									#"reactome", 
									"uniprotswissprot"))
			if (save == TRUE){
			write.table(bm.extra.annot, paste(bm.out, '/biomart.extra.tab', sep=''), 
			row.names=T, col.names=T, sep='\t')
			}
		} else {
			bm.extra.annot <- read.table(paste(bm.out, '/biomart.extra.tab', sep=''), 
			header=T, sep='\t')
		}

		# Now that we have all the pieces, merge them all together
		# into a single annotation variable
		#bm.annot <- merge(bm.ensembl.annot, bm.entrez.annot, by="ensembl_gene_id")
		#bm.annot <- merge(bm.annot, bm.go.annot, by="ensembl_gene_id")
		#bm.annot <- merge(bm.annot, bm.goslim.annot, by="ensembl_gene_id")
		#bm.annot <- merge(bm.annot, bm.extra.annot,  by="ensembl_gene_id")

		# SAVE ANNOTATION
		# ---------------
		#write.table(bm.annot, file=paste(bm.out, '/biomaRt.annotation.txt', sep=''), 
		#    sep='\t', row.names=T, col.names=T)

	} else {
		bm.ensembl.annot <- read.table(paste(bm.out, '/biomart.ensembl.tab', sep=''), 
			header=T, sep='\t')
		bm.entrez.annot <- read.table(paste(bm.out, '/biomart.entrez.tab', sep=''), 
			header=T, sep='\t')
		bm.go.annot <- read.table(paste(bm.out, '/biomart.go.tab', sep=''), 
			header=T, sep='\t')
		bm.goslim.annot <- read.table(paste(bm.out, '/biomart.goslim.tab', sep=''), 
			header=T, sep='\t')
		bm.fam.annot <- read.table(paste(bm.out, '/biomart.fam.tab', sep=''), 
			header=T, sep='\t')
		bm.prosite.annot <- read.table(paste(bm.out, '/biomart.prosite.tab', sep=''), 
			header=T, sep='\t')
		bm.sfam.annot <- read.table(paste(bm.out, '/biomart.sfam.tab', sep=''), 
			header=T, sep='\t')
		bm.extra.annot <- read.table(paste(bm.out, '/biomart.extra.tab', sep=''), 
			header=T, sep='\t')

	# this is too big to use
	#bm.annot <- read.table(file=paste(bm.out, '/biomaRt.annotation.txt', sep=''), 
	#    sep='\t', header=T)
	}


	if ( ! file.exists(paste(bm.out, 'biomart.annotation.1st.tab', sep='/')) ) {
	    
    # One possible way to do it would be to filter the queries above
    # to retrieve the annotation matching ensembl_ids
    #
    # We can set a field to use to filter the output data
    # Set the filter type and values
    #ourFilterType <- "ensembl_gene_id"
    # and the values to select from that field
    #filterValues <- rownames(fit)
    #
    # and then obtain the specified annotation from records that match the values
    # specified in the filter field
    #fit.bm.extra.annot <- getBM(attributes=c(
    #                       "ensembl_gene_id", 
    #                       "pdb",
    #                       "reactome", 
    #                       "uniprotswissprot"), 
    #                   mart=mart.db,
    #                   filters=ourFilterType,
    #                   values=filterValues)
    #                   
    # deduplicate selecting the first annotation
    # fit.bm.extra.annot.1 <- fit.bm.extra.annot[ ! duplicated(fit.bm.extra.annot$ensembl_gene_id), ]
    # then we would repeat this for each annotation subset and merge all of them
    # at the end...

		# We will eduplicate everything first and match aftwerards
		bm.ensembl.annot.1 <- bm.ensembl.annot[ ! duplicated(bm.ensembl.annot$ensembl_gene_id), ]
		bm.entrez.annot.1 <- bm.entrez.annot[ ! duplicated(bm.entrez.annot$ensembl_gene_id), ]
		bm.go.annot.1 <- bm.go.annot[ ! duplicated(bm.go.annot$ensembl_gene_id), ]
		bm.goslim.annot.1 <- bm.goslim.annot[ ! duplicated(bm.goslim.annot$ensembl_gene_id), ]
		bm.extra.annot.1 <- bm.extra.annot[ ! duplicated(bm.extra.annot$ensembl_gene_id), ]

		bm.annot.1 <- merge(bm.ensembl.annot.1, bm.entrez.annot.1, by = "ensembl_gene_id")
		bm.annot.1 <- merge(bm.annot.1, bm.go.annot.1, by = "ensembl_gene_id")
		bm.annot.1 <- merge(bm.annot.1, bm.goslim.annot.1, by = "ensembl_gene_id")
		bm.annot.1 <- merge(bm.annot.1, bm.extra.annot.1,  by = "ensembl_gene_id")

		# we save it to avoid repeating this in the future
		if (save == TRUE){	
			write.table(bm.annot.1, file = paste(bm.out, '/biomart.annotation.1st.tab', sep = ''), 
				sep = '\t', row.names = T, col.names = T)
			write.table(bm.annot.1, file = paste(annotation.dir, '/biomart.annotation.1st.tab', sep = ''), 
				sep = '\t', row.names = T, col.names = T)
		}
	    # or even do it all at once?
	    ##bm.annot.1 <- bm.annot[ ! duplicated(bm.annot$ensembl_gene_id), ]
	    ##write.table(bm.annot.1, 
	    ##            file=paste(folder, '/biomart.annotation.1st.txt', sep=''), 
	    ##            sep='\t', row.names=T, col.names=T)

	} else {
	    bm.annot.1 <- read.table(
	    		file=paste(bm.out, '/biomart.annotation.1st.tab', sep=''), 
				sep='\t', 
		        header=T)
	}

	return(bm.annot.1)
}

biomart.ann <- getBiomaRtAnnotation(biomart.ds.name = biomart.ds.name,
									annotation.dir = annotation.dir,
									save = TRUE)
head(biomart.ann)


# as a byproduct we know at this point annotation files have been saved

# And now we are ready with our reference info at hand...


#################################################################################
#
# Get the alignments and feature counts
#
#################################################################################

# Note: since we will be looking at infected cells with aberrant RNA viral genomes
# which might present recombinations of virus-host reads, we will restrict ourselves
# to matches in both ends.

if ( ALIGN == TRUE) {

	align.fastq(fastq.dir = fastq.dir, reference = reference, aln.dir = aln.dir, save.dir = aln.dir)

	compute.feature.counts(	aln.dir = aln.dir,
							annotation = annotation,
							paired.end = paired.end,
							count.dir = 'analysis.out/counts')
}

# Alignment files
bam.files <- list.files(path = aln.dir, pattern = '.bam$', full.names = TRUE, ignore.case = TRUE)


### ADRIAN ###

if(HTSEQ == FALSE){

	#### If Rsubread is used:
	# One single file will be generated with the counts for all samples

	fc <- readRDS(file=paste(out.dir, '/featureCounts.rds', sep=''))

	#   get the feature counts
	#	we use EnsEMBL annotation from release 6a
	#	this will tag the position counts and match them against the 
	#	genes annotated in the GFF3 file for the reference, obtaining
	#	a list of genes and the number of reads that mapped to them
	#	as an indicator of their expression level

}else{

	#### If HTSeq is used:
	# One count file will be generated per sample. We read the files into a list
	# and convert it into a single data frame.

	fc <- merge.HTSeq.counts(count.dir = 'analysis.out/counts')
}


#################################################################
#                                                               #
#################################################################

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     E D G E R
# ---------------------------------------------------------------

## LOAD FEATURE COUNTS

if (HTSEQ == FALSE){
	countDataFile <- paste('analysis.out/counts',
					'featureCounts.csv', sep = '/')
	countData <- read.csv(countDataFile,
					row.names=1) %>%
					as.matrix()
} else {
	countDataFile <- paste('analysis.out/counts',
					'htseq_counts_merged.tab', sep = '/')	

	countData <- read.table(countDataFile,
				row.names=1) %>%
				as.matrix()
}

# Discard genes with no counts in any of the samples

counts <- countData[rowSums(countData)>1, ]
head(counts)
dim(counts)

## LOAD SAMPLE INFORMATION### ADRIAN ### --- ALREADY READ BEFORE INTO VARIABLE TARGET

sampleData <- read.table("SampleInfo.tab", header = TRUE)
sampleData$Sample <- relevel(as.factor(sampleData$Sample), "DF1")
sampleData$CellStatus <- relevel(as.factor(sampleData$CellStatus), "Parental")
sampleData$VirusStatus <- relevel(as.factor(sampleData$VirusStatus), "Mock")

#target <- read.table('SampleInfo.tab', header=TRUE, row.names=1)
target <- sampleData
print(target)


# Then we proceed to the analysis.. for this, we will need other packages:
library(edgeR)
library(limma)
library(RColorBrewer)
library(gplots)

	########################################
	## FIRST FILTERING BY COUNT THRESHOLD ##
	########################################

# The next steps are to distinguish the genes whose expression is significant
# from  the ones that have an 'insignificant' expression that could be by
# chance or irrelevant to the  situation of the cell. So we set up the
# threshold of expression for each gene to a minimum 10-15 counts. Since
# we will work with normalized CPM (counts per million) data, we need
# to check which number of CPM corresponds to  to this amount of counts in
# order to filter the data

# we'll plot the correlation of cpm and counts to see which is the number of
# cpm that corresponds to 10-15 counts minimum. we see this information
# graphically we use for example column 1... we could check every file
# independently but since they are all similar in number of reads there 
# should be no need
# convert counts to CPM

my.cpm <- cpm(counts)

out.png <- paste(out.dir, 'edgeR/img/edgeR_CPMcounts.png', sep = '/')

# Plot counts-counts per million (cpm) correlation in order to define a count threshold.
counts.cpm.plot <- function(counts, cpm, out.png = NULL) {

	plot(x = as.matrix(counts), y = as.matrix(my.cpm), xlim = c(0,20), ylim = c(0,5),
			pch = 1, cex = 0.8, xlab = 'Counts', ylab = 'Counts per Million (CPM)')
	abline(v = 10, col = 2, lwd = 2)
	abline(v = 15, col = 2, lwd = 2)
	#abline(h = 0.3, col = 2)
	arrows(10, 3, 15, 3, angle = 20, code = 3, col = 1, length = 0.2, lwd = 3, lty = 1)
	text(x=12.5, y=3.2, 'Expression\nthreshold', cex= 1)

	#if (! is.null(out.png)) {
		as.png({	plot(x = as.matrix(counts), y = as.matrix(my.cpm), xlim = c(0,20), ylim = c(0,5),
							pch = 1, cex = 0.8, xlab = 'Counts', ylab = 'Counts per Million (CPM)')
					abline(v = 10, col = 2, lwd = 2)
					abline(v = 15, col = 2, lwd = 2)
					#abline(h = 0.3, col = 2)
					arrows(10, 3, 15, 3, angle = 20, code = 3, col = 1, length = 0.2, lwd = 3, lty = 1)
					text(x=12.5, y=3.2, 'Expression\nthreshold', cex= 1)
		}, out.png, overwrite = TRUE)
	#}
}

counts.cpm.plot(counts = counts, cpm = my.cpm, out.png = out.png)

# we will set up the threshold to 0.5 according on the plot we have drawn
# this command will return a table of trues and falses, then, we want to keep
# only the genes or features that exceed the threshold at least in three 
# different samples (experiments .bam files)

thres <- my.cpm > 0.5
keep <- rowSums(thres) >= 3
table(keep)


# Note that we could as well have used the counts instead of cpm to filter,
# yet it would be less approximate and yield a different selection

# then we store in a different variable the genes whose counts exceed the
# threshold and visualise the content of the new variable to see the amount of
# remaining genes

counts.keep <- counts[keep,]
dim(counts.keep)

##################################################################################
if ( FALSE) {
    # some paranoid manual checks
    # we are using counts instead of cpm
    cpmavg <- data.frame(vd000.0=apply(counts.keep[,13:15],1,mean), 
                         vd000.1=apply(counts.keep[,1:3],1,mean),
                         vd001.0=apply(counts.keep[,10:12],1,mean), 
                         vd010.0=apply(counts.keep[,7:9],1,mean), 
                         vd100.0=apply(counts.keep[,4:6],1,mean) 
                         )

    f.c000.0 <- data.frame(vd000.0_vd000.1=(cpmavg$vd000.0 / cpmavg$vd000.1),
                           vd000.0_vd001.0=(cpmavg$vd000.0 / cpmavg$vd001.0),
                           vd000.0_vd010.0=(cpmavg$vd000.0 / cpmavg$vd010.0),
                           vd000.0_vd100.0=(cpmavg$vd000.0 / cpmavg$vd100.0)
                          ) 
    row.names(f.c000.0) <- row.names(cpmavg)
    l.f.c000.0 <- log2(f.c000.0)
    write.table(f.c000.0, file=paste(out.dir, '/edgeR/hand.fc_000.0.tab', sep=''), sep='\t')
    write.table(l.f.c000.0, file=paste(out.dir, '/edgeR/hand.lfc_000.0.tab', sep=''), sep='\t')
    write.table(cpmavg, file=paste(out.dir, '/edgeR/hand.cpmavg.tab', sep=''), sep='\t')
}
##################################################################################


	############################################
	## STATISTICS USING DELIST OBJECT (EDGER) ##
	############################################

# We have manipulated the data discarding whatever is not of high interest
# and now we need to see the differencial expression and highlight the
# differences among the cells

# Convert the counts.keep to a DGEList and define sample groups.

dge <- DGEList(counts.keep)
dge$samples$group <- target$Sample

# Perform TMM normalization
dge <- calcNormFactors(dge)
dge$samples

# Plot the library size if the different samples.
out.png <- paste(out.dir, '/edgeR/img/edgeR_sample_lib_size.png', sep='')
as.png(barplot(dge$samples$lib.size, cex.names= 1, main = "Library Size", col = dge$samples$group, names.arg=dge$samples$group, ylab = "Reads"), out.png, overwrite=TRUE)

# Now, do some quality control plots, barplots and boxplots.
# We need normalized data counts so we take the logarithm
# we need to group together all the experiment data that correspond to the
# same  cell type (target) and asign a different color to each one so as to
# distinguish them graphically; we know we have 4 groups so we have to specify
# 4 different colors

logcpm <- cpm(dge$counts, log=TRUE)

### ADRIAN ### Target$PFU as factor.

#group.col <- c('blue', 'green', 'yellow', 'purple', 'red')[factor(target$PFU[1:15])] # 1:15 Do not consider Mock

group.col <- c('blue', 'green', 'yellow', 'purple', 'red', 'brown')[dge$samples$group]

out.png <- paste(out.dir, 'edgeR/img/edgeR_log2_cpm.png', sep='/')
as.png( {
        par(mfrow=c(1,1))
        boxplot(logcpm, xlab='', ylab='Log2 counts per million', 
		col = group.col, las=2, outline=FALSE)
        abline(h=median(logcpm), col='red')
        title('Boxplots of logCPMs unnormalised')
    	}, out.png, overwrite=TRUE
)

# then we check the MDS plot to see any significant difference between the
# groups
out.png <- paste(out.dir, '/edgeR/img/edgeR_mds_plot.png', sep='')
as.png( {
        par(mfrow= c(1,1))
        plotMDS(dge, col=group.col)
    }, out.png, overwrite=FALSE)


# now is time to see the differences in the expression (variance). We have
# to apply a funcion that calculates the variance by rows (genes) and then
# retrieve the n.genes most DE genes


### ADRIAN ### --- Why we calculate log counts again. See logcpm.
logcounts <- cpm(dge, log =TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:n.genes]
highly_var <- logcounts[select_var,]
dim(highly_var)

# we now set up the colours we will use for the heatmap plot 
# and then we create all the colors in between in the palette

mypalette <- brewer.pal(11, 'RdYlBu')
morecolors <- colorRampPalette(mypalette)

# we plot the heatmap.2 (gplots) without a line (trace), scale by row
# (difference in color) margins (something about the labels used), also we
# reverse the colors because by default the red is associated with low
# expression and we are not familiar with it
out.png <- paste(out.dir, '/edgeR/img/edgeR_heatmap.png', sep='')
as.png( {
        #margins <- par("mar")
        #par(mar=c(25, 5, 5, 10))
        heatmap.2(highly_var, 
                  col= rev(morecolors(50)), 
                  trace='none', 
                  ColSideColors=group.col,
                  scale='row', 
                  margins= c(15,6))
        #par(mar=margins)
    }, out.png)

# This plot is of limited use. We'd better have other names for rows and 
# columns and plot the n first (most variable) to see them well

### ---- ADRIAN ---- ### 	This is not working

n <- 50
high_var <- highly_var
colnames(high_var) <- gsub("_R1.bam", "", colnames(highly_var))
name <- ensembldb::select(ens.db, keys=rownames(high_var), 
                   column='GENEID', keytype='GENEID', 
                   columns=c('GENENAME'))
out.png <- paste(out.dir, '/edgeR/img/edgeR_heatmap.', n,'.png', sep='')
as.png( {
        margins <- par("mar")
        par(mar=c(25, 5, 5, 10))
        heatmap.2(high_var[1:n,1:15], 
                      col=rev(morecolors(50)), 
                      trace='none', 
                      ColSideColors=group.col,
                      scale='row', 
                      margins= c(15,5),
                      labRow=name[1:n, 2])
        par(mar=margins)
    }, out.png )


# We can automate the estimation of the dispersion and add it to the dge object
dge <- estimateCommonDisp(dge)

# And now we can estimate gene-wise dispersion estimates allowing for a
# possible trend with average count size. These will allow us to use a GLM
# instead of a plain LM.
# This gives us the BCV (Biological Coefficient of Variation) between samples
dge <- estimateGLMTrendedDisp(dge)
dge <- estimateTagwiseDisp(dge)

out.png <- paste(out.dir, '/edgeR/img/edgeR_BCV_dispersions.png', sep='')
as.png(plotBCV(dge), out.png)

############## JR #########################
# 
# Chapter 3 
# 
# Specific experimental designs 3.1 Introduction In this chapter, we outline
# the principles for setting up the design matrix and forming contrasts for
# some typical experimental designs.
# 
# Throughout this chapter we will assume that the read alignment, normalization
# and dispersion estimation steps described in the previous chapter have
# already been completed. We will assume that a DGEList object y has been
# created containing the read counts, library sizes, normalization factors and
# dispersion estimates.
# 
# 3.2 Two or more groups
# 
# 3.2.1 Introduction
# 
# The simplest and most common type of experimental design is that in which a
# number of experimental conditions are compared on the basis of independent
# biological replicates of each condition. Suppose that there are three
# experimental conditions to be compared, treatments A, B and C, say. The
# samples component of the DGEList data object might look like:
# 
# > y$samples
# group lib.size norm.factors
# Sample1 A 100001 1
# Sample2 A 100002 1
# Sample3 B 100003 1
# Sample4 B 100004 1
# Sample5 C 100005 1
# 
# Note that it is not necessary to have multiple replicates for all the
# conditions, although it is usually desirable to do so. By default, the
# conditions will be listed in alphabetical order, regardless of the order that
# the data were read:
# 
# > levels(y$samples$group)
# [1] "A" "B" "C"
# 29
# 
# edgeR User’s Guide
# 
# 3.2.2 Classic approach
# 
# The classic edgeR approach is to make pairwise comparisons between the
# groups. For example,
# 
# > et <- exactTest(y, pair=c("A","B"))
# > topTags(et)
# 
# will find genes differentially expressed (DE) in B vs A. Similarly
# 
# > et <- exactTest(y, pair=c("A","C"))
# 
# for C vs A, or
# 
# > et <- exactTest(y, pair=c("C","B"))
# 
# for B vs C.
# 
# Alternatively, the conditions to be compared can be specified by number, so
# that
# 
# > et <- exactTest(y, pair=c(3,2))
# 
# is equivalent to pair=c("C","B"), given that the second and third levels of
# group are B and C respectively.
# 
# Note that the levels of group are in alphabetical order by default, but can
# be easily changed.
# 
# Suppose for example that C is a control or reference level to which
# conditions A and B are to be compared. Then one might redefine the group
# levels, in a new data object, so that C is the first level:
# 
# > y2 <- y
# > y2$samples$group <- relevel(y2$samples$group, ref="C")
# > levels(y2$samples$group)
# [1] "C" "A" "B"
# 
# Now
# 
# > et <- exactTest(y2, pair=c("A","B"))
# 
# would still compare B to A, but
# 
# > et <- exactTest(y2, pair=c(1,2))
# 
# would now compare A to C.
# 
# When pair is not specified, the default is to compare the first two group
# levels, so
# 
# > et <- exactTest(y)
# 
# compares B to A, whereas
# 
# > et <- exactTest(y2)
# 
# compares A to C.
# 
# 
# 
# 3.2.3 GLM approach
# 
# The glm approach to multiple groups is similar to the classic approach, but
# permits more general comparisons to be made. The glm approach requires a
# design matrix to describe the treatment conditions. We will usually use the
# model.matrix function to construct the design matrix, although it could be
# constructed manually. There are always many equivalent ways to define this
# matrix. Perhaps the simplest way is to define a coefficient for the
# expression level of each group:
# 
# > design <- model.matrix(~0+group, data=y$samples)
# > colnames(design) <- levels(y$samples$group)
# > design
# A B C
# Sample1 1 0 0
# Sample2 1 0 0
# Sample3 0 1 0
# Sample4 0 1 0
# Sample5 0 0 1
# attr(,"assign")
# [1] 1 1 1
# attr(,"contrasts")
# attr(,"contrasts")$group
# [1] "contr.treatment"
# 
# Here, the 0+ in the model formula is an instruction not to include an
# intercept column and instead to include a column for each group.
# 
# One can compare any of the treatment groups using the contrast argument of
# the glmQLFTest or glmLRT function. For example,
# 
# > fit <- glmQLFit(y, design)
# > qlf <- glmQLFTest(fit, contrast=c(-1,1,0))
# > topTags(qlf)
# 
# will compare B to A. The meaning of the contrast is to make the comparison
# -1*A + 1*B + 0*C, which is of course is simply B-A.
# 
# The contrast vector can be constructed using makeContrasts if that is
# convenient. The above comparison could have been made by
# 
# > BvsA <- makeContrasts(B-A, levels=design)
# > qlf <- glmQLFTest(fit, contrast=BvsA)
# 
# One could make three pairwise comparisons between the groups by
# 
# > my.contrasts <- makeContrasts(BvsA=B-A, CvsB=C-B, CvsA=C-A, levels=design)
# > qlf.BvsA <- glmQLFTest(fit, contrast=my.contrasts[,"BvsA"])
# > topTags(qlf.BvsA)
# > qlf.CvsB <- glmQLFTest(fit, contrast=my.contrasts[,"CvsB"])
# > topTags(qlf.CvsB)
# > qlf.CvsA <- glmQLFTest(fit, contrast=my.contrasts[,"CvsA"])
# > topTags(qlf.CvsA)
# 
# which would compare B to A, C to B and C to A respectively.
# 
# 
# Any comparison can be made. For example,
# 
# > qlf <- glmQLFTest(fit, contrast=c(-0.5,-0.5,1))
# 
# would compare C to the average of A and B. Alternatively, this same contrast
# could have been specified by
# 
# > my.contrast <- makeContrasts(C-(A+B)/2, levels=design)
# > qlf <- glmQLFTest(fit, contrast=my.contrast)
# 
# with the same results.
# 
############## JR #########################

# we know that the variability we see in the expression depends on infection so
# we have to take this into account and create a model
# 0+ forces the design to include all groups and not have
# an intercep (reference) column

#design.column <- "Src"
#design.column <- "Infected"
#design.column <- "viral.dose"
design.column <- "Sample"

#design <- model.matrix(~0 + target[ , design.column])	# no reference
design <- model.matrix(~ target[ , design.column])
colnames(design) <- levels(as.factor(target[ , design.column]))
rownames(design) <- rownames(target)

# This is intended to be used by the Gallus gallus experiment
#design.persist.inf <- model.matrix(~ 0 + target$Persistent * target$Infected)
#design.persist.inf <- model.matrix(~ target$Persistent * target$Infected)
#design.persist.inf

# IF WE DO NOT INCLUDE THE "~ 0 + " IN THE FORMULA, MODEL.MATRIX()
# WILL USE AS REFERENCE THE FIRST ALPHABETICAL ORDER LEVEL !!!
# WHICH HERE WOULD BE 0.1 !!!
# Another problem is that we are limited to the comparisons defined by
# coefficents, if we want more control we need to use limma::makeContrasts
#
# If we use as formula "~ table[ , design.column ]" then we are limited to 
# comparisons defined by the fitting coefficients (ref vs. variable-in-coeff).
# If we use "~ 0 + table[ , design.column]" we would have diffs for
# all levels, but no reference.


# Let us test for differential expression with the DGE data using a GLM:
# first, fit genewise GLMs
gfit <- glmFit(dge, design)
names(gfit)
head(coef(gfit))

# We can now conduct Likelihood Rato tests and show the top genes
# for the selected comparison coefficients (reference vs. coeff)
lrt <- glmLRT(gfit, coef=1)	# coef = 1... length(gfit$coefficients)
topTags(lrt)


# let's do a VOOM analysis with the formula employed
# in the current design

DO_VOOM=TRUE
if (DO_VOOM) {
    #voom transformation of the data 
    v <- voom(dge, design, plot=FALSE)
    out.png <- paste(out.dir, '/edgeR/img/edgeR_voom.png', sep='')
    as.png( {
            par(mfrow= c(1,1))
            voom(dge, design, plot=TRUE)
        }, out.png )



    # Carry on a summary variation analysis using the VOOM transformed data
    #
    # we fit the results of the voom depending on the model we created with our
    # design, we overwrite the fit with the use of eBayes (stat model), and
    # depending on that we decide abut which tests fit the best and summarise the
    # results: topTable gives us helpful information about everything

    v.fit <- lmFit(v, design)
    v.fit <- eBayes(v.fit)
    results <- decideTests(v.fit)
	cat("\nResults:\n")    
	summary(results)
    #topTable(fit, coef= 1, sort.by='p')

    # SAVE TOP 'N' GENES
    # ------------------
    # save the table of the n.genes top expressed genes sorted in various orders
    eR.save.top <- function(fit, out.dir, n.genes=500, sort.by='p', coef=1, p.value=0.01) {
        # default adjustment is BH
        tt <- topTable(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(out.dir, '/edgeR/cmp_coef=', coef, '_top_', n.genes, '_by', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top(v.fit, out.dir, n.genes, 'p')
    eR.save.top(v.fit, out.dir, n.genes, 'B')
    eR.save.top(v.fit, out.dir, n.genes, 'logFC')
    eR.save.top(v.fit, out.dir, n.genes, 'AveExpr')


    # the problem here is that we are limited to the comparisons defined by
    # coefficents, if we want more control we need to use makeContrasts


    #---------------------------------------------------------------------------

    #Now is time to connect all the results we have with the existing
    #information  we know from the literature, so we will retrieve infos from
    #ENSEMBL and connect with the genes  we have kept. We are interested only in
    #genes and transcriptomes (non-characterized genes) from the organism
    #Coturnix japonica. there is no database available for this organism so we
    #have to create it by ourselves

    # we created two different databases one from the gtf file (small archive)
    # and one from the data we extracted from ENSEMBL and stored it in a sqlite 
    # file

    # we'll use ens.db from above

    # As long as it works perfectly, we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma") which is a
    # list. I.e. we assign the annotation to a new element named 'genes'  
    # of this list.
    #	This works if they are in the same order
    #fit$genes <- ens.ann.1
    #	This checks they go in the same order)

############################-- ADRIAN --##############################################
### ANOTATE GENES ###
#####################

#	ALREADY DONE BEFOREHAND
#    ens.db <- prepareEnsemblDb(	target.organism = target.organism,
#							    annotation.dir = annotation.dir,
#							    use.online.annotation = TRUE,
#							    annotation = annotation,
#							    ensembl.version = ensembl.version)
#    
#    biomart.ann <- getBiomaRtAnnotation(biomart.ds.name = biomart.ds.name,
#										annotation.dir = annotation.dir,
#										save = TRUE)

    ens.ann <- annotateGenesFromENSEMBL(	ensembl.db = ens.db,
				                    	    gene.names = rownames(dge),
				                    	    save = TRUE,
				                    	    out.dir = annotation.dir)

    ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENENAME), ]

    bm.annot.1 <- biomart.ann[ ! duplicated(biomart.ann$entrezgene_accession), ]

    #v.fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1$GENEID), ]
    v.fit$genes <- ens.ann.1[ match(rownames(dge), ens.ann.1$GENEID), ]
        
    # SAVE THE ANNOTATED FIT
    # ----------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(v.fit, file=paste(out.dir, '/edgeR/annotatedVOOMfit.rds', sep=''))
    save(v.fit, file=paste(out.dir, '/edgeR/annotatedVOOMfit.RData', sep=''))

    # SAVE TOP 'N' GENES ANNOTATED
    # ----------------------------
    # we can now run again the topTable and we will have all the annotation 
    # information linked 
    #	n.genes <- 500 ALREADY DEFINED ABOVE
    eR.save.top.annotated <- function( 	fit,
										out.dir,
										n.genes=500,
										sort.by='p',
										coef=1,
										p.value=0.01) {
        # default adjustment is BH
        tt <- topTable(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(out.dir, '/edgeR/cmp=', coef, '_top_', n.genes, '_annotated_by_', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')
    }

    eR.save.top.annotated(v.fit, out.dir, n.genes, sort.by='p')
    eR.save.top.annotated(v.fit, out.dir, n.genes, sort.by='logFC')
    eR.save.top.annotated(v.fit, out.dir, n.genes, sort.by='AveExpr')


    #now create the volcano plot for only the top 1/2 genes
    out.png <- paste(out.dir, '/edgeR/img/edgeR_volcanoplot.png', sep='')
    as.png(volcanoplot(v.fit, highlight=n.genes/2, coef=1, names=fit$genes$SYMBOL),
        out.png)

    #Testing relative to a threshold: 1 means a 2x fold change
    threshold=1
    v.fit.thres <- treat(v.fit, lfc=threshold)
    res.thres <- decideTests(v.fit.thres)
    summary(res.thres)
    topTreat(v.fit.thres, coef=1, sort.by='p')

    # SAVE TOP RESULTS RELATIVE TO THRESHOLD
    # --------------------------------------
    #	n.genes <- 500 ALREADY DEFINED ABOVE
    eR.save.top.ann.thresh <- function(fit,
										out.dir,
										n.genes=500,
										sort.by='p',
										coef=1,
										p.value=0.01,
										threshold=1) {

        # default adjustment is BH
        tt <- topTreat(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(out.dir, '/edgeR/cmp=', coef, '_top_', n.genes, '_annotated_lfc>=', threshold, '_by_', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top.ann.thresh(v.fit.thres, out.dir, n.genes, 'p')
    eR.save.top.ann.thresh(v.fit.thres, out.dir, n.genes, 'logFC')
    eR.save.top.ann.thresh(v.fit.thres, out.dir, n.genes, 'AveExpr')


    # and now we need to annotate 'fit$genes' with the biomaRt data
    v.fit$genes <- merge(v.fit$genes, bm.annot.1, by.x="GENEID", by.y="ensembl_gene_id")
    head(v.fit$genes)

    v.fit.thres$genes <- merge(v.fit.thres$genes, bm.annot.1, by.x="GENEID", by.y="ensembl_gene_id")
    head(v.fit.thres$genes)


    # save fit data with the extended version containing an extra
    # annotated bucket:

    # save all the contents of 'fit' in an RDS file
    saveRDS(v.fit, file=paste(out.dir, '/edgeR/annotatedVOOMfit+.rds', sep=''))
	save(v.fit, v.fit.thres, file=paste(out.dir, '/edgeR/annotatedVOOMfit+.RData', sep=''))

    # save top N genes, annotated-
    eR.save.top.bm.ann.thresh <- function(	fit,
											out.dir,
											n.genes=500,
											sort.by='p',
											coef=1,
											p.value=0.01,
											threshold=1) {

        # default adjustment is BH
        tt <- topTable(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(out.dir, '/edgeR/cmp=', coef, '_top_', n.genes, '_bm.annotated_lfc>=', threshold, '_by_', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top.bm.ann.thresh(v.fit.thres, out.dir, n.genes, 'p')
    eR.save.top.bm.ann.thresh(v.fit.thres, out.dir, n.genes, 'B')
    eR.save.top.bm.ann.thresh(v.fit.thres, out.dir, n.genes, 'logFC')
    eR.save.top.bm.ann.thresh(v.fit.thres, out.dir, n.genes, 'AveExpr')

}  # if (VOOM == TRUE)


# -----------------------------------------------------------------
# Getting beyond here is easier with DESeq2
# -----------------------------------------------------------------
# 
# # redefine the design so we can make any kind of comparison by
# # using a null reference ("0 + ")
# design.column <- "viral.dose"
# design <- model.matrix(~0 + target[ , design.column])
# design
# 
# # fit now using a negative binomial model with Quasi-Likelihood fit
# # which is similar to glmFit but also estimates QL dispersion values
# # calling limma squeezeVar().
# qlfit <- glmQLFit(dge, design, robust=TRUE, abundance.trend=TRUE)
# 
# # And now we can conduct likelihood ratio tests and show the top genes
# lrt <- glmLRT(qlfit, coef=1)
# topTags(lrt)
# 
# exps <- levels(as.factor(target[ , design.column]))
# n.exps <- length(levels(as.factor(target[ , design.column])))
# exps
# # Pairwise comparison using contrasts
# #	This requires that we know the number of experiments and
# # design the contrasts vector
# # PFU
# #"0.1" "1"   "10"  "100" "wt" 
# pfu_0.1_wt <- glmQLFTest(qlfit, contrast=c(1,0,0,0,-1))
# # viral.dose
# # "vd_000.0" "vd_000.1" "vd_001.0" "vd_010.0" "vd_100.0"
# vd_000.1_000.0 <- glmQLFTest(qlfit, contrast=c(-1,1,0,0,0))
# 
# # pairwise comparison using formula
# # this requires we loop over all the 'exps' against all other 'exps'
# # but requires that all the levels in the design.column are valid
# # R identifiers.
# # for viral.dose
# cmp_0.1_wt <- makeContrasts(vd_000.1-vd_000.0, levels=design)
# topTags(glmQLFTest(qlfit, contrast=cmp_0.1_wt))
# 

# ---------------------------------------------------------------------
# Do all comparisons at once

eR.annotate.fit <- function(fit, end.db, biomart) {

    ens.ann <- annotateGenesFromENSEMBL(	ensembl.db = ens.db,
				                gene.names = rownames(fit),
				                save = TRUE,
				                out.dir = annotation.dir)

    # use only first annotation
    ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENENAME), ]

    # As long as it works perfectly we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma"), DGELRT or
    # DGELM (package edgeR) which is a list. I.e. we assign the annotation 
    # to a new element named 'genes' of this list.
    #	This works if they are in the same order
    #fit$genes <- ens.ann.1
    #	This works by matching rownames and gene IDs
    fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1$GENENAME), ]
    rownames(fit$genes) <- rownames(fit)
    #
    # add biomart annotation
    fit$genes <- merge(fit$genes, biomart, by.x="GENEID", by.y="ensembl_gene_id")
    #head(fit$genes)


    return (fit)
}

eR.save.fit <- function(fit, name) {
    # SAVE THE ANNOTATED FIT
    # ----------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(fit, file=paste(name,'.rds', sep=''))
    save(fit, file=paste(name, '.RData', sep=''))
}

eR.save.top.fit <- function(fit, file, n.genes=500, sort.by='PValue', p.value=0.01) {    # default adjust.method is BH
    
    # default p.value is 1 (all genes)
    tt <- topTags(fit, n=n.genes, sort.by=sort.by, p.value=p.value)
    n <- dim(tt)[1]

    # cap at the maximum number of genes requested
    if (n > n.genes) n <- n.genes
    file.name <- paste(file, '_top_', n, '_by_', sort.by, '.tab', sep='')
    write.table(tt$table[1:n, ], file.name, 
	row.names=F, col.names=T, sep='\t')
    # if fit is annotated the annotation will also be saved
}

# for Coturnix we use 'viral.dose' as basis for comparison
#design.column <- 'viral.dose'
# For Coturnix we will use 'Sample', 'CellStatus' or 'VirusStatus'

# For Gallus we will use ''

design.column <- 'Sample'
design <- model.matrix(~0 + target[ , design.column])

grps <- levels(as.factor(target[ , design.column]))
n.grps <- length(grps)
colnames(design) <- grps

qlfit <- glmQLFit(dge, design, robust=TRUE, abundance.trend=TRUE)
png.file <- paste(out.dir, '/edgeR/img/edgeR_QLFit_', design.column, '.png', sep='')
as.png(plotQLDisp(qlfit), png.file, overwrite=T)

# if we wanted to apply a log-fold-change theshold, we could do
# the calculations using qlfit here before testing for contrasts 
# using e.g.
# qlf.cmp <- glmTreat(glfit, cef=1..ncol(glfit$design), lfc=threshold)
# or
# qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)
#
cat("\n                         ")
cat("\n\tPerforming Pairwise Comparison by:", design.column)
cat("\n\t==================================================\n")

eR.data <- list()
for (i in grps) {
    for (j in grps) {
        if (i == j) next	# ignore self-comparisons
        
	comp.name <- paste(i, '_vs_', j, sep = "")        
	cat("\nComputing DGE:", comp.name,       '\n')
	cat(  "===================================\n")

	formula <- paste(i,"-",j)		
        cmp <- makeContrasts(formula, levels = design)
        
		# glmQLFTest is similar to glmLRT except it uses Bayes quasi-likelihood
        # the P-values are always >= those produced by glmLRT
        qlf.cmp <- glmQLFTest(qlfit, contrast = cmp)
        # or if a threshold has been defined
        #qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)
        
        # get summary of up/down regulated genes
        print(summary(decideTests(qlf.cmp)))
        png.file <- paste(out.dir, "/edgeR/img/edgeR_QLF_MD_", 
                          design.column, '_',comp.name, '.png', sep='')
        as.png( { 
                plotMD(qlf.cmp)
                abline(h = c(-1, 1), col="darkgreen")
                }, png.file, overwrite = TRUE)
        #continue.on.enter("Press [RETURN] to continue: ")
        
        # ANNOTATE the results
        qlf.cmp <- eR.annotate.fit(qlf.cmp, ens.db, bm.annot.1)#Got.1)

        # qlf.cmp is a list of tables, if we want to save it,
        # we'll need to save the whole object
        # will add .rds and .RData to the files created
        
		name <- paste(out.dir, '/edgeR/fit_', design.column, '_', comp.name, '_annot', sep='')
        if ( ! file.exists(paste(name, '.RData', sep='')) || ! file.exists(paste(name, '.RData', sep=''))) {
			eR.save.fit(qlf.cmp, name)
		}        
		# defaults: n.genes=500, sort.by='PValue', p.valu=0.01
        # we'll save all (<=100.000) significant genes
        
		# qlf.cmp$table is a table with logFC, logCPM, F and PValue
        # that is what weill be saved when using topTags and write.table
        # will add "_top_" n "_by_" sort.by
        # defaults: n.genes=500, sort.by='PValue', p.value=0.01
        
		name <- paste(out.dir, '/edgeR/comp_', design.column, '_', i, '_-_', j, '_annot', sep='')
		if ( ! file.exists(paste(name, '_top_', n, '_by_', 'PValue', '.tab', sep=''))){
			eR.save.top.fit(qlf.cmp, file = name, n.genes = 100000,
							sort.by = 'PValue', p.value=0.01)
        }

		# we can use limma to test for over-representation of gene
		# ontology (GO) terms or KEGG pathways with goana() or kegga()
		# using the entrez.gene.ids of DE genes, to an FDR of 0.05 (default)
		# we can use species.KEGG="gga" or "cjo"
		# ( see https://www.kegg.jp/kegg/catalog/org_list.html )
		#
		    # eR.go <- goana(qlf.cmp, species="Cj")
		# eR.go <- goana(qlf.cmp, species.="Gg")
		# eR.kegg <- kegga(qlf.cmp, species.KEGG="cjo")
		# eR.kegg <- kegga(qlf.cmp, species.KEGG="gga")
		# topGO(go, sort="up", number=n.genes)
		# topKEGG(keg, sort="up", number=n.genes)
		    
        qlf.result <- list(
							eR.cmp=qlf.cmp
							#, eR.go=eR.go
							#, eR.kegg=eR.kegg
							)
        eR.data[[formula]] <- qlf.result 
		print(topTags(qlf.cmp))
    }
}

# now eR.data is a list where each element is a comparison A-B (A minus B),
# i.e. each A-B is a list of tables, one of them named "table" and containing
# logFC, logCPM, F and PValue
#
# We would like to have FDR-corrected p-values as well, which can be got by
# defining a threshold. But that implies we know of a meaningful one, which
# we don't here.





#################################################################
#
#################################################################
                                                                #
#################################################################
#
#################################################################
                                                                #
#################################################################

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     D E S E Q 2
# ---------------------------------------------------------------

### FUNCTIONS FOR DESEQ2 ANALYSIS

source("./scripts/DESeq_functions.R")

# 1) ANNOTATE RESULTS

annotateDESeqResults <- function(results,
				ensembl.db,
				biobiomart.ds.name,
				annotation.dir,		# For saving annotation dataset
				save = TRUE,
				out.file.base,
				out.dir ){
	
	cat("\n	Annotating genes ...\n")
	
	## Get gene annotation from ENSEMBL databsase and BioMart database

	if ( ! file.exists(paste(annotation.dir, 'ensembl.annotation.txt', sep = '/'))) {
		ensembl.ann <- annotateGenesFromENSEMBL(ensembl.db = ensembl.db,
							gene.names = rownames(results),
							save = save,
							out.dir = annotation.dir)
	} else {
		ensembl.ann <- read.table(paste(annotation.dir, 'ensembl.annotation.txt', sep = '/'))
	}

	if ( ! file.exists(paste(annotation.dir, 'biomart.annotation.1st.tab', sep='/'))) {
		biomart.ann <- getBiomaRtAnnotation(biomart.ds.name = biobiomart.ds.name,
							annotation.dir = annotation.dir,
							save = save)
	} else {
		biomart.ann <- read.table(paste(annotation.dir, 'biomart.annotation.1st.tab', sep='/'))
	}

	# Remove duplicated records preserving only first entry

	ensembl.ann <- ensembl.ann[ ! duplicated(ensembl.ann$GENENAME), ]
	biomart.ann <- biomart.ann[ ! duplicated(biomart.ann$ensembl_gene_id), ]


	## Add annotation to DESeq results data.frame to keep everything in one place

	# Add row names (ensembl.gene.id) as an additional column
	res.ann <- results
	res.ann$ensembl.gene.id <- rownames(res.ann)
		
	# ENSEMB GENENAMES match the GENE.IDs from our results dataframe
	#res.ann <- cbind(res.ann, ensembl.ann[ match(res.ann$ensembl.gene.id, ensembl.ann$GENENAME), ])
	for (i in colnames(ensembl.ann)){
		res.ann[i] <- ensembl.ann[match(res.ann$ensembl.gene.id, ensembl.ann$GENENAME), i]
	}	
	
	# BIOMART ENTREZ ACCESSION match the GENE.IDs from our results dataframe
	#res.ann <- cbind(res.ann, biomart.ann[ match(res.ann$ensembl.gene.id, biomart.ann$entrezgene_accession), ])
	for (i in colnames(biomart.ann)){
		res.ann[i] <- biomart.ann[match(res.ann$ensembl.gene.id, biomart.ann$entrezgene_accession), i]
	}
		
	## Save annotated results as data frame (table)
	if (save){
		write.table(data.frame(res.ann),
					file = paste(out.dir, '/DESeq2/', out.file.base, '_ann_results.tab', sep = ''),
					row.names = TRUE, col.names = TRUE)

		saveRDS(res.ann, paste(out.dir, "/DESeq2/", out.file.base, "_ann_results.rds", sep=""))
	}
	return(res.ann)
}

# 2) COMPARE AND ANNOTATE RESULTS FROM DESEQ2 OBJECT
ddsCompareAnnotate <- function(	dds, 
								contrast,			# c(factor, numerator, denominator) 
								filterFun = ihw,	# IHW increases statistical power
								alpha = 0.01,
								annotate = TRUE, 
								ensembl.db,
								biobiomart.ds.name,
								annotation.dir,		# For saving annotation tables
								save = TRUE,
								out.file.base,
								out.dir ) {

	results <- results(dds,
					contrast = contrast, 
					pAdjustMethod="BH",
					filterFun = filterFun,
					alpha = alpha)
	
	print(head(results))
	print(mcols(results)$description)	# Description of the Results fields

    ## Provide default output base name
    if (missing(out.file.base)) {
		out.file.base <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
	}

	if (save){
	
		## Save comparison results
		write.table(data.frame(results),
				file = paste(out.dir, '/DESeq2/', out.file.base, '_raw_results.tab', sep = ''),
				row.names = TRUE, col.names = TRUE)

		saveRDS(results, paste(out.dir, "/DESeq2/", out.file.base, "_raw_results.rds", sep=""))
	
		## Save summary
		sink(paste(out.dir, '/DESeq2/', out.file.base, '_summary.txt', sep=''), split=T)
		summary(results)
		sink()
	}
	
	if (annotate){
		results <- annotateDESeqResults(results = results,
										ensembl.db = ensembl.db,
										biobiomart.ds.name = biobiomart.ds.name,
										annotation.dir = annotation.dir,
										save = save,
										out.file.base = out.file.base,
										out.dir = out.dir )
	}
	return(results)
}



# 3) IDENTIFY TOP N DIFFERENTIALLY EXPRESSED GENES
top_n <- function(ds, n, col = "padj", decreasing = F) {
	
	df <- data.frame(ds)
	ordered.idx <- order(abs(df[col]), decreasing = decreasing)
	ordered.df <- df[ordered.idx, ]
	top.n <- ordered.df[1:n, ]
	
	return(top.n)
}


# 4) PERFORM DESEQ2 COMPARISON, ANNOTATE AND PLOT THE RESULTS

analysePlotDESeq <- function(	#results,
								dds, 
								contrast,			# c(factor, numerator, denominator) 
								filterFun = ihw,	# IHW increases statistical power
								alpha = 0.01,
								shrnk.type,
								annotate = TRUE,
								ensembl.db,
								biobiomart.ds.name,
								annotation.dir,		# For saving annotation data
								save = TRUE,
								overwrite = FALSE,
								out.file.base,
								out.dir ) {

    ## Provide default output base name
    if (missing(out.file.base)) {
		out.file.base <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
	}

    ## Check if final comparison results already exist. If so 
	if ( ! file.exists(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep = "")) || (overwrite == TRUE)){

		## Check if DESeqResults object is provided
		
		results <- ddsCompareAnnotate(	
							dds = dds, 
							contrast = contrast,
							filterFun = ihw,
							alpha = alpha,
							annotate = annotate, 
							ensembl.db = ensembl.db,
							biobiomart.ds.name = biobiomart.ds.name,
							annotation.dir = annotation.dir,	
							save = save,
							out.file.base = out.file.base,
							out.dir = out.dir)

		## Histogram plot p-value distribution
		if (save) {
			out.png <- paste(out.dir, '/DESeq2/img/DESeq2_', out.file.base, '_padj_hist.png', sep='')
		} else out.png <- NULL

		as.png( {
				margins <- par("mar")
				par(mar = c(5, 5, 5, 5))
				hist(results$padj, main=paste("Adj. p-value distribution", "\n", out.file.base, sep = ''),
				breaks = 1/alpha, xlab = "Adj. p-value")
				par(mar = margins)
				},
		out.png, overwrite = TRUE)
		
		
		## Shrinkage of data to improve visualization and ranking
		cat("\n\tShrinking Log 2 Fold Change ...\n")    
		coef <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
		shrunk.lfc <- lfcShrink(dds = dds, coef = coef, type = shrnk.type)
		
		if (save) {
		    out.png <- paste(out.dir, "/DESeq2/img/DESeq2_", out.file.base, "_shrunk_", shrnk.type, "_MA.png", sep='') 
		} else out.png <- NULL
		
		as.png( {

			dim <- par("mfrow")
			par(mfrow = c(2,1))
		   
			# plotMA shows the log2 fold changes attributable to a given variable
			# over the mean of normalized counts for all the samples in the dataset
			# Points above alpha are colored, outliers are shown as directed triangles    
			DESeq2::plotMA(results, alpha = alpha, main = "Log2 Fold Change")
			
			# it is useful to look at the shrunk l2fc values, which removes the noise
			# associated with l2fc changes from low-count genes
			DESeq2::plotMA(shrunk.lfc, alpha = alpha, main = paste("Shrunken LFC (", shrnk.type, ")", sep = ""))
			
			# after plotMA, one may identify interesting points interactively using
			# identify() and clicking on them:
			# idx <- identify(res$baseMean, res$log2FoldChange)
			# rownames(res[idx, ])
			#
			# Alternatively, looking at the plot and deciding which coordinates are
			# of interest, one can use, e.g. 
			# res_wt_vs_PC[ (res_wt_vs_PC$log2FoldChange) > 4) 
			#               & (res_wt_vs_PC$baseMean > 1000), ]
			# or
			# res_wt_vs_PC[  (abs(res_wt_vs_PC$log2FoldChange) > 4) 
			#              & (res_wt_vs_PC$baseMean > 1000), ]
			par(mfrow = dim)
			},
			
		out.png, overwrite = TRUE)
		
		## Save shrunken results
		if (save) {
		
		    out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_raw_results.rds", sep='')
		    cat("\n\tSaving", out.file, '\n')
		    saveRDS(shrunk.lfc, file=out.file)
		
			out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_raw_results.tab", sep='')
			cat("\tSaving", out.file, '\n')
			write.table(shrunk.lfc, 
				        file=out.file,
				        row.names=T, col.names=T, sep='\t')
		}
		
		## Annotate shrunken data for visualization in MA plot	
		if (annotate) {
			
			ann.shrunk <- annotateDESeqResults(results = shrunk.lfc,
									ensembl.db = ensembl.db,
									biobiomart.ds.name = biobiomart.ds.name,
									annotation.dir = annotation.dir,
									save = FALSE)	
			if (save) {

			## Save shrunken annotated results		
				out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_ann_results.rds", sep='')
				cat("\n\tSaving", out.file, '\n')
				saveRDS(ann.shrunk, file=out.file)
			
				out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_ann_results.tab", sep='')
				cat("\tSaving", out.file, '\n')
				write.table(ann.shrunk, 
					        file=out.file,
					        row.names=T, col.names=T, sep='\t')
			}
		# Else, unnannotated shrunk.lfc will be represented in 
		} else ann.shrunk <- shrunk.lfc
		
		if (FALSE) {
		    if (save) {
		        out.file <- paste(out.dir, "/DESeq2/img/DESeq2_", out.file.base, "_shrunk_", shrnk.type, "_l2FC.png", sep='')
		        cat("\n\tPlotting", out.file, '\n')
		    } else out.file <- NULL
		    
		    as.png( {
					ggplot2::ggplot(data.frame(ann.shrunk), 
					ggplot2::aes(x = log2(baseMean), y=log2FoldChange),
					environment=environment()) + 	# this is supposed to make it work in a local env
					
					ggplot2::geom_point(ggplot2::aes(colour = padj < alpha), shape = 10, size = 1) +
					
					#ggplot2::geom_text(data=~top_n(.x, 10, wt=-padj),
					ggplot2::geom_text(data = top_n(ann.shrunk, 10, "padj"),
					ggplot2::aes(label = rownames(ann.shrunk))) +
					ggplot2::labs(x="Log2 Mean of Normalised Counts", y="Log2 Fold Change")
		    
		    }, out.file , overwrite = TRUE)
		}
		
		# Find statistically significant changes. Ignore genes with NAN p-value.
		signif <- results[ (results$padj < alpha) & ! is.na(results$padj) , ]
		signif$abs_lfc <- abs(signif$log2FoldChange)
		
		if (save) {
		    out.file <- paste(out.dir, '/DESeq2/signif/signif_', out.file.base, "_α<", alpha, ".tab", sep='')
		    cat("\n\tSaving", out.file, '\n')
		    write.table(data.frame(signif), file = out.file, sep = '\t', row.names = T, col.names = T)
		}
		
		# Order the data, firstly by the decreasing absolute 
		# value of log2FoldChange and secondly by the increasing pvalue....
		srt.signif <- signif[ order( signif$abs_lfc,
								signif$padj,
								decreasing=c(T, F)), ]
		if (save) {
		    out.file <- paste(out.dir, "/DESeq2/signif/signif_sorted_", out.file.base, ".tab", sep='')
		    cat("\tSaving", out.file,'\n')
		    write.table(srt.signif, file = out.file, sep = '\t', row.names = T, col.names = T)
		}
		
		## If shrunken data is unnanotated, return null.
		if ( ! annotate) ann.shrunk <- NULL
		
		# Generate output list containing results    
		res.list <- list(result = results, 
		            shrunk = shrunk.lfc, 
		            shrunk.annot = ann.shrunk, 
		            signif = srt.signif)
		
		# Save results as an object to avoid redundant analyses.
		cat("\n\tSaving final results ...")		
		saveRDS(res.list, file=paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep=''))
		save(res.list, file=paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.RData", sep=''))

	} else {
	
		res.list <- load(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.RData", sep = ""))
	}
		
	return(res.list)

		cat("\n =========================== \n")
		cat(  "| A N A L Y S I S   D O N E |\n")
		cat(  " =========================== \n")
}

#############################################################################################
############################ --- START DESEQ2 ANALYSIS --- ##################################
#############################################################################################

##################
###   ADRIAN   ###	
##################

## Once the necessary functions have been loaded we can start the Differential
## Expression Analysis. We first need to load the count data files.

if (HTSEQ == FALSE){
	countDataFile <- paste('analysis.out/counts',
					'featureCounts.csv', sep = '/')
	countData <- read.csv(countDataFile,
				 row.names=1) %>%
				 as.matrix()
} else {
	countDataFile <- paste('analysis.out/counts',
					'htseq_counts_merged.tab', sep = '/')	

	countData <- read.table(countDataFile,
				 row.names=1) %>%
				 as.matrix()
}

## Discard genes with no counts in any of the samples
countData <- countData[rowSums(countData)>1, ]
head(countData)
dim(countData)

## Read table containing the sample information
sampleData <- read.table("SampleInfo.tab", header = TRUE)
sampleData$Sample <- relevel(as.factor(sampleData$Sample), "DF1")
sampleData$CellStatus <- relevel(as.factor(sampleData$CellStatus), "Parental")
sampleData$VirusStatus <- relevel(as.factor(sampleData$VirusStatus), "Mock")


cat("\n\t ==================================================\n")
cat(  "\t|   DIFFERENTIAL EXPRESSION ANALYSIS WITH DESEQ2   |\n")
cat(  "\t ==================================================\n")

	###################################################
	## COMPARE BY SAMPLE (PARENTAL-PERSISTANT-CURED) ##
	###################################################

# Do DESeq analysis comparing samples by Cell Status and save results
if ( ! file.exists(paste(out.dir, "DESeq2/dds.Sample.rds", sep='/'))) {

    dds.sample <- DESeqDataSetFromMatrix(countData, sampleData,
				design=~Sample, tidy=F)

    dds.sample <- DESeq(dds.sample)

    # SAVE DESEQ ANALYSIS
    # -------------------
    save(dds.sample, file=paste(out.dir, "DESeq2/dds.Sample.RData", sep='/'))
    saveRDS(dds.sample, file=paste(out.dir, "DESeq2/dds.Sample.rds", sep='/'))
} else {
    dds.sample <- readRDS(file=paste(out.dir, "DESeq2/dds.Sample.rds", sep='/'))
	#dds.sample <- load(file=paste(out.dir, "DESeq2/dds.Sample.RData", sep='/'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.sample))

### COMPARE BY SAMPLE AND ANNOTATE THE RESULTS

dds <- dds.sample
comparisons <- resultsNames(dds)[2:length(resultsNames(dds))]

for (comparison in comparisons){
	if (comparison == "Intercept") next
#for (comp in levels(sampleData$Sample)){
#	for (j in levels(sampleData$Sample)){
#		if (i == j) next
		cat("\n                                    ")						
		cat("\nComputing DGE:",	  comparison,	"\n")
		cat(  "==================================\n")
		#contrast <- c("Sample", "DF1.PC", "DF1")
		contrast <- str_split(comparison, pattern = "_")[[1]][c(1, 2, 4)]	#Ex: c("Sample", "DF1.PC", "DF1")
		#print(contrast)

		#out.file.base <- paste("dds.Sample.DF1_PC.vs.DF1")
		out.file.base <- paste("dds", comparison, sep = ".")

		dds.sample.results <- analysePlotDESeq(
							dds = dds,
							contrast = contrast,
							filterFun = ihw,
							alpha = 0.01,
							shrnk.type = 'apeglm',
							annotate = TRUE, 
							ensembl.db = ensembl.db,
							biobiomart.ds.name = biobiomart.ds.name,
							annotation.dir = annotation.dir,
							save = TRUE,
							out.file.base = out.file.base,
							out.dir = out.dir)
#	}
}

	########################################################
	## COMPARE BY CELL STATUS (PARENTAL-PERSISTANT-CURED) ##
	########################################################

# Do DESeq analysis comparing samples by Cell Status and save results
if ( ! file.exists(paste(out.dir, "DESeq2/dds.CellStatus.rds", sep='/'))) {

    dds.cell <- DESeqDataSetFromMatrix(countData, sampleData,
				design=~VirusStatus + CellStatus, tidy=F)

# We need to add VirusStatus to consider the differences within the same
# Cell type due to the infection(batch effects)

    dds.cell <- DESeq(dds.cell)

    # SAVE DESEQ ANALYSIS
    # -------------------
    save(dds.cell, file=paste(out.dir, "DESeq2/dds.CellStatus.RData", sep='/'))
    saveRDS(dds.cell, file=paste(out.dir, "DESeq2/dds.CellStatus.rds", sep='/'))
} else {
    dds.cell <- readRDS(file=paste(out.dir, "DESeq2/dds.CellStatus.rds", sep='/'))
	#dds.cell <- load(file=paste(out.dir, "DESeq2/dds.CellStatus.RData", sep='/'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.cell))

### COMPARE BY CELL STATUS AND ANNOTATE THE RESULTS

dds <- dds.cell
comparisons <- resultsNames(dds)[2:length(resultsNames(dds))]

for (comparison in comparisons){
	if (comparison == "Intercept") next
		cat("\n                                    ")				
		cat("\nComputing DGE:",	  comparison,	"\n")
		cat(  "==================================\n")
		contrast <- str_split(comparison, pattern = "_")[[1]][c(1, 2, 4)]	#Ex: c("Sample", "DF1.PC", "DF1")
		#print(contrast)

		out.file.base <- paste("dds", comparison, sep = ".")

		dds.cell.results <- analysePlotDESeq(
							dds = dds,
							contrast = contrast,
							filterFun = ihw,
							alpha = 0.01,
							shrnk.type = 'apeglm',
							annotate = TRUE, 
							ensembl.db = ensembl.db,
							biobiomart.ds.name = biobiomart.ds.name,
							annotation.dir = annotation.dir,
							save = TRUE,
							out.file.base = out.file.base,
							out.dir = out.dir)

}

	#############################################
	## COMPARE BY VIRUS STATUS (MOCK-INFECTED) ##
	#############################################

# Do DESeq analysis comparing samples by Virus Status and save results
if ( ! file.exists(paste(out.dir, "DESeq2/dds.VirusStatus.rds", sep='/'))) {

    dds.virus <- DESeqDataSetFromMatrix(countData, sampleData,
				design=~CellStatus + VirusStatus, tidy=F)

# We need to add CellStatus to consider the differences within the same
# Infection due to the cell type(batch effects)

    dds.virus <- DESeq(dds.virus)

    # SAVE DESEQ ANALYSIS
    # -------------------
    save(dds.virus, file=paste(out.dir, "DESeq2/dds.VirusStatus.RData", sep='/'))
    saveRDS(dds.virus, file=paste(out.dir, "DESeq2/dds.VirusStatus.rds", sep='/'))
} else {
    dds.virus <- readRDS(file=paste(out.dir, "DESeq2/dds.VirusStatus.rds", sep='/'))
	#dds.virus <- load(file=paste(out.dir, "DESeq2/dds.VirusStatus.RData", sep='/'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.virus))

### COMPARE BY INFECTION STATUS AND ANNOTATE THE RESULTS

dds <- dds.virus
comparisons <- resultsNames(dds)[2:length(resultsNames(dds))]

for (comparison in comparisons){
	if (comparison == "Intercept") next
		cat("\n                                            ")		
		cat("\nComputing DGE:",	  comparison,	        "\n")
		cat(  "==========================================\n")
		contrast <- str_split(comparison, pattern = "_")[[1]][c(1, 2, 4)]	#Ex: c("Sample", "DF1.PC", "DF1")
		#print(contrast)

		out.file.base <- paste("dds", comparison, sep = ".")

		dds.virus.results <- analysePlotDESeq(
							dds = dds,
							contrast = contrast,
							filterFun = ihw,
							alpha = 0.01,
							shrnk.type = 'apeglm',
							annotate = TRUE, 
							ensembl.db = ensembl.db,
							biobiomart.ds.name = biobiomart.ds.name,
							annotation.dir = annotation.dir,
							save = TRUE,
							out.file.base = out.file.base,
							out.dir = out.dir)

}

	###########################################
	## COMPARE BY BOTH CELL AND VIRUS STATUS ##
	###########################################

# Do DESeq analysis comparing samples by Virus Status and save results
if ( ! file.exists(paste(out.dir, "DESeq2/dds.FullCompare.rds", sep='/'))) {

    dds.both <- DESeqDataSetFromMatrix(countData, sampleData,
				design=~CellStatus + VirusStatus, tidy=F)
                                
    dds.both <- DESeq(dds.both)

    # SAVE DESEQ ANALYSIS
    # -------------------
    save(dds.both, file=paste(out.dir, "DESeq2/dds.FullCompare.RData", sep='/'))
    saveRDS(dds.both, file=paste(out.dir, "DESeq2/dds.FullCompare.rds", sep='/'))
} else {
    dds.both <- readRDS(file=paste(out.dir, "DESeq2/dds.FullCompare.rds", sep='/'))
	#dds.both <- load(file=paste(out.dir, DESeq2/dds.FullCompare.RData", sep='/'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.both))


### COMPARE BY CELL AND INFECTION STATUS AND ANNOTATE THE RESULTS

dds <- dds.both
comparisons <- resultsNames(dds)[2:length(resultsNames(dds))]

for (comparison in comparisons){
	if (comparison == "Intercept") next
		cat("\n                                    ")
		cat("\nComputing DGE:",	  comparison,	"\n")
		cat(  "==================================\n")
		contrast <- str_split(comparison, pattern = "_")[[1]][c(1, 2, 4)]	#Ex: c("Sample", "DF1.PC", "DF1")
		#print(contrast)

		out.file.base <- paste("dds", comparison, sep = ".")

		dds.both.results <- analysePlotDESeq(
							dds = dds,
							contrast = contrast,
							filterFun = ihw,
							alpha = 0.01,
							shrnk.type = 'apeglm',
							annotate = TRUE, 
							ensembl.db = ensembl.db,
							biobiomart.ds.name = biobiomart.ds.name,
							annotation.dir = annotation.dir,
							save = TRUE,
							out.file.base = out.file.base,
							out.dir = out.dir)

}


#
# Extract annotation using Org object if avaiable
#

if (FALSE) 
if ( ! is.null(org.db) ) {

	# And now we should be able to use the Org package if we successfully built it
	# at the beginning.

	geneSymbols <- mapIds(org.db, 
			  keys = as.character(res.ann$ENTREZID), 
			  columns = "ENTREZID", 
			  keytype = "ENTREZID", 
			  multiVals = "first")

	# Check the information we can retrieve from the database
	keytypes(org.db)
	
	# Retrieve Gene Ontology IDs
	GO.id <- AnnotationDbi::select(org.db, 
				    keys=as.character(res.ann$ENTREZID), 
				    columns=c("ENTREZID", "GO", "GOALL", "ONTOLOGY","ONTOLOGYALL"), 
				    keytype="ENTREZID", 
				    multiVals="CharacterList")

	# Retrieve Gene Ontology descriptions
	GO.description <- AnnotationDbi::select(GO.db,
					keys=GO.id$GO, 
					columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"),
					keytype= "GOID")

}



# Generating this takes long. We should save it and load from file
# when we run the script
ds.data <- list()
#x <- 0
contr <- "Sample"
grps <- levels(sampleData[ , contr ])
for (a in grps) {
    for (b in grps) {
        print(paste(a, b))
        if (a == b) next
        res <- plot_and_save( 
                        dds, contr,
                        x=a, y=b,
                        filterFun=ihw, alpha=0.01,
                        ensembl.ann=ens.ann.1,
                        biomart.ann=bm.annot.1,
                        outDir=out.dir
                        )
        print(names(res))
	name <- paste(contr, "_", a, "_", b, sep='')
        #x <- x + 1
        #ds.data[[x]] <- res
        ds.data[[name]] <- res
        #names(ds.data)[x] <- name
        #stop()
    }
}
print(names(ds.data))




# -------------------------------
# Do Gene Set Enrichment Analysis
# -------------------------------

GO_fgsea <- function (ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(out.dir, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      top.n=20,
                      top.biblio=5,
                      verbose=FALSE) {
    # Do GSEA on GO terms using fgsea

    # Rank all genes on their fold change.
    #	Here we exclude genes for which we have no EntrezID and
    #	use shrunk LFC values
#    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    gseaDat <- filter(ann.shrunk.lfc, !is.na(GENEID))

    ranks <- gseaDat$lfc
    #names(ranks) <- gseaDat$ENTREZID
    names(ranks) <- gseaDat$GENEID
    head(ranks)
    #uranks <- ranks[!duplicated(sort(names(ranks)))]

    # plot all the ranked fold changes
    out.png <- paste(out.dir, '/', out.name, '.barplot.png', sep='')
    as.png(barplot(sort(ranks, decreasing=T)), out.png, overwrite = TRUE)
    
    # load pathways
    #pathways <- ann.shrunk.lfc$ENTREZID
    #pathways <- ann.shrunk.lfc$go_id
    #names(pathways) <- gseaDat$ENTREZID
    #names(pathways) <- gseaDat$GENEID
    #upathways <- pathways[!duplicated(sort(names(pathway)))]
    
    # create a list of go_id terms where each term contains a vector
    # of emsembl_gene_id in that term
    #   first recover the annotation (in case we are re-run and do not
    #   have it
    if ( ! exists(substitute(bm.go.annot)) ) {
        bm.go.annot <- read.table(paste(out.dir, '/biomart.go.tab', sep=''), 
	            header=T, sep='\t')
    }
    #if ( ! exists(substitute(bm.goslim.annot)) ) {
    #    bm.goslim.annot <- read.table(paste(out.dir, '/biomart.goslim.tab', sep=''), 
    #	            header=T, sep='\t')
    #}

    if (use.description == FALSE) {
        # split() will divide the gene-ids by their go-id
        pathways <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$go_id)
    } else {
        # Do the same but with large gene ontology names
        # split() will divide the gene-ids by their go-name
        pathways.go <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$name_1006)
    }
    # do analysis
    # the resulting table contains enrichment scores and p-values
    out.file <- sprintf("%s/go_fgsea_10-%d.RData", out.dir, max.size)
    if (file.exists(out.file)) {
        load(file=out.file)
    } else {
        fgseaRes <- fgsea(pathways=pathways.go, 
                          stats=ranks, 
		          minSize=10, 
		          maxSize=max.size, 
                          nPermSimple=100000 
                          )
        save(fgseaRes, file=out.file)
    }
    if (verbose == TRUE)
        head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    if (verbose == TRUE) {
        # plot enrichment score
        sorted.fgsea.res <- fgseaRes[order(padj, -abs(NES)), ]
        sfr.names <- sorted.fgsea.res$pathway
        for (i in 1:top.n) {
            if (use.description == FALSE)
                descr <- bm.go.annot[bm.go.annot$go_id == sfr.names[i], "name_1006"]
            else
                descr <- sfr.names[i]
            print(
                plotEnrichment(pathways.go[[ sfr.names[i] ]], ranks) +
                    labs(title=descr)
                )
            ans <- readline("Press RETURN to continue: ")
            if (ans == "q") break
        }
    }
    # gsea table plot of top.n gene families
    #   top_n() is now deprecated
    topUp <- fgseaRes %>%
        filter(ES > 0) %>%
        top_n(top.n, wt=-padj)

    topDown <- fgseaRes %>%
        filter(ES < 0) %>%
        top_n(top.n, wt=-padj)

    topPathways <- bind_rows(topUp, topDown) %>%
        arrange(-ES)
    # last resort when "pos" is used    
    #topPathways <- sorted.fgsea.res[1:top.n, ]

    # do the plots and save descriptions
    out.file <- paste(out.dir, "/topUp.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topUp[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topUp', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topUp.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }
    out.file <- paste(out.dir, "/topDn.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topDown[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topDown', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topDn.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }

    out.png <- paste(out.dir, '/', out.name, '.GSEAtable.png', sep='')
    as.png(
        plotGseaTable(pathways.go[topPathways$pathway], 
                      ranks, 
                      fgseaRes, 
                      gseaParam = 0.5)
    , out.png, width=1024, height=100*top.n, overwrite=TRUE)
    
    
    # and now do a plot of the interest in citations during
    # the last ten years for the top 5 sets
    out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
    cur.year <- as.integer(format(Sys.Date(), "%Y"))
    terms <- topDown[1:top.biblio]$pathway
    as.png(
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE),
        out.png, overwrite = TRUE)

    
    # finally, return fgseaRes
    return(fgseaRes)
}



GO_KEGG_clusterProfiler <- function(ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(out.dir, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      OrgDb = org.Cjaponica.eg.db,
                      kegg_organism = "cjo",	# (cjo = coturnix japonica)
                                             # (gga = gallus gallus)
                      top.n=10,
                      top.biblio=5,
                      verbose=FALSE) {
        
    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    ranks <- gseaDat$lfc
    names(ranks) <- gseaDat$ENTREZID
    ranks<-na.omit(ranks)

    s.ranks <- sort(ranks, decreasing=T)
    
    ### G O annotation
    #
    gse <- gseGO(geneList=s.ranks, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = OrgDb, 
                 pAdjustMethod = "fdr")
                 #pAdjustMethod = "none")
    if (dim(gse)[1] == 0) {
        # Try without correction issuing a warning
        gse <- gseGO(geneList=s.ranks, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = OrgDb, 
                 #pAdjustMethod = "fdr")
                 pAdjustMethod = "none")
        out.name <- paste(out.name, '.raw_p', sep='')
    }
    if (dim(gse)[1] == 0) {
        # give up
        cat('', file=paste(out.dir, '/',
                out.name, 'NO_ENRICHED_GO', sep=''))
#        return()
    } else {
        out.png <- paste(out.dir, '/', out.name, '.GOdotplot.png', sep='')
        as.png(
            dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
            , out.png, overwrite = TRUE)

        out.png <- paste(out.dir, '/', out.name, '.GOemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(gse), showCategory = top.n)
            , out.png, overwrite = TRUE)

        out.png <- paste(out.dir, '/', out.name, '.GOridgeplot.png', sep='')
        as.png(
            ridgeplot(gse) + labs(x = "enrichment distribution")
            , out.png, overwrite = TRUE)

        out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
        cur.year <- as.integer(format(Sys.Date(), "%Y"))
        terms <- gse$Description[1:top.biblio]
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE)

        for (i in 1:dim(gse)[1] ) {
        out.png <- paste(out.dir, '/', 
                out.name, '.GOcnetplot.', i, '.', gse$ID[i], '.png', sep='')
        # categorySize can be either 'pvalue' or 'geneNum'
        as.png(
            cnetplot(gse, categorySize="pvalue", foldChange=s.ranks, 
                     showCategory=i)
            , out.png, overwrite = TRUE)

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                    out.name, '.GOgseaplot.', i, '.', gse$ID[i], '.png', sep='')
            gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
        }
        out.file <- paste(out.dir, '/', out.name, '.topGO.tab', sep='')
        write.table(gse, file=out.file, row.names=T, col.names=T, sep='\t')
    }
    
    
    ### K E G G annotation
    # let's try with KEGG (from the ENTREZID which is the same a ncbi-genid)
    kse <- gseKEGG(geneList     = s.ranks,
               organism     = kegg_organism,
               #nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid",
               nPermSimple = 100000)
               
    if (dim(kse)[1] == 0) {
        # Try without correction issuing a warning
        kse <- gseKEGG(geneList     = s.ranks,
                   organism     = kegg_organism,
                   #nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "none",
                   keyType       = "ncbi-geneid",
                   nPermSimple = 100000)
        out.name <- paste(out.name, 'raw_p', sep='')
    }
    if (dim(kse)[1] == 0) {
        # give up
        cat('', file=paste(out.dir, '/',
                out.name, 'NO_ENRICHED_KEGG', sep=''))
#        return()
    } else {
        out.png <- paste(out.dir, '/', out.name, '.KEGGdotplot.png', sep='')
        as.png( 
            dotplot(kse, showCategory = 10, title = "Enriched Pathways" , 
                    split=".sign") + facet_grid(.~.sign)
            , out.png, overwrite = TRUE)

        out.png <- paste(out.dir, '/', out.name, '.KEGGemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(kse))
            , out.png, overwrite = TRUE)


        # categorySize can be either 'pvalue' or 'geneNum'
        out.png <- paste(out.dir, '/', out.name, '.KEGGcnetplot.png', sep='')
        as.png(
            cnetplot(kse, categorySize="pvalue", foldChange=s.ranks)
            , out.png, overwrite = TRUE)

        out.png <- paste(out.dir, '/', out.name, '.KEGGridgeplot.png', sep='')
        as.png(
            ridgeplot(kse) + labs(x = "enrichment distribution")
            , out.png, overwrite = TRUE)

        cur.dir <- getwd()
        setwd(out.dir)
        for (i in 1:dim(kse)[1]) {
            # for each of the pathways in kse

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                out.name, '.GSEAplot.', i, '.', kse$ID[i], '.png', sep='')
            gseaplot(kse, by = "all", title = kse$Description[i], geneSetID = i)

            # Produce the native KEGG plot (PNG)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism)

            # Produce a different plot (PDF) (different from the previous one)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism, kegg.native = F)
        }
        setwd(cur.dir)

        out.file <- paste(out.dir, '/', out.name, '.topKEGG.tab', sep='')
        write.table(kse, file=out.file, row.names=T, col.names=T, sep='\t')
    }
}


# we do not know what is the best upper limit, so we'll try several
ms <- 500 # default value in GSEA, should work as well as the others
for (ms in c(50, 100, 250, 500)) {
    for (cmp.name in names(ds.data)) {
        # for each comparison name
        cmp.data <- ds.data[[ cmp.name ]]
        cat('Doing fGSEA on', cmp.name, '\n') 
        ann.shrunk.lfc <- cmp.data[[ "shrunk.annot" ]]
        #out.dir <- paste(out.dir, "DESeq2/GO_fgsea", cmp.name, sep='/')
        out.dir <- sprintf("%s/DESeq2/GO_fgsea/max.size=%03d/%s",
                out.dir, ms, cmp.name)
        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
        gogsea <- GO_fgsea(ann.shrunk.lfc, 
                   max.size=ms,
                   out.dir=out.dir, 
                   out.name='fgsea',
                   use.description=TRUE)
        # try dotplot(gogsea), emapplot(gogsea)... etc from clusterProfiler
        # we should add fgsea results to the cmp.data list
        ds.data[[cmp.name]][['go.fgsea']] <- gogsea

        out.cprof <- sprintf("%s/DESeq2/GO+KEGG_cProf/max.size=%03d/%s",
                out.dir, ms, cmp.name)
        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
        GO_KEGG_clusterProfiler(
                      ann.shrunk.lfc, 
		      max.size=ms,
                      out.dir=out.cprof, 
                      out.name='cProf',
                      use.description=TRUE,
                      verbose=FALSE)
    }
}



continue.on.enter("You may want to maximize your terminal before continuing ")
options(width=200)
n <- 10
for (i in names(ds.data)) {
    cat("\nMost significant', n, 'genes for", i, '\n')
    print(head(ds.data[[i]]$signif.annot[ , c("log2FoldChange", "entrezgene_description")] ), n)
    continue.on.enter("Press [ENTER] to continue ")
}
options(width=80)
continue.on.enter("Done, you can restore your terminal now ")


# Plot counts of the gene with maximal l2FC (up or down)
threshold <- 0
for (n in names(ds.data)) {
    #res <- data[["PFU_wt_0.1"]]$signif
    res <- ds.data[[n]]$signif.annot
    up <- res[ res$log2FoldChange > threshold, ]
    down <- res[ res$log2FoldChange < -threshold, ]

    most.up <- rownames(up)[which.max(up$log2FoldChange)]
    cat("\ncounts of most overexpressed gene in", n, ":\n",
        most.up, up$log2FoldChange[which.max(up$log2FoldChange)], up$entrezgene_description[which.max(up$log2FoldChange)], 
        "\n")
    #plotCounts(dds, gene=rownames(res)[which.max(res$log2FoldChange)], intgroup="PFU")
    print(plotCounts(dds, gene=most.up, intgroup="PFU"))
    k <- continue.on.key()
    if (k == "q") break
    
    most.down <- rownames(down)[which.min(down$log2FoldChange)]
    cat("\ncounts of most underexpressed gene in", n, ":\n",
        most.down, down$log2FoldChange[which.min(down$log2FoldChange)], down$entrezgene_description[which.min(down$log2FoldChange)],
        "\n")
    #plotCounts(dds, gene=rownames(res)[which.min(res$log2FoldChange)], intgroup="PFU")
    print(plotCounts(dds, gene=most.down, intgroup="PFU"))
    continue.on.key()
    if (k == "q") break
}
cat("\n")




# Analyze GO representation
# -------------------------
#
# we use the biomaRt database from above
#	mart.db might have been already assigned above, but not
#	necessarily
#ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
mart.db <- useMart("ensembl", mart.name)
#listAttributes(mart.db)

for (n in names(ds.data)) {
    cat("Processing GOs for", n, '\n')
    sig <- as.data.frame(ds.data[[n]]$signif)
    sig.a <- as.data.frame(ds.data[[n]]$signif.annot)
    sig$ensembl_gene_id <- rownames(sig)
    ourFilterType <- "ensembl_gene_id"
    filterValues <- sig$ensembl_gene_id
    # this will give us ALL GO annotations
    # we should use the saved searches from above to save bandwidth.
#    gos <- getBM(attributes=c("ensembl_gene_id", 
#                              "go_id", "name_1006", 
#                              "definition_1006" ), 
#                 mart=mart.db,
#                 filters=ourFilterType,
#                 values=filterValues)
    gos <- bm.go.annot[ bm.go.annot$ensembl_gene_id %in% filterValues,
                        c("ensembl_gene_id", 
                          "go_id", "name_1006", 
                          "definition_1006" ) ]
    sig.a <- merge(sig, gos, by="ensembl_gene_id")
    table(sig.a$go_id)	# this gives counts, but we'd like to multiply those
                            # counts by log2FC
    #sig.a[ , c(1, 3, 10)]
    # aggregate log2FC by go_id and sum the values
    print(
    aggregate(sig.a$log2FoldChange, 
              by=list(Category=sig.a$go_id),
              FUN=sum)
    )
    # or with formula interface
    gosums <- aggregate(log2FoldChange ~ go_id, sig.a, sum)
    gosums <- gosums[ order(gosums$log2FoldChange), ]
    # annotate gosums
    gosums.a <- cbind(gosums, bm.annot.1[ match(gosums$go_id, bm.annot.1$go_id), ])
    
    out.dir <- paste(out.dir, '/DESeq2/go', sep='')
    dir.create(out.dir, showWarning=FALSE)
    out.file <- paste(out.dir, '/GO_', n, "_sum.tab", sep='')
    write.table(gosums.a, out.file, sep='\t')
    
    goavgs <- aggregate(log2FoldChange ~ go_id, sig.a, mean)
    goavgs <- gosums[ order(gosums$log2FoldChange), ]
    # annotate gosums
    goavgs.a <- cbind(goavgs, bm.annot.1[ match(goavgs$go_id, bm.annot.1$go_id), ])
    out.file <- paste(out.dir, '/GO_', n, "_average.tab", sep='')
    write.table(goavgs.a, out.file, sep='\t')
}


# Reannotate GO with ancestry
# library(GO.db)

annotate.go.ancestry <- function (annot.data, l2fc.threshold, outDir) {
    res <- annot.data
    
    goBPanc <- as.list(GOBPANCESTOR)
    # remove GO terms that do not have any ancestor
    goBPanc <- goBPanc[ ! is.na(goBPanc) ]

    goCCanc <- as.list(GOCCANCESTOR)
    # remove GO terms that do not have any ancestor
    goCCanc <- goCCanc[ ! is.na(goCCanc) ]

    goMFanc <- as.list(GOMFANCESTOR)
    # remove GO terms that do not have any ancestor
    goMFanc <- goMFanc[ ! is.na(goMFanc) ]
    
    res <- res[ abs(res$log2FoldChange) > l2fc.threshold, ]
    go.l2fc <- res[ , c("ensembl_gene_id", "log2FoldChange", "GENENAME", "go_id") ]
    
    for (i in 1:nrow(res)) {
        if ((i %% 100) == 0) cat(".")
        if (is.null(res[i, "go_id"])) next
        if (is.na(res[i, "go_id"])) next
        if (res[i, "go_id"] == '') next
        anc <- goBPanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) & ! length(anc) == 0) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
        anc <- goCCanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) ) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
        anc <- goMFanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) ) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
    }
    # aggregate data by GOID (will result in c(Category, x) columns
    go.l2fc.sum <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=sum)
    go.l2fc.avg <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=mean)
    
    # rename Category, x to go_id, l2fc
    colnames(go.l2fc.sum) <- c('go_id', 'sum.l2fc')
    colnames(go.l2fc.avg) <- c('go_id', 'avg.l2fc')
  
    
    # annotate
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc <- cbind(go.l2fc, godesc)
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc.sum$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc.sum <- cbind(go.l2fc.sum, godesc)
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc.avg$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc.avg <- cbind(go.l2fc.avg, godesc)


    # sort by l2FC and save
    go.l2fc <- go.l2fc[ order(go.l2fc$log2FoldChange, decreasing=T), ]
    go.l2fc.sum <- go.l2fc.sum[ order(go.l2fc.sum$sum.l2fc, decreasing=T), ]
    go.l2fc.avg <- go.l2fc.avg[ order(go.l2fc.avg$avg.l2fc, decreasing=T), ]
    # sort by abs(log2FC)
    go.l2fc.sum.abs <- go.l2fc.sum[ order(abs(go.l2fc.sum$sum.l2fc), decreasing=T), ]
    go.l2fc.avg.abs <- go.l2fc.avg[ order(abs(go.l2fc.avg$avg.l2fc), decreasing=T), ]
    
    # save
    out.file <- paste(outDir, '/GOANC_', n, ".tab", sep='')
    write.table(go.l2fc, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_sum.tab", sep='')
    write.table(go.l2fc.sum, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_average.tab", sep='')
    write.table(go.l2fc.avg, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_sum_abs.tab", sep='')
    write.table(go.l2fc.sum.abs, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_average_abs.tab", sep='')
    write.table(go.l2fc.avg.abs, out.file, sep='\t')
    
}
 
 
l2fc.threshold <- 0
out.goanc <- paste(out.dir, '/Deseq2/goanc')
dir.create(out.goanc, showWarning=FALSE)

for (n in names(ds.data)) {
    cat("\nTracing GO ancestry for", n, '\n')
    res <- ds.data[[n]]$signif.annot
        
    annotate.go.ancestry(res, l2fc.threshold, out.goanc)
    
    #ans <- readline("Press RETURN to continue: ")
    #if (ans == "q") break
}


# Analyze PFAM representation
# ---------------------------
#
# for PFAM, we can use
#
# Get PFAM database and use AC -> DE mapping
#library(PFAM.db)
db <- PFAMDE
mk <- PFAMDE[mappedkeys(PFAMDE)]
xx <- as.list(mk)
pfam.table <- toTable(PFAMDE)

out.pfam <- paste(out.dir, '/DESeq2/pfam', sep='')
dir.create(out.pfam, showWarning=FALSE)


#for (i in pfam.ids) print(xx[[i]])

# Get PFAM families and sort them by their representation in the dataset
options(width=200)
for (n in names(ds.data)) {
    sig <- as.data.frame(ds.data[[n]]$signif) 
    names(sig)
    sig$ensembl_gene_id <- rownames(sig)
    sig.a <- cbind(sig, 
           bm.fam.annot[ match(sig$ensembl_gene_id, bm.fam.annot$ensembl_gene_id), ])
    pfam.ids <- sig.a$pfam[ ! is.na(sig.a$pfam) ]
    aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)

    sig.a <- merge(sig, bm.fam.annot, by="ensembl_gene_id")
    head(sig.a, 10)
    
    sig.sum <- aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)
    sig.sum <- sig.sum[ order(sig.sum$x, decreasing=T), ]
    sig.sum <- sig.sum[ ! is.na(sig.sum$Category), ]
    sig.sum <- sig.sum[ sig.sum$Category != '', ]
    names(sig.sum) <- c("pfam", "sum.l2FC")
    # annotate
    for (i in 1:length(sig.sum$pfam)) { 
        pf <- sig.sum$pfam[i] ; 
        if (! is.null(xx[[pf]])) { 
            sig.sum$pfam.de[i] <- xx[[pf]] 
        } else { 
            sig.sum$pfam.de[i] <- '' 
        }
    }
    cat("\nMost represented PFAM families in", n, "\n")
    print(head(sig.sum, 20)) ; print(tail(sig.sum, 20))
    out.file <- paste(out.dir, '/PFAM_', n, "_sum.tab", sep='')
    write.table(sig.sum, out.file, sep='\t')

    sig.avg <- aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=mean)
    sig.avg <- sig.avg[ order(sig.avg$x, decreasing=T), ]
    sig.avg <- sig.avg[ ! is.na(sig.avg$Category), ]
    sig.avg <- sig.avg[ sig.avg$Category != '', ]
    names(sig.avg) <- c("pfam", "avg.l2FC")
    # annotate
    for (i in 1:length(sig.avg$pfam)) { 
        pf <- sig.avg$pfam[i] ; 
        if (! is.null(xx[[pf]])) { 
            sig.avg$pfam.de[i] <- xx[[pf]] 
        } else { 
            sig.avg$pfam.de[i] <- '' 
        }
    }
    out.file <- paste(fout.dir, '/PFAM_', n, "_average.tab", sep='')
    write.table(sig.avg, out.file, sep='\t')
    
    # save also the files sorted by abs(l2fc)
    sig.sum.abs <- sig.sum[ order(abs(sig.sum$sum.l2FC), decreasing=T), ]
    sig.avg.abs <- sig.avg[ order(abs(sig.avg$avg.l2FC), decreasing=T), ]
    out.file <- paste(out.dir, '/PFAM_', n, "_sum_abs.tab", sep='')
    write.table(sig.sum.abs, out.file, sep='\t')
    out.file <- paste(out.dir, '/PFAM_', n, "_average_abs.tab", sep='')
    write.table(sig.avg.abs, out.file, sep='\t')


    #continue.on.key()
    ans <- continue.on.enter(prompt="Press return to continue ")
    if (ans == "q") break
}
cat('\n')
options(width=80)




# DO PCA AND CLUSTERING ANALYSES
# ------------------------------

## CREATE THE COLUMNS FOR THE PCA ANALYSIS 
# we need to retrieve from the dataframes only the log2Fold change and the
# gene names
# we will take the data from the unsorted tables signif_* 
# we save the rownames as a distinct column and we pass it as a column 
# we delete column1 (gene-names) before clustering because it is not needed

# now that we have the data we need to keep only the common rows to all of them
# we use the function intersect by pairs and then all together
# the common genes are rows from common2
# finally we create a dataframe with everything we have


# This should go up above all 
references <- c('wt')		# Cj
#references <- c('wt', 'PC')	# Gg
#contrasts.column <- "Src"
contrasts.column <- "PFU"

for (ref in references) {
    # Find genes that change w.r.t. the reference strain

    # create a convenience text variable to simplify/unify filenames below
    ccol_ref <- paste(contrasts.column, ref, sep='_')

    name <- paste(contrasts.column, ref, levels(as.factor(target[ , contrasts.column]))[1], sep='_')
    common <- rownames(ds.data[[name]]$signif)
    for (i in levels(as.factor(target[ , contrasts.column]))) {
        if (i == ref) next
        name <- paste(ccol_ref, i, sep='_')
        print(name)
        common <- intersect(common, rownames(ds.data[[name]]$signif))
    }
    length(common)	# 681 in Coturnix, 1670 in Gallus

    data.table <- data.frame(genes=common)
    #for (i in c("wt", "PC", "P", "1.0", "0.1")) {
    # compare ref (wt) against the different PFU-infected samples
    for (i in levels(as.factor(target[ , contrasts.column]))) {
    # we do not want to include PC this time
    #for (i in c("P", "1.0", "0.1")) {
        if (i == ref) next
        name <- paste(ccol_ref, i, sep='_')
        print(name)
        data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
    }
    rownames(data.table) <- common

    data.annot <- ds.data[[name]]$signif.annot[common, ]

    par(mfrow=c(1,1))

    set.seed(1963)

    # First have a general look at the methods to get a feeling for the
    # best number of clusters


    dif <- data.table[ , -1]
    by.row <- 1
    by.col <- 2
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)
    # Do a scatterplot matrix
    car::scatterplotMatrix(dif)
    out.png <- paste(out.dir, "/DESeq2/img/", 
               'scatterplot_matrix_', ccol_ref, ".png", 
               sep='')
    as.png( {
            print(car::scatterplotMatrix(dif))
            }, out.png, , overwrite = TRUE)


    # Try to guess the optimum number of K-means clusters
    # NBClust
    out.log <- paste(out.dir, "/DESeq2/cluster/NBClust_", ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    library("NbClust")
    # predict best number of clusters for hierarchical and k-means clustering
    nbc <- NbClust(dif, diss=NULL, 
            distance="euclidean", method="complete", 
            min.nc=3, max.nc=10, 
            index="all", alphaBeale=0.1)
    # 4 for the C. japonica wt vs others data
    # 5 for G. gallus wt vs infected.data
    sink()


    # Let the user see the various clusters and make a decision
    #
    out.log <- paste(out.dir, "/DESeq2/cluster/kmeans/DESeq2_kmeans_all_", 
                     ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    for (i in 2:10) {
        cat("Clustering with K-means (", i, " clusters)\n")
        cl <- kmeans(nor, i)				# NOTE: nor
        print(table(cl$cluster))
        print(fviz_cluster(cl, geom = "point", data=nor))
        out.png <- sprintf("%s/DESeq2/cluster/kmeans/DESeq2_kmeans_%s_nc=%03d.png", 
        	out.dir, ccol_ref, i)
        as.png(fviz_cluster(cl, geom = "point", data=nor,
                    main=paste("K-means clsutering nc =", i) ),
                out.png, overwrite = TRUE)
        #continue.on.enter("Press [ENTER] to continue ")
    }
    sink()
    
    cat("The log and plots for K-means clustering have been saved to\n", 
        "\t", out.dir, "/DESeq2/cluster/kmeans/\n", 
        "please, inspect them and select the best number of clusters\n",
        sep='')
    
    #kmeans.nc <- readline("How many clusters should I use for k-means? ")
    #kmeans.nc <- as.numeric(kmeans.nc)
    # for gg
    #kmeans.nc <- 5
    # for Cj
    kmeans.nc <- 4


    out.log <- paste(out.dir, "/DESeq2/cluster/pam/DESeq2_pam_all_", 
                     ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    for (i in 2:10) {
        cat("\nClustering with PAM (", i, " clusters )\n")
        cl <- pam(nor, i)					# NOTE:nor
        print(table(cl$cluster))
        print(fviz_cluster(cl, geom = "point"))
        out.png <- sprintf("%s/DESeq2/cluster/pam/DESeq2_pam_%s_nc=%03d.png", 
            out.dir, ccol_ref, i)
        as.png(fviz_cluster(cl, geom = "point",
                    main=paste("Partitioning Around Medoids clsutering nc =", i) ), 
                out.png, overwrite = TRUE)
        #continue.on.enter("Press [ENTER] to continue ")
    }
    sink()
    
    cat("The log and plots for PAM clustering have been saved to\n", 
        "\t", out.dir, "/DESeq2/cluster/pam/\n", 
        "please, inspect them and select the best number of clusters\n",
        sep='')
    
    #pam.nc <- readline("How many clusters should I use for PAM? ")
    #pam.nc <- as.numeric(pam.nc)
    # for gg
    pam.nc <- 5

    out.log <- paste(out.dir, "/DESeq2/cluster/dbscan/DESeq2_dbscan_all_", 
                     ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    for (e in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.0)) {
        cat("\nClustering with DBScan ( eps =", e, " )\n")
        cl <- dbscan(dif, eps=e, MinPts=10, showplot=1)	# NOTE: dif
        print(table(cl$cluster))
        print(fviz_cluster(cl, data=dif, 
                     show.clust.cent=FALSE, labelsize=4,
                     ellipse=TRUE, ellipse.type="convex"))
        out.png <- sprintf("%s/DESeq2/cluster/dbscan/DESeq2_dbscan_%s_eps=%03.2f.png", 
                     out.dir, ccol_ref, e)
        as.png(fviz_cluster(cl, data=dif, 
                      show.clust.cent=FALSE, 
                      geom="point", #labelsize=4,
                      ellipse=TRUE, ellipse.type="convex",
                      main=paste("Density based clustering eps =", e) ),
                out.png, overwrite = TRUE)
        #continue.on.enter("Press [ENTER] to continue ")
    }
    sink()
    
    cat("The log and plots for DBscan clustering have been saved to\n", 
        "\t", out.dir, "/DESeq2/cluster/dbscan/\n", 
        "please, inspect them and select the best epsilon value\n",
        sep='')
    
    #dbscan.eps <- readline("Which epsilon should I use for DBscan? ")
    #dbscan.eps <- as.numeric(dbscan.eps)
    # for gg
    dbscan.eps <- 1.5
    # for Cj
    dbscan.eps <- 1.0 # 0.8


    hclust.nc <- 4	# we'll set it by hand for now

    # Now, proceed in detail, method by method, with a deeper and more 
    # detailed analysis

    # possible distances are c("euclidean", "maximum", "manhattan", "canberra",
    #	"binary", "minkowski", "pearson", "spearman", "kendall")
    # possible methods are c("ward.D", "ward.D2", "single", "complete",
    #	"average" (UPGMA), "mcquitty" (WPGMA), "median" (WPGMC), "centroid" (UPGMC))
    out.log <- paste(out.dir, "/DESeq2/cluster/hcluster/DESeq2_hcluster_", contrasts.column, "_", ref, ".txt", sep='')
    sink(out.log, split=T)
    h.cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=ds.data[[name]]$signif.annot[common, ],
        clusters=hclust.nc,
        data.table=dif,
        annotated.data=data.annot
    )
    sink()


    boot <- 100
 
    hclust_out <- paste(out.dir, "/DESeq2/cluster/hclust_", ccol_ref, sep='')
    dir.create(hclust_out, showWarnings=FALSE)
    
    hcut_out <- paste(out.dir, "/DESeq2/cluster/hcut_", ccol_ref, sep='')
    dir.create(hcut_out, showWarnings=FALSE)
    
    out.log <- paste(hcut_out, "/DESeq2_hcut_", 
                     ccol_ref, ".txt", sep='')
    sink(out.log, split=T)
    hc.cl <- cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=sa.data[[name]]$signif.annot[common, ],
        data.table=dif,
        annotated.data=data.annot,
        FUN=hcut,
        clusters=hclust.nc,
        estimate=TRUE,
        gap_bootstrap=boot,
        output.folder=hcut_out
        )
    sink()


    kmeans_out <- paste(out.dir, "/DESeq2/cluster/kmeans_", ccol_ref, sep='')
    dir.create(kmeans_out, showWarnings=FALSE)

    out.log <- paste(kmeans_out, "/DESeq2_kmeans_", 
                     ccol_ref, ".txt", sep='')
    sink(out.log, split=T)
    km.cl <- cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=ds.data[[name]]$signif.annot[common, ],
        data.table=dif,
        annotated.data=data.annot,
        FUN=kmeans,
        algorithm="Hartigan-Wong",
        clusters=kmeans.nc,
        nstart=kmeans.nc*10,
        estimate=T,
        gap_bootstrap=boot,
        output.folder=kmeans_out
        )
    sink()


    pam_out <- paste(out.dir, "/DESeq2/cluster/pam_", ccol_ref, sep='')
    dir.create(pam_out, showWarnings=FALSE)

    out.log <- paste(pam_out, "/DESeq2_pam_", 
                     ccol_ref, ".txt", sep='')
    sink(out.log, split=T)
    pam.cl <- cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=ds.data[[name]]$signif.annot[common, ],
        data.table=dif,
        annotated.data=data.annot,
        FUN=pam,
        distance="euclidean", # see dist()
        clusters=pam.nc,	# n. of clusters
        nstart=pam.nc*10,	# n. of random start sets to choose
        estimate=T,
        gap_bootstrap=boot,
        output.folder=pam_out
        )
    sink()


    dbscan_out <- paste(out.dir, "/DESeq2/cluster/dbscan_", ccol_ref, sep='')
    dir.create(dbscan_out, showWarnings=FALSE)

    out.log <- paste(dbscan_out, "/DESeq2_dbscan_", 
                     ccol_ref, ".txt", sep='')
    sink(out.log, split=T)
    dbs.cl <- cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=ds.data[[name]]$signif.annot[common, ],
        data.table=dif,
        annotated.data=data.annot,
        FUN=dbscan,
        distance="euclidean", # see dist()
        eps=dbscan.eps,
        normalize=F,
        estimate=T,
        gap_bootstrap=boot,
        output.folder=dbscan_out
        )
    sink()



    # Now cluster by experiment
    # -------------------------

## ALREADY DONE ABOVE
##    # Find genes that change w.r.t. the reference strain
##
##     name <- paste(contrasts.column, ref, levels(as.factor(target[ , contrasts.column]))[1], sep='_')
##     common <- rownames(ds.data[[name]]$signif)
##     for (i in levels(as.factor(target[ , contrasts.column]))) {
##         if (i == ref) next
##         name <- paste(contrasts.column, ref, i, sep='_')
##         print(name)
##         common <- intersect(common, rownames(ds.data[[name]]$signif))
##     }
##     length(common)	# 681 in Coturnix, 1670 in Gallus
## 
##     data.table <- data.frame(genes=common)
##     #for (i in c("wt", "PC", "P", "1.0", "0.1")) {
##     for (i in levels(as.factor(target[ , contrasts.column]))) {
##     # we do not want to include PC this time
##     #for (i in c("P", "1.0", "0.1")) {
##         if (i == ref) next
##         name <- paste(contrasts.column, ref, i, sep='_')
##         print(name)
##         data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
##     }
##     rownames(data.table) <- common
## 

    fid <- t(dif)
    means <- apply(fid, by.col, mean)
    sds <- apply(fid, by.col, sd)
    ron <- scale(fid,center=means,scale=sds)

    # here we have a small number of rows and can set a maximum number of clusters
    maxclust <- nrow(fid) - 1
    
    # hclust
    distan = dist(ron, method="euclidean")
    hcl <- hclust(distan)
    plot(hcl,labels=rownames(fid),main='Default from hclust')
    out.png <- sprintf(
            "%s/DESeq2_hcluster_grps_%s.png",
	    hclust_out, ccol_ref)
    as.png(plot(hcl,labels=rownames(fid),main='Default from hclust'), out.png, overwrite = TRUE)

    # kmeans
    fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust)
    out.png <- sprintf(
            "%s/DESeq2_kmeans_grps_silhouette_%s.png",
	    kmeans_out, ccol_ref)
    as.png(fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust), out.png, overwrite = TRUE)
    kcl <- kmeans(ron, centers=3, nstart=100)
    fviz_cluster(kcl, data=ron)
    out.png <- sprintf(
            "%s/DESeq2_kmeans_grps_%s_nc=%03d.png",
	    kmeans_out, ccol_ref, 3)
    as.png(fviz_cluster(kcl, data=ron), out.png, overwrite = TRUE)

    # pam
    fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust)
    out.png <- sprintf(
            "%s/DESeq2_pam_grps_silhouette_%s.png",
	    pam_out, ccol_ref)
    as.png(fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust), out.png, overwrite = TRUE)
    pcl <- pam(ron, k=3, diss=F)
    fviz_cluster(pcl, data=ron)
    out.png <- sprintf(
            "%s/DESeq2_pam_grps_%s_nc=%03d.png",
	    pam_out, ccol_ref, 3)
    as.png(fviz_cluster(pcl, data=ron) , out.png, overwrite = TRUE)

    # dbscan
    # this results in the same two PCs but one cluster
    fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust)
    out.png <- sprintf(
            "%s/DESeq2_dbscan_grps_silhouette_%s.png",
	    dbscan_out, ccol_ref)
    as.png(fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust), out.png, overwrite = TRUE)
    ### JR ### eps should be tuned for each experiment
    dcl <- dbscan(ron, eps=0.2, MinPts=2, showplot=1)
    fviz_cluster(dcl, data=ron)
    out.png <- sprintf(
            "%s/DESeq2_dbscan_grps_%s_eps=%03.2f.png",
	    dbscan_out, ccol_ref, 0.2)
    as.png(fviz_cluster(dcl, data=ron), out.png, overwrite = TRUE)

    # this fails
    tryCatch( {
        NbClust(ron, diss=NULL, 
                distance="euclidean", method="complete", 
                min.nc=3, max.nc=10, 
                index="all", alphaBeale=0.1)
    } )
}

while (sink.number() > 0) sink()

