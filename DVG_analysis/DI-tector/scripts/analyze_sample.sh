#!/bin/bash
#
# usage: ../analyze.sh (called from within a subdirectory)

# If $1 argument is provided will be set as sample. Otherwise, sample = none.
sample=${1:-none}

if [ ! -d $sample ] ; then
    exit
fi

# Enter working directory (see end of sctipt)
cd $sample
echo "Sample $sample"    

sampledir=`realpath .`
base=`basename $sampledir`

# format is:
# =================================
# = 3' cb/sb DVG
# =================================
# DVG's type      Length  BP_Pos  RI_Pos  Delta_Positions Ref     Counts  %_to_Virus
#.... data ....

# columns in counts
#DVG's type	Length	BP_Pos	RI_Pos	Delta_Positions	Ref	Counts	%_to_Virus
# columns in output sorted
#DVG's type	Length	BP_Pos	RI_Pos	Delta_Positions	Segmentation	MAPQ_F	MAPQ_L	RNAME_F	RNAME_L	CIGAR_F	CIGAR_L	MD_CIGAR_F	MD_CIGAR_L	POS_F	POS_L	QNAME_F	SEQ_FL_ori

# BP is breaking point
# RI is reinitiation site
# Count is counted DI entities
# Type is one of
# 5cb	5' snap back/copy back		(5' extreme hairpin (sb) or panhandle (cb))
# 3cb	3' snap back/copy back		(3' extreme hairpin (sb) or panhandle (db))
# Insertion
# Deletion
#
# F	first read segment
# L	last read segment
#
# QNAME	query name
# RNAME reference name


r1_counts=`realpath ${base}_R1_counts.txt`
r1_outsrt=`realpath ${base}_R1_output_sorted.txt`

r2_counts=`realpath ${base}_R2_counts.txt`
r2_outsrt=`realpath ${base}_R2_output_sorted.txt`


mkdir -p out

# obtain list of DVG types considered
cat $r1_outsrt \
| grep -v ^DVG \
| cut -f1 \
| sort \
| uniq \
> out/DVG_types

cat $r1_counts | grep -v -E "^=|^$|^DVG" \
| sed -e 's/None.*/None	0	0	0	0	X	0	0/' \
> out/R1counts.txt

echo "DVG's type	Length	BP_Pos	RI_Pos	Delta_Positions	Ref	Counts	%_to_Virus" \
> out/R1counts+header.txt
cat out/R1counts.txt >> out/R1counts+header.txt


cat $r2_counts | grep -v -E "^=|^$|^DVG" \
| sed -e 's/None.*/None	0	0	0	0	X	0	0/' \
> out/R2counts.txt

echo "DVG's type	Length	BP_Pos	RI_Pos	Delta_Positions	Ref	Counts	%_to_Virus" \
> out/R2counts+header.txt
cat out/R2counts.txt >> out/R2counts+header.txt

cat $r1_outsrt | grep -v '^DVG' > out/R1outsrt.txt
cat $r2_outsrt | grep -v '^DVG' > out/R2outsrt.txt

# counts
grep 'A|A' $r1_counts > out/R1counts_AA.txt
grep 'A|B' $r1_counts > out/R1counts_AB.txt
grep 'B|A' $r1_counts > out/R1counts_BA.txt
grep 'B|B' $r1_counts > out/R1counts_BB.txt

grep 'A|A' $r2_counts > out/R2counts_AA.txt
grep 'A|B' $r2_counts > out/R2counts_AB.txt
grep 'B|A' $r2_counts > out/R2counts_BA.txt
grep 'B|B' $r2_counts > out/R2counts_BB.txt

# output sorted
grep '	A	A' $r1_outsrt > out/R1outsrt_AA.txt
grep '	A	B' $r1_outsrt > out/R1outsrt_AB.txt
grep '	B	A' $r1_outsrt > out/R1outsrt_BA.txt
grep '	B	B' $r1_outsrt > out/R1outsrt_BB.txt

grep '	A	A' $r2_outsrt > out/R2outsrt_AA.txt
grep '	A	B' $r2_outsrt > out/R2outsrt_AB.txt
grep '	B	A' $r2_outsrt > out/R2outsrt_BA.txt
grep '	B	B' $r2_outsrt > out/R2outsrt_BB.txt

# get count of detected DVG's
cat $r1_outsrt \
| grep -v ^DVG \
| cut -f1 \
| sort \
| uniq -c \
> out/R1_DVGs.cnt

cat $r2_outsrt \
| grep -v ^DVG \
| cut -f1 \
| sort \
| uniq -c \
> out/R2_DVGs.cnt


for i in out/R[12]counts_[AB][AB].txt ; do 
    echo -n "$i	"
    cat $i  \
    | cut -f7 \
    | while read val ; do sum=$((sum + val)) ; echo $sum ; done \
    | tail -1
done > out/count_totals_by_segment.txt

wc -l out/R[12]outsrt_[AB][AB].txt | awk '{print $2 "\t" $1}' > out/outsrt_totals_by_segment.txt

for i in out/R[12]counts_[AB][AB].txt ; do 
    echo "$i"
    while read dvgtype ; do
        echo -ne "\t$dvgtype\t"
        grep "^$dvgtype" $i \
        | cut -f7 \
        | while read val ; do sum=$((sum + val)) ; echo $sum ; done \
        | tail -1
    done < out/DVG_types
done > out/count_by_segment_and_DVG.txt

for i in out/R[12]outsrt_[AB][AB].txt ; do 
    echo "$i" ; 
    cut -f1 $i | sort | uniq -c \
    | while read count type ; do echo -e "\t${type}\t$count" ; done
done > out/outsrt_by_segment_and_DVG.txt

while read dvgtype ; do
    typnam=`echo $dvgtype | sed -e 's/ /_/g' -e 's-/-|-g' | tr -d "'"`
    grep "^$dvgtype" out/R1counts.txt > out/R1_counts_${typnam}.txt
    grep "^$dvgtype" out/R2counts.txt > out/R2_counts_${typnam}.txt

    grep "^$dvgtype" out/R1outsrt.txt > out/R1_outsrt_${typnam}.txt
    grep "^$dvgtype" out/R2outsrt.txt > out/R2_outsrt_${typnam}.txt

done < out/DVG_types


wc -l out/R[12]_counts_*.txt | awk '{print $2 "\t" $1}' > out/count_totals_by_DVGtype.txt
wc -l out/R[12]_outsrt_*.txt | awk '{print $2 "\t" $1}' > out/outsrt_totals_by_DVGtype.txt


# prepare plots
#---------------

cat out/R1counts+header.txt R2counts.txt | tr "'" '.' > out/zcounts12
cat out/R2counts+header.txt | tr "'" '.' > out/zcounts2
cat out/R1counts+header.txt | tr "'" '.' > out/zcounts1

R --vanilla <<END
r1 <- read.table("out/zcounts1", header=TRUE, 
  stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
r2 <- read.table("out/zcounts2", header=TRUE, 
  stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
r12 <- read.table("out/zcounts12", header=TRUE, 
  stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(lattice)
png("out/plot_R1_byDVG.png", width=1024, height=1024)
xyplot(Counts ~ BP_Pos | DVG.s.type, groups=Ref, 
  type="p", 
  pch=16, 
  auto.key=list(border=TRUE), 
  par.settings=simpleTheme(pch=16), 
  scales=list(x=list(relation='free'), y=list(relation='free')), 
  data=r1)
dev.off()
png("out/plot_R2_byDVG.png", width=1024, height=1024)
xyplot(Counts ~ BP_Pos | DVG.s.type, groups=Ref, 
  type="p", 
  pch=16, 
  auto.key=list(border=TRUE), 
  par.settings=simpleTheme(pch=16), 
  scales=list(x=list(relation='free'), y=list(relation='free')), 
  data=r2)
dev.off()
png("out/plot_R12_byDVG.png", width=1024, height=1024)
xyplot(Counts ~ BP_Pos | DVG.s.type, groups=Ref, 
  type="p", pch=16, 
  auto.key=list(border=TRUE), 
  par.settings=simpleTheme(pch=16), 
  scales=list(x=list(relation='free'), y=list(relation='free')), 
  data=r12)
dev.off()

png("out/plot_R1_bySegment.png", width=1024, height=1024)
xyplot(Counts ~ BP_Pos | Ref, groups=DVG.s.type, 
  type="p", 
  pch=16, 
  auto.key=list(border=TRUE), 
  par.settings=simpleTheme(pch=16), 
  scales=list(x=list(relation='free'), y=list(relation='free')), 
  data=r1)
dev.off()
png("out/plot_R2_bySegment.png", width=1024, height=1024)
xyplot(Counts ~ BP_Pos | Ref, groups=DVG.s.type,
  type="p", 
  pch=16, 
  auto.key=list(border=TRUE), 
  par.settings=simpleTheme(pch=16), 
  scales=list(x=list(relation='free'), 
  y=list(relation='free')), data=r2)
dev.off()
png("out/plot_R12_bySegment.png", width=1024, height=1024)
xyplot(Counts ~ BP_Pos | Ref, groups=DVG.s.type,
  type="p", pch=16, 
  auto.key=list(border=TRUE), 
  par.settings=simpleTheme(pch=16), 
  scales=list(x=list(relation='free'), 
  y=list(relation='free')), data=r12)
dev.off()

END

rm out/zcounts1 out/zcounts2 out/zcounts12

# exit working directory (see start of script)
cd -
