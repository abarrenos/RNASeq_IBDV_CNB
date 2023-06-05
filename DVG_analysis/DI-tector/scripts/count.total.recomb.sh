for i in [01m]*[123] ; do
    cnt=0
    # compute total recombinant reads (sum of all types)
    # tail -n 4 will produce
    # 5' cb/sb DVGs         : 406 reads
    # 3' cb/sb DVGs         : 1,235 reads
    # DVGs with deletion    : 8,207 reads
    # DVGs with insertion   : 37,378 reads
    # we remove everything up to the colon, the comma and reads$
    # then we add the numbers and outut the name and total count
    tail -q -n 4 $i/*_R1_DI-tector_summary.txt \
    | sed -e 's/^.*://g' -e 's/reads$//g' \
    | tr -d ',' \
    | while read n ; do
        cnt=$((cnt + n))
        echo "${i}_R1	$cnt" > $i/out/tot.rec.reads.R1
    done
    # same for R2
    tail -q -n 4 $i/*_R2_DI-tector_summary.txt \
    | sed -e 's/^.*://g' -e 's/reads$//g' \
    | tr -d ',' \
    | while read n ; do
        cnt=$((cnt + n))
        echo "4{i}_R2	$cnt" > $i/out/tot.rec.reads.R2
    done
    # and now for both together
    tail -q -n 4 $i/*_R[12]_DI-tector_summary.txt \
    | sed -e 's/^.*://g' -e 's/reads$//g' \
    | tr -d ',' \
    | while read n ; do
        cnt=$((cnt + n))
        echo "$i	$cnt" > $i/out/tot.rec.reads.R12
    done
done
# put all together
echo "ID	tot.count" > rec.cnt
cat [01m]*/out/tot.rec.reads.R12>> rec.cnt
