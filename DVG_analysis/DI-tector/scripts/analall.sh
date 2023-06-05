
for sample in D_* ; do
    ./analyze_sample.sh $sample
done \
|& tee log/analyze_sample.log

continue


for i in DF*-1 ; do
    exp=${i%-1}		# remove _sample-1 to get experiment name    
    #exp=`echo $i | sed -e 's/_sample.*//g'`
    ./analyze_group.sh $exp


done \
|& tee log/analyze_group.log
