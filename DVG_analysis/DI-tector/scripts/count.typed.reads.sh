
echo "ID	type	count"	> counts/rec.typed.cnt
for i in 000.0pfu 000.1pfu 001.0pfu 010.0pfu 100.0pfu mock ; do
  tot_cb5=0
  tot_cb3=0
  tot_del=0
  tot_ins=0
  while read cb5 cb3 del ins err ; do
        #echo $cb5 $cb3 $del $ins
        tot_cb5=$((tot_cb5 + cb5))
        #echo "$i	5' cb/sb	$tot_cb5"
        tot_cb3=$((tot_cb3 + cb3))
        #echo "$i	3' cb/sb	$tot_cb3"
        tot_del=$((tot_del + del))
        #echo "$i	deletion	$tot_del"
        tot_ins=$((tot_ins + ins))
        #echo "$i	insertion	$tot_ins"
  done <<< $(
  for j in 1 2 3 ; do
    replica=${i}-${j}
    # compute total recombinant reads --sum of each type--
    # tail -n 5 will produce
    # 5' cb/sb DVGs         : 406 reads
    # 3' cb/sb DVGs         : 1,235 reads
    # DVGs with deletion    : 8,207 reads
    # DVGs with insertion   : 37,378 reads
    # Errors                : 0
    # we remove everything up to the colon, the comma and reads$
    # then we put the four numbers in a single line sent to stdout
    tail -q -n 5 $replica/*_R1_DI-tector_summary.txt \
    | sed -e 's/^.*://g' -e 's/reads$//g' \
    | tr -d ',' \
    | tr -d '\n'
    echo ''
    tail -q -n 5 $replica/*_R1_DI-tector_summary.txt \
    | sed -e 's/^.*://g' -e 's/reads$//g' \
    | tr -d ',' \
    | tr -d '\n'
    echo ''
  done )
    
  echo "$i	\"5' cb/sb\"	$tot_cb5"
  echo "$i	\"3' cb/sb\"	$tot_cb3"
  echo "$i	\"deletion\"	$tot_del"
  echo "$i	\"insertion\"	$tot_ins"
done \
>> counts/rec.typed.cnt
