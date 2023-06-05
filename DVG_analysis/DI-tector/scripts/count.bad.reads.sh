for i in [01m]*[123] ; do
    R1_A_AB=`grep -c '	A	[AB]	' $i/*_R1_DI-tector_output_sorted.txt`
    R1_A_A=`grep -c '	A	A	' $i/*_R1_DI-tector_output_sorted.txt`
    R1_A_B=`grep -c '	A	B	' $i/*_R1_DI-tector_output_sorted.txt`
    R1_B_A=`grep -c '	B	A	' $i/*_R1_DI-tector_output_sorted.txt`
    R1_B_B=`grep -c '	B	B	' $i/*_R1_DI-tector_output_sorted.txt`
    R1_B_AB=`grep -c '	B	[AB]	' $i/*_R1_DI-tector_output_sorted.txt`
    echo "$i	$R1_A_A	$R1_A_B	$R1_A_AB	$R1_B_A	$R1_B_B	$R1_B_AB" \
    	> $i/out/tot.bad.R1

    R2_A_AB=`grep -c '	A	[AB]	' $i/*_R2_DI-tector_output_sorted.txt`
    R2_A_A=`grep -c '	A	A	' $i/*_R2_DI-tector_output_sorted.txt`
    R2_A_B=`grep -c '	A	B	' $i/*_R2_DI-tector_output_sorted.txt`
    R2_B_A=`grep -c '	B	A	' $i/*_R2_DI-tector_output_sorted.txt`
    R2_B_B=`grep -c '	B	B	' $i/*_R2_DI-tector_output_sorted.txt`
    R2_B_AB=`grep -c '	B	[AB]	' $i/*_R2_DI-tector_output_sorted.txt`
    echo "$i	$R2_A_A	$R2_A_B	$R2_A_AB	$R2_B_A	$R2_B_B	$R2_B_AB" \
    	> $i/out/tot.bad.R2
        
    RA_A=$(( R1_A_A + R2_A_A ))
    RA_B=$(( R1_A_B + R2_A_B ))
    RA_AB=$(( R1_A_AB + R2_A_AB ))
    RB_A=$(( R1_B_A + R2_B_A ))
    RB_B=$(( R1_B_B + R2_B_B ))
    RB_AB=$(( R1_B_AB + R2_B_AB ))
    echo "$i	$RA_A	$RA_B	$RA_AB	$RB_A	$RB_B	$RB_AB" \
    	> $i/out/tot.bad.R12
done

cat > rec.bad.cnt <<END
# ID	experiment (PFU-replica)					
# A-A	bad reads starting in A and ending in A					
# A-B	bad reads starting in A and ending in B					
# A-AB	bad reads starting in A (and ending in either A or B)					
# B-B	bad reads starting in B and ending in B					
# B-A	bad reads starting in B and ending in A					
# B-AB	bad reads starting in B (and ending in either A or B)					
# NOTE: needs further processing with a spreadsheet to calculate totals
# for the three replicas together
# 						
ID	A-A	A-B	A-AB	B-A	B-B	B-AB
END
cat [01m]*/out/tot.bad.R12 >> rec.bad.cnt
