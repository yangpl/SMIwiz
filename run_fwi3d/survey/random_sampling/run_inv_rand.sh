batchsize=50 #each group has less than or equal to 128 shots
ntotal=400 #total number of shots used in the inversion
j=0
k=0
for i in `shuf -i 1-$ntotal -n $ntotal`; do
    if [ `expr $j % $batchsize` -eq 0 ]; then 
        k=$[ $k + 1 ]
	echo 'batch' $k
    	echo $i > batch$k #write into a new file
    else
    	echo $i >> batch$k #append into an existing file
    fi
    j=$[$j + 1]
    echo $i
done
nbatch=$k
echo 'nbatch=' $nbatch

