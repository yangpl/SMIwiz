batchsize=4 #each group has less than or equal to 128 shots
ntotal=23 #total number of shots used in the inversion
j=0
k=0
for i in `seq 1 $ntotal`; do
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

for k in `seq 1 $nbatch`; do
    echo "#input parameters
===============================================================================
mode=0 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion

acquifile=acqui.txt
vpfile=vp _init
rhofile=rho _init
stffile=fricker

freesurf=1 //free surface boundary condition
fm=25
nt=2000 //number of time steps
dt=0.0005 //temporal sampling
nb=20  //Sponge ABC
n1=141 //size of input FD model
n2=141 //size of input FD model
d1=5 //grid spacing of the input FD model
d2=5 //grid spacing of the input FD model

" >inputpar.txt
    shots=`echo $(echo $(cat batch$k)) | tr ' ' ','`
    echo "shots="$shots >> inputpar.txt

    j=0 #count the actual number of shots in this batch
    for i in $(cat batch$k); do
	j=$[$j + 1]
    done
    echo batch$k, 'nshot='$j, 'shot_idx='$shots
    mpirun -n $j ../bin/SMIwiz $(cat inputpar.txt)
done
