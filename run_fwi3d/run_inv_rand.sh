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

for k in `seq 1 $nbatch`; do
    echo "#input parameters
===============================================================================
mode=1 //0=modeling; 1=FWI; 2=RTM; 3=LSRTM; 4=FWI gradient; 5=source inversion

acquifile=acqui.txt
rhofile=rho
stffile=fricker

freesurf=1 //free surface boundary condition
fm=6
nt=2000 //number of time steps
dt=0.004 //temporal sampling
nb=20  //number of layers for absorbing boundaries
n1=61 //size of input FD model
n2=201 //size of input FD model
n3=201 //size of input FD model
d1=50 //grid spacing of the input FD model
d2=50 //grid spacing of the input FD model
d3=50  //grid spacing of the input FD model

===============================================================================
niter=3 //number of iterations using l-BFGS
nls=10 //number of line search per iteration
npair=5 //memory length in l-BFGS
preco=1 //precondition

family=1
npar=1
idxpar=1 //only 1 parameter - velocity
bound=1 //bound the inversion parameters
vpmin=2000
vpmax=6000
rhomin=1000
rhomax=2800

===============================================================================
dxwdat=200 //dx for data weighting
xwdat=0,0.3,0.7,1,1 //weights for dx increment

" >inputpar$k.txt
    shots=`echo $(echo $(cat batch$k)) | tr ' ' ','`
    echo "shots="$shots >> inputpar$k.txt
    if [ $k -eq 1 ]; then
	echo "vpfile=vp_init">>inputpar$k.txt
    else #use updated vp as initial model
	echo "vpfile=param_final">>inputpar$k.txt 
    fi

    j=0 #count the actual number of shots in this batch
    for i in $(cat batch$k); do
	j=$[$j + 1]
    done
    echo batch$k, 'nshot='$j, 'shot_idx='$shots
    mpirun -n $j ../bin/SMIwiz $(cat inputpar$k.txt)
done
