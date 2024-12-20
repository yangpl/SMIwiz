#!/bin/bash
#SBATCH --job-name=SMIwiz
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load mpi/2021.6.0
module load mkl/2022.1.0

ulimit -s unlimited
ulimit -l unlimited

#==========================================================
#suwaveform type=ricker1 dt=0.0005 ns=2000 fpeak=25 |sustrip > fricker
#makevel nz=141 nx=141 v000=1000 > rho 

echo "#input parameters
===============================================================================
mode=1 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion

acquifile=acqui.txt
vpfile=vp_init
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


===============================================================================
niter=50 //number of iterations using l-BFGS
nls=10 //number of line search per iteration
npair=5 //memory length in l-BFGS
preco=2 //precondition

family=1
npar=1
idxpar=1 ,2 //only 1 parameter - velocity
bound=1 //bound the inversion parameters
vpmin=1500
vpmax=3700
rhomin=1800
rhomax=2100

===============================================================================
dxwdat=100 //dx for data weighting
xwdat=0,0.3,0.7,1,1 //weights for dx increment

===============================================================================
muteopt=0 //0, no mute; 1, front mute; 2, tail mute; 3, front and tail mute
ntaper=20 //number of points for taper
xmute1=50,775.8,2227.4 
tmute1=0.1,0.4,1 
xmute2=50,603.3,1954.3 
tmute2=0.3,0.64,1

r1=2
r2=2
repeat=2

" >inputpar.txt

export OMP_NUM_THREADS=1
mpirun -n 1 ../bin/SMIwiz $(cat inputpar.txt)

