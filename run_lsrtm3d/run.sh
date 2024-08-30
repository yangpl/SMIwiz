#!/bin/bash
#SBATCH --job-name=SMIwiz
#SBATCH --partition=cpu
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=56
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load mpi/2021.6.0
module load mkl/2022.1.0

ulimit -s unlimited
ulimit -l unlimited

#suwaveform type=ricker1 dt=0.0018 ns=1500 fpeak=12 |sustrip > fricker
#makevel nz=201 nx=201 v000=40 > fbathy

echo "#input parameters
===============================================================================
mode=3 //0=modeling; 1=FWI; 2=RTM; 3=LSRTM; 4=FWI gradient; 5=source inversion

acquifile=acqui.txt
vpfile=vp_init true
rhofile=rho_init true
stffile=fricker
bathyfile=fbathy

freesurf=0 //free surface boundary condition
order=8 //FD order in space
nt=1500 //number of time steps
dt=0.0018 //temporal sampling
nb=20  //number of layers for absorbing boundaries
n1=101 //size of input FD model
n2=201 //size of input FD model
n3=201 //size of input FD model
d1=20 //grid spacing of the input FD model
d2=25 //grid spacing of the input FD model
d3=25 //grid spacing of the input FD model

muteopt=1 //0, no mute; 1, front mute; 2, tail mute; 3, front and tail mute
ntaper=10 //number of points for taper
xmute1=0,2000,4000,8000
tmute1=0.2,0.85,1.65,3.3

niter=10
preco=0
family=2
npar=2
idxpar=1,2

" >inputpar.txt

export OMP_NUM_THREADS=1
mpirun -n 256 ../bin/SMIwiz $(cat inputpar.txt)

