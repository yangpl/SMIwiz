suwaveform type=ricker1 dt=0.004 ns=2000 fpeak=6 |sustrip > fricker
makevel nz=201 nx=201 v000=100 > fbathy

echo "#input parameters
===============================================================================
mode=0 //0=modeling; 1=FWI; 2=RTM; 3=LSRTM; 4=FWI gradient; 5=source inversion

acquifile=acqui.txt
vpfile=vp
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

" >inputpar.txt

mpirun -n 1 ../bin/SMIwiz $(cat inputpar.txt)

