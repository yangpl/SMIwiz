#suwaveform type=ricker1 dt=0.004 ns=500 fpeak=5 |sustrip > fricker

echo "#input parameters
===============================================================================
mode=4 //0=modeling; 1=FWI; 2=RTM; 3=LSRTM; 4=FWI gradient; 5=source inversion

check=1 //check backward reconstructed wavefield is the same as forward field
itcheck=250 //check snapshot at it=itcheck

acquifile=acqui.txt
vpfile=vp
rhofile=rho
stffile=fricker

freesurf=1 //free surface boundary condition
fm=5 //centre frequency of wavelet
nt=500 //number of time steps
dt=0.008 //temporal sampling
nb=20  //number of layers for absorbing boundaries
n1=81 //size of input FD model
n2=181 //size of input FD model
n3=181
d1=50 //grid spacing of the input FD model
d2=50 //grid spacing of the input FD model
d3=50

===============================================================================
niter=50 //number of iterations using l-BFGS
nls=20 //number of line search per iteration
npair=5 //memory length in l-BFGS
preco=1 //precondition

npar=1 //only 1 parameter - velocity
idxpar=1
bound=1 //bound the inversion parameters
vpmin=1200
vpmax=3700
rhomin=1000
rhomax=1000

=============================================================================
dxwdat=200 //dx for data weighting
xwdat=0,0.3,0.7,1,1 //weights for dx increment


" >inputpar.txt

mpirun -n 1 ../bin/SMIwiz $(cat inputpar.txt)

