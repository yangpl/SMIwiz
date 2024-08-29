suwaveform type=ricker1 dt=0.001 ns=4000 fpeak=12 |sustrip > fricker
makevel nz=1 nx=767 v000=40 > fbathy

echo " 
===============================================================================
mode=3 0 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion

acquifile=acqui.txt
vpfile=vp_init true
rhofile=rho_init true
stffile=fricker
bathyfile=fbathy

freesurf=0 //free surface boundary condition
order=8
nt=4000 //number of time steps
dt=0.001 //temporal sampling
nb=20  //Sponge ABC
n1=251 //size of input FD model
n2=767 //size of input FD model
d1=12 //grid spacing of the input FD model
d2=12 //grid spacing of the input FD model

===============================================================================
dxwdat=200 //dx for data weighting
xwdat=1,1,1 //weights for dx increment

muteopt=1 //0, no mute; 1, front mute; 2, tail mute; 3, front and tail mute
ntaper=20 //number of points for taper
xmute1=0,1800,4000,8000
tmute1=0.25,1.37,2.85,5.7

===============================================================================
niter=15 //number of iterations using l-BFGS
preco=0 //precondition

family=2
npar=2 //only 1 parameter - velocity
idxpar=1,2

" >inputpar.txt

export OMP_NUM_THREADS=1
mpirun -n 55 ../bin/SMIwiz $(cat inputpar.txt)
