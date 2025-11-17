#=====================================================                                                                                                
makevel nz=1 nx=1568 v000=300 > fbathy

echo " 
===============================================================================
mode=5 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion
suopt=1

vpfile=vp_init
rhofile=rho_init
stffile=stf
bathyfile=fbathy

freesurf=0 //free surface boundary condition
order=8
nt=2500 //number of time steps
dt=0.0024 //temporal sampling
nb=20  //Sponge ABC
n1=201 //size of input FD model
n2=1568 //size of input FD model
d1=20 //grid spacing of the input FD model
d2=20 //grid spacing of the input FD model

===============================================================================
dxwdat=200 //dx for data weighting
xwdat=1,1,1 //weights for dx increment

muteopt=1 //0, no mute; 1, front mute; 2, tail mute; 3, front and tail mute
ntaper=20 //number of points for taper
xmute1=262,3237
tmute1=0.9,3

===============================================================================
niter=15 //number of iterations using l-BFGS
preco=1 //precondition

family=2
npar=2 //only 1 parameter - velocity
idxpar=1,2

" >inputpar.txt

export OMP_NUM_THREADS=1
mpirun -np 85 ../bin/SMIwiz $(cat inputpar.txt)
