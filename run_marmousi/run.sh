suwaveform type=ricker1 dt=0.002 ns=2500 fpeak=5 |sustrip > fricker
makevel nz=1 nx=461 v000=40 > fbathy
# makevel nz=151 nx=461 v000=1000 > rho
# makevel nz=151 nx=461 v000=1100 > rho_init

echo " 
===============================================================================
mode=1 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion

acquifile=acqui.txt
vpfile=vp_marm_init true
rhofile=rho_init true
stffile=fricker
bathyfile=fbathy

freesurf=1 //free surface boundary condition
fm=5
nt=2500 //number of time steps
dt=0.002 //temporal sampling
nb=20  //Sponge ABC
n1=151 //size of input FD model
n2=461 //size of input FD model
d1=20 //grid spacing of the input FD model
d2=20 //grid spacing of the input FD model
dr=1 //decimation ratio from CFL to Nyquist

===============================================================================
dxwdat=200 //dx for data weighting
xwdat=0,0.3,0.7,1,1,1 //weights for dx increment

===============================================================================
muteopt=0 //0, no mute; 1, front mute; 2, tail mute; 3, front and tail mute
ntaper=20 //number of points for taper
xmute1=50,775.8,2227.4 
tmute1=0.1,0.4,1 
xmute2=50,603.3,1954.3 
tmute2=0.3,0.64,1

===============================================================================
niter=50 //number of iterations using l-BFGS
nls=10 //number of line search per iteration
npair=5 //memory length in l-BFGS
preco=1 //precondition

family=1
npar=2 //only 1 parameter - velocity
idxpar=1,2
bound=1 //bound the inversion parameters
vpmin=1500 //lower bound 
vpmax=5500 //upper bound 
rhomin=1000 //lower bound 
rhomax=3000 //upper bound 
scaleopt=1 //0, no scaling; 1=log; 2=unity-based normalization

gamma1=0 // penalty for Tikhonov regularization 
gamma2=0 // penalty for TV regularization

" >inputpar.txt

mpirun -n 24 ../bin/SMIwiz $(cat inputpar.txt)
