suwaveform type=ricker1 dt=0.002 ns=2500 fpeak=10 |sustrip > fricker
makevel nz=1 nx=461 v000=40 > fbathy

echo " 
===============================================================================
mode=3 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion

acquifile=acqui.txt
vpfile=vp_marm_init true
rhofile=rho_init true
stffile=fricker
bathyfile=fbathy

freesurf=1 //free surface boundary condition
fm=10
order=8
nt=2500 //number of time steps
dt=0.002 //temporal sampling
nb=20  //Sponge ABC
n1=151 //size of input FD model
n2=461 //size of input FD model
d1=20 //grid spacing of the input FD model
d2=20 //grid spacing of the input FD model


===============================================================================
dxwdat=200 //dx for data weighting
xwdat=1,1,1 //weights for dx increment

muteopt=1 //0, no mute; 1, front mute; 2, tail mute; 3, front and tail mute
ntaper=20 //number of points for taper
xmute1=0,3000,6781
tmute1=0.3,2.3,5.
xmute2=0,3000,6781
tmute2=0.3,2.3,5.


===============================================================================
niter=50 //number of iterations using l-BFGS
nls=10 //number of line search per iteration
npair=5 //memory length in l-BFGS
preco=2 //precondition

family=2
npar=2 //only 1 parameter - velocity
idxpar=1,2
bound=1 //bound the inversion parameters
vpmin=1500 //lower bound 
vpmax=5500 //upper bound 
rhomin=1000 //lower bound 
rhomax=3000 //upper bound 


gamma1=0 // penalty for Tikhonov regularization 
gamma2=0 // penalty for TV regularization

" >inputpar.txt

mpirun -n 24 ../bin/SMIwiz $(cat inputpar.txt)
