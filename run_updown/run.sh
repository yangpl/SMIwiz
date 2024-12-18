#suwaveform type=ricker1 dt=0.0005 ns=2000 fpeak=25 |sustrip > fricker
#makevel nz=141 nx=141 v000=1000 > rho 

echo "#input parameters
===============================================================================
mode=10 //0, forward modeling; 1, FWI; 2, RTM; 3, LSRTM; 4, FWI gradient; 5, source inversion
check=1
itcheck=300
ntaper=10

acquifile=acqui.txt
vpfile=vp 
rhofile=rho 
stffile=fricker

freesurf=0 //free surface boundary condition
fm=25
nt=600 //number of time steps
dt=0.0005 //temporal sampling
nb=20  //Sponge ABC
n1=141 //size of input FD model
n2=141 //size of input FD model
d1=5 //grid spacing of the input FD model
d2=5 //grid spacing of the input FD model
order=8


" >inputpar.txt

mpirun -n 2 ../bin/SMIwiz $(cat inputpar.txt)

