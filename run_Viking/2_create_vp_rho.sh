#Our initial model has CMP# = 119:50:2019 (total=39),and CMP spacing = 12.5m
fx2=3093.5 #1618.5+(119-1)*12.5 #=3093.5m
dx2=625 #50*12.5 #=625m
nx2=39
#length of X=23.75km (=(NX-1)*DX)
#in depth
fx1=0
dx1=6.25
nx1=640
#length of Z=4km
#We estimate
vmin=1500 #m/s
vmax=4410 #m/s
#We want to resample to dx=dz=h=lambda_min/5=20m (fine)
#dispersion relation requires: fmax = vmin/lambda_min = 15Hz
#stability condition: dt=0.6*h/vmax=0.00272s (fine)
                              
#Resample the model vintz (obtained from velocity analysis)
#output model size
h=20
nz=201 #=ceiling(4000/20)
nx=1568 #=ceiling((3093.5+28250)/20)
unisam2 < vintz  fx1=$fx1 dx1=$dx1 nx1=$nx1 fx2=$fx2 dx2=$dx2 nx2=$nx2 n1=$nz n2=$nx> vp_init
#ximage < vinit.30m n1=$nz d1=$h d2=$h

gfortran gardner.f90
./a.out<<EOF
vp_init rho_init
201 1568
20 0 1500
EOF
