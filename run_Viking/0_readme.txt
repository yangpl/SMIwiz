
LSRTM workflow for 2D Viking Graben Line 12

This folder contains a 2D processing flow for Viking Graben Line 12, including:

    data preprocessing

    RTM

    data-domain LSRTM

    image-domain LSRTM

### Step 1: Data conversion and shot selection (every 300 m, 85 shots)

Convert the original data, select 1 shot every 300 m , and apply filtering:

bash 1_split_shots.sh

to view the first shot:

suximage < dat_0001.su

### Step 2: Build background velocity and density models

Build background vp and rho models using Gardner law and resampling:

bash 2_creat_vp_rho.sh

### Step 3: Source estimation from direct arrivals

Use all selected shots, pick direct waves to estimate Source:

mode = 5
muteopt= 2 (tail mute)
xmute2 = 262, 3237
tmute2 = 0.90, 2.76

Output file: stf

to view sourc wavelet:

xwigb < stf n1=2500

### Step 4: RTM using reflected waves

Use reflected waves for RTM imaging:

mode = 2
stffile = stf
muteopt = 1
xmute1 = 262, 3237
tmute1 = 0.90, 2.76

Output file: param_final_rtm

to view RTM image:

ximage < param_final_rtm n1=201

### Step 5: Vp-Ip Data-domain LSRTM

Run Vp-Ip data-domain LSRTM:

mode = 3
stffile = stf
niter = 15
family = 2
npar = 2
idxpar = 1, 2
muteopt = 1

Output file: param_final

### Step 6: PSF-based image-domain LSRTM

(1) Compute Hessian:

mode = 7
mdopt = 1
muteopt= 1

Output: param_final_m1, param_final_m2

(2) Run PSF-based image-domain LSRTM (single core):

mode = 8

mpirun -np 1 ...

Output: param_final_decon

### Step 7: FFT-based image-domain LSRTM

(1) Compute Hessian:

mode = 7
mdopt = 2
muteopt= 1

Output: param_final_m1, param_final_m2

(2) Run FFT-based image-domain LSRTM (single core):

mode = 9

mpirun -np 1 ...

Output: param_final_decon


