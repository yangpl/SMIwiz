# SMIwiz
Seismic Modelling and Imaging wiz: An integrated toolbox

Author: Pengliang Yang, Harbin Institute of Technology, China

Email: ypl.2100@gmail.com

Programming language: C, Shell, Fortran

Operating System: Linux

Software dependencies: MPI, FFTW

Solution method: High-order finite-dfference time-domain (FDTD) for modelling
on staggered grid; Quasi-Newton LBFGS algorithm for nonlinear optimization;
line search to estimate step length based on Wolfe condition

Governing equation: 1st order acoustic wave equation

## Credit

* Pengliang Yang, SMIwiz: An integrated toolbox for multidimensional seismic modelling and imaging, Computer Physics Communications 2024, volume 295, 109011 [doi:10.1016/j.cpc.2023.109011](https://doi.org/10.1016/j.cpc.2023.109011)

* Zhengyu Ji and Pengliang Yang: SMIwiz-2.0: Extended functionalities for wav­field decomposition, linearized and nonlinear inversion, Computer Physics Communications 2025, 109503 [doi:10.1016/j.cpc.2025.109503](https://doi.org/10.1016/j.cpc.2025.109503) 

* Pengliang Yang and Zhengyu Ji. A comparative study of data- and image- domain LSRTM under velocity-impedance parametrization. Computers & Geosciences, 208:106091, 2026 [doi:10.1016/j.cageo.2025.106091]{https://doi.org/10.1016/j.cageo.2025.106091}

A vedio has been recorded to explain the design of SMIwiz, thanks to the seminar invitation from the editor of Computer Physics Communictions: https://cassyni.com/events/P4W1QfiGXffuSf6Rv5VzJZ

## Code structure

* src: the source code in .c 

* include: the header files in .h

* doc: documents for theoretic background

* bin: the folder to store executable after compilation

* run_fwi2d: a quick FWI example in layered medium, special acquisition geometry allows you to do FWI using only 1 shot, completing 50 iterations within 2 min on your laptop.

* run_marmousi: Example for 2D FWI on Marmousi model

* run_rtm: Example for 2D RTM

* run_fbrec3d: Example for reproducing 3D wavefield reconstruction via deimation and interpolation

  Note that you must first generate acquisition file acqui.txt in run_fbrec3d/survey by running: 
  gfortran generate_acqui.f90; ./a.out. 
  Please copy it to /run_fbrec3d before numerical test. (This is because acqui.txt for this test is too large (400 shots * 10000 receivers) to be uploaded)

* run_lsrtm2d: Example to do two-parameter LSRTM in data-domain

* run_lsrtm3d: Example of a 3D two-parameter LSRTM in data-domain on Overthrust model

* run_decon_1par: Example for migration deconvolution using PSF and FFT-based Wiener filtering
  Note that you have to first run in mode=7, then run in mode=8 (with mdopt=1 for PSF) or run mode=9 (with mdopt=2 for FFT-Wiener).
  
* run_updown: Example for up-going and down-going wavefield decomposition

* run_adcig: Example for extracting ADCIG using Poynting vector.
git checkout d27813a8d3752639df2aea199e534738f75d8c7f

* run_Viking: A real data exampling using Vking Graben Line 12 for linearized inversion (data-domain LSRTM, image-domain LSRTM using PSF and nonstationary deblurring filter)

  Please follow the workflow listed in  /run_Viking/0_readme.txt to reproduce this example. Version: git checkout f3d2d693ae42d0689348cbaf527ba90a99943644 (2025/11/20 updated for VTI modelling)


## Instructions to run

1. go to /src and compile:

cd /src;

make
	
2. go to running template and test:

cd ../run_fwi2d

bash run.sh

To run RTM, you need:

cd ../run_rtm

bash run.sh


To do FWI/RTM, you first need to generate observed data from true models using mode=0.

Then, you start FWI in mode=1 with initial models (or RTM in mode=2).

	
