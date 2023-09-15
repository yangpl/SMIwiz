# SMIwiz
Seismic modelling and imaging wiz

Author: Pengliang Yang, Harbin Institute of Technology, China

Email: ypl.2100@gmail.com

Programming language: C, Shell

Operating System: Linux

Software dependencies: MPI, FFTW

Solution method: High-order finite-dfference time-domain (FDTD) on staggered non-uniform grid

Governing equation: 1st order acoustic wave equation

Code structure:
===============

* src: the source code in .c and .cu/.cuh.

* include: the header files in .h

* doc: documents for theoretic background

* bin: the folderto store executable after compilation

* run_fwi2d: a quick FWI example in layered medium, special acquisition geometry allows you to do FWI using only 1 shot, completing 50 iterations within 2 min.

* run_marmousi: Example for 2D FWI on Marmousi model

* run_rtm: Example for 2D RTM

* run_fbrec3d: Example for reproducing 3D wavefield reconstruction via deimation and interpolation
  Note that you must first generate acquisition file acqui.txt in run_fbrec3d/survey by running: 
  gfortran generate_acqui.f90; ./a.out. 
  Please copy it to /run_fbrec3d before numerical test.

Instructions to run
===================
