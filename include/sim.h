#ifndef _sim_h_
#define _sim_h_

typedef struct {
  int mode;//0=modelling;1=FWI, 2=RTM;3= LSRTM, 4=FWI gradient
  int order;//order of FD scheme
  int nt;//number of time steps
  int n1, n2, n3, nb;//number grid points in 1st, 2nd and 3rd coorindate, number of ABC layers
  int n1pad, n2pad, n3pad;//model dimensions 
  int n123, n123pad;
  float d1, d2, d3;//grid spacing
  float volume;//cell volume
  float dt;//temporal sampling
  float freq;//dominant frequency for PML
  float *stf;//source time function
  int eachopt;//stf for each shot
  int aniso;//0=isotropic, 1=VTI, 2=TTI anisotropy
  
  int freesurf;//1, stress-free surface condition;0, no free surface
  int sign_dt;//sign of dt used for back propagation
  float cfl;//CFL number for stability
  float vmax, vmin;//maximum and minimum velocity
  float rhomax, rhomin;//maximum and minimum density

  int ri;//radius for Kaiser windowed sinc interpolation, 2*ri+1 points in total    
  float *pmla, *pmlb;//CPML damping factor
  float *pmla_mh, *pmlb_mh;//CPML damping factor with half grid shift
  float *pmla_ph, *pmlb_ph;//CPML damping factor with half grid shift
  
  float ***vp, ***rho;//model of original size
  float ***ip, ***dm;//ip and dm=dln(ip) in RWI
  float ***epsil, ***delta;//TI medium parameters
  float ***azimul, ***dip;//TTI Euler angles
  float ***kappa, ***buz, ***bux, ***buy;//bulk modulus and buoyancy=1/rho 

  //computing box for forward and adjoint wavefields
  int ibox;//1=add computing box;0=no computing box;
  int *i1min_fwd, *i1max_fwd, *i2min_fwd, *i2max_fwd, *i3min_fwd, *i3max_fwd;
  int *i1min_adj, *i1max_adj, *i2min_adj, *i2max_adj, *i3min_adj, *i3max_adj;

  //store boundary for reverse propagation
  float **face1, **face2, **face3;
  float **face1_, **face2_, **face3_;

  //incident field and PML variables
  float ***vz1, ***vx1, ***vy1, ***p1;
  float ***memD1p1, ***memD2p1, ***memD3p1;
  float ***memD1vz1, ***memD2vx1, ***memD3vy1;
  float ***ph1, ***pv1;//sigma_H, sigma_V
  
  //adjoint field and PML variables
  float ***vz2, ***vx2, ***vy2, ***p2;
  float ***memD1p2, ***memD2p2, ***memD3p2;
  float ***memD1vz2, ***memD2vx2, ***memD3vy2;
    
  //incident field and PML variables (used in RWI and LSRTM for scattering field)
  float ***vz0, ***vx0, ***vy0, ***p0;
  float ***memD1p0, ***memD2p0, ***memD3p0;
  float ***memD1vz0, ***memD2vx0, ***memD3vy0;

  //adjoint field and PML varialbes (used in RWI)
  float ***vz3, ***vx3, ***vy3, ***p3;
  float ***memD1p3, ***memD2p3, ***memD3p3;
  float ***memD1vz3, ***memD2vx3, ***memD3vy3;

  float ***divv, ***dvzdt, ***dvxdt, ***dvydt;//backup divergence of v
  float ***divv0, ***dvzdt0, ***dvxdt0, ***dvydt0;//used only in RWI
  float **dobs, **dcal, **dres;//observed, calculated and residual data
  int dr;//decimation ratio
  int mt;//mt=nt/dr, number of time steps after decimation
  
  
  //data weighting and muting option
  int muteopt;//mute options
  int na;//number of discrete angles
  float da;//angle interval
  int awh;//half window length of angles

  int itcheck;//at step itcheck, check wavefield snapshot

  int nw1, nw2, nw3;//window length
  int cw1, cw2, cw3;//window center=nw/2

  //Primary-secondary (PS) decomposition
  int nr;//ratio between coarse and fine grid
  int i1start, i2start, i3start;
  int i1end, i2end, i3end; 
} sim_t;//simulator for forward and adjoint equations

#endif
