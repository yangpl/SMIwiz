#ifndef _sim_h_
#define _sim_h_

typedef struct {
  int mode; //0=modelling; 1=FWI, 2=RTM; 3= LSRTM, 4=FWI gradient
  int order; //order of FD scheme
  int nt; /* number of time steps */
  int nt_verb; //verbose display every nt_verb time steps
  int n1, n2, n3, nb;
  int n1pad, n2pad, n3pad;/* model dimensions */
  int n123, n123pad;
  float d1, d2, d3; /* grid spacing */
  float volume;/* cell volume */
  float dt; /* temporal sampling */
  float freq; /* dominant frequency for PML */
  float *stf; /* source time function */

  int ps;//1=PS decomposition; 0=no PS decomposition
  int freesurf;//1, stress-free surface condition; 0, no free surface
  int sign_dt;
  float cfl;
  float vmax, vmin;/* maximum and minimum velocity */
  float rhomax, rhomin;/* maximum and minimum density */

  int ri;//radius for Kaiser windowed sinc interpolation, 2*ri+1 points in total    
  float *pmla, *pmlb;  /* CPML profile */
  
  float ***vp, ***rho;/* model of original size */
  float ***ip_smooth;/* smooth Ip model for RWI */
  float ***kappa, ***buz, ***bux, ***buy;/* buoyancy=1/rho and bulk modulus */

  //computing box for forward and adjoint wavefields
  int ibox;//1=add computing box; 0=no computing box;
  int *i1min_fwd, *i1max_fwd, *i2min_fwd, *i2max_fwd, *i3min_fwd, *i3max_fwd;
  int *i1min_adj, *i1max_adj, *i2min_adj, *i2max_adj, *i3min_adj, *i3max_adj;

  //store boundary for reverse propagation
  float **face1, **face2, **face3;
  float **face1_, **face2_, **face3_;

  float ***vz1, ***vx1, ***vy1, ***p1;/* incident field */
  float ***memD1p1, ***memD2p1, ***memD3p1;/* memory variables for incident field */
  float ***memD1vz1, ***memD2vx1, ***memD3vy1;
  
  float ***vz2, ***vx2, ***vy2, ***p2;/* adjoint field */
  float ***memD1p2, ***memD2p2, ***memD3p2;/* memory variables for adjoint field */
  float ***memD1vz2, ***memD2vx2, ***memD3vy2;

  float ***vz0, ***vx0, ***vy0, ***p0;/* scattering field */
  float ***memD1p0, ***memD2p0, ***memD3p0;  /* memory variables for scattering field */
  float ***memD1vz0, ***memD2vx0, ***memD3vy0;

  float ***divv, ***dvzdt, ***dvxdt, ***dvydt;/* backup divergence of v */
  float **dcal, **dobs, **dres;
  int dr; /* decimation ratio */
  int mt; /* mt=nt/dr */

  //data weighting and muting option
  int muteopt;
  int na;//number of discrete angles
  float da;//angle interval
  int awh;//half window length of angles

  //check snapshot for radiation pattern and source field reconstruction
  int check, itcheck;

  int nw1, nw2, nw3;//window length
  int cw1, cw2, cw3;//window center=nw/2

  
  int nr;//ratio between coarse and fine grid
  int i1start, i2start, i3start;
  int i1end, i2end, i3end;  
} sim_t; //simulator for forward and adjoint equations


#endif
