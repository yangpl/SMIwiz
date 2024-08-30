#ifndef _fwi_h_
#define _fwi_h_

typedef struct {
  int n; //total number of unknowns in FWI
  int family;//1=vp-rho; 2=vp-ip
  int npar;//number of parameters to invert by FWI
  int *idxpar;//index of the inversion parameters
  float rhomin, rhomax, vpmin, vpmax;
  float *minpar, *maxpar;

  float **bathy; //bathymetry to prescribe water bottom
  int **ibathy;//index of bathymetry at depth z
  float transition; //depth of transition zone
  int itransition; //index of transition zone

  int niter; //maximum number of iterations
  int iter; //iteration index
  int restart;//restart iterations
  float fcost, fcost_dat, fcost_mod;
  float fcost_mod1, fcost_mod2;
  
  float alpha;//scaling factor for misfit function
  int firstgrad;//if first gradient computation, firstgrad=1, otherwise firstgrad=0
  float gamma1;//penalty parameter for Tikhonov regularization
  float gamma2;//penalty parameter for TV regularization
  int r1, r2, r3, repeat;//smoothing parameters

  int isrcpershot;// estimate source per shot or not
  float *hess;//pseudo-Hessian

  int rwi;//1=RWI; 0=FWI
  int preco;//0=no precondition; 1=depth precondition; 2=pseudo-Hessian precondition
  int mdopt;//options for migration deconvolution 
} fwi_t;

#endif
