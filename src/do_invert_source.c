/* Source wavelet estimation according to Pratt 1999 Geophysics
 *---------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *-------------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include <fftw3.h>
#include "sim.h"
#include "acq.h"
#include "fwi.h"
 
void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_free(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim, float ***vp, float ***rho, float ***kappa, float ***buz, float ***bux, float ***buy);
void extend_model_free(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_free(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_free(sim_t *sim);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);

void write_data(sim_t *sim, acq_t *acq);

void do_invert_source(sim_t *sim, acq_t *acq)
/*< invert source time function in FWI >*/
{
  FILE *fp;
  int it,irec,ntpow2;
  char *stffile;

  if(!getparstring("stffile",&stffile)) err("must give stffile= ");
  //0=one source for all shots; 1=one source for each shot
  
  //-----------------------------------------------------------------------
  //1. simulate synthetic data (Green's function) using Dirac delta function
  //-----------------------------------------------------------------------
  memset(sim->stf, 0, sim->nt*sizeof(float)); sim->stf[0] = 1;
  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  extend_model(sim, sim->vp, sim->rho, sim->kappa, sim->buz, sim->bux, sim->buy);
  computing_box_init(acq, sim, 0);
  fdtd_init(sim, 1);
  fdtd_null(sim, 1);
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);
  }
  extend_model_free(sim);
  fdtd_free(sim, 1);
  cpml_free(sim);
  computing_box_free(sim, 0);

  //----------------------------------------------------------------------
  //2. estimate wavelet in frequency domain, see eqn 17 Pratt Geophysics 1999
  //----------------------------------------------------------------------
  ntpow2=1;
  while(ntpow2<2*sim->nt) { ntpow2 *= 2; } //ntpow2=2^power, ntpow2>=nt

  float *den=alloc1float(ntpow2);
  fftw_complex *num=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_complex *ft_dcal=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_complex *ft_dobs=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_complex *tmp=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_plan fft=fftw_plan_dft_1d(ntpow2, tmp, tmp, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan ifft=fftw_plan_dft_1d(ntpow2, tmp, tmp, FFTW_BACKWARD, FFTW_ESTIMATE);

  for(it=0; it<ntpow2; it++) { 
    num[it]=0.; 
    den[it]=0.;
  }
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++) tmp[it]=sim->dcal[irec][it]*acq->xweight[irec];
    for(it=sim->nt; it<ntpow2; it++) tmp[it]=0.;
    fftw_execute(fft);
    memcpy(ft_dcal, tmp, ntpow2*sizeof(fftw_complex));

    for(it=0; it<sim->nt; it++) tmp[it]=sim->dobs[irec][it]*acq->xweight[irec];
    for(it=sim->nt; it<ntpow2; it++) tmp[it]=0.;
    fftw_execute(fft);
    memcpy(ft_dobs, tmp, ntpow2*sizeof(fftw_complex));

    for(it=0; it<ntpow2; it++){
      num[it]+=conj(ft_dcal[it])*ft_dobs[it];
      den[it]+=creal(conj(ft_dcal[it])*ft_dcal[it]);
    }
  }

  float *num_real=alloc1float(ntpow2);
  float *num_imag=alloc1float(ntpow2);
  float *num_sum_real=alloc1float(ntpow2);
  float *num_sum_imag=alloc1float(ntpow2);
  float *den_sum=alloc1float(ntpow2);
  for(it=0; it<ntpow2; it++) {
    num_real[it]=creal(num[it]);
    num_imag[it]=cimag(num[it]);
  }
  if(sim->eachopt){
    memcpy(num_sum_real, num_real, ntpow2*sizeof(float));
    memcpy(num_sum_imag, num_imag, ntpow2*sizeof(float));
    memcpy(den_sum, den, ntpow2*sizeof(float));
  }else{
    MPI_Allreduce(num_real,num_sum_real,ntpow2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(num_imag,num_sum_imag,ntpow2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(den,den_sum,ntpow2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  }
  for(it=0; it<ntpow2; it++){
    num[it]=num_sum_real[it]+I*num_sum_imag[it];
    den[it]=den_sum[it];
  }

  free1float(num_real);
  free1float(num_imag);
  free1float(num_sum_real);
  free1float(num_sum_imag);
  free1float(den_sum);
  MPI_Barrier(MPI_COMM_WORLD);
  //num[] and den[] sums over all traces, then do the division
  float maxval = fabs(den[0]);
  for(it=0; it<ntpow2; it++){ 
    if(maxval<fabs(den[it])) maxval=fabs(den[it]);
  }
  float eps=1.e-4*maxval;
  for(it=0; it<ntpow2; it++) tmp[it]=num[it]/(den[it]+eps);//freq domain Wiener filter
  //Inverse Fourier transform back to time domain
  fftw_execute(ifft);//keep in mind there is missing factor 1/N after this step
  //truncate real part of src to be of length nt
  for(it=0; it<sim->nt; it++) sim->stf[it]=creal(tmp[it])/ntpow2;

  if(sim->eachopt){
    char number[sizeof("0000")];
    char fname[10];
    sprintf(number, "%04d", acq->shot_idx[iproc]);
    snprintf(fname, sizeof(fname), "%s_%s", stffile, number);
    
    fp=fopen(fname,"wb");
    fwrite(sim->stf, sim->nt*sizeof(float), 1, fp);
    fclose(fp);
  }else{
    if(iproc==0){
      fp=fopen(stffile,"wb");
      fwrite(sim->stf, sim->nt*sizeof(float), 1, fp);
      fclose(fp);
    }
  }


  free1float(den);
  fftw_free(num);
  fftw_free(tmp);
  fftw_free(ft_dobs);
  fftw_free(ft_dcal);
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);

}
