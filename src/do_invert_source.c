/*  Source wavelet estimation according to Pratt 1999 Geophysics paper
 *-----------------------------------------------------------------------
 *
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
#include "mpi_info.h"

void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_close(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim);
void extend_model_close(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_close(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_close(sim_t *sim);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);

void read_data(sim_t *sim, acq_t *acq);
void write_data(sim_t *sim, acq_t *acq);
void setup_data_weight(acq_t *acq, sim_t *sim);


void do_invert_source(sim_t *sim, acq_t *acq)
/*< invert source time function in FWI >*/
{
  FILE *fp;
  int it,irec,ntpow2;
  char *fname;

  sim->dcal = alloc2float(sim->nt, acq->nrec);
  sim->dobs = alloc2float(sim->nt, acq->nrec);
  read_data(sim, acq);
  setup_data_weight(acq, sim);
  
  //-----------------------------------------------------------------------
  //1. simulate synthetic data (Green's function) using Dirac delta function
  //-----------------------------------------------------------------------
  memset(sim->stf, 0, sim->nt*sizeof(float)); sim->stf[0] = 1;
  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  extend_model(sim);
  computing_box_init(acq, sim, 0);
  fdtd_init(sim, 1);
  fdtd_null(sim, 1);
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    fdtd_update_v(sim, 1, it, 0);
    fdtd_update_p(sim, 1, it, 0);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);
  }
  extend_model_close(sim);
  fdtd_close(sim, 1);
  cpml_close(sim);
  computing_box_close(sim, 0);


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
  MPI_Allreduce(num_real,num_sum_real,ntpow2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(num_imag,num_sum_imag,ntpow2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(den,den_sum,ntpow2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
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
  float maxval=fabs(den[0]);
  for(it=0; it<ntpow2; it++){ 
    if(maxval<fabs(den[it])) maxval=fabs(den[it]);
  }
  float eps=1.e-4*maxval;
  for(it=0; it<ntpow2; it++) tmp[it]=num[it]/(den[it]+eps);//freq domain Wiener filter
  //Inverse Fourier transform back to time domain
  fftw_execute(ifft);//keep in mind there is missing factor 1/N after this step
  //truncate real part of src to be of length nt
  for(it=0; it<sim->nt; it++) sim->stf[it]=creal(tmp[it])/ntpow2;

  fname="src_inverted";
  if(iproc==0){
    fp=fopen(fname,"wb");
    fwrite(sim->stf, sim->nt*sizeof(float), 1, fp);
    fclose(fp);
  }


  free1float(den);
  fftw_free(num);
  fftw_free(tmp);
  fftw_free(ft_dobs);
  fftw_free(ft_dcal);
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);

  free2float(sim->dcal);
}
