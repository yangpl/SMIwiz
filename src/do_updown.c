/* Up-down wavefield separation using Hilbert transform
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
 
#include <mpi.h>
#include <fftw3.h>

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

void do_updown(sim_t *sim, acq_t *acq)
{
  int i1, i2, i3, i1_, i2_, i3_, it, j;
  int ntaper;
  FILE *fp;

  if(!getparint("itcheck", &sim->itcheck)) sim->itcheck = sim->nt/2;
  if(!getparint("ntaper",&ntaper)) ntaper = 10;

  //---------------------------------------------------------------
  float *hstf = alloc1float(sim->nt);//hilbert transform of source time function

  int ntpow2 = 1;
  while(ntpow2<=sim->nt) { ntpow2 *= 2; } //ntpow2=2^power, ntpow2>=nt

  fftw_complex *tmp_time = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_plan fft_time = fftw_plan_dft_1d(ntpow2, tmp_time, tmp_time, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan ifft_time = fftw_plan_dft_1d(ntpow2, tmp_time, tmp_time, FFTW_BACKWARD, FFTW_ESTIMATE);
  memset(tmp_time, 0, ntpow2*sizeof(fftw_complex));
  for(it=0; it<sim->nt; it++) tmp_time[it] = sim->stf[it];
  fftw_execute(fft_time);
  for(it=0; it<ntpow2; it++){
    if(it==0) tmp_time[it] = tmp_time[it];
    else if(it<=ntpow2/2) tmp_time[it] *= 2;
    else tmp_time[it] = 0;
  }
  fftw_execute(ifft_time);
  for(it=0; it<sim->nt; it++) hstf[it] = cimag(tmp_time[it])/ntpow2;

  //-----------------------------------------------------------
  int nzpow2 = 1;
  while(nzpow2<=sim->n1pad) { nzpow2 *= 2; } //nzpow2=2^power, nzpow2>=nz
  int n[] = {nzpow2};
  int howmany = sim->n2pad*sim->n3pad;
  int len = nzpow2*howmany;
  fftw_complex *tmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*len);
  fftw_plan fft_zaxis = fftw_plan_many_dft(1, n, howmany,
					   tmp, n,  1, nzpow2,
					   tmp, n,  1, nzpow2,
					   FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan ifft_zaxis = fftw_plan_many_dft(1, n, howmany,
					    tmp, n,  1, nzpow2,
					    tmp, n,  1, nzpow2,
					    FFTW_BACKWARD, FFTW_ESTIMATE);
  
  float ***pu = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);

  sim->sign_dt = 1;
  check_cfl(sim);
  cpml_init(sim);
  extend_model_init(sim);
  extend_model(sim, sim->vp, sim->rho, sim->kappa, sim->buz, sim->bux, sim->buy);
  computing_box_init(acq, sim, 0);
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_null(sim, 1);//flag=1, incident field
  fdtd_init(sim, 0);//flag=0, incident field
  fdtd_null(sim, 0);//flag=0, incident field

  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);//flag=1
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);//flag=1
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);

    fdtd_update_v(sim, 0, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);//flag=1
    fdtd_update_p(sim, 0, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);//flag=1
    inject_source(sim, acq, sim->p0, hstf[it]);

    memset(tmp, 0, len*sizeof(fftw_complex));
    for(i3=0; i3<sim->n3pad; i3++){
      for(i2=0; i2<sim->n2pad; i2++){
	for(i1=0; i1<sim->n1pad; i1++){
	  j = i1 + nzpow2*(i2 + sim->n2pad*i3);
	  tmp[j] = sim->p1[i3][i2][i1] + I*sim->p0[i3][i2][i1];//analytic wavefield
	}
      }
    }
    fftw_execute(fft_zaxis);
    for(i3=0; i3<sim->n3pad; i3++){
      for(i2=0; i2<sim->n2pad; i2++){
	for(i1=0; i1<nzpow2; i1++){
	  j = i1 + nzpow2*(i2 + sim->n2pad*i3);

	  if(i1>nzpow2/2) tmp[j] = 0.;//set negative kz to 0, keep only positive kz
	  else{
	    if(i1<ntaper) tmp[j] *= sin(0.5*PI*i1/ntaper);//apply a taper
	    else if(i1>nzpow2/2-ntaper) tmp[j] *= sin(0.5*PI*(nzpow2/2-i1)/ntaper);//apply a taper
	  }
	}
      }
    }
    fftw_execute(ifft_zaxis);
    for(i3=0; i3<sim->n3pad; i3++){
      for(i2=0; i2<sim->n2pad; i2++){
	for(i1=0; i1<sim->n1pad; i1++){
	  j = i1 + nzpow2*(i2 + sim->n2pad*i3);
	  pu[i3][i2][i1] = creal(tmp[j])/nzpow2;
	}
      }
    }
    
    
    if(iproc==0 && it==sim->itcheck){
      fp = fopen("wave_up.bin", "wb");
      for(i3=0; i3<sim->n3; i3++){
      	i3_ = (sim->n3>1)?i3 + sim->nb:0;
      	for(i2=0; i2<sim->n2; i2++){
      	  i2_ = i2 + sim->nb;
      	  for(i1=0; i1<sim->n1; i1++){
      	    i1_ = i1 + sim->nb;

	    float val = pu[i3_][i2_][i1_]; 
      	    fwrite(&val, sizeof(float), 1, fp);
      	  }
      	}
      }
      fclose(fp);

      fp = fopen("wave_down.bin", "wb");
      for(i3=0; i3<sim->n3; i3++){
      	i3_ = (sim->n3>1)?i3 + sim->nb:0;
      	for(i2=0; i2<sim->n2; i2++){
      	  i2_ = i2 + sim->nb;
      	  for(i1=0; i1<sim->n1; i1++){
      	    i1_ = i1 + sim->nb;

	    float val = sim->p1[i3_][i2_][i1_] - pu[i3_][i2_][i1_]; 
      	    fwrite(&val, sizeof(float), 1, fp);
      	  }
      	}
      }
      fclose(fp);
    }//end if
  }
  write_data(sim, acq);

  extend_model_free(sim);
  fdtd_free(sim, 1);//flag=1, incident field
  fdtd_free(sim, 0);//flag=0, incident field
  cpml_free(sim);
  computing_box_free(sim, 0);

  free1float(hstf);
  free3float(pu);

}
