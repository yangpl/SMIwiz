/* Adjoint source using AWI (adaptive waveform inversion)
 *--------------------------------------------------------------------------
 * Reference:
 *   Warner, Michael, and Lluís Guasch. "Adaptive waveform inversion: Theory."
 *   Geophysics 81.6 (2016): R429-R445.
 *---------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *--------------------------------------------------------------------------*/
#include "cstd.h"
#include "acq.h"
#include "sim.h"
#include "fwi.h"
#include <fftw3.h>

//there might be a missing minus sign in the AWI adjoint source, not tested yet!
float awi_adjoint_source(acq_t *acqui, sim_t *sim, fwi_t *fwi)
{
  float fcost;
  float s1, s2, tmp, maxval, eps;
  int irec,it;
  
  int ntpow2=1;
  while(ntpow2<2*sim->nt) { ntpow2 *= 2; } // then ntpow2 is the 2^power, >=2*nt

  float **v= alloc2float(sim->nt, acqui->nrec);
  float *den = alloc1float(ntpow2);
  fftw_complex *num = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_complex *ft_dcal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2*acqui->nrec);
  fftw_complex *ft_dobs = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2*acqui->nrec);
  fftw_complex *ffttmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_plan fft = fftw_plan_dft_1d(ntpow2, ffttmp, ffttmp, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan ifft = fftw_plan_dft_1d(ntpow2, ffttmp, ffttmp, FFTW_BACKWARD, FFTW_ESTIMATE);

  //----------------------------------------------------------------------
  s1 = 0;
  s2 = 0;
  for(irec=0; irec<acqui->nrec; irec++){
    memset(ffttmp, 0, ntpow2*sizeof(fftw_complex));
    for(it=0; it<sim->nt; it++) ffttmp[it]=sim->dcal[irec][it]*acqui->xweight[irec];
    fftw_execute(fft);
    memcpy(&ft_dcal[ntpow2*irec], ffttmp, ntpow2*sizeof(fftw_complex));

    memset(ffttmp, 0, ntpow2*sizeof(fftw_complex));
    for(it=0; it<sim->nt; it++) ffttmp[it]=sim->dobs[irec][it]*acqui->xweight[irec];
    fftw_execute(fft);
    memcpy(&ft_dobs[ntpow2*irec], ffttmp, ntpow2*sizeof(fftw_complex));

    maxval=fabs(den[0]);
    for(it=0; it<ntpow2; it++){
      num[it] = conj(ft_dobs[it+ntpow2*irec])*ft_dcal[it+ntpow2*irec];
      den[it] = conj(ft_dobs[it+ntpow2*irec])*ft_dobs[it+ntpow2*irec];
      if(maxval<fabs(den[it])) maxval = fabs(den[it]);
    }
    eps=1.e-5*maxval;
    for(it=0; it<ntpow2; it++) ffttmp[it]=num[it]/(den[it]+eps);//freq domain Wiener filter
    //Inverse Fourier transform back to time domain
    fftw_execute(ifft);//keep in mind there is missing factor 1/N after this step
    //truncate real part of src to be of length nt
    for(it=0; it<sim->nt; it++) {
      v[irec][it] = creal(ffttmp[it])/ntpow2;
      tmp = fabs(it)*sim->dt*v[irec][it];
      s1 += tmp*tmp;
      s2 += v[irec][it]*v[irec][it];
    }
  }
  fcost = 0.5*s1/s2;
  for(irec=0; irec<acqui->nrec; irec++){
    for(it=0; it<sim->nt; it++) {
      tmp = fabs(it)*sim->nt;
      ffttmp[it] = -(tmp*tmp - 2.*fcost)/s2*v[irec][it];
    }

    for(it=sim->nt; it<ntpow2; it++) ffttmp[it] = 0.;
    fftw_execute(fft);
    memcpy(&ft_dcal[ntpow2*irec], ffttmp, ntpow2*sizeof(fftw_complex));
    maxval = 0;
    for(it=0; it<ntpow2; it++){
      num[it] = ft_dcal[it+ntpow2*irec]*ft_dobs[it+ntpow2*irec];
      den[it] = conj(ft_dobs[it+ntpow2*irec])*ft_dobs[it+ntpow2*irec];
      if(maxval<fabs(den[it])) maxval = fabs(den[it]);
    }
    eps=1.e-5*maxval;
    for(it=0; it<ntpow2; it++) ffttmp[it]=num[it]/(den[it]+eps);//freq domain Wiener filter
    fftw_execute(ifft);//keep in mind there is missing factor 1/N after this step
    for(it=0; it<sim->nt; it++) sim->dres[irec][it] = creal(ffttmp[it])/ntpow2;
  }

  free2float(v);  
  free1float(den);
  fftw_free(num);
  fftw_free(ffttmp);
  fftw_free(ft_dobs);
  fftw_free(ft_dcal);
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);

  return fcost;
}

