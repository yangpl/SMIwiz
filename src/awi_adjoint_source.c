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

static int awi_lag(int iw, int nw)
{
  return (iw <= nw/2) ? iw : nw - iw;
}

float awi_adjoint_source(acq_t *acqui, sim_t *sim, fwi_t *fwi)
{
  (void)fwi;
  float fcost = 0.;
  double s1, s2, tmp, maxval, eps;
  double lag_time, trace_weight;
  int irec,it,lag;
  
  int ntpow2=1;
  while(ntpow2<2*sim->nt) { ntpow2 *= 2; } // then ntpow2 is the 2^power, >=2*nt

  float **v= alloc2float(ntpow2, acqui->nrec);
  float *ftrace = alloc1float(acqui->nrec);
  float *s2trace = alloc1float(acqui->nrec);
  float *epstrace = alloc1float(acqui->nrec);
  double *den = alloc1double(ntpow2);
  fftw_complex *num = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_complex *ft_dcal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2*acqui->nrec);
  fftw_complex *ft_dobs = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2*acqui->nrec);
  fftw_complex *ffttmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntpow2);
  fftw_plan fft = fftw_plan_dft_1d(ntpow2, ffttmp, ffttmp, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan ifft = fftw_plan_dft_1d(ntpow2, ffttmp, ffttmp, FFTW_BACKWARD, FFTW_ESTIMATE);

  //----------------------------------------------------------------------
  for(irec=0; irec<acqui->nrec; irec++){
    trace_weight = acqui->xweight[irec]*acqui->xweight[irec];
    ftrace[irec] = 0.;
    s2trace[irec] = 0.;
    epstrace[irec] = 0.;

    memset(ffttmp, 0, ntpow2*sizeof(fftw_complex));
    for(it=0; it<sim->nt; it++) ffttmp[it]=sim->dcal[irec][it]*acqui->wdat[irec][it];
    fftw_execute(fft);
    memcpy(&ft_dcal[ntpow2*irec], ffttmp, ntpow2*sizeof(fftw_complex));

    memset(ffttmp, 0, ntpow2*sizeof(fftw_complex));
    for(it=0; it<sim->nt; it++) ffttmp[it]=sim->dobs[irec][it]*acqui->wdat[irec][it];
    fftw_execute(fft);
    memcpy(&ft_dobs[ntpow2*irec], ffttmp, ntpow2*sizeof(fftw_complex));

    maxval = 0.;
    for(it=0; it<ntpow2; it++){
      num[it] = conj(ft_dobs[it+ntpow2*irec])*ft_dcal[it+ntpow2*irec];
      den[it] = creal(conj(ft_dobs[it+ntpow2*irec])*ft_dobs[it+ntpow2*irec]);
      if(maxval<den[it]) maxval = den[it];
    }
    if(!(maxval > 0.) || !(trace_weight > 0.)){
      for(it=0; it<sim->nt; it++) sim->dres[irec][it] = 0.;
      continue;
    }
    eps=1.e-5*maxval;
    epstrace[irec] = eps;
    for(it=0; it<ntpow2; it++) ffttmp[it]=num[it]/(den[it]+eps);//freq domain Wiener filter
    //Inverse Fourier transform back to time domain
    fftw_execute(ifft);//keep in mind there is missing factor 1/N after this step
    s1 = 0.;
    s2 = 0.;
    for(it=0; it<ntpow2; it++) {
      lag = awi_lag(it, ntpow2);
      v[irec][it] = creal(ffttmp[it])/ntpow2;
      if(lag < sim->nt){
	lag_time = lag*sim->dt;
	tmp = lag_time*v[irec][it];
	s1 += tmp*tmp;
	s2 += v[irec][it]*v[irec][it];
      }else{
	v[irec][it] = 0.;
      }
    }
    if(!(s2 > 0.)){
      for(it=0; it<sim->nt; it++) sim->dres[irec][it] = 0.;
      continue;
    }
    ftrace[irec] = 0.5*s1/s2;
    s2trace[irec] = s2;
    fcost += trace_weight*ftrace[irec];
  }

  for(irec=0; irec<acqui->nrec; irec++){
    trace_weight = acqui->xweight[irec]*acqui->xweight[irec];
    if(!(s2trace[irec] > 0.) || !(epstrace[irec] > 0.) || !(trace_weight > 0.)){
      for(it=0; it<sim->nt; it++) sim->dres[irec][it] = 0.;
      continue;
    }
    for(it=0; it<ntpow2; it++) {
      lag = awi_lag(it, ntpow2);
      if(lag < sim->nt){
	lag_time = lag*sim->dt;
	/* inject_adjoint_source follows the dobs-dcal convention used by the
	 * L2 objective, so this is the negative of paper equation 14. */
	ffttmp[it] = (2.*ftrace[irec] - lag_time*lag_time)/s2trace[irec]*v[irec][it];
      }else{
	ffttmp[it] = 0.;
      }
    }
    fftw_execute(fft);
    memcpy(&ft_dcal[ntpow2*irec], ffttmp, ntpow2*sizeof(fftw_complex));
    for(it=0; it<ntpow2; it++){
      num[it] = ft_dcal[it+ntpow2*irec]*ft_dobs[it+ntpow2*irec];
      den[it] = creal(conj(ft_dobs[it+ntpow2*irec])*ft_dobs[it+ntpow2*irec]);
    }
    eps=epstrace[irec];
    for(it=0; it<ntpow2; it++) ffttmp[it]=num[it]/(den[it]+eps);//freq domain Wiener filter
    fftw_execute(ifft);//keep in mind there is missing factor 1/N after this step
    for(it=0; it<sim->nt; it++){
      sim->dres[irec][it] = trace_weight*acqui->wdat[irec][it]*
	creal(ffttmp[it])/ntpow2;
    }
  }

  free2float(v);  
  free1float(ftrace);
  free1float(s2trace);
  free1float(epstrace);
  free1double(den);
  fftw_free(num);
  fftw_free(ffttmp);
  fftw_free(ft_dobs);
  fftw_free(ft_dcal);
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);

  return fcost;
}
