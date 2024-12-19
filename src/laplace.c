/* Laplace filtering, output= -laplace(input)
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *---------------------------------------------------------------------*/
 #include "cstd.h"
#include <fftw3.h>

fftw_complex *tmpfft;
fftw_plan fft_plan, ifft_plan;

int n1fft, n2fft, n3fft;
int n1, n2, n3, npar;
float *k1, *k2, *k3;

//return an integer m=2^k, where m>=n
int fft_next_fast_size(int n)
{
  int m;
  
  m = 1;
  while(m<n) m *= 2;
  return m;
}

void laplace_init(int n1_, int n2_, int n3_, int npar_, float d1_, float d2_, float d3_)
{
  int i1, i2, i3;
  float d1, d2, d3;
  float dk1, dk2, dk3;
  static float PI = 3.14159265359;

  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  npar = npar_;
  d1 = d1_;
  d2 = d2_;
  d3 = d3_;

  n1fft = fft_next_fast_size(n1);
  n2fft = fft_next_fast_size(n2);
  n3fft = (n3>1)?fft_next_fast_size(n3):1;
  
  k1 = alloc1float(n1fft);
  k2 = alloc1float(n2fft);

  /* pre-compute the discrete wavenumber - k1 */
  dk1 = 2.0*PI/(d1*n1fft);
  k1[0] = 0;
  for(i1 = 1; i1<(n1fft+1)/2; i1++) {
    k1[i1] = i1*dk1;
    k1[n1fft-i1] = -i1*dk1;
  }
  if(n1fft%2==0) k1[n1fft/2] = (n1fft/2)*dk1;/* Nyquist freq*/
  /* pre-compute the discrete wavenumber - k2 */
  dk2 = 2.0*PI/(d2*n2fft);
  k2[0] = 0;
  for(i2 = 1; i2<(n2fft+1)/2; i2++) {
    k2[i2] = i2*dk2;
    k2[n2fft-i2] = -i2*dk2;
  }
  if(n2fft%2==0) k2[n2fft/2] = (n2fft/2)*dk2;/* Nyquist freq*/

  if(n3>1){/* pre-compute the discrete wavenumber - k3 */
    k3 = alloc1float(n3fft);
    dk3 = 2.0*PI/(d3*n3fft);
    for(i3 = 1; i3<(n3fft+1)/2; i3++) {
      k3[i3] = i3*dk3;
      k3[n3fft-i3] = -i3*dk3;
    }
    if(n3fft%2==0) k3[n3fft/2] = (n3fft/2)*dk3;/* Nyquist freq*/
  }

  tmpfft = fftw_malloc(sizeof(fftw_complex)*n1fft*n2fft*n3fft);
  fft_plan = fftw_plan_dft_3d(n1fft, n2fft, n3fft, tmpfft, tmpfft, FFTW_FORWARD, FFTW_ESTIMATE);
  ifft_plan = fftw_plan_dft_3d(n1fft, n2fft, n3fft, tmpfft, tmpfft, FFTW_BACKWARD, FFTW_ESTIMATE);  
}

void laplace_close()
{
  fftw_free(tmpfft);
  fftw_destroy_plan(fft_plan);
  fftw_destroy_plan(ifft_plan);  

  free1float(k1);
  free1float(k2);
  if(n3>1) free1float(k3);
}

void laplace_apply(double *in, double *out)
{
  int i1, i2, i3, ipar;
  int i, j, n123fft;
  float kk;

  n123fft = n1fft*n2fft*n3fft;
  for(ipar=0; ipar<npar; ipar++){
    //------------------------------------------
    for(i3=0; i3<n3fft; i3++){
      for(i2=0; i2<n2fft; i2++){
	for(i1=0; i1<n1fft; i1++){
	  i = i1 + n1fft*(i2 + n2fft*i3);
	  j = i1 + n1*(i2 + n2*(i3 + n3*ipar));
	  tmpfft[i] = (i1<n1 && i2<n2 && i3<n3)?in[j]:0;	
	}
      }
    }
    fftw_execute(fft_plan);
    for(i3=0; i3<n3fft; i3++){
      for(i2=0; i2<n2fft; i2++){
	for(i1=0; i1<n1fft; i1++){
	  if(n3>1) kk = k1[i1]*k1[i1] + k2[i2]*k2[i2] + k3[i3]*k3[i3];
	  else 	 kk = k1[i1]*k1[i1] + k2[i2]*k2[i2];
	
	  i = i1 + n1fft*(i2 + n2fft*i3);
	  tmpfft[i] *= kk;
	}
      }
    }
    fftw_execute(ifft_plan);
    for(i3=0; i3<n3; i3++){
      for(i2=0; i2<n2; i2++){
	for(i1=0; i1<n1; i1++){
	  i = i1 + n1fft*(i2 + n2fft*i3);
	  j = i1 + n1*(i2 + n2*(i3 + n3*ipar));
	  out[j] = -creal(tmpfft[i]/n123fft);
	}
      }
    }
  }//end for ipar
}

/*
int main()
{
  int i;
  int n1_ = 151;
  int n2_ = 461;
  int n3_ = 1;
  int npar_ = 1;
  float d1_ = 20;
  float d2_ = 20;
  float d3_ = 20;
  int n = n1_*n2_*n3_*npar_;

  float *in;
  double *f;
  FILE *fp;

  in = alloc1float(n);
  f = alloc1double(n);
  
  fp = fopen("param_final_rtm", "rb");
  fread(in, n*sizeof(float), 1, fp);
  fclose(fp);
  printf("read completed\n");
  
  for(i=0; i<n; i++) f[i] = in[i];
  laplace_init(n1_, n2_, n3_, npar_, d1_, d2_, d3_);

  laplace_apply(f, f);

  laplace_close();


  for(i=0; i<n; i++) in[i] = f[i];
  fp = fopen("param_final_rtm_filtered", "wb");
  fwrite(in, n*sizeof(float), 1, fp);
  fclose(fp);
  printf("write completed\n");

  free1float(in);
  free1double(f);
}

*/
