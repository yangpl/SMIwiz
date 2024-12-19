/* model regularization by Tikhonov or Total variation (TV)
 *--------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *--------------------------------------------------------------------------*/
#include "cstd.h"

float regularization_tikhonov(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3)
{
  int i1, i2, i3, k0, kp1, km1;
  float s1, s2, s3, t1, t2, t3, tmp, fcost_mod;
  float _d1 = 1./(d1*d1);
  float _d2 = 1./(d2*d2);
  float _d3 = 1./(d3*d3);
  
  memset(g, 0, n1*n2*n3*sizeof(float));

  if(n3>1){
    s1 = s2 = s3 = 0.;
    for(i3=1; i3<n3-1; i3++){
      for(i2=1; i2<n2-1; i2++){
	for(i1=1; i1<n1-1; i1++){
	  k0 = i1 + n1*(i2 + n2*i3);

	  km1 = k0-1;
	  kp1 = k0+1;
	  t1 = -(x[km1] -2.0*x[k0]+x[kp1])*_d1;
	  tmp = x[kp1] - x[k0];
	  s1 += tmp*tmp;

	  km1 = k0-n1;
	  kp1 = k0+n1;
	  t2 = -(x[km1] -2.0*x[k0]+x[kp1])*_d2;
	  tmp = x[kp1] - x[k0];
	  s2 += tmp*tmp;

	  km1 = k0-n1*n2;
	  kp1 = k0+n1*n2;
	  t3 = -(x[km1] -2.0*x[k0]+x[kp1])*_d3;
	  tmp = x[kp1] - x[k0];
	  s3 += tmp*tmp;

	  g[k0] = t1+t2+t3;
	}
      }
    }
    s1 *= _d1;
    s2 *= _d2;
    s3 *= _d3;
  
    fcost_mod = s1 + s2 + s3;
  }else{
    i3 = 0;

    s1 = s2 = 0;
    for(i2=1; i2<n2-1; i2++){
      for(i1=1; i1<n1-1; i1++){
	k0 = i1 + n1*(i2 + n2*i3);

	km1 = k0-1;
	kp1 = k0+1;
	t1 = -(x[km1] -2.0*x[k0]+x[kp1])*_d1;
	tmp = x[kp1] - x[k0];
	s1 += tmp*tmp;

	km1 = k0-n1;
	kp1 = k0+n1;
	t2 = -(x[km1] -2.0*x[k0]+x[kp1])*_d2;
	g[k0] = t1 + t2;
	tmp = x[kp1] - x[k0];
	s2 += tmp*tmp;
      }
    }
    s1 *= _d1;
    s2 *= _d2;

    fcost_mod = s1 + s2;
  }
  
  return fcost_mod;
  
}


float regularization_tv(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3)
{
  int i1, i2, i3, k0, kp1, km1;
  float s1, s2, s3, tmp, fcost_mod;
  float beta = 1e-3/MAX(d1,MAX(d2,d3));
  float beta2 = beta*beta;

  float *g1, *g2, *g3;
  g1 = alloc1float(n1*n2*n3);
  g2 = alloc1float(n1*n2*n3);
  g3 = alloc1float(n1*n2*n3);

  memset(g, 0, n1*n2*n3*sizeof(float));
  memset(g1, 0, n1*n2*n3*sizeof(float));
  memset(g2, 0, n1*n2*n3*sizeof(float));
  memset(g3, 0, n1*n2*n3*sizeof(float));
  
  fcost_mod = 0;
  if(n3>1){
    for(i3=1; i3<n3-1; i3++){
      for(i2=1; i2<n2-1; i2++){
	for(i1=1; i1<n1-1; i1++){
	  k0 = i1 + n1*(i2 + n2*i3);

	  kp1 = k0+1;
	  s1 = (x[kp1] - x[k0])/d1;

	  kp1 = k0 + n1;
	  s2 = (x[kp1] - x[k0])/d2;

	  kp1 = k0 + n1*n2;
	  s3 = (x[kp1] - x[k0])/d3;

	  tmp = sqrtf(s1*s1 + s2*s2 + s3*s3 + beta2);
	  fcost_mod += tmp;

	  g1[k0] = s1/tmp;
	  g2[k0] = s2/tmp;
	  g3[k0] = s3/tmp;
	}
      }
    }
  
    for(i3=1; i3<n3-1; i3++){
      for(i2=1; i2<n2-1; i2++){
	for(i1=1; i1<n1-1; i1++){
	  k0 = i1 + n1*(i2 + n2*i3);

	  km1 = k0 - 1;
	  s1 = (g1[k0] - g1[km1])/d1;

	  km1 = k0 - n1;
	  s2 = (g2[k0] - g2[km1])/d2;

	  km1 = k0 - n1*n2;
	  s3 = (g3[k0] - g3[km1])/d3;
	  g[k0] = -(s1 + s2 + s3);
	}
      }
    }

  }else{
    i3 = 0;

    for(i2=1; i2<n2-1; i2++){
      for(i1=1; i1<n1-1; i1++){
	k0 = i1 + n1*i2 ;

	kp1 = k0+1;
	s1 = (x[kp1] - x[k0])/d1;

	kp1 = k0 + n1;
	s2 = (x[kp1] - x[k0])/d2;

	tmp = sqrtf(s1*s1 + s2*s2 + beta2);
	fcost_mod += tmp;
	g1[k0] = s1/tmp;
	g2[k0] = s2/tmp;
      }
    }

    for(i2=1; i2<n2-1; i2++){
      for(i1=1; i1<n1-1; i1++){
	k0 = i1 + n1*(i2 + n2*i3);

	km1 = k0 - 1;
	s1 = (g1[k0] - g1[km1])/d1;

	km1 = k0 - n1;
	s2 = (g2[k0] - g2[km1])/d2;

	g[k0] = -(s1 + s2 );
      }
    }

  }
  
  free1float(g1);
  free1float(g2);
  free1float(g3);
  
  return fcost_mod;
}
  
