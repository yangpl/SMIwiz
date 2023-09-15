/*
  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
  Homepage: https://yangpl.wordpress.com
  E-mail: ypl.2100@gmail.com
*/
#include <math.h>

//=============================================================
double sinc(double x)
{
  const float PI = 3.141592653589793238462643;
  static float eps=1e-15;
  float pix;

  if(fabs(x)>eps) {
    pix=PI*x;
    return sin(pix)/pix;
  } return 1;

}

//=============================================================
double bessi0(double x)
/*< Evaluate modified Bessel function In(x) and n=0.  >*/
{
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

double kaiser_windowed_sinc(double x, double dx, int r)
{
  double xr, w;  
  static float b1[10]={1.24,2.94,4.53,6.31,7.91,9.42,10.95,12.53,14.09,14.18};//monopole
  /* static float b2[10]={0.00,3.33,4.96,6.42,7.77,9.52,11.11,12.52,14.25,16.09};//dipole */
  //float b =b2[r-1];
  //w =(cos(PI*x)-sinc(x))*bessi0(b2_*sx)/(x*i0b2);
  float b = b1[r-1];//r=4
  
  w = 0.;
  xr = fabs(x/(dx*r));
  if(xr<1){
    if(xr>1e-7) w = sinc(x)*bessi0(b*sqrt(1.-xr*xr) )/bessi0(b);
    else w = 1.;
  }
  
  return w;
}
