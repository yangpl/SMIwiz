#include <stdlib.h>
#include <stdio.h>

int main()
{
  int n1 = 81;
  int n2 = 181;
  int n3 = 181;
  float *rho, *vp, *vpinit;
  FILE *fp;

  rho = malloc(n1*n2*n3*sizeof(float));
  vp = malloc(n1*n2*n3*sizeof(float));
  vpinit = malloc(n1*n2*n3*sizeof(float));

  for(int i=0; i<n1*n2*n3; i++){
    rho[i] = 1000;
    vp[i] = 2100;
    vpinit[i] = 2000;
  }

  fp = fopen("rho", "wb");
  fwrite(rho, sizeof(float), n1*n2*n3, fp);
  fclose(fp);

  fp = fopen("vp", "wb");
  fwrite(vp, sizeof(float), n1*n2*n3, fp);
  fclose(fp);

  fp = fopen("vp_init", "wb");
  fwrite(vpinit, sizeof(float), n1*n2*n3, fp);
  fclose(fp);
  
  free(rho);
  free(vp);
  free(vpinit);

}
