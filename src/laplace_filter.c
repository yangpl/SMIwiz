#include "cstd.h"
#include "sim.h"

void laplace_filter(sim_t *sim, float ***in, float ***out)
{
  int i1, i2 ,i3;
  float diff1, diff2, diff3;
  int i1min = 1;
  int i1max = sim->n1-2;
  int i2min = 1;
  int i2max = sim->n2-2;
  int i3min = (sim->n3>1)?1:0;
  int i3max = (sim->n3>1)?sim->n3-2:0;

  memset(&out[0][0][0], 0, sim->n123*sizeof(float));
  for(i3=i3min; i3<=i3max; i3++){
    for(i2=i2min; i2<=i2max; i2++){
      for(i1=i1min; i1<=i1max; i1++){
	diff1 = (in[i3][i2][i1+1] - 2.*in[i3][i2][i1] + in[i3][i2][i1-1])/(sim->d1*sim->d1);
	diff2 = (in[i3][i2+1][i1] - 2.*in[i3][i2][i1] + in[i3][i2-1][i1])/(sim->d2*sim->d2);
	diff3 = (sim->n3>1)?(in[i3+1][i2][i1] - 2.*in[i3][i2][i1] + in[i3-1][i2][i1])/(sim->d3*sim->d3):0.;
	out[i3][i2][i1] = diff1 + diff2 + diff3;
      }
    }
  }


}
