/*---------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *-------------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "mpi_info.h"

void read_data(sim_t *sim, acq_t *acq)
/*< read observed data according to shot/process index >*/
{
  char fname[sizeof("dat_0000")];
  FILE *fp;

  sprintf(fname, "dat_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"rb");
  if(fp==NULL) { fprintf(stderr,"read_data, error opening file\n"); exit(1);}
  fread(&sim->dobs[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
}

void write_data(sim_t *sim, acq_t *acq)
/*< write synthetic data for each shot/process >*/
{
  char fname[sizeof "dat_0000"];
  FILE *fp;

  sprintf(fname, "dat_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"write_data, error opening file\n"); exit(1);}
  fwrite(&sim->dcal[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
}
