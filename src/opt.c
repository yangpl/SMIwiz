/* interface to implement l-BFGS algorithm 

   The Wolfe conditions can be found in Nocedal, Numerical Optimization,
   2nd edition, p.33

   The linesearch method in this code is based first on a bracketing
   strategy, then on a dichotomy algorithm. See a full description in:
   Numerical Optimizationn Theoretical and  Practical Aspects,
   J.F.Bonnans, J.C.Gilbert,  C. Lemaréchal, C.A. Sagastizábal,
   Springer-Verlag

   Reference:
   [1] Numerical Optimization, Nocedal, 2nd edition, 2006
   Algorithm 7.4 p. 178, Algorithm 7.5 p. 179
   [2] https://en.wikipedia.org/wiki/Limited-memory_BFGS

   Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
   Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
   Homepage: https://yangpl.wordpress.com
   E-mail: ypl.2100@gmail.com
*/
#include "cstd.h"
#include "opt.h"

float l2norm(int n, float *a)
/*< L2 norm of a vector >*/
{
  int i;
  double sum=0.;

  for(i=0; i<n; i++)	sum += a[i]*a[i];
  sum=sqrtf(sum);
  
  return sum;
}

float dotprod(int n, float *a, float *b)
/*< dot product of two vectors >*/
{
  int i;
  double sum=0.;

  for (i=0; i<n; i++)	sum += a[i]*b[i];
  return sum;
}

void flipsign(int n, float *a, float *b)
/*< reverse the sign of the vector >*/
{
  int i;
  for (i=0; i<n; i++)	b[i]=-a[i];
}


void lbfgs_save(int n, float *x, float *g, float **sk, float **yk, opt_t *opt)
/*< save current model and gient >*/
{
  int i;
  if(opt->kpair < opt->npair){
    //save x and g if number of stored pairs does not exceeds maximum number of buffers
    memcpy(sk[opt->kpair], x, n*sizeof(float));//sk[opt->kpair]<--x
    memcpy(yk[opt->kpair], g, n*sizeof(float));//yk[opt->kpair]<--g
    opt->kpair += 1; // always we have: kpair<=npair
  }else{
    // otherwise erase the oldest pair and save the new one by shift
    for(i=0; i < opt->npair-1; i++){
      memcpy(sk[i], sk[i+1], n*sizeof(float)); //sk[i+1]<--sk[i]
      memcpy(yk[i], yk[i+1], n*sizeof(float)); //yk[i+1]<--yk[i]
    }
    memcpy(sk[opt->npair-1], x, n*sizeof(float));
    memcpy(yk[opt->npair-1], g, n*sizeof(float));
  }
}

void lbfgs_update(int n, float *x, float *g, float **sk, float **yk, opt_t *opt)
/*< update current sk and yk >*/
{
  int i,j;
  j=opt->kpair-1;// kpair has been stored, index ranges between 0 and kpair-1
  for(i=0; i<n; i++){
    sk[j][i]=x[i]-sk[j][i];
    yk[j][i]=g[i]-yk[j][i];
  }
}

void lbfgs_descent(int n, float *g, float *d, float **sk, float **yk, opt_t *opt)
/*< calculate search direction (two-loop recursion) >*/
{
  int i, j;
  float *rho, *q, *alpha, tmp0, tmp1, gamma, beta;

  // safeguard
  tmp0=l2norm(n, sk[opt->kpair-1]);
  tmp1=l2norm(n, yk[opt->kpair-1]);
  if(!( tmp0>0. && tmp1>0.)){
    flipsign(n, g, d); //descent direction= -gradient
    return;
  }

  q=alloc1float(n);
  rho=alloc1float(opt->kpair);
  alpha=alloc1float(opt->kpair);

  //store the gradient in vector q
  memcpy(q, g, n*sizeof(float));

  // first loop
  for(i=opt->kpair-1; i>=0; i--){
    // calculate rho
    tmp0=dotprod(n, yk[i], sk[i]);
    tmp1=dotprod(n, sk[i], q);
    rho[i]=1./tmp0;
    alpha[i]=rho[i]*tmp1;
    for(j=0; j<n; j++) q[j] -= alpha[i]*yk[i][j];
  }

  // tmp0=dotprod(n, yk[opt->kpair-1], sk[opt->kpair-1]); //same as 1./rho[opt->kpair-1]
  tmp0=1./rho[opt->kpair-1];
  tmp1=dotprod(n, yk[opt->kpair-1], yk[opt->kpair-1]);
  gamma=tmp0/tmp1;     // initial Hessian=gamma* I
  for (j=0; j<n; j++)	d[j]=gamma*q[j];

  // second loop
  for(i=0; i<opt->kpair; i++){
    tmp0=dotprod(n, yk[i], d);
    beta=rho[i]*tmp0;
    tmp1=alpha[i]-beta;
    for(j=0; j<n; j++) d[j] += tmp1*sk[i][j];
  }

  // descent direction = - H^(-1)*g
  for(j=0; j<n; j++)	d[j]=-d[j];

  // deallocate variables
  free1float(q);
  free1float(alpha);
  free1float(rho);
}

int lbfgs_descent1(int n, float *g, float *q, float *rho, float *alpha,
		    float **sk, float **yk, opt_t *opt)
/*< calculate search direction (1st loop of the two-loop recursion) >*/
{
  int loop1 = 0;
  int i, j;
  float tmp0, tmp1;

  // store gradient in vector q
  memcpy(q, g, n*sizeof(float));

  // safeguard a descent direction from negative gradient
  // d=-q=-g
  tmp0=l2norm(n, opt->sk[opt->kpair-1]);
  tmp1=l2norm(n, opt->yk[opt->kpair-1]);
  if(!( tmp0>0. && tmp1>0.))	return loop1;

  for(i=opt->kpair-1; i>=0; i--){
    // calculate rho
    tmp0=dotprod(n, yk[i], sk[i]);
    tmp1=dotprod(n, sk[i], q);
    rho[i]=1./tmp0;
    alpha[i]=rho[i]*tmp1;
    for(j=0; j<n; j++) q[j] -= alpha[i]*yk[i][j];
  }
  loop1 = 1; //done the 1st loop
  return loop1;
}

void lbfgs_descent2(int n, float *g, float *q, float *rho, float *alpha,
		    float **sk, float **yk, opt_t *opt)
/*< calculate search direction (2nd loop of the two-loop recursion) >*/
{
  int i, j;
  float tmp0, tmp1, gamma, beta;

  // tmp0=dotprod(n, yk[opt->kpair-1], sk[opt->kpair-1]); //same as 1./rho[opt->kpair-1]
  tmp0=1./rho[opt->kpair-1];
  tmp1=dotprod(n, yk[opt->kpair-1], yk[opt->kpair-1]);
  gamma=tmp0/tmp1; // initial Hessian=gamma* I
  for (j=0; j<n; j++)	q[j]=gamma*q[j];

  // second loop
  for(i=0; i<opt->kpair; i++){
    tmp0=dotprod(n, yk[i], q);
    beta=rho[i]*tmp0;
    tmp1=alpha[i]-beta;
    for(j=0; j<n; j++) q[j] += tmp1*sk[i][j];
  }
}


void boundx(float *x, int n, float *xmin, float *xmax)
/*< clip x based on lower and upper bounds >*/
{
  int i;

  for(i=0; i<n; i++){
    if(x[i]<xmin[i]) x[i] = xmin[i];
    if(x[i]>xmax[i]) x[i] = xmax[i];
  }
}

void line_search(int n, //dimension of x
		 float *x, //input vector x
		 float *g, //gradient of misfit function
		 float *d, //descent direction
		 opt_fg fg, //subroutine to evaluation function and gradient
		 opt_t *opt) //pointer of l-BFGS optimization parameters
/*< bisection line search  based on Wolfe condition >*/
{
  int j;
  float gxd, c1_gxd, c2_gxd, fcost, fxx, alpha1, alpha2;
  float *xk;
  float inf = 0x3f3f3f3f;//this is actual infinity in computer
  
  //use estimated stepsize from previous iteration when iter>1
  opt->alpha = 1.;
  alpha1 = 0;
  alpha2 = inf;
  
  xk=alloc1float(n);  // allocate memory for store current x
  memcpy(xk, x, n*sizeof(float)); // store x at k-th iteration
  //m3=the slope of the function of alpha along search d
  gxd = dotprod(n, g, d);//<G[f(x)]|d>
  c1_gxd = opt->c1*gxd;//c1*<G[f(x)]|d>
  c2_gxd = opt->c2*gxd;//c2*<G[f(x)]|d>
  for(opt->ils=0; opt->ils<opt->nls; opt->ils++){
    for(j=0; j<n; j++) x[j] = xk[j] + opt->alpha*d[j];//update x

    //clip x by lower+upper bounds (l-BFGS-B, bounded l-BFGS)
    if(opt->bound==1) boundx(x, n, opt->xmin, opt->xmax);
    fcost = fg(x, g);// function and gradient evaluation
    opt->igrad++; //update the counter for function and gradient evaluation in total

    //m3=the slope of the function of alpha along search d
    gxd = dotprod(n, g, d);//<G[f(x+alp*d)]|d>
    fxx = opt->fk + opt->alpha * c1_gxd;

    //check Wolfe condition for current step length, see Nocedal book, eqn 3.6
    //eqn 3.6a: f(x + alp*d) <= f(x) + m1*alpha*<G[f(x)]|d>
    //eqn 3.6b: <G[f(x+alp*d)]|d>  >= m2*<G[f(x)]|d>
    if(fcost > fxx){
      if(opt->verb) printf("Wolfe condition 1 fails: insufficient misfit decrease!\n");
      alpha2 = opt->alpha;
      opt->alpha = 0.5*(alpha1+alpha2);//shrink the search interval
    }else if(gxd < c2_gxd){
      if(opt->verb) printf("Wolfe condition 2 fails: stepsize is too small!\n");
      alpha1 = opt->alpha;
      if(alpha2<inf) opt->alpha = 0.5*(alpha1 + alpha2);//shrink the search interval
      else                opt->alpha *= 5.; //extend the search interval if alpha2==0
    }else{//conditions satisfied, terminate line search
      break;
    }

    if(opt->verb){
      printf("#line search %d, alp1=%f alp2=%f alp=%f\n", opt->ils, alpha1, alpha2, opt->alpha);
      printf("------------------------------------------------\n");
    }
  }

  if(fcost <= opt->fk) {
    opt->ls_fail = 0;
    opt->fk = fcost;//fcost was not increased with this stepsize, accept it
  }else{
    opt->ls_fail = 1; //line search fails, exit
  }
  
  free1float(xk);
}


void cg_solve(int n, //dimension of x
	      float *x, //input vector x
	      float *g, //gradient of misfit function
	      float *d, //descent direction
	      opt_Hv Hv, //subroutine to evaluation function and gradient
	      opt_t *opt) //pointer of l-BFGS optimization parameters
{
  int i,k;
  float rsold, rsnew, rs0, pAp, alp, beta;
  float tol = 1e-3;

  float *r = alloc1float(n);
  float *p = alloc1float(n);
  float *Ap = alloc1float(n);
  //assume x0=d=0; then r0=b-Ax0=b=-g
  memset(d, 0, n*sizeof(float));
  rsold = 0;
  for(i=0; i<n; i++) {
    r[i] = -g[i];  
    p[i] = r[i];
    rsold += r[i]*r[i];
  }
  rs0 = rsold;
  for(k=0; k<opt->ncg; k++){
    Hv(x, p, Ap);//compute Ap

    pAp = 0;
    for(i=0; i<n; i++) pAp += p[i]*Ap[i];
    alp = rsold/pAp;

    rsnew = 0;
    for(i=0; i<n; i++){
      d[i] += alp*p[i];
      r[i] -= alp*Ap[i];
      rsnew += r[i]*r[i];
    }

    if(rsnew<tol*rs0) break;

    beta = rsnew/rsold;
    for(i=0; i<n; i++) p[i] = r[i] + beta*p[i];
    rsold = rsnew;
  }

  free1float(r);
  free1float(p);
  free1float(Ap);
}
