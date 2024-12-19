#ifndef _opt_h_
#define _opt_h_

typedef struct {
  int niter;//total number of iterations
  float tol;//convergence tolerance
  
  int npair;//maximum number pair allowed, l-BFGS memory length 
  int kpair;//actual number of pair already stored
  int nls;//maximum number of line search
  int ils;//actual number of line search
  int igrad;//number of function and gradient evaluation
  
  float c1;//Wolfe condition, Nocedal value=1e-4
  float c2;//Wolfe condition, Nocedal value=0.9
  float alpha;//true step length, initial value=1
  
  float f0;//initial misfit
  float fk;//misfit at k-th iteration
  float gk_norm;
  
  float *x;//unknown parameter
  float *g;//gradient
  float *pg;//preconditioned gradient
  float *d;//descent direction
  float **sk, **yk;//storing vectors in two-loop recursion
  float *q, *alp, *rho;//parameters in two-loop recursion
  
  int bound;//1=clip x using upper and lower bounds,0=not
  float *xmin;//lower bound for x
  float *xmax;//upper bound for x
  
  int verb;//verbose information to print or not
  int preco;//precondition or not
  int loop1;//flag to check 1st loop in two-loop recursion is done or not
  int ls_fail; //0=misfit decreases, accept it; 1=line search fails
    
  int ncg;//number of CG to solve normal equation for Newton method
  int method;//method=1, lbfgs; method=2, Newton-CG
} opt_t;

typedef float (*opt_fg)(float*,float*); 
typedef void (*opt_Hv)(float*, float *, float*);

#endif
