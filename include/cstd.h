#ifndef cstd_h
#define cstd_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <limits.h>

#ifndef NULL
#define NULL	((void *)0)
#endif
#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef NINT
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#endif

/* allocate and free multi-dimensional arrays */
void *alloc1 (size_t n1, size_t size);
void *realloc1 (void *v, size_t n1, size_t size);
void **alloc2 (size_t n1, size_t n2, size_t size);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6, 
                   size_t size);
void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free4 (void ****p);
void free5 (void *****p);
void free6 (void ******p);

int *alloc1int (size_t n1);
int *realloc1int (int *v, size_t n1);
int **alloc2int (size_t n1, size_t n2);
int ***alloc3int (size_t n1, size_t n2, size_t n3);
int ****alloc4int (size_t n1, size_t n2, size_t n3, size_t n4);
int *****alloc5int (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free1int (int *p);
void free2int (int **p);
void free3int (int ***p);
void free4int (int ****p);
void free5int (int *****p);

float *alloc1float (size_t n1);
float *realloc1float (float *v, size_t n1);
float **alloc2float (size_t n1, size_t n2);
float ***alloc3float (size_t n1, size_t n2, size_t n3);
float ****alloc4float (size_t n1, size_t n2, size_t n3, size_t n4);
float *****alloc5float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ******alloc6float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6);
void free1float (float *p);
void free2float (float **p);
void free3float (float ***p);
void free4float (float ****p);
void free5float (float *****p);
void free6float (float ******p);

double *alloc1double (size_t n1);
double *realloc1double (double *v, size_t n1);
double **alloc2double (size_t n1, size_t n2);
double ***alloc3double (size_t n1, size_t n2, size_t n3);
void free1double (double *p);
void free2double (double **p);
void free3double (double ***p);


float _Complex *alloc1complexf(size_t n1);
float _Complex *realloc1complexf(float _Complex *v, size_t n1);
float _Complex **alloc2complexf(size_t n1, size_t n2);
float _Complex ***alloc3complexf(size_t n1, size_t n2, size_t n3);
float _Complex ****alloc4complexf(size_t n1, size_t n2, size_t n3, size_t n4);
void free1complexf(float _Complex *p);
void free2complexf(float _Complex **p);
void free3complexf(float _Complex ***p);
void free4complexf(float _Complex ****p);

double _Complex *alloc1complex(size_t n1);
double _Complex *realloc1complex(double _Complex *v, size_t n1);
double _Complex **alloc2complex(size_t n1, size_t n2);
double _Complex ***alloc3complex(size_t n1, size_t n2, size_t n3);
double _Complex ****alloc4complex(size_t n1, size_t n2, size_t n3, size_t n4);
void free1complex(double _Complex *p);
void free2complex(double _Complex **p);
void free3complex(double _Complex ***p);
void free4complex(double _Complex ****p);

char *alloc1char(size_t n1);
char *realloc1char(char *v, size_t n1);
void free1char(char *p);

/* GLOBAL DECLARATIONS */
extern int xargc; 
extern char **xargv;

/* FUNCTION PROTOTYPES */

#ifdef __cplusplus  /* if C++, specify external C linkage */
extern "C" {
#endif

/* getpar parameter parsing */
void initargs(int argc, char **argv);
int getparint(char *name, int *p);
int getparuint(char *name, unsigned int *p);
int getparshort(char *name, short *p);
int getparushort(char *name, unsigned short *p);
int getparlong(char *name, long *p);
int getparulong(char *name, unsigned long *p);
int getparfloat(char *name, float *p);
int getpardouble(char *name, double *p);
int getparstring(char *name, char **p);
int getparstringarray(char *name, char **p);
int getnparint(int n, char *name, int *p);
int getnparuint(int n, char *name, unsigned int *p);
int getnparshort(int n, char *name, short *p);
int getnparushort(int n, char *name, unsigned short *p);
int getnparlong(int n, char *name, long *p);
int getnparulong(int n, char *name, unsigned long *p);
int getnparfloat(int n, char *name, float *p);
int getnpardouble(int n, char *name, double *p);
int getnparstring(int n, char *name, char **p);
int getnparstringarray(int n, char *name, char **p);
int getnpar(int n, char *name, char *type, void *ptr);
int countparname(char *name);
int countparval(char *name);
int countnparval(int n, char *name);
void checkpars( void );

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

#ifdef __cplusplus  /* if C++ (external C linkage is being specified) */
}
#endif

void err(char *fmt, ...);
void warn(char *fmt, ...);

#endif
