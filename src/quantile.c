/* this is Hoare's quantile algorithm implemented in Madagascar!
 *--------------------------------------------------------------*/
double quantile(int q    /* quantile */, 
		int n    /* array length */, 
		double* a /* array [n] */) 
/*< find quantile (caution: a is changed) >*/ 
{
  double *i, *j, ak, *low, *hi, buf, *k;

  low = a;
  hi = a+n-1;
  k = a+q; 
  while (low<hi) {
    ak = *k;
    i = low; j = hi;
    do {
      while (*i < ak) i++;     
      while (*j > ak) j--;     
      if (i<=j) {
	buf = *i;
	*i++ = *j;
	*j-- = buf;
      }
    } while (i<=j);
    if (j<k) low = i; 
    if (k<i) hi = j;
  }
  return (*k);
}
