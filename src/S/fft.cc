#include "fft.h"
#include <complex>
using namespace std;

/* treats inp as a numbits number and bitreverses it. 
 * inp < 2^(numbits) for meaningful bit-reversal
 */ 
int bitrev(int inp, int numbits)
{
  int i, rev=0;
  for (i=0; i < numbits; i++)
  {
    rev = (rev << 1) | (inp & 1);
    inp >>= 1;
  }
  return rev;
}
/* returns log n (to the base 2), if n is positive and power of 2 */ 
int log_2(int n) 
{
  int res; 
  for (res=0; n >= 2; res++) 
    n = n >> 1; 
  return res; 
}
 
/* W will contain roots of unity so that W[bitrev(i,log2n-1)] = e^(2*pi*i/n)
 * n should be a power of 2
 * Note: W is bit-reversal permuted because fft(..) goes faster if this is done.
 *       see that function for more details on why we treat 'i' as a (log2n-1) bit number.
 */
void compute_w(int n, complex<double> *w)
{
  int i, br;
  int log2n = log_2(n);

  for (i=0; i<(n/2); i++)
  {
    br = bitrev(i,log2n-1); 
    w[br] = complex<double>(cos(((double)i*2.0*M_PI)/((double)n)),
                      -sin(((double)i*2.0*M_PI)/((double)n)));  
  }
}
/* permutes the array using a bit-reversal permutation */ 
void permute_bitrev(int n, complex<double> *in) 
{ 
  int i, bri, log2n;
  complex<double> t;

  log2n = log_2(n); 
  
  for (i=0; i<n; i++)
  {
      bri = bitrev(i, log2n);

      /* skip already swapped elements */
      if (bri <= i) continue;
      t = in[i];
      in[i] = in[bri];
      in[bri] = t;
  }  
} 
/*
 1D fft algorithm
 */
void fft(int n, //length of signal, should be a power of 2
	 complex<double> * out, //input - output
	 complex<double> * _w  //roots of unity
       )
{
 int m, g, b;
 int i, mt, k;
 complex<double> w, t, u;
 /* for each stage */
 for (m = n; m >= 2; m >>=1)
 {
  mt = m>>1;
  /* for each group of butterfly */ 
  /*each group use only one root of unity */
  for (g=0,k=0; g<n; g+=m,k++) 
  {
   w = _w[k];
   /* for each butterfly */ 
   for (b=g; b<(g+mt); b++) 
   {
    u = out[b];
    t = w*out[b+mt];
    out[b] = u + t;
    out[b+mt] = u -t;
   }
  }
 }
 permute_bitrev(n, out);
}
/*
 2D fft
*/
void fft2(int n, 
          complex<double> ** out,
          complex<double> * _w)
{
 int i, j;
 complex<double> *x = new complex<double>[n];
 for(i = 0; i < n ; i++)
 {
  fft(n, out[i], _w);
 }
 for (i = 0; i < n; i++)
 {
  for(j = 0; j < n; j++) x[j] = out[j][i];
  fft(n, x, _w); 
  for(j = 0; j < n; j++) out[j][i] = x[j];
 }
 delete[] x;
}
/*
 3D fft
 */
void fft3(int n, 
          complex<double> *** out,
          complex<double> *_w)
{
 int i, j, k;
 complex<double> *x = new complex<double> [n];
 for (k = 0; k < n; k++)
 {
  for(i = 0; i < n ; i++)
   {
    fft(n, out[k][i], _w);
   }
  for (i = 0; i < n; i++)
   {
    for(j = 0; j < n; j++) x[j] = out[k][j][i];
    fft(n, x, _w); 
    for(j = 0; j < n; j++) out[k][j][i] = x[j];
   }
 }
 for(i = 0; i < n; i++)
 {
  for(j = 0; j < n; j++)
  {
   for(k = 0; k < n; k++) x[k] = out[k][i][j];
   fft(n, x, _w);
   for(k = 0; k < n; k++) out[k][i][j]  = x[k];
  }
 }
 delete[] x;
}
