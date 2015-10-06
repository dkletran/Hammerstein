#include "general.h"
#include "fft.h"

void EstimationSpectre3(TMoment3_1D<complex<double> > *bisp,//bispectrum
			TMoment3_1D<double> *cum//cumulant
                        )
{
 int i,j, cord[2];
 int N = cum->size();
 complex<double> **cum3;
 complex<double> *w = new complex<double>[N];
 cum3 = new complex<double> *[N];
 for (i = 0; i < N; i++) cum3[i] = new complex<double> [N];
 for (i = 0; i < N; i++)
  for(j = 0; j < N; j++)
  {
   cord[0] = (i<N/2)?i: i-N;
   cord[1] = (j<N/2)?j: j-N;
   cum->SymetrieMoment(3, cord);
   if(cord[0] > N/2){ 
    cum3[i][j] = 0;
   }else
   {
   cum3[i][j] = ( *cum)(cord[0], cord[1]);
   }
  }
 compute_w(N, w);
 fft2(N, cum3, w);
 for(i = 0; i <=N/2; i++)
  for(j = 0; j <=i; j++)
   (*bisp)(i, j) = cum3[i][j];
 delete[] w;
 for (i = 0; i < N; i++) delete[] cum3[i];
 delete [] cum3;
}
void EstimationSpectre4(TMoment4_1D<complex<double> > *trisp,//trispectrum
			TMoment4_1D<double> *cum//cumulant
                        )
{
 int i, j, k, cord[3];
 int N = cum->size();
 complex<double> ***cum4;
 complex<double> *w = new complex<double>[N];
 cum4 = new complex<double> **[N];
 for (i = 0; i < N; i++) cum4[i] = new complex<double> *[N];
 for (i = 0; i < N; i++)
  for(j = 0; j < N; j++)
   cum4[i][j] = new complex<double> [N];
 for (i = 0; i < N; i++)
  for (j = 0; j < N; j++)
   for (k = 0; k < N; k++)
  {
   cord[0] = (i<N/2)?i:i-N;
   cord[1] = (j<N/2)?j:j-N;
   cord[2] = (k<N/2)?k:k-N;
   cum->SymetrieMoment(4, cord);
   if(cord[0] > N/2){ 
    cum4[i][j][k] = 0;
   }else
   {
   cum4[i][j][k] = ( *cum)(cord[0], cord[1], cord[2]);
   }
  }
 compute_w(N, w);
 fft3(N, cum4, w);
 for(i = 0; i <=N/2; i++)
  for(j = 0; j <=i; j++)
   for(k = 0; k<=j; k++)
   (*trisp)(i, j, k) = cum4[i][j][k];
 delete[] w;
 for (i = 0; i < N; i++)
 {
  for (j = 0; j< N; j++) delete[] cum4[i][j];
  delete[] cum4[i];
 }
 delete [] cum4;
}
