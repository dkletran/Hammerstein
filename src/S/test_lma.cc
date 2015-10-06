/*
 Programme pour le test de l'algorithme LMA
*/

#include <iostream>
#include <iomanip>
#include <complex>
#include "Vecteur.h"
#include "steepest_descent.h"
#include "general.h"
#include "fft.h"
//#include "ConjugateGradient.h"
#include "LM_algo.h"
using namespace std;
#define _USE_MATH_DEFINES
#define pi M_PI
#define NUMPOINTS 	16
#define DIM 		2*NUMPOINTS	
complex<double> S3[NUMPOINTS][NUMPOINTS]; //order-3 cumulant
complex<double> S4[NUMPOINTS][NUMPOINTS][NUMPOINTS];//order-4 cumulant

void estSpectra()
{ 
  int i, j, k, l;
  complex<double> H1[NUMPOINTS], H2[NUMPOINTS];
  complex<double> iu = complex<double>(0, 1);
  for(i = 0; i<NUMPOINTS; i++)
  {
    H1[i] = i*1.0;//(1.0 - 2.0*0.5*cos(pi/4)*exp(-iu*(i*pi/NUMPOINTS)) + 0.5*0.5*exp(-iu*(2*pi*i/NUMPOINTS)));
    H2[i] = 1.0;//(1.0 - 2.0*0.9*cos(pi/6)*exp(-iu*(i*pi/NUMPOINTS)) + 0.9*0.9*exp(-iu*(2*pi*i/NUMPOINTS)));
  }
  //for (i = 0; i<NUMPOINTS; i++)
   //     cout<<setprecision(5)<<abs(H1[i])<<"<"<<arg(H1[i])<<endl;
  for (i = 0; i<NUMPOINTS; i++)
   for(j = i; j < NUMPOINTS-i; j++)
   {
     k = i+j;
     S3[i][j] =       2.0*H1[i]*H1[j]*conj(H2[k])
                    + 2.0*H1[i]*H2[j]*conj(H1[k])
		    + 2.0*H2[i]*H1[j]*conj(H1[k])
		    + 8.0*H2[i]*H2[j]*conj(H2[k]);
   }
  for (i = 0; i<NUMPOINTS; i++)
   for (j = i; j<NUMPOINTS-i; j++)
    for (k = j; k<NUMPOINTS-i-j; k++)
    {
      l = i+j+k;
      S4[i][j][k] =       8.0*H1[i]*H1[j]*H2[k]*conj(H2[l])
			+ 8.0*H2[i]*H2[j]*H1[k]*conj(H1[l])
			+ 8.0*H1[i]*H2[j]*H1[k]*conj(H2[l])
			+ 8.0*H1[i]*H2[j]*H2[k]*conj(H1[l])
			+ 8.0*H2[i]*H1[j]*H1[k]*conj(H2[l])
			+ 8.0*H2[i]*H1[j]*H2[k]*conj(H1[l])
			+ 48.0*H2[i]*H2[j]*H2[k]*conj(H2[l]);


   }
}



void ffx(double *fx, TVecteur<double> *jr, double **J, TVecteur<double> &x)
{
 int i, j, k, l, s, i1, j1, k1, l1, s1;
 complex<double> tmp1, tmp2, H1_i, H1_j, H1_k, H1_l, H2_i, H2_j, H2_k, H2_l;
 double *H1, *H2, *jr1, *jr2;
 H1 = x.addr();
 H2 = &H1[2*NUMPOINTS];
 jr1 = jr->addr();
 jr2 = &jr1[2*NUMPOINTS];
 for (i = 0; i < 2*NUMPOINTS; i++) jr1[i] =  jr2[i] = 0.0;
 for (*fx = 0, s = 0, i = 0; i < 2*NUMPOINTS; i++, i++)
  for (j = i; j < 2*(NUMPOINTS-i); j++, j++, s++, s++)
  {
   k = i+j;
   i1 = i+1; j1 = j+1; k1 = k+1; s1 = s+1;
   J[s][i] = J[s][i+2*NUMPOINTS] = J[s][i1] = J[s][i1+2*NUMPOINTS]
   = J[s][j] = J[s][j+2*NUMPOINTS] = J[s][j1] = J[s][j1+2*NUMPOINTS]
   = J[s][k] = J[s][k+2*NUMPOINTS] = J[s][k1] = J[s][k1+2*NUMPOINTS]
   = J[s1][i] = J[s1][i+2*NUMPOINTS] = J[s1][i1] = J[s1][i1+2*NUMPOINTS]
   = J[s1][j] = J[s1][j+2*NUMPOINTS] = J[s1][j1] = J[s1][j1+2*NUMPOINTS]
   = J[s1][k] = J[s1][k+2*NUMPOINTS] = J[s1][k1] = J[s1][k1+2*NUMPOINTS] = 0.0;   
   H1_i = complex<double>(H1[i], H1[i1]);
   H1_j = complex<double>(H1[j], H1[j1]);
   H1_k = complex<double>(H1[k], -H1[k1]);
   H2_i = complex<double>(H2[i], H2[i1]);
   H2_j = complex<double>(H2[j], H2[j1]);
   H2_k = complex<double>(H2[k], -H2[k1]);
   tmp1 = S3[i/2][j/2] - 2.0*H1_i*H1_j*H2_k
			  - 2.0*H1_i*H2_j*H1_k
			  - 2.0*H2_i*H1_j*H1_k
			  - 8.0*H2_i*H2_j*H2_k;
 
   *fx += abs(tmp1*tmp1);
   tmp2 = 2.0*(H1_j*H2_k + H2_j*H1_k);
   J[s][i] -= tmp2.real();
   J[s][i1] += tmp2.imag();
   J[s1][i] -= tmp2.imag();
   J[s1][i1] -= tmp2.real();
   jr1[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr1[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*H1_j*H1_k+ 8.0*H2_j*H2_k;
   J[s][i+2*NUMPOINTS] -= tmp2.real();
   J[s][i1+2*NUMPOINTS] += tmp2.imag();
   J[s1][i+2*NUMPOINTS] -= tmp2.imag();
   J[s1][i1+2*NUMPOINTS] -= tmp2.real();   
   jr2[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr2[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*(H1_i*H2_k + H2_i*H1_k);
   J[s][j] -= tmp2.real();
   J[s][j1] += tmp2.imag();
   J[s1][j] -= tmp2.imag();
   J[s1][j1] -= tmp2.real();
   jr1[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr1[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*H1_i*H1_k+ 8.0*H2_i*H2_k;
   J[s][j+2*NUMPOINTS] -= tmp2.real();
   J[s][j1+2*NUMPOINTS] += tmp2.imag();
   J[s1][j+2*NUMPOINTS] -= tmp2.imag();
   J[s1][j1+2*NUMPOINTS] -= tmp2.real();   
   jr2[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr2[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();
   
   tmp2 = 2.0*(H1_i*H2_j + H2_i*H1_j);
   J[s][k] -= tmp2.real();
   J[s][k1] -= tmp2.imag();
   J[s1][k] -= tmp2.imag();
   J[s1][k1] += tmp2.real();
   jr1[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr1[k1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*H1_i*H1_j+ 8.0*H2_i*H2_j;
   
   J[s][k+2*NUMPOINTS] -= tmp2.real();
   J[s][k1+2*NUMPOINTS] -= tmp2.imag();
   J[s1][k+2*NUMPOINTS] -= tmp2.imag();
   J[s1][k1+2*NUMPOINTS] += tmp2.real();   
   jr2[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr2[k1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();
  }
 for (i = 0; i < 2*NUMPOINTS; i++, i++)
  for (j = i; j < 2*(NUMPOINTS-i); j++, j++)
    for (k = j; k < 2*(NUMPOINTS-i-j);k++, k++,s++,s++)
    {
     l = i+j+k;
     i1 = i+1; j1 = j+1; k1 = k+1; l1 = l+1; s1 = s+1;
     J[s][i] = J[s][i+2*NUMPOINTS] = J[s][i1] = J[s][i1+2*NUMPOINTS]
     = J[s][j] = J[s][j+2*NUMPOINTS] = J[s][j1] = J[s][j1+2*NUMPOINTS]
     = J[s][k] = J[s][k+2*NUMPOINTS] = J[s][k1] = J[s][k1+2*NUMPOINTS]
     = J[s][l] = J[s][l+2*NUMPOINTS] = J[s][l1] = J[s][l1+2*NUMPOINTS]
     = J[s1][i] = J[s1][i+2*NUMPOINTS] = J[s1][i1] = J[s1][i1+2*NUMPOINTS]
     = J[s1][j] = J[s1][j+2*NUMPOINTS] = J[s1][j1] = J[s1][j1+2*NUMPOINTS]
     = J[s1][k] = J[s1][k+2*NUMPOINTS] = J[s1][k1] = J[s1][k1+2*NUMPOINTS] 
     = J[s1][l] = J[s1][l+2*NUMPOINTS] = J[s1][l1] = J[s1][l1+2*NUMPOINTS] = 0.0;   
     H1_i = complex<double>(H1[i], H1[i1]);
     H1_j = complex<double>(H1[j], H1[j1]);
     H1_k = complex<double>(H1[k], H1[k1]);
     H1_l = complex<double>(H1[l], -H1[l1]);
     H2_i = complex<double>(H2[i], H2[i1]);
     H2_j = complex<double>(H2[j], H2[j1]);
     H2_k = complex<double>(H2[k], H2[k1]);     
     H2_l = complex<double>(H2[l], -H2[l1]);
     tmp1 =  S4[i/2][j/2][k/2]- 8.0*H1_i*H1_j*H2_k*H2_l
	                - 8.0*H2_i*H2_j*H1_k*H1_l
		        - 8.0*H1_i*H2_j*H1_k*H2_l
			- 8.0*H1_i*H2_j*H2_k*H1_l
			- 8.0*H2_i*H1_j*H1_k*H2_l
			- 8.0*H2_i*H1_j*H2_k*H1_l
			- 48.0*H2_i*H2_j*H2_k*H2_l;
 
    *fx += abs(tmp1*tmp1);
     
     tmp2 = 8.0*(H1_j*H2_k*H2_l+H2_j*H1_k*H2_l+H2_j*H2_k*H1_l);
     J[s][i] -= tmp2.real();
     J[s][i1] += tmp2.imag();
     J[s1][i] -= tmp2.imag();
     J[s1][i1] -= tmp2.real();
     jr1[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_j*H1_k*H1_l+H1_j*H1_k*H2_l+H1_j*H2_k*H1_l+6.0*H2_j*H2_k*H2_l);
     J[s][i+2*NUMPOINTS] -= tmp2.real();
     J[s][i1+2*NUMPOINTS] += tmp2.imag();
     J[s1][i+2*NUMPOINTS] -= tmp2.imag();
     J[s1][i1+2*NUMPOINTS] -= tmp2.real();   
     jr2[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H1_i*H2_k*H2_l+H2_i*H1_k*H2_l+H2_i*H2_k*H1_l);
     J[s][j] -= tmp2.real();
     J[s][j1] += tmp2.imag();
     J[s1][j] -= tmp2.imag();
     J[s1][j1] -= tmp2.real();
     jr1[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_i*H1_k*H1_l+H1_i*H1_k*H2_l+H1_i*H2_k*H1_l+6.0*H2_i*H2_k*H2_l);
     J[s][j+2*NUMPOINTS] -= tmp2.real();
     J[s][j1+2*NUMPOINTS] += tmp2.imag();
     J[s1][j+2*NUMPOINTS] -= tmp2.imag();
     J[s1][j1+2*NUMPOINTS] -= tmp2.real();   
     jr2[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real(); 

     tmp2 = 8.0*(H1_i*H2_j*H2_l+H2_i*H1_j*H2_l+H2_i*H2_j*H1_l);
     J[s][k] -= tmp2.real();
     J[s][k1] += tmp2.imag();
     J[s1][k] -= tmp2.imag();
     J[s1][k1] -= tmp2.real();
     jr1[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[k1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_i*H1_j*H1_l+H1_i*H1_j*H2_l+H1_i*H2_j*H1_l+6.0*H2_i*H2_j*H2_l);
     J[s][k+2*NUMPOINTS] -= tmp2.real();
     J[s][k1+2*NUMPOINTS] += tmp2.imag();
     J[s1][k+2*NUMPOINTS] -= tmp2.imag();
     J[s1][k1+2*NUMPOINTS] -= tmp2.real();   
     jr2[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[k1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H1_i*H2_j*H2_k+H2_i*H1_j*H2_k+H2_i*H2_j*H1_k);
     J[s][l] -= tmp2.real();
     J[s][l1] -= tmp2.imag();
     J[s1][l] -= tmp2.imag();
     J[s1][l1] += tmp2.real();
     jr1[l] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[l1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_i*H1_j*H1_k+H1_i*H1_j*H2_k+H1_i*H2_j*H1_k+6.0*H2_i*H2_j*H2_k);
   
     J[s][l+2*NUMPOINTS] -= tmp2.real();
     J[s][l1+2*NUMPOINTS] -= tmp2.imag();
     J[s1][l+2*NUMPOINTS] -= tmp2.imag();
     J[s1][l1+2*NUMPOINTS] += tmp2.real();   
     jr2[l] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[l1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();
    }
}
int main(int argc, const char* argv[])
{
 int i, j, k, l, s;
 complex<double> tmp;
 for (s=0, i = 0; i < NUMPOINTS; i++)
  for (j = i; j < NUMPOINTS-i; j++, s++)
   S3[i][j] = 14.0;
 for (i = 0; i < NUMPOINTS; i++)
  for (j = i; j < NUMPOINTS-i; j++)
    for (k = j; k < NUMPOINTS-i-j;k++,s++)
     S4[i][j][k] = 96;
 
 LM_Minimiser<double> lm(2*DIM, 2*s);
 TVecteur<double> x(2*DIM, 0);
 for (i = 0; i < NUMPOINTS; i++) {x(2*i) = 1.50*i;x(2*i+1) = 0.50;}
 for (; i < 2*NUMPOINTS; i++) {x(2*i) = 1.50;x(2*i+1) = 0.1;}
 estSpectra();
 lm.set_f(ffx);
 lm.set_x0(x);
 lm.execute(0.0001, 0.000000001, 100000000, 10, 100000, &x);
 cout <<"Solution:"<<endl;
 for(i = 0; i<DIM; i++)
 {  
  tmp = complex<double>(x(2*i), x(2*i+1));
  cout << abs(tmp)<<"<"<<arg(tmp) <<endl;
  //cout<<x(i)<<endl;
 }

 return 0;
}

