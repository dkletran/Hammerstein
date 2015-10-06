
/*
 Programme pour le test de l'algorithme steepest descent
*/

#include <iostream>
#include <iomanip>
#include <complex>
#include "Vecteur.h"
#include "steepest_descent.h"
#include "general.h"
#include "fft.h"


using namespace std;
#define _USE_MATH_DEFINES
#define pi M_PI
#define NUMPOINTS 	8
#define NUMPOINTS1      (NUMPOINTS+1)
TMoment3_1D<complex<double> > bisp(2*NUMPOINTS);
TMoment4_1D<complex<double> >trisp(2*NUMPOINTS);

/*
 La function objective 
*/

double targetc(TVecteur<complex<double> > &x)
{
   int i, j, k, l;
   double result;
   complex<double> tmp;
   complex<double> *H1, *H2;
   H1 = x.addr();
   H2 = &H1[NUMPOINTS1];
   result = 0;
   for (i = 0; i<NUMPOINTS; i++)
      for(j = i; j <= NUMPOINTS-i; j++)
	{
	   k = i+j;
	   tmp = bisp(j,i) - 2.0*H1[i]*H1[j]*conj(H2[k])
			  - 2.0*H1[i]*H2[j]*conj(H1[k])
			  - 2.0*H2[i]*H1[j]*conj(H1[k])
			  - 8.0*H2[i]*H2[j]*conj(H2[k]);
           
	  result += abs(tmp*tmp);
	}
   for (i = 0; i<NUMPOINTS; i++)
      for(j = i; j<NUMPOINTS-i; j++)
	for (k = j; k<=NUMPOINTS-i-j; k++)
	  {
	     l = i+j+k;
	     tmp =  trisp(k,j,i) - 8.0*H1[i]*H1[j]*H2[k]*conj(H2[l])
				- 8.0*H2[i]*H2[j]*H1[k]*conj(H1[l])
				- 8.0*H1[i]*H2[j]*H1[k]*conj(H2[l])
				- 8.0*H1[i]*H2[j]*H2[k]*conj(H1[l])
				- 8.0*H2[i]*H1[j]*H1[k]*conj(H2[l])
				- 8.0*H2[i]*H1[j]*H2[k]*conj(H1[l])
				- 48.0*H2[i]*H2[j]*H2[k]*conj(H2[l]);
           result += abs(tmp*tmp); 
	  }
   return result;
}
/*
 Le gradient
*/
TVecteur<complex<double> > gradientc(TVecteur<complex<double> > &x)
{
 int i, j, k, l;
 complex<double> *H1, *H2, *grad1, *grad2, tmp, tmp2;
 TVecteur<complex<double> > grad(2*NUMPOINTS1);
 H1 = x.addr();
 H2 = &H1[NUMPOINTS1];
 grad1 = grad.addr();
 grad2 = &grad1[NUMPOINTS1];
 for (i = 0; i < NUMPOINTS1; i++)
   {
    grad1[i] = 0.0;
    grad2[i] = 0.0;
   }
  for (i = 0; i<NUMPOINTS; i++)
     for(j = i; j<=NUMPOINTS-i; j++)
	{
	   k = i+j;
	   tmp = bisp(j,i) - 2.0*H1[i]*H1[j]*conj(H2[k])
			   - 2.0*H1[i]*H2[j]*conj(H1[k])
			   - 2.0*H2[i]*H1[j]*conj(H1[k])
			   - 8.0*H2[i]*H2[j]*conj(H2[k]);

           grad1[i] -= 2.0*(conj(H1[j])*H2[k] + conj(H2[j])*H1[k])*tmp;
	   grad2[i] -= (2.0*conj(H1[j])*H1[k]+ 8.0*conj(H2[j])*H2[k])*tmp;
	   grad1[j] -= 2.0*(conj(H1[i])*H2[k] + conj(H2[i])*H1[k])*tmp;
	   grad2[j] -= (2.0*conj(H1[i])*H1[k]+ 8.0*conj(H2[i])*H2[k])*tmp;
	   grad1[k] -= 2.0*(H1[i]*H2[j] + H2[i]*H1[j])*conj(tmp);
	   grad2[k] -= (2.0*H1[i]*H1[j]+ 8.0*H2[i]*H2[j])*conj(tmp);
	} 
    for (i = 0; i<NUMPOINTS; i++)
      for(j = i; j<NUMPOINTS-i; j++)
	for (k = j; k<=NUMPOINTS-i-j; k++)
	  {
	     l = i+j+k;
	     tmp =  trisp(k,j,i) - 8.0*H1[i]*H1[j]*H2[k]*conj(H2[l])
				- 8.0*H2[i]*H2[j]*H1[k]*conj(H1[l])
				- 8.0*H1[i]*H2[j]*H1[k]*conj(H2[l])
				- 8.0*H1[i]*H2[j]*H2[k]*conj(H1[l])
				- 8.0*H2[i]*H1[j]*H1[k]*conj(H2[l])
				- 8.0*H2[i]*H1[j]*H2[k]*conj(H1[l])
				- 48.0*H2[i]*H2[j]*H2[k]*conj(H2[l]);

             grad1[i] -= 8.0*(conj(H1[j]*H2[k])*H2[l]+conj(H2[j]*H1[k])*H2[l]+conj(H2[j]*H2[k])*H1[l])*tmp;
	     grad2[i] -= 8.0*(conj(H2[j]*H1[k])*H1[l]+conj(H1[j]*H1[k])*H2[l]+conj(H1[j]*H2[k])*H1[l]+6.0*conj(H2[j]*H2[k])*H2[l])*tmp;
	     grad1[j] -= 8.0*(conj(H1[i]*H2[k])*H2[l]+conj(H2[i]*H1[k])*H2[l]+conj(H2[i]*H2[k])*H1[l])*tmp;
	     grad2[j] -= 8.0*(conj(H2[i]*H1[k])*H1[l]+conj(H1[i]*H1[k])*H2[l]+conj(H1[i]*H2[k])*H1[l]+6.0*conj(H2[i]*H2[k])*H2[l])*tmp;
	     grad1[k] -= 8.0*(conj(H1[i]*H2[j])*H2[l]+conj(H2[i]*H1[j])*H2[l]+conj(H2[i]*H2[j])*H1[l])*tmp;
	     grad2[k] -= 8.0*(conj(H2[i]*H1[j])*H1[l]+conj(H1[i]*H1[j])*H2[l]+conj(H1[i]*H2[j])*H1[l]+6.0*conj(H2[i]*H2[j])*H2[l])*tmp;
	     grad1[l] -= 8.0*(H1[i]*H2[j]*H2[k]+H2[i]*H1[j]*H2[k]+H2[i]*H2[j]*H1[k])*conj(tmp);
	     grad2[l] -= 8.0*(H2[i]*H1[j]*H1[k]+H1[i]*H1[j]*H2[k]+H1[i]*H2[j]*H1[k]+6.0*H2[i]*H2[j]*H2[k])*conj(tmp);
	  }
 
 return grad;
}
int main(int argc, const char* argv[])
{
 int i, j, k, l, s;
 
 TMoment4_1D<double> Moment4(2*NUMPOINTS), Cumulant4(2*NUMPOINTS);
 TMoment3_1D<double> Moment3(2*NUMPOINTS), Cumulant3(2*NUMPOINTS);
 TVecteur<double> Moment2(2*NUMPOINTS), Cumulant2(2*NUMPOINTS);
 unsigned short int seed16v[3];
 complex<double> tmp;
 double mean;
 seed16v[0]=time(NULL);
 struct drand48_data *tampon=new  struct drand48_data [1];


 
 TVecteur<double> lineaire(3);
 TVecteur<double> quadratique(3);
 TVecteur<double> x(1024*1024*8);
 TVecteur<double> y(1024*1024*8);


 /*
 Les noyaux du filtre Hammerstein 
 */
 
 lineaire(0) = 2.0; lineaire(1) = 0.0; lineaire(2) = 0.0;
 quadratique(0) = 1.0;// quadratique(1) = 1; quadratique(2) = 1;
 /* 
 Init le solver
 */

 SteepestDescent<complex<double> > std(2*NUMPOINTS1);
 TVecteur<complex<double> > x0(2*NUMPOINTS1);
 

 std.set_f(targetc);
 std.set_g(gradientc);
 for (i = 0; i < NUMPOINTS1; i++) {x0(i) = 2.0;}
 for (; i < 2*NUMPOINTS1; i++) {x0(i) = 1.0;}
 seed48_r(seed16v,tampon);
 for(k=0; k<1000; k++)
 {
  /*
  Générer du signal gaussien
 */


 GenereVecteurInitialGaussien(x, 1.0, tampon); 
 
 /*
  Filtrage Hammerstein
 */
 
 FiltrageHammerstein(x, y, lineaire, quadratique);
 /*
  Ajouter du bruit (blanc, gaussien)
 */

 GenereVecteurInitialGaussien(x, 1.0, tampon);
 y=y+x;
 

 /*
  Estimer le bispectre et le trispectre
 */


 cout<< "Estimation moments ..."<<endl;
 EstimationMoment(&Moment4, &Moment3, &Moment2, &y, 0);
 mean = 0.0;
 for(i = 0; i<y.size(); i++) mean += y(i);
 mean /= y.size();
 cout <<"Estimation cumulants..."<<endl;
 EstimationCumulant(&Cumulant4, &Cumulant3, &Cumulant2, &Moment4, &Moment3, &Moment2, mean);
 cout <<"Estimation polyspectre..."<<endl;
 EstimationSpectre3(&bisp, &Cumulant3);
 EstimationSpectre4(&trisp, &Cumulant4);
 

 
 /*
 Optimisation
 */
 cout <<"Solving..."<<endl;
 std.set_x0(x0);
 cout << "Return code: "<<
 std.execute(0.00001, 0.00000000001, 100000, &x0);
 cout <<endl;
 cout <<"Solution:"<<endl;
 for(i = 0; i<2*NUMPOINTS1; i++)
 {  
  tmp = x0(i);
  cout << abs(tmp)<<"<"<<arg(tmp) <<endl;

 }
 }

 return 0;
}

