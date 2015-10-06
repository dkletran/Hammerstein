/*
 Programme pour tester l'estimation des polyspectres
*/
#include <iostream>
#include <iomanip>
#include <complex>
#include "Vecteur.h"
#include "steepest_descent.h"
#include "general.h"
#include "fft.h"

int main(int argc, const char* argv[])
{
 TVecteur<double> lineaire(3);
 TVecteur<double> quadratique(3);
 TVecteur<double> x(128*128);
 TVecteur<double> y(128*128);
 //TVecteur<complex<double> >  w(1024);
 //TVecteur<complex<double> > ts(128, (complex<double>)1.0);
 int i;
 lineaire(0) = 1.0; //lineaire(1) = -0.5043; lineaire(2) = 0.25;
 quadratique(0) = 1.0;//quadratique(1) = -0.8776; quadratique(2) = 0.25;
 unsigned short int seed16v[3];
 seed16v[0]=time(NULL);
 struct drand48_data *tampon=new  struct drand48_data [1];
 seed48_r(seed16v,tampon);
 GenereVecteurInitialGaussien(x, 1.0, tampon);
 
 FiltrageHammerstein(x, y, lineaire, quadratique); 
 //GenereVecteurInitialGaussien(x, 1.0, tampon);
 //y=y+x;
 TMoment4_1D<double> Moment4(32), Cumulant4(32);
 TMoment3_1D<double> Moment3(32), Cumulant3(32);
 TVecteur<double> Moment2(32), Cumulant2(32);
 TMoment3_1D<complex<double> > bisp(32);
 TMoment4_1D<complex<double> >trisp(32);
 cout<< "Estimation moments ..."<<endl;
 EstimationMoment(&Moment4, &Moment3, &Moment2, &y, 0);
 double mean = 0.0;
 for(i = 0; i<y.size(); i++) mean += y(i);
 mean /= y.size();
 cout <<"Estimation cumulants..."<<endl;
 EstimationCumulant(&Cumulant4, &Cumulant3, &Cumulant2, &Moment4, &Moment3, &Moment2, mean);
 cout <<"Estimation poly spectre..."<<endl;
 EstimationSpectre3(&bisp, &Cumulant3);
 EstimationSpectre4(&trisp, &Cumulant4);

 cout<<"S4(0, 0, 0) = "<<trisp(0,0,0)<<endl;
 
 

 return 0;
}
