#include "general.h"

void GenereVecteurInitialGaussien(TVecteur<double> & VecteurEntree,
				  float coef,
				  struct drand48_data *tampon)
{
 double resultat,v1,v2,fac, rsq;
 int i;
 int N = VecteurEntree.size();
 for(i = 0; i<N; i+=2)
 {
  do{
    drand48_r(tampon, &resultat);
    v2=2.0*((resultat))-1.0; 
    drand48_r(tampon, &resultat);
    v1=2.0*((resultat))-1.0; 	
    rsq=v1*v1+v2*v2;
  }while (rsq>= 1.0 || rsq==0.0);
  fac=sqrt(-2.0*log(rsq)/rsq);
  VecteurEntree(i) = v1*fac*coef;
  VecteurEntree((i+1)%N) = v2*fac*coef;
 }
}
