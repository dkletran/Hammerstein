
#include "general.h"

void EstimationCumulant(TMoment4_1D<double> *Cumulant4,
			TMoment3_1D<double> *Cumulant3,
			TVecteur<double> *Cumulant2,
			TMoment4_1D<double> *Moment4,
			TMoment3_1D<double> *Moment3,
			TVecteur<double> *Moment2,
			double Moyenne)
{

  int I1=(Moment4->size())/2;


  for(int i=0; i<=I1;i++)
    for(int j=0; j<=i;j++)
      for(int k=0; k<=j;k++)
	(*Cumulant4)(i,j,k)=(*Moment4)(i,j,k)
	  -(*Moment2)(i)*(*Moment2)(j-k)
	  -(*Moment2)(j)*(*Moment2)(i-k)
	  -(*Moment2)(k)*(*Moment2)(i-j)
	  -Moyenne*
	  ((*Moment3)(i-k,j-k)+(*Moment3)(i,j)
	   +(*Moment3)(i,k)+(*Moment3)(j,k))
	  +2*pow(Moyenne,2)*
	  ((*Moment2)(i)+(*Moment2)(j)+(*Moment2)(k)
	   +(*Moment2)(i-j)+(*Moment2)(i-k)
	   +(*Moment2)(j-k))
	  -6*pow(Moyenne,4);
  
  
  for(int i=0; i<=I1;i++)
    for(int j=0; j<=i;j++)
      {
	(*Cumulant3)(i,j)=(*Moment3)(i,j)
	  -Moyenne*(*Moment2)(i)
	  -Moyenne*(*Moment2)(j)-Moyenne*(*Moment2)(i-j)+
	  2*pow(Moyenne,3);
		
      }


  for(int i=0; i<2*I1;i++)
    (*Cumulant2)(i)=(*Moment2)(i)-pow(Moyenne,2);

}






