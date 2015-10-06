//  



#ifndef _Moment_1D_
#define _Moment_1D_

#include <iostream>    
#include <stdio.h>  
using namespace std;
template <class T>


class TMoment_1D
{ 
  protected :

  public:
  void SymetrieMoment(int ordre, int *cord);    
  void SymetrieMoment(int ordre, int *cord) const; 

};  //* End of TMoment_1D<T> Class *//

   
template <class T>
void TMoment_1D<T>::SymetrieMoment(int ordre,int *cord)
{	
  int i,j,ValeurNegative,PlusGrand,IndicePlusGrand,ok;
  
  ValeurNegative=0;

  ok=0;

  int *cord_tam= new int [6];

  
  ///////////////////////////////////////////////
  //
  // Recherche de la plus petite valeur negative.
  //
  ///////////////////////////////////////////////

  for(i=0;i<ordre-1;i++)  
    {
      if(cord[i]<ValeurNegative)	
	  ValeurNegative=cord[i];
    }
    
  for(i=0;i<ordre-1;i++)
    {
      if((cord[i]==ValeurNegative)&&(ok==0))
	{
	  cord_tam[i]=(-ValeurNegative);
	  ok=1;
	}
      else	
	  cord_tam[i]=cord[i]-ValeurNegative;
	
    }

 
  // Un probleme dans la partie ci dessus peut se poser si 
  // toutes les valeurs sont égales à 0
  // Point 0



  for(i=0;i<ordre-1;i++)
    {
      cord[i]=0;
      PlusGrand=0;
      // Initialisation de IndicePlusGrand n'atant utile que dans le cas
      // du point 0. C'est le seul cas ou IndicePlusGrand peut ne pas etre
      // Initialise correctement car on ne passe jamais dans la boucle

      // IndicePlusGrand=i;
      for(j=0;j<ordre-1;j++)
	{
	if(cord_tam[j]>PlusGrand)
	  {
	    PlusGrand=cord_tam[j];
	    IndicePlusGrand=j;
	  }
	}
      cord[i]=PlusGrand;

 //     cout << "  " << IndicePlusGrand <<"\n";
      cord_tam[IndicePlusGrand]=0;	
    }

  delete [] cord_tam;
}
   


template <class T>
void TMoment_1D<T>::SymetrieMoment(int ordre,int cord[]) const
{	
  int i,j,ValeurNegative,PlusGrand,IndicePlusGrand,ok;

  int *cord_tam= new int [6];

  ValeurNegative=0;
  ok=0;

/*   cout <<"On rentre dans symetrie" <<endl << flush; */
  

  // Recherche de la plus petite valeur negative.
  for(i=0;i<ordre-1;i++)  
      if(cord[i]<ValeurNegative)	
	  ValeurNegative=cord[i];
	
/*   for(i=0;i<ordre-1;i++) */
/*     cerr << "---"<<cord[i]<<"\n\n"; */
    
  for(i=0;i<ordre-1;i++)
    {
      if((cord[i]==ValeurNegative)&&(ok==0))
	{
	  cord_tam[i]=(-ValeurNegative);
	  ok=1;
	}
      else	
	  cord_tam[i]=cord[i]-ValeurNegative;
	
    }

  for(i=0;i<ordre-1;i++)
    {
      cord[i]=0;
      PlusGrand=0; 
      // Initialisation de IndicePlusGrand n'atant utile que dans le cas
      // du point 0. C'est le seul cas ou IndicePlusGrand peut ne pas etre
      // Initialise correctement car on ne passe jamais dans la boucle

      IndicePlusGrand=i;
      for(j=0;j<ordre-1;j++)
	if(cord_tam[j]>PlusGrand)
	  {
	    PlusGrand=cord_tam[j];
	    IndicePlusGrand=j;
	  }
      cord[i]=PlusGrand;
      cord_tam[IndicePlusGrand]=0;	
    } 

  delete [] cord_tam;

}


#endif 

