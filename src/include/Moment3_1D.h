//  



#ifndef _Moment3_1D_
#define _Moment3_1D_

#include <iostream>       
#include "Moment_1D.h"
/* using namespace std; */

template <class T>
class TMoment3_1D : public TMoment_1D<T>
{ 

  friend class TMoment_1D<T>;

  protected :
    T     **v;
    int   fenetre;  // size of original matrix, not submatrix
    
 public:


        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


    TMoment3_1D();
    TMoment3_1D(int);
    virtual ~TMoment3_1D();


        /*::::::::::::::::::::::::::::::::*/

        /*  Indices and access operations */

        /*::::::::::::::::::::::::::::::::*/

    inline int size() const;   // submatrix size  
    inline T operator()(int i, int j) const;  
    inline T & operator()(int i, int j);
    TMoment3_1D<T>& operator=(T s);
    


};  //* End of TMoment3_1D<T> Class *//

   


template <class T>
inline int TMoment3_1D<T>::size() const
{
  /////////////////////////////////////////
  // 
  // On renvoie toujours la taille entier du 
  // support
  //
  //////////////////////////////////////////
   
  return 2*(this->fenetre);
    
}

template <class T>
inline T TMoment3_1D<T>::operator()(int i, int j) const
{

  /////////////////////////////////////////////////
  // 
  // Procedure de lecture seule, faisant les symétries ...
  //
  ///////////////////////////////////////////////////
  int *cord= new int [6];
  cord[0]=i;
  cord[1]=j;
  
  TMoment_1D<T>::SymetrieMoment(3,cord); 
  

  if ((cord[0] > this->fenetre) || (cord[1] > this->fenetre))
    {
      cerr <<"Erreur dans le calcul de symomment";
      return((T)0);
    }

  int n1=cord[0];
  int n2=cord[1];

  delete [] cord;

  return v[n1][n2];

 
}

template <class T>
inline T& TMoment3_1D<T>::operator()(int i, int j) 
{ 

/*   cerr << "On est dans la procedure d'ecriture"<<endl; */
  /////////////////////////////////////////////////
  // 
  // Procedure de lecture-ecriture ne travaillant
  // que sur le support non redondant
  //
  ///////////////////////////////////////////////////
   //cerr <<i <<"  "<<j<<endl; 
 return v[i][j];
}



#endif 

