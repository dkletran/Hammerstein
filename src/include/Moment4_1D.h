//  



#ifndef _TMoment4_1D_
#define _TMoment4_1D_

#include <iostream>       // for formatted printing of matrices
#include "Moment_1D.h"
using namespace std;
template <class T>

//TMoment4_1D<T> operator*(T,TMoment4_1D<T>);

class TMoment4_1D : public  TMoment_1D<T>
{ 
  protected :
    T     ***v;
    int   fenetre;  // size of original matrix, not submatrix
    
 public:



        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


    TMoment4_1D();
    TMoment4_1D(int);
    virtual ~TMoment4_1D();


        /*::::::::::::::::::::::::::::::::*/

        /*  Indices and access operations */

        /*::::::::::::::::::::::::::::::::*/

    inline int size() const;   // submatrix size  
    inline T operator()(int i, int j, int k) const;  
    inline T& operator()(int i, int j, int k);
    
      

};  //* End of TMoment4_1D<T> Class *//

   


template <class T>
inline int TMoment4_1D<T>::size() const
{
   
  return 2*(this->fenetre);
    
}

template <class T>
inline T TMoment4_1D<T>::operator()(int i, int j, int k) const
{
  int cord[6];
  cord[0]=i;
  cord[1]=j;
  cord[2]=k;


  TMoment_1D<T>::SymetrieMoment(4,cord);




  return v[cord[0]][cord[1]][cord[2]];
}

template <class T>
inline T& TMoment4_1D<T>::operator()(int i, int j, int k) 
{
  
  
  return v[i][j][k];
  
}



#endif 

