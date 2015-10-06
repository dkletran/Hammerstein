#include <complex>
#include <iostream>
#include "Moment4_1D.h"


template <class T>
TMoment4_1D<T>::~TMoment4_1D()
{
  int i, j;
 
  for (i = 0; i <= fenetre; i++ )
    for(j = 0; j <=i; j++)
    	delete[] v[i][j];
  for (i = 0; i <= fenetre; i++)
    delete[] v[i];
    
  delete[] v;

  
  v=NULL;
  fenetre = 0;
  

}

template <class T>
TMoment4_1D<T>::TMoment4_1D() 
{
  
  fenetre = 0;
  v=NULL;
}


template <class T>
TMoment4_1D<T>::TMoment4_1D(int n) 
{
  int i,j,k;

  ////////////////////////////////////////////////
  // La Taille passe est la taille entiere
  // et non plus la taille du support non redondant
  // Cependant le parametre sauvegarde (Private)
  // est la taille du support.
  ////////////////////////////////////////////////

  fenetre = n/2;
 

  v=new T** [fenetre+1];

  for ( i = 0; i <=fenetre; i++ )
    v[i]= new T* [ i+1 ];


  for ( i = 0; i <=fenetre; i++ )
    for ( j = 0; j <=i; j++ )
      v[i][j]=new T [j+1];
  for ( i = 0; i <=fenetre; i++ )
    for ( j = 0; j <=i; j++ )
      for (k = 0; k <= j; k++)
        v[i][j][k] = ((T)0);
}
/*
template <class T>
TMoment4_1D<T>& TMoment4_1D<T>::operator=(T s)
{
  int i,j; 

  for ( i = 0; i <=fenetre; i++ )
    for ( j = 0; j <=i; j++ )
     for ( k = 0; k <=j; k++ )
      v[i][j][k]=s;

    return *this;

}
*/


template class TMoment4_1D<short>;
template class TMoment4_1D<float>;
template class TMoment4_1D<double>;
template class TMoment4_1D<complex<double> >;

