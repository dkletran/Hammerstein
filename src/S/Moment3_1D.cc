
#include <complex>
#include <iostream>
#include "Moment3_1D.h"


template <class T>
TMoment3_1D<T>::~TMoment3_1D()
{
  

 
  for (int i = 0; i <= fenetre; i++ )
    delete[] v[i];

    
  delete[] v;

  
  v=NULL;
  fenetre = 0;
  

}

template <class T>
TMoment3_1D<T>::TMoment3_1D() 
{
  
  fenetre = 0;
  v=NULL;
}


template <class T>
TMoment3_1D<T>::TMoment3_1D(int n) 
{
  int i,j;

  ////////////////////////////////////////////////
  // La Taille passe est la taille entiere
  // et non plus la taille du support non redondant
  // Cependant le parametre sauvegarde (Private)
  // est la taille du support.
  ////////////////////////////////////////////////

  fenetre = n/2;
 

  v=new T* [fenetre+1];

  for ( i = 0; i <=fenetre; i++ )
    v[i]= new T [ i+1 ];


  for ( i = 0; i <=fenetre; i++ )
    for ( j = 0; j <=i; j++ )
      v[i][j]=((T)(0));
}

template <class T>
TMoment3_1D<T>& TMoment3_1D<T>::operator=(T s)
{
  int i,j; 

  for ( i = 0; i <=fenetre; i++ )
    for ( j = 0; j <=i; j++ )
      v[i][j]=s;

    return *this;

}


// template <class T>
// TMoment3_1D<T>&  TMoment3_1D<T>::resize(int m)
// {


//     TMoment3_1D<T> tmp1(0);
//     ref(tmp1);

//     cerr <<&tmp1<<"  tmp1" << endl;

//     // now, reference an MxN matrix
//     TMoment3_1D<T> tmp(m);
//     ref(tmp);
//     cerr <<&tmp<<"  tmp" << tmp.size()<< " " <<(this)<< endl;

//     return *this;

// }



// template <class T>
// TMoment3_1D<T>& TMoment3_1D<T>::ref(const TMoment3_1D<T>& s)
// {
//   int i,j;
//         // handle trivial M.ref(M) case
//         if (this == &s) return *this;
//         else
//         {
// 	  fenetre=s.fenetre;
	 	
// 	  v=new (T *)[fenetre+1];

// 	  for(i = 0; i <=fenetre; i++ )
// 	    v[i]= new (T) [i+1];

// 	  for(i = 0; i <=fenetre; i++ )
// 	    for ( j = 0; j <=i; j++ )
// 	      v[i][j]=((T)(0));

	  
// 	  return *this;
//         }
// }

template class TMoment3_1D<short>;
template class TMoment3_1D<float>;
template class TMoment3_1D<double>;
template class TMoment3_1D<complex<double> >;
