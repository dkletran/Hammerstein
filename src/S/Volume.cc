
#include <complex>
#include <iostream>
#include "Volume.h"

/* =========================== TMatrice ============================== */

template <class T>
TVolume<T>::~TVolume()
{

  
  for (int i = 0; i < dimx; i++ )
    for (int j = 0; j < dimy; j++ )
      delete[] v[i][j];

  for (int i = 0; i < dimx; i++ )
    delete[] v[i];

    
  delete[] v;

  
  v=NULL;
  dimx = 0;
  dimy = 0;
  dimz = 0;
}

template <class T>
TVolume<T>::TVolume() 
{
  
  dimx = dimy = dimz= 0;
  v=NULL;
}


template <class T>
TVolume<T>::TVolume(int m, int n, int p) 
{
  dimx = m;
  dimy = n;
  dimz = p;
  int i,j,k;

  v=new T** [dimx];

  for (i = 0; i < dimx; i++ )
   v[i]= new T* [ dimy ];

  for (i = 0; i < dimx; i++ )
    for (j = 0; j < dimy; j++ )
   v[i][j]= new T  [ dimz ];

  for (i = 0; i < dimx; i++ )
    for (j = 0; j < dimy; j++ )
      for (k = 0; k < dimz; k++ )
	v[i][j][k]=0;

}


template <class T>
TVolume<T>& TVolume<T>::ref(const TVolume<T>& s)
{
  int i,j;
        // handle trivial M.ref(M) case
        if (this == &s) return *this;
        else
        {
	  dimx=s.dimx;  // size of original matrix, not submatrix
	  dimy=s.dimy;
	  dimz=s.dimz;
	
	  v=new T** [dimx];

	  for (i = 0; i < dimx; i++ )
	    v[i]= new T* [ dimy ];

	  for (i = 0; i < dimx; i++ ) 
	    for (j = 0; j < dimy; j++ )
	      v[i][j]= new T  [ dimz ];
	  
	  return *this;
        }
}


template <class T>
TVolume<T>&  TVolume<T>::resize(int m, int n, int p)
{

    // first, reference 0x0x0 matrix, potentially freeing memory
    // this allows one to resize a matrix > 1/2 of the available
    // memory

    TVolume<T> tmp1(0,0,0);
    ref(tmp1);

    // now, reference an MxNxP matrix
    TVolume<T> tmp(m,n,p);
    ref(tmp);


    return *this;

}

template <class T>
TVolume<T>& TVolume<T>::copy(const TVolume& X) 
{
    
  
  // current scheme in copy() is to detach the left-hand-side
  // from whatever it was pointing to.
  //
  resize(X);
  
  int i,j,k, M = X.size(0),  N = X.size(1), P= X.size(2);
  for (i=0; i<M; i++)
    for (j=0; j<N; j++)
      for (k=0; k<P; k++)
      (*this)(i,j,k) = X(i,j,k);
  


  return *this;
}

template <class T>
TVolume<T>&  TVolume<T>::resize(const TVolume &A)
{
  

   resize(A.size(0), A.size(1), A.size(2));

 }
template class TVolume<bool>;

template class TVolume<short>;
template class TVolume<int>;
template class TVolume<long>;
template class TVolume<float>;
template class TVolume<double>;
//template class TVolume< complex<float> >;
//template class TVolume< complex<double> >;
