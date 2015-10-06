

#include <complex>
#include <iostream>
#include "Image.h"

/* =========================== TMatrice ============================== */

template <class T>
TImage<T>::~TImage()
{

//   cerr << this << "deb des im "<<endl;
  for (int i = 0; i < ligne; i++ )
    delete[] v[i];
    

    
  delete[] v;

  
  v=NULL;
  ligne = 0;
  colonne = 0;

//   cerr << this << "fin des im "<<endl;


}

template <class T>
TImage<T>::TImage() 
{
  
  ligne = colonne = 0;
  v=NULL;
}


template <class T>
TImage<T>::TImage(int m, int n) 
{

//   cerr << this << endl;
  int i,j;
  ligne = m;
  colonne = n;

  v=new T* [ligne];

  for ( i = 0; i < ligne; i++ )
      v[i]= new T [ colonne ];

  for ( i = 0; i < ligne; i++ )
    for ( j = 0; j < colonne; j++ )
      v[i][j]=((T)(0));
}


template <class T>
TImage<T>::TImage(T *d, int m, int n) 
{
  ligne = m;
  colonne = n;
  int i,j;

  v=new T* [ligne];
  for (i = 0; i < ligne; i++ )
    v[i]= new T [ colonne ];

  for (i = 0; i < ligne; i++ )
    for (j = 0; j < colonne; j++ )
      v[i][j]=(*d);
  
}



////////////////////////////
template <class T>
TImage<T> TImage<T>::operator/(T reel) const
{

  
  
  TImage<T> newmatrice(this->size(0),this->size(1));
  
  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)/reel;
      
  
  return newmatrice;
}

////////////////////////////
template <class T>
TImage<T> TImage<T>::operator*(T reel) const
{
  
  TImage<T> newmatrice(this->size(0),this->size(1));

  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)*reel;
   
  return newmatrice;
}

////////////////////////////
template <class T>
TImage<T>& TImage<T>::operator=(T s)
{
  int i,j;
  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      (*this)(i,j)=s;
  return *this;

}


////////////////////////////
template <class T>
TImage<T> TImage<T>::operator+(T reel) const
{
    
  TImage<T> newmatrice(this->size(0),this->size(1));

  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)+reel;
   
  return(newmatrice);
}

////////////////////////////
template <class T>
TImage<T> TImage<T>::operator-(T reel) const
{
    
  TImage<T> newmatrice(this->size(0),this->size(1));

  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)-reel;
   
  return(newmatrice);
}

////////////////////////////
template <class T>
TImage<T> TImage<T>::operator*(TImage postmatrice) const
{
    

  if ((this->size(0) != postmatrice.size(0)) ||
      (this->size(1) != postmatrice.size(1)))
  {

    cerr <<
      "Erreur dans les dimensions de matrices pour la multiplication point a points \n";
    exit(1);
  }


  TImage newmatrice(this->size(0),this->size(1));

  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)*postmatrice(i,j);
   
  return (newmatrice);
}

////////////////////////////
template <class T>
TImage<T> TImage<T>::operator/(TImage postmatrice) const
{
    

  if ((this->size(0) != postmatrice.size(0)) ||
      (this->size(1) != postmatrice.size(1)))
  {

    cerr <<
      "Erreur dans les dimensions de matrices pour la division point a points \n";
    exit(1);
  }


  TImage newmatrice(this->size(0),this->size(1));

  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)/postmatrice(i,j);
   
  return (newmatrice);
}

template <class T>
TImage<T> TImage<T>::operator+(TImage postmatrice) const
{
    

  if ((this->size(0) != postmatrice.size(0)) ||
      (this->size(1) != postmatrice.size(1)))
  {

    cerr <<
      "Erreur dans les dimensions de matrices pour l'addition \n";
    exit(1);
  }


  TImage newmatrice(this->size(0),this->size(1));
  
  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)+postmatrice(i,j);
   
  return (newmatrice);
}

template <class T>
TImage<T> TImage<T>::operator-(TImage postmatrice) const
{
    

  if ((this->size(0) != postmatrice.size(0)) ||
      (this->size(1) != postmatrice.size(1)))
  {

    cerr <<
      "Erreur dans les dimensions de matrices pour la soustraction \n";
    exit(1);
  }


  TImage<T> newmatrice(this->size(0),this->size(1));

  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      newmatrice(i,j)=(*this)(i,j)-postmatrice(i,j);
   
  return (newmatrice);
}


template <class T>
TImage<T>&  TImage<T>::resize(int m, int n)
{

    // first, reference 0x0 matrix, potentially freeing memory
    // this allows one to resize a matrix > 1/2 of the available
    // memory

    TImage<T> tmp1(0,0);
    ref(tmp1);

    
    // now, reference an MxN matrix
    TImage<T> tmp(m,n);
    ref(tmp);


    return *this;

}

template <class T>
TImage<T> TImage<T>::operator*=(TImage<T> postmatrice) const
{
  if ((this->size(0) != postmatrice.size(0)) ||
      (this->size(1) != postmatrice.size(1)))
  {
    cerr <<
      "Erreur dans les dimensions de matrices pour la multiplication point a points \n";
    exit(1);
  }


  for (int i=0;i<this->size(0);i++)
    for (int j=0;j<this->size(1);j++)
      (*this)(i,j)=(*this)(i,j)*postmatrice(i,j);
   
  
}

template <class T>
TImage<T>&  TImage<T>::resize(const TImage &A)
{
  

   resize(A.size(0), A.size(1));

 }

//
// Constructeur par copy
template <class T>
TImage<T>::TImage(const TImage<T>& X) 
{
  int i,j;
  
  
  ligne = X.size(0);
  colonne = X.size(1);
  v=new T* [ligne];

  
  for ( i = 0; i < ligne; i++ )
    v[i]= new T [ colonne ];


  for (i = 0; i < ligne; i++ )
    for (j = 0; j < colonne; j++ )
      v[i][j] = X.v[i][j];

 
}




template <class T>
TImage<T>& TImage<T>::copy(const TImage& X) 
{
    
  
  // current scheme in copy() is to detach the left-hand-side
  // from whatever it was pointing to.
  //
  resize(X);
  
  int i,j, M = X.size(0),  N = X.size(1);
  for (i=0; i<M; i++)
    for (j=0; j<N; j++)
      (*this)(i,j) = X(i,j);
  


  return *this;
}

template <class T>
TImage<T>& TImage<T>::inject(const TImage& s)
{
 
  int i, j,  M=size(0), N=size(1);
  
  for (j=0;j<N;j++)
    for (i=0;i<M; i++)
      //   operator()(i,j) = s(i,j);
  
  return *this;
}




template <class T>
TImage<T>& TImage<T>::ref(const TImage<T>& s)
{
  int i,j;
        // handle trivial M.ref(M) case
        if (this == &s) return *this;
        else
        {
	  ligne=s.ligne;  // size of original matrix, not submatrix
	  colonne=s.colonne;
	
	  v=new T* [ligne];

	  for (i = 0; i < ligne; i++ )
	    v[i]= new T [ colonne ];
	  
	  return *this;
        }
}

template <class T>
TImage<T>&  TImage<T>::free()
{

  for (int i = 0; i < this->size(0); i++ )
    delete[] (this->v)[i];
    
  delete[] this->v;

  
  this->v=NULL;
  this->ligne = 0;
  this->colonne = 0;
  
}

template class TImage<unsigned char>;
template class TImage<char>;
template class TImage<short>;
template class TImage<int>;
template class TImage<long>;
template class TImage<float>;
template class TImage<double>;
//template class TImage<complex<float> >;
//template class TImage<complex<double> >;
