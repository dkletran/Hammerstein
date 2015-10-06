


#include <complex>
#include <iostream>
#include "Vecteur.h"



                                                                       
template <class T>
TVecteur<T>::~TVecteur()
{

  // cerr << this << "  debut des"<<endl;
  delete[] p;
  p=NULL;
  flength=0;
//       cerr << this << "  fin des"<<endl;
   
}




template <class T>
TVecteur<T>::TVecteur(long n)
{                                                                      
  p = new T [n];
  flength= n; 

 for (long i=0; i<n; i++)                                            
        p[i] = ((T)(0)); 
  
                                                          
  
}   

template <class T>
TVecteur<T>::TVecteur(int n)
{                                                                      
  p = new T [n];
  flength= n; 

 for (long i=0; i<n; i++)                                            
        p[i] = ((T)(0)); 
  
                                                          
  
}                                                                 

template <class T>
TVecteur<T>::TVecteur()
{                                                                      
  p = NULL;
  flength= 0;                                                          
  
}     

template <class T>
 TVecteur<T>::TVecteur(T *d, long n)
{                                                                      
    p =d;
    flength= n;                                                       
}                                                                      
  

template <class T>
 TVecteur<T>::TVecteur(T *d, int n)
{                                                                      
    p =d;
    flength= n;                                                       
}                                                                      
 
                                                                     


template <class T>
TVecteur<T>::TVecteur(long n, T scalar)
{                                                                      
    flength= n;
    p = new T [n];    
    for (long i=0; i<n; i++)                                            
        p[i] = scalar;                                            
}                                                                      
        

template <class T>
TVecteur<T>::TVecteur(int n, T scalar)
{                                                                      
    flength= n;
    p = new T [n];    
    for (long i=0; i<n; i++)                                            
        p[i] = scalar;                                            
} 
                                                               


// this actually frees memory first, then resizes it.  it reduces
// internal fragmentation of memory pool, and the resizing of
// matrices > 1/2 available memory.
template <class T>
long TVecteur<T>::resize(long d)
{

    if (d<0)                // do nothing if invalid size
    {
        return size();
    }
    else
    {
        delete[] p;
        flength = d;
        if (d>0)
            p = new T[d];
        else
            p = NULL;
    }
    return size();
}

template <class T>
long TVecteur<T>::resize(int d)
{

    if (d<0)                // do nothing if invalid size
    {
        return size();
    }
    else
    {
        delete[] p;
        flength = d;
        if (d>0)
            p = new T[d];
        else
            p = NULL;
    }
    return size();
}

////////////////////////////
template <class T>
TVecteur<T>& TVecteur<T>::operator=(T s)
{
  long i,j;
  for (long i=0;i<flength;i++)
    p[i]=s;

}



template <class T>
TVecteur<T>& TVecteur<T>::inject( TVecteur<T>& m)
{
    if (m.size() != size())
    {
          cerr << "TVecteur<T>::inject(): vector sizes do not match.\n";
      return *this;
    }
    long N = size();
    for (long i=0; i<N; i++)
        (*this)(i) = m(i);

    return *this;
}

template <class T>
TVecteur<T> TVecteur<T>::operator/(const T& reel) const
{
    
  TVecteur<T> newvecteur(this->size());

  for (long i=0;i<this->size();i++)
    newvecteur(i)=this->p[i]/reel;
   
  return (newvecteur);
}

template <class T>
TVecteur<T> TVecteur<T>::operator * (const T& reel) const
{
    
  TVecteur<T> newvecteur(this->size());

  for (long i=0;i<this->size();i++)
    newvecteur(i)=this->p[i]*reel;
   
  return (newvecteur);
}

template <class T>
TVecteur<T> TVecteur<T>::operator+(const T& reel) const
{
    
  TVecteur<T> newvecteur(this->size());

  for (long i=0;i<this->size();i++)
    newvecteur(i)=this->p[i]+reel;
  
  return (newvecteur);
}

template <class T>
TVecteur<T> TVecteur<T>::operator+(const TVecteur<T>& vect) const
{
    
  TVecteur<T> newvecteur(this->size());
  
  if (this->size() != vect.size())
    {
      cerr << "Dimension non conforme pour l'additon des vecteurs \n";
      newvecteur=(*this);
    }
  else 
    for (long i=0;i<this->size();i++)
      newvecteur(i)=this->p[i]+vect(i);
   
  
  return (newvecteur);
}


template <class T>
TVecteur<T> TVecteur<T>::operator-(const TVecteur<T>& vect) const
{
    
  TVecteur<T> newvecteur(this->size());
  
  if (this->size() != vect.size())
    {
      cerr << "Dimension non conforme pour la soustraction des vecteurs \n";
      newvecteur=(*this);
    }
  else 
    for (long i=0;i<this->size();i++)
      newvecteur(i)=this->p[i]-vect(i);
   
  
  return (newvecteur);
}

template <class T>
TVecteur<T> TVecteur<T>::operator/(const TVecteur<T>& vect) const
{
    
  TVecteur<T> newvecteur(this->size());
  
  if (this->size() != vect.size())
    {
      cerr << "Dimension non conforme pour la division des vecteurs \n";
      newvecteur=(*this);
    }
  else 
    for (long i=0;i<this->size();i++)
      newvecteur(i)=this->p[i]/vect(i);
   
  
  return (newvecteur);
}

template <class T>
TVecteur<T> TVecteur<T>::operator*(const TVecteur<T>& vect) const
{
    
  TVecteur<T> newvecteur(this->size());
  
  if (this->size() != vect.size())
    {
      cerr << "Dimension non conforme pour la division des vecteurs \n";
      newvecteur=(*this);
    }
  else 
    for (long i=0;i<this->size();i++)
      newvecteur(i)=this->p[i]*vect(i);
   
  
  return (newvecteur);
}

template <class T>
TVecteur<T> TVecteur<T>::operator-(const T& reel) const
{
    
  TVecteur<T> newvecteur(this->size());

  for (long i=0;i<this->size();i++)
    newvecteur(i)=this->p[i]-reel;
  
  return (newvecteur);
}



template <class T>
TVecteur<T>& TVecteur<T>::copy(const TVecteur<T> &m)
{

        if (null()) resize(m.size());

        if (size() != m.size())
           cerr << "TVecteur<T>::copy(TVecteur<T> &): incompatible vector \
                sizes : "<< size() << " vs. " << m.size() << ".\n";
        else

    resize(0);                  // free up destination

    long N = m.size();
    TVecteur<T> tmp(N);

    for (long i=0; i<N; i++)     // should use memcpy() here...
       (*this)(i) = m(i);

    ref(tmp);
    return *this;
}




template class TVecteur<bool>;
template class TVecteur<short>;
template class TVecteur<int>;
template class TVecteur<long>;
template class TVecteur<float>;
template class TVecteur<double>;
template class TVecteur<complex<float> >;
template class TVecteur<complex<double> >;
