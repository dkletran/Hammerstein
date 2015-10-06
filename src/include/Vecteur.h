//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector Float Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_FLOAT_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to float*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _TVecteur_
#define _TVecteur_    

using namespace std;

#include <iostream>      
#include <complex>       
#include <cstdlib>      


template <class T>
                    
class TVecteur
{      

/*   friend class  iostream;                                                                */
 private:
  long flength;
  T *p;
  
 public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
  TVecteur();
  TVecteur(long);  
  TVecteur(int);  
             
  TVecteur(long, T);   
  TVecteur(int, T);   
  
  // statement.
  TVecteur(T*, long);
  TVecteur(T*, int);

  virtual ~TVecteur() ;                              
  
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
  
  inline T&       operator[](long); 
  inline T&       operator[](long) const;  //read only
  inline T&       operator()(long); 
  inline T&       operator()(long) const; // read only

  inline T& operator[](int); 
  inline T& operator[](int) const; 
  inline T&       operator()(int); 
  inline T&       operator()(int) const; // read only

  inline  operator T*();
  inline long          size() const;
  inline long          null() const;
  long          resize(long d);  
  long          resize(int d);

  void centrage();
  inline T* addr() const;
  inline double norm() const;
  void normalise();                                                                    
  /*::::::::::::::*/                                             
  /*  Assignment  */                                             
  /*::::::::::::::*/                                             
  
  inline  TVecteur<T>& operator=(const TVecteur<T>&);
  TVecteur<T>& operator=(T);
  inline  TVecteur<T>& ref(const TVecteur<T> &);
  
      
  TVecteur<T>& inject(TVecteur<T>&);
  TVecteur<T>& copy(const TVecteur<T>&);
  
 
  /* I/O */  
  TVecteur<T> operator/(const T&) const;  
  TVecteur<T> operator*(const T&) const ; 
  TVecteur<T> operator+(const T&) const;
  TVecteur<T> operator-(const T&) const;
  TVecteur<T> operator+(const TVecteur<T>&) const;
  TVecteur<T> operator/(const TVecteur<T>&) const;
  TVecteur<T> operator*(const TVecteur<T>&) const; 
  TVecteur<T> operator-(const TVecteur<T>&) const;
/*    friend iostream& operator<<<> (std::iostream&,const TVecteur<T>&);              */
};         /* friend ostream& operator << <> (ostream&,const TVecteur<T>&);   */       
                                                              


    // operators and member functions

template <class T>
inline long TVecteur<T>::null()  const
{
    return (size() == 0) ;
}

template <class T>
inline long TVecteur<T>::size() const
{
    return   flength;
}


template <class T>
inline T* TVecteur<T>::addr() const
{
    return p;
}

template <class T>
inline TVecteur<T>::operator T*() 
{
    return p;
}

template <class T>
inline T& TVecteur<T>::operator()(long i)
{
 
   
    if(i>=flength)
      {

	cerr<< "Erreur dans un vecteur"<<endl;
	exit(0);
      } 

    return p[i];

}

template <class T>
inline T& TVecteur<T>::operator()(int i)
{
 
   
    if(i>=flength)
      {

	cerr<< "Erreur dans un vecteur"<<endl;
	exit(0);
      } 

    return p[i];

}


template <class T>
inline T& TVecteur<T>::operator()(long i) const
{

   if(i>=flength)
      {

	cerr<< "Erreur dans un vecetur"<<endl;
	exit(0);
      } 
    return p[i];
}


template <class T>
inline T& TVecteur<T>::operator()(int i) const
{

   if(i>=flength)
      {

	cerr<< "Erreur dans un vecetur"<<endl;
	exit(0);
      } 
    return p[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()

template <class T>
inline T& TVecteur<T>::operator[](long i)
{

    return p[i];
}



//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()

template <class T>
inline T& TVecteur<T>::operator[](long i) const
{
    return p[i];
}


template <class T>
inline T& TVecteur<T>::operator[](int i)
{

    return p[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()

template <class T>
inline T& TVecteur<T>::operator[](int i) const
{
    return p[i];
}



template <class T>
inline TVecteur<T>& TVecteur<T>::operator=(const TVecteur<T>& m)
{

    delete[] p;
    flength = m.size();
    p = new T[flength];
    for(int i = 0; i<flength; i++)
      p[i] = m[i];
    return  *this;
}

template <class T>
inline TVecteur<T>& TVecteur<T>::ref(const TVecteur<T>& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorFloat...
    
                   
            
  p = m.p;
  flength = m.flength;
       
  return *this;
}

///////////////////////////////////////////////
// Fonction non membre pour la premultiplication

template <class T>
TVecteur<T> operator*(const T& reel, const TVecteur<T>& vecteur)
{
  return(vecteur*reel);
}


////////////////////
template <class T>
TVecteur<T> operator/(const T& reel,const TVecteur<T>& vecteur)
{
  TVecteur<T> newvecteur(vecteur.size());
  long i;
  for (i=0;i<vecteur.size();i++)
    newvecteur(i)=reel/vecteur(i);

  return(newvecteur);
}

////////////////////
template <class T>
void TVecteur<T>::centrage() 
{
    
  T somme;

  for (long i=0;i<this->size();i++)
    somme+=(*this)(i);

  somme/=((T)(this->size()));

  for (long i=0;i<this->size();i++)
    (*this)(i)-=somme;

}
////////////////////
template <class T>
inline double TVecteur<T>::norm() const
{
    
  double s = 0;

  for (long i=0;i<this->size();i++)
    s+=abs((*this)(i)*(*this)(i));
  return sqrt(s);
}
////////////////////
template <class T>
void TVecteur<T>::normalise() 
{
    
  double s = 0;

  for (long i=0;i<this->size();i++)
    s += abs((*this)(i)*(*this)(i));
  s = sqrt(s);
 
  if (s != 0.0)
   for (long i=0;i<this->size();i++)
     (*this)(i)/=s;

}
////////////////////
//template <class T> 
/* ostream& operator<<(ostream &s, const TVecteur<T> &m) */
/* { */
/*         if (m.p) */
/*         { */
/*                 int n = m.size(); */
/*                 for (int i=0; i<n; i++) */
/* 		  s << m(i) << "  "; */
/*                  s << "\n"; */
/*         } */
/*         else  s << "NULL TVecteur<T>.\n"; */
/*     return  s  ; */
/* } */


#endif 


