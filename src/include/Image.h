  

#ifndef _TImage_
#define _TImage_

#include <stdlib.h>       // for formatted printing of matrices
#include <iostream> 
using namespace std;
template <class T>



class TImage
{ 
  protected :
    T     **v;
    int   ligne;  // size of original matrix, not submatrix
    int   colonne;   // size of this submatrix
   

 public:



        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


    TImage();
    TImage(int, int);
    TImage(T*, int, int);
    TImage(const TImage<T>&);
    virtual ~TImage();


        /*::::::::::::::::::::::::::::::::*/

        /*  Indices and access operations */

        /*::::::::::::::::::::::::::::::::*/

    inline int size(int d) const;   // submatrix size    
    void centrage();
  
    inline T& operator()(int i, int j);
    inline T& operator()(int i, int j) const;
    inline T* operator[](int i);

    TImage<T> operator/(T reel) const;
    TImage<T> operator*(T reel) const; 
    TImage<T> operator+(T reel) const; 
    TImage<T> operator-(T reel) const;
    TImage<T> operator/(TImage<T> postmatrice) const;
    TImage<T> operator*(TImage<T> postmatrice) const;
    TImage<T> operator+(TImage<T> postmatrice) const; 
    TImage<T> operator-(TImage<T> postmatrice) const;
    TImage<T> operator*=(TImage<T> postmatrice) const;

            TImage<T>& operator=(T s);
    inline  TImage<T>& operator=(const TImage<T>& s); //copy

    TImage<T>& resize(int m, int n);
    TImage<T>& resize(const TImage<T>& s);
    TImage<T>& ref(const TImage<T>& s);
    TImage<T>& inject(const TImage<T>& s);
    TImage<T>& copy(const TImage<T>& s);
    TImage<T>& free();
    //friend TImage<T> operator*(T,TImage<T>);

    

    //* I/O *//
    //friend ostream& operator<<<>(ostream&, const TImage<T>&);
    ostream & Info(ostream & s)
      {
        s << "Size: (" << size(0) << "x" << size(1) << ") " ;
	return s;
      };

};  //* End of TImage<T> Class *//

   


template <class T>
inline int TImage<T>::size(int d) const
{
    if(d==0)
      return(this->ligne);
    else
      return(this->colonne);
}

template <class T>
inline T& TImage<T>::operator()(int i, int j) const
{
   if((i>=ligne) || (j>=colonne))
    {
      cerr <<"erreur dans image" << "   " << this << endl;
      exit(0);
    }

    return v[i][j];
}

template <class T>
inline T* TImage<T>::operator[](int i)
{
    return v[i];
}

template <class T>
inline T& TImage<T>::operator()(int i, int j) 
{

 /*  if((i>=ligne) || (j>=colonne)) */
/*     { */
/*       cerr <<"erreur dans image" << endl; */
/*       exit(0); */
/*     } */
    return v[i][j];
}

template <class T>
inline TImage<T>& TImage<T>::operator=(const TImage<T>& s)
{

    return copy(s);
}

template <class T>
TImage<T> operator*(T reel,TImage<T> matrice)
{
  
  return(matrice*reel);
}

template <class T>
void TImage<T>::centrage() 
{
    
  T somme;

  for (int i=0;i<this->size(0);i++)  
    for (int j=0;j<this->size(1);j++)
      somme+=(*this)(i,j);

  somme/=(((T)(this->size(0)))*((T)(this->size(1))));

 for (int i=0;i<this->size(0);i++)  
    for (int j=0;j<this->size(1);j++)
      (*this)(i,j)-=somme;

}



template <class T>
ostream& operator<<(ostream& s, const TImage<T>& G)
{
  
      
  int i,j;
  
  for (i=0; i<G.size(0); i++)
    {
      for (j=0; j<G.size(1); j++)
	s << G(i,j) << "  ";
      s << "\n";
    }
    
  return s;
}


#endif 

