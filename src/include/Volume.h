//


#ifndef _TVolume_
#define _TVolume_

#include <iostream>       // for formatted printing of volume


template <class T>



class TVolume
{ 
  protected :
    T     ***v;
    int   dimx;  // size of original matrix, not submatrix
    int   dimy;   // size of this submatrix
    int   dimz;

 public:



        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


    TVolume();
    TVolume(int, int, int);
    TVolume(T*, int, int, int );
    virtual ~TVolume();


        /*::::::::::::::::::::::::::::::::*/

        /*  Indices and access operations */

        /*::::::::::::::::::::::::::::::::*/

    inline int size(int d) const;   // submatrix size    
    inline T& operator()(int i, int j, int k);
    inline T& operator()(int i, int j, int k) const;
    inline T* operator()(int i, int j); 

    inline  TVolume<T>& operator=(const TVolume<T>& s); //copy

    TVolume<T>& resize(int m, int n,int p);
    TVolume<T>& resize(const TVolume<T>& s);
    TVolume<T>& ref(const TVolume<T>& s);
    TVolume<T>& copy(const TVolume<T>& s);

};  //* End of TVolume<T> Class *//

   


template <class T>
inline int TVolume<T>::size(int d) const
{
    if(d==0)
      return(this->dimx);
    else if (d==1)
      return(this->dimy);
    else
      return(this->dimz);
}

template <class T>
inline T& TVolume<T>::operator()(int i, int j, int k) const
{
    return v[i][j][k];
}


template <class T>
inline T* TVolume<T>::operator()(int i, int j)
{
    return v[i][j];
}

template <class T>
inline T& TVolume<T>::operator()(int i, int j, int k) 
{
    return v[i][j][k];
}

template <class T>
inline TVolume<T>& TVolume<T>::operator=(const TVolume<T>& s)
{

    return copy(s);
}

#endif 
