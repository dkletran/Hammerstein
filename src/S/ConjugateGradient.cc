/*
 Conjugate gradient (méthode du gradient conjugué) est 
 un algorithme pour résoudre des systèmes d'équations linéaires dont la matrice 
 est symétrique et définie positive.
  
*/
#include "ConjugateGradient.h"

template <class T>
MatriceSymetrique<T>::MatriceSymetrique()
{
 n = 0;
 v = NULL;
}
template <class T>
MatriceSymetrique<T>::MatriceSymetrique(int nn)
{
 int i;
 n = nn;
 v = new T*[n];
 for(i = 0; i<n; i++)
  v[i] = new T[i+1];
}
template <class T>
MatriceSymetrique<T>::MatriceSymetrique(int nn, T x)
{
 int i, j;
 n = nn;
 v = new T*[n];
 for(i = 0; i<n; i++)
  v[i] = new T[i+1];
 for(i = 0; i<n; i++)
  for(j = 0; j<=i; j++)
   v[i][j] = x;
}
template <class T>
MatriceSymetrique<T>::~MatriceSymetrique()
{
 int i;
 for (i = 0; i<n; i++)
  delete[] v[i];
 delete[] v;
}

template <class T> 
inline int MatriceSymetrique<T>::size() const
{
 return n;
}
template <class T>
T MatriceSymetrique<T>::conjugate(const TVecteur<T> &v1, const TVecteur<T> &v2) const
{
 int i, j;
 T result;
 result = (T)0;
 for(i = 0; i < n; i++)
 {
  for(j = 0; j <= i; j++)
    result += v1[i]*v[i][j]*v2[j];
  for(; j <n; j++)
   result += v1[i]*v[j][i]*v2[j];
 }
 return result;
}
template <class T>
int MatriceSymetrique<T>::conjugate_gradient(TVecteur<T> & b, TVecteur<T> & x) const
/*
 Gradient conjugué algorithme: résoudre le système Ax = b
  A: matrice symétrique (l'objet this)
  b, x: paramètre de la function
*/
{
 int i,j;
 T alpha, rsold, rsnew;
 TVecteur<T> r(n), p(n), Ap(n);

 r = b - (*this)*x;
 p = r;
 rsold = r^r;
 j = 0;
 
 for(i = 0; i<100000000; i++)
 {
  Ap = (*this)*p;
  alpha = rsold/(p^Ap);
  
  x = x + alpha*p;
  r = r - alpha*Ap;
  rsnew = r^r;
  if (abs(rsnew) < 0.000000000001)
  { 
   j =1;
   break;
  }
  p = r + (rsnew/rsold)*p;
  rsold = rsnew;
 
 }
 return j;
}
template class MatriceSymetrique<double>;
template class MatriceSymetrique<complex<double> >;
