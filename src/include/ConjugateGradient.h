#include "Vecteur.h"
using namespace std;

template <class T>
class MatriceSymetrique
{
 private:
  int n;
  T **v;
 public:
  /*Constructors and Destructors*/
  MatriceSymetrique();
  MatriceSymetrique(int);
  MatriceSymetrique(int, T);
  virtual ~MatriceSymetrique();
  /*Operator*/
  inline T& operator()(int, int);
  inline TVecteur<T> operator*(const TVecteur<T> &) const;
  inline int size() const;
  /*Conjugate with 2 vectors */
  T conjugate(const TVecteur<T> &, const TVecteur<T>&) const;
  /*Conjugate gradient algorithm*/
  int conjugate_gradient(TVecteur<T> &, TVecteur<T> &) const;
};
template <class T>
inline T & MatriceSymetrique<T>::operator()(int i, int j)
{
 if(i>j)
 {
  return v[i][j];
 }
 else 
 {
  return v[j][i];
 }
}
/*Vector and matrices operators*/
template <class T>
inline T operator^(const TVecteur<T> &v1, const TVecteur<T> &v2)
{
 T result;
 result = (T)0;
 for (int i= 0; i< v1.size(); i++) result += v1(i)*v2(i);
 return result;
}
template <class T>
inline TVecteur<T> MatriceSymetrique<T>::operator*(const TVecteur<T> &vv) const
{
 int i,j;
 TVecteur<T> newvecteur(n);
 for(i = 0; i< n; i++)
 {
  for(j = 0; j <=i ; j++)
   newvecteur(i) += v[i][j]*vv(j);
  for(; j < n; j++)
   newvecteur(i) += v[j][i]*vv(j);
 }
 return (newvecteur);
}

