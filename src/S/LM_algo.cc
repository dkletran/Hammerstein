/*
  algorithme LMA
*/
#include <complex>
#include "Vecteur.h"
#include "LM_algo.h"
#include <iostream>
#include <limits>
using namespace std;


template <class T>
LM_Minimiser<T>::LM_Minimiser()
{

 f = 0;

}

template <class T>
LM_Minimiser<T>::LM_Minimiser(TVecteur<T> &x)
{

 m = d = x.size();
 xx0.resize(d);
 xx0 = x;
 f = 0;
}

template <class T>
LM_Minimiser<T>::LM_Minimiser(int dd)
{

 m = d = dd;
 xx0.resize(dd);

}

template <class T>
LM_Minimiser<T>::LM_Minimiser(int dd, int mm)
{

 m = mm;
 d = dd;
 xx0.resize(dd);

}

template <class T>
LM_Minimiser<T>::LM_Minimiser(TVecteur<T> &x, void (*ff)(double *, TVecteur<T> *, T **, TVecteur<T> &))
{

 d = x.size();
 xx0.resize(d);
 xx0 = x;
 f = ff;
}
template <class T>
LM_Minimiser<T>::~LM_Minimiser()
{


}
//Set functions
template <class T>
void LM_Minimiser<T>::set_d(int dd)
{
 d = dd;
 xx0.resize(dd);
}
template <class T>
void LM_Minimiser<T>::set_x0(TVecteur<T> &x)
{
 xx0 = x;
}
template <class T>
void LM_Minimiser<T>::set_f(void (*ff)(double *, TVecteur<T> *, T **, TVecteur<T> &))
{
 f = ff;
}
template <class T>
void LM_Minimiser<T>::set_m(int mm)
{
 m = mm;
}

//dimension
template <class T>
int LM_Minimiser<T>::dim()
{
 return d;
}


template <class T>
int LM_Minimiser<T>::execute(double e, double b, T lamda0, T sigma, int loops, TVecteur<T> *x)
 /*
 Algorithme de  LMA
  e: Critère d'arrêt pour la fonction objectif (arrêter si f(x) <e) (critère 1)
  b: critère d'arrêt pour le step (arrêter si abs(x(k+1) - x(k))<b) (critère 2)
  lamda, sigma: paramètre pour la mis à jour le parametre lamda
  loops: nombre maximum d'itéraction	(critère 3)
  x: le résultat
 Code de retour: 
  0: retour par critère 1
  1: retour par critère 3
  2: retour par critère 2
  -1: erreur
*/
{
 int i, j, k, l, flag;
 T **J;
 TVecteur<T> x0(d);
 if ((f==0)||(x0.size() == 0)||(m < d))
 {
   return -1;
 }
 J = new T * [m];
 for (i = 0; i < m; i++)
  J[i] = new T[d];
 for (i = 0; i < m; i++)
  for (j = 0; j < d; j++)
   J[i][j] = (T)0;
 TVecteur<T> jr(d), delta(d);
 
 MatriceSymetrique<T > A(d);
 double fx0new, fx0old;
 T lamda = lamda0;
 int dir = 0;
 flag = 1;
 *x = xx0;
 fx0old = numeric_limits<double>::max();
 for(i = 0; i < loops; i++)
 {
   f(&fx0new, &jr, J, *x);
   cout << "Iteractions: " <<i<< "   Objective function:f(x) ="<< fx0new<<endl;
   if (fx0new <e)
   {
      flag = 0;
      break;
   }      
   if (fx0new <= fx0old)
   {  
    lamda /= sigma;
    x0 = *x;
   }
   else
   {
    lamda =lamda0;
    *x = x0; 
   }
   fx0old = fx0new;
   for (j = 0; j < d; j++)
   { 
    for (k = 0; k <= j; k++)
    {
     A(j, k) = (T)0;
     for (l = 0; l < m; l++)
       A(j, k) += J[l][j]*J[l][k];
    }
    A(j, j) +=  lamda;
   }
   A.conjugate_gradient(jr, delta);
   *x = *x - delta; 
   
   if (delta.norm() < b)
   {
    flag = 2;
     break;
   }

  }
 
 for (i = 0; i < m; i++)
  delete[] J[i];
 delete[] J;
 
 return flag;
}
template class LM_Minimiser<double>;
template class LM_Minimiser<complex<double> >;

