/*
 Algorithme steepest_descent
*/
#include <complex>
#include "Vecteur.h"
#include "steepest_descent.h"

using namespace std;

//Constructor & Destructor functions
template <class T>
SteepestDescent<T>::SteepestDescent()
{
// cout << "construction an instance of steepest descent"<<endl;
 f = 0;
 g = 0;
}

template <class T>
SteepestDescent<T>::SteepestDescent(TVecteur<T> &x)
{
// cout << "construction an instance of steepest descent"<<endl;
 d = x.size();
 x0.resize(d);
 x0 = x;
 f = 0;
 g = 0;
}

template <class T>
SteepestDescent<T>::SteepestDescent(int dd)
{
// cout << "construction an instance of steepest descent"<<endl;
 d = dd;
 x0.resize(dd);

}

template <class T>
SteepestDescent<T>::SteepestDescent(TVecteur<T> &x, double (*ff)(TVecteur<T>&) ,TVecteur<T> (*gg)(TVecteur<T>&))
{
// cout << "construction an instance of steepest descent"<<endl;
 d = x.size();
 x0.resize(d);
 x0 = x;
 f = ff;
 g = gg;
}
template <class T>
SteepestDescent<T>::~SteepestDescent()
{
// cout << "destruction an instance of steepest descent"<<endl;
}
//Set functions
template <class T>
void SteepestDescent<T>::set_d(int dd)
{
 d = dd;
 x0.resize(dd);
}
template <class T>
void SteepestDescent<T>::set_x0(TVecteur<T> &x)
{
 x0 = x;
}
template <class T>
void SteepestDescent<T>::set_f(double (*ff)(TVecteur<T>&))
{
 f = ff;
}
template <class T>
void SteepestDescent<T>::set_g(TVecteur<T> (*gg)(TVecteur<T>&))
{
 g = gg;
}
//Get dimension
template <class T>
int SteepestDescent<T>::dim()
{
 return d;
}

//Algorithm steepest descent
template <class T>
int SteepestDescent<T>::execute(double e, double b, int loops, TVecteur<T> *x)
 /*
 Algorithme steepest descent (gradient descent)
  e: Critère d'arrêt pour la fonction objectif (arrêter si f(x) <e) (critère 1)
  b: critère d'arrêt pour le step (arrêter si abs(x(k+1) - x(k))<b) (critère 2)
  loops: nombre maximum d'itéraction	(critère 3)
  x: le résultat
 Code de retour: 
  0: retour par critère 1
  1: retour par critère 3
  2: retour par critère 2
  -1: erreur
*/
{
 if ((f==0)||(g==0)||(x0.size() == 0))
 {
   return -1;
 }
 
 TVecteur<T> gr(d);
 TVecteur<T> xx0(d);
 TVecteur<T> xx1(d);
 TVecteur<T> xx2(d); 
 double fx0, fx1, fx2;
 double hmin;
 double h0, h1, h2;
 double h = 1;
 int i, j;
 int flag;
 xx0 = x0;
 *x = x0;
 for(i = 0; i < loops; i++)
 {
   fx0 = f(xx0);
   cout<<"Iteractions: "<<i<<"   Objective function:f(x) ="<< fx0<<"  Step: h = "<<h<<endl;
   if (fx0 <e)
   {
     *x = xx0; 
      return 0;
   }      
   gr = g(xx0);
   gr.normalise();
   xx1 = xx0 - ((T)h)*gr;
   xx2 = xx0 - ((T)2*h)*gr;
   fx1 = f(xx1);
   fx2 = f(xx2);
   j = 0;
   flag = 1;
   while(flag)
   {
      if (fx0 <= fx1)
      {
        xx2 = xx1;
        fx2 = fx1;
        h = h/2;
        xx1 = xx0 -((T)h)*gr;
        fx1 = f(xx1);
      }
      else if (fx2 <= fx1)
      {
       xx1 = xx2;
       fx1 = fx2;
       h = 2*h;
       xx2 = xx0 -((T)2*h)*gr;
       fx2 = f(xx2);
      } else
      {
       flag = 0;
      }
      j = j + 1;
   }
   hmin = h/2*(4*fx1 - 3*fx0 -fx2)/(2*fx1 - fx0 - fx2); 
   xx0 =  xx0 -(T)hmin*gr;
   h0 = abs(hmin);
   h1 = abs(hmin -h);
   h2 = abs(hmin -2*h);
   if(h0 < h) h = h0;
   if(h1 < h) h = h1;
   if(h2 < h) h = h2;
   if (h == 0) h = hmin;
   if (h < b)
   {
    *x = xx0;
     return 2;
   }
  }
  *x = xx0;
 return 1;
}
template class SteepestDescent<double>;
template class SteepestDescent<complex<double> >;

