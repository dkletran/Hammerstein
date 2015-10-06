//////////////////////
////////////////////
#include "Vecteur.h"
using namespace std;

template <class T>
class SteepestDescent
{
 private:
  int d;
  TVecteur<T> x0;//initial point
  double (*f)(TVecteur<T> &x);//target function
  TVecteur<T> (*g)(TVecteur<T> &x);//function that return the "steepest direction"
 public:
  SteepestDescent();
  SteepestDescent(int dd);
  SteepestDescent(TVecteur<T> &x);
  SteepestDescent(TVecteur<T> &x, double (*ff)(TVecteur<T> &x), TVecteur<T> (*gg)(TVecteur<T> &x));
  ~SteepestDescent();
  void set_x0(TVecteur<T> &x);
  void set_f(double (*ff)(TVecteur<T> &x));
  void set_g(TVecteur<T> (*gg)(TVecteur<T> &x));
  void set_d(int dd);
  int dim();
  int execute(double e, double b, int loop, TVecteur<T> *x);
};
