#include "Vecteur.h"
#include "ConjugateGradient.h"
using namespace std;

template <class T>
class LM_Minimiser
{
 private:
  int d;
  int m;
  TVecteur<T> xx0;//initial point
  void (*f)(double *, TVecteur<T> *, T **, TVecteur<T> &);//target function, jacobian...

 public:
  LM_Minimiser();
  LM_Minimiser(int dd);
  LM_Minimiser(int dd, int mm);
  LM_Minimiser(TVecteur<T> &x);
  LM_Minimiser(TVecteur<T> &x, void (*ff)(double *, TVecteur<T> *, T **, TVecteur<T> &));
  ~LM_Minimiser();
  void set_x0(TVecteur<T> &x);
  void set_f(void (*ff)(double *, TVecteur<T> *, T **, TVecteur<T> &));
  void set_d(int dd);
  void set_m(int mm);
  int dim();
  int execute(double e, double b, T lamda0, T sigma, int loop, TVecteur<T> *x);
};
