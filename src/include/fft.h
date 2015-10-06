#include <complex>

using namespace std;
void fft(int n, complex<double> * _out,  complex<double> *w);
void fft2(int n, complex<double> ** _out, complex<double> * _w);
void fft3(int n, complex<double> *** _out, complex<double> *_w);
void compute_w(int n, complex<double> *_w);

