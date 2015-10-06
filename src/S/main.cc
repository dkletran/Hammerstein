/*
Simulations
Le nombre de réalisations est defini par NUMTEST. Le résolution (dans le domaine fréquentiel)
de l'estimation est defini par NUMPOINTS (doit être une puissance de 2)
A chaque itéraction:
  + Générer du signal (i.i.d Gaussien)
  + Filtrage Hammerstein
  + Ajouter un bruit blanc, Gaussien
  + Estimer le bispectre et trispectre (méthode de corrélogram)
  + Optimisation (en utilisant LMA) pour estimer les deux noyaux du modèle Hammerstein
  + Corriger la phase
Après tout: calculer le biais et la variance de cette estimation 
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Vecteur.h"
#include "general.h"
#include "fft.h"
#include "LM_algo.h"
#define _USE_MATH_DEFINES
#define pi M_PI
#define NUMPOINTS 	8 //nombre de points à estimer dans l'interval [0, pi]
#define NUMPOINTS1 	(NUMPOINTS+1)  //il faut estimer NUMPOINTS+1 points de 0 à pi (y compris pi)
#define NUMTEST		1000    //nombre de réalisation du signal
#define LENGTHSIGNAL    1024*1024*16
TMoment3_1D<complex<double> > bisp(2*NUMPOINTS);
TMoment4_1D<complex<double> >trisp(2*NUMPOINTS);
/*
 La fonction qui calculer la  function objective et la matrice jacobienne
 ce qui sert à l'agorithme d'optimisation Levenberg-Marquardt (LM_algo)
*/
void ffx(double *fx, TVecteur<double> *jr, double **J, TVecteur<double> &x)
{
 int i, j, k, l, s, i1, j1, k1, l1, s1;
 complex<double> tmp1, tmp2, H1_i, H1_j, H1_k, H1_l, H2_i, H2_j, H2_k, H2_l;
 double *H1, *H2, *jr1, *jr2;
 H1 = x.addr();
 H2 = &H1[2*NUMPOINTS1];
 jr1 = jr->addr();
 jr2 = &jr1[2*NUMPOINTS1];
 for (i = 0; i < 2*NUMPOINTS; i++) jr1[i] =  jr2[i] = 0.0;
 for (*fx = 0, s = 0, i = 0; i < 2*NUMPOINTS; i++, i++)
  for (j = i; j <= 2*(NUMPOINTS-i); j++, j++, s++, s++)
  {
   k = i+j;
   i1 = i+1; j1 = j+1; k1 = k+1; s1 = s+1;
   J[s][i] = J[s][i+2*NUMPOINTS1] = J[s][i1] = J[s][i1+2*NUMPOINTS1]
   = J[s][j] = J[s][j+2*NUMPOINTS1] = J[s][j1] = J[s][j1+2*NUMPOINTS1]
   = J[s][k] = J[s][k+2*NUMPOINTS1] = J[s][k1] = J[s][k1+2*NUMPOINTS1]
   = J[s1][i] = J[s1][i+2*NUMPOINTS1] = J[s1][i1] = J[s1][i1+2*NUMPOINTS1]
   = J[s1][j] = J[s1][j+2*NUMPOINTS1] = J[s1][j1] = J[s1][j1+2*NUMPOINTS1]
   = J[s1][k] = J[s1][k+2*NUMPOINTS1] = J[s1][k1] = J[s1][k1+2*NUMPOINTS1] = 0.0;   
   H1_i = complex<double>(H1[i], H1[i1]);
   H1_j = complex<double>(H1[j], H1[j1]);
   H1_k = complex<double>(H1[k], -H1[k1]);
   H2_i = complex<double>(H2[i], H2[i1]);
   H2_j = complex<double>(H2[j], H2[j1]);
   H2_k = complex<double>(H2[k], -H2[k1]);
   tmp1 = bisp(j/2,i/2) - 2.0*H1_i*H1_j*H2_k
			  - 2.0*H1_i*H2_j*H1_k
			  - 2.0*H2_i*H1_j*H1_k
			  - 8.0*H2_i*H2_j*H2_k;
 
   *fx += abs(tmp1*tmp1);
   tmp2 = 2.0*(H1_j*H2_k + H2_j*H1_k);
   J[s][i] -= tmp2.real();
   J[s][i1] += tmp2.imag();
   J[s1][i] -= tmp2.imag();
   J[s1][i1] -= tmp2.real();
   jr1[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr1[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*H1_j*H1_k+ 8.0*H2_j*H2_k;
   J[s][i+2*NUMPOINTS1] -= tmp2.real();
   J[s][i1+2*NUMPOINTS1] += tmp2.imag();
   J[s1][i+2*NUMPOINTS1] -= tmp2.imag();
   J[s1][i1+2*NUMPOINTS1] -= tmp2.real();   
   jr2[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr2[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*(H1_i*H2_k + H2_i*H1_k);
   J[s][j] -= tmp2.real();
   J[s][j1] += tmp2.imag();
   J[s1][j] -= tmp2.imag();
   J[s1][j1] -= tmp2.real();
   jr1[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr1[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*H1_i*H1_k+ 8.0*H2_i*H2_k;
   J[s][j+2*NUMPOINTS1] -= tmp2.real();
   J[s][j1+2*NUMPOINTS1] += tmp2.imag();
   J[s1][j+2*NUMPOINTS1] -= tmp2.imag();
   J[s1][j1+2*NUMPOINTS1] -= tmp2.real();   
   jr2[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr2[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();
   
   tmp2 = 2.0*(H1_i*H2_j + H2_i*H1_j);
   J[s][k] -= tmp2.real();
   J[s][k1] -= tmp2.imag();
   J[s1][k] -= tmp2.imag();
   J[s1][k1] += tmp2.real();
   jr1[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr1[k1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

   tmp2 = 2.0*H1_i*H1_j+ 8.0*H2_i*H2_j;
   
   J[s][k+2*NUMPOINTS1] -= tmp2.real();
   J[s][k1+2*NUMPOINTS1] -= tmp2.imag();
   J[s1][k+2*NUMPOINTS1] -= tmp2.imag();
   J[s1][k1+2*NUMPOINTS1] += tmp2.real();   
   jr2[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
   jr2[k1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();
  }
 for (i = 0; i < 2*NUMPOINTS; i++, i++)
  for (j = i; j < 2*(NUMPOINTS-i); j++, j++)
    for (k = j; k <= 2*(NUMPOINTS-i-j);k++, k++,s++,s++)
    {
     l = i+j+k;
     i1 = i+1; j1 = j+1; k1 = k+1; l1 = l+1; s1 = s+1;
     J[s][i] = J[s][i+2*NUMPOINTS1] = J[s][i1] = J[s][i1+2*NUMPOINTS1]
     = J[s][j] = J[s][j+2*NUMPOINTS1] = J[s][j1] = J[s][j1+2*NUMPOINTS1]
     = J[s][k] = J[s][k+2*NUMPOINTS1] = J[s][k1] = J[s][k1+2*NUMPOINTS1]
     = J[s][l] = J[s][l+2*NUMPOINTS1] = J[s][l1] = J[s][l1+2*NUMPOINTS1]
     = J[s1][i] = J[s1][i+2*NUMPOINTS1] = J[s1][i1] = J[s1][i1+2*NUMPOINTS1]
     = J[s1][j] = J[s1][j+2*NUMPOINTS1] = J[s1][j1] = J[s1][j1+2*NUMPOINTS1]
     = J[s1][k] = J[s1][k+2*NUMPOINTS1] = J[s1][k1] = J[s1][k1+2*NUMPOINTS1] 
     = J[s1][l] = J[s1][l+2*NUMPOINTS1] = J[s1][l1] = J[s1][l1+2*NUMPOINTS1] = 0.0;   
     H1_i = complex<double>(H1[i], H1[i1]);
     H1_j = complex<double>(H1[j], H1[j1]);
     H1_k = complex<double>(H1[k], H1[k1]);
     H1_l = complex<double>(H1[l], -H1[l1]);
     H2_i = complex<double>(H2[i], H2[i1]);
     H2_j = complex<double>(H2[j], H2[j1]);
     H2_k = complex<double>(H2[k], H2[k1]);     
     H2_l = complex<double>(H2[l], -H2[l1]);
     tmp1 =  trisp(k/2,j/2,i/2)- 8.0*H1_i*H1_j*H2_k*H2_l
	                - 8.0*H2_i*H2_j*H1_k*H1_l
		        - 8.0*H1_i*H2_j*H1_k*H2_l
			- 8.0*H1_i*H2_j*H2_k*H1_l
			- 8.0*H2_i*H1_j*H1_k*H2_l
			- 8.0*H2_i*H1_j*H2_k*H1_l
			- 48.0*H2_i*H2_j*H2_k*H2_l;
 
    *fx += abs(tmp1*tmp1);
     
     tmp2 = 8.0*(H1_j*H2_k*H2_l+H2_j*H1_k*H2_l+H2_j*H2_k*H1_l);
     J[s][i] -= tmp2.real();
     J[s][i1] += tmp2.imag();
     J[s1][i] -= tmp2.imag();
     J[s1][i1] -= tmp2.real();
     jr1[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_j*H1_k*H1_l+H1_j*H1_k*H2_l+H1_j*H2_k*H1_l+6.0*H2_j*H2_k*H2_l);
     J[s][i+2*NUMPOINTS1] -= tmp2.real();
     J[s][i1+2*NUMPOINTS1] += tmp2.imag();
     J[s1][i+2*NUMPOINTS1] -= tmp2.imag();
     J[s1][i1+2*NUMPOINTS1] -= tmp2.real();   
     jr2[i] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[i1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H1_i*H2_k*H2_l+H2_i*H1_k*H2_l+H2_i*H2_k*H1_l);
     J[s][j] -= tmp2.real();
     J[s][j1] += tmp2.imag();
     J[s1][j] -= tmp2.imag();
     J[s1][j1] -= tmp2.real();
     jr1[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_i*H1_k*H1_l+H1_i*H1_k*H2_l+H1_i*H2_k*H1_l+6.0*H2_i*H2_k*H2_l);
     J[s][j+2*NUMPOINTS1] -= tmp2.real();
     J[s][j1+2*NUMPOINTS1] += tmp2.imag();
     J[s1][j+2*NUMPOINTS1] -= tmp2.imag();
     J[s1][j1+2*NUMPOINTS1] -= tmp2.real();   
     jr2[j] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[j1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real(); 

     tmp2 = 8.0*(H1_i*H2_j*H2_l+H2_i*H1_j*H2_l+H2_i*H2_j*H1_l);
     J[s][k] -= tmp2.real();
     J[s][k1] += tmp2.imag();
     J[s1][k] -= tmp2.imag();
     J[s1][k1] -= tmp2.real();
     jr1[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[k1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_i*H1_j*H1_l+H1_i*H1_j*H2_l+H1_i*H2_j*H1_l+6.0*H2_i*H2_j*H2_l);
     J[s][k+2*NUMPOINTS1] -= tmp2.real();
     J[s][k1+2*NUMPOINTS1] += tmp2.imag();
     J[s1][k+2*NUMPOINTS1] -= tmp2.imag();
     J[s1][k1+2*NUMPOINTS1] -= tmp2.real();   
     jr2[k] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[k1] -= tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H1_i*H2_j*H2_k+H2_i*H1_j*H2_k+H2_i*H2_j*H1_k);
     J[s][l] -= tmp2.real();
     J[s][l1] -= tmp2.imag();
     J[s1][l] -= tmp2.imag();
     J[s1][l1] += tmp2.real();
     jr1[l] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr1[l1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();

     tmp2 = 8.0*(H2_i*H1_j*H1_k+H1_i*H1_j*H2_k+H1_i*H2_j*H1_k+6.0*H2_i*H2_j*H2_k);
   
     J[s][l+2*NUMPOINTS1] -= tmp2.real();
     J[s][l1+2*NUMPOINTS1] -= tmp2.imag();
     J[s1][l+2*NUMPOINTS1] -= tmp2.imag();
     J[s1][l1+2*NUMPOINTS1] += tmp2.real();   
     jr2[l] -= tmp2.real()*tmp1.real()+tmp2.imag()*tmp1.imag();
     jr2[l1] += tmp2.real()*tmp1.imag() - tmp2.imag()*tmp1.real();
    }
}
int main(int argc, const char* argv[])
{
 int i, j, k, l, s;
 
 TMoment4_1D<double> Moment4(2*NUMPOINTS), Cumulant4(2*NUMPOINTS);
 TMoment3_1D<double> Moment3(2*NUMPOINTS), Cumulant3(2*NUMPOINTS);
 TVecteur<double> Moment2(2*NUMPOINTS), Cumulant2(2*NUMPOINTS);
 TVecteur<complex<double> > H1(NUMPOINTS1), H1bar(NUMPOINTS1), H2(NUMPOINTS1), H2bar(NUMPOINTS1);
 TVecteur<double> var1(NUMPOINTS1), var2(NUMPOINTS1), bias1(NUMPOINTS1), bias2(NUMPOINTS1);
 complex<double> I = complex<double>(0, 1.0);
 unsigned short int seed16v[3];
 complex<double> tmp;
 double mean, k1, k2;
 seed16v[0]=time(NULL);
 struct drand48_data *tampon=new  struct drand48_data [1];
 seed48_r(seed16v,tampon);
 /*
 Ouvrir un fichier (texte) pour sauvegarder le resultat
 */
 fstream myfile; 
 myfile.open("data.txt",fstream::in | fstream::out | fstream::app); 
 if (myfile.is_open())
  {
  cout<< "data.txt open" <<endl;
  }
 else cout << "Unable to open file";

 
 TVecteur<double> lineaire(3);
 TVecteur<double> quadratique(3);
 TVecteur<double> x(LENGTHSIGNAL);
 TVecteur<double> y(LENGTHSIGNAL);


 /*
 Les noyaux du filtre Hammerstein 
 */
 
 lineaire(0) = 1.0; lineaire(1) = 1.0; lineaire(2) = 1.0;
 quadratique(0) = 1.0; quadratique(1) = 1.0; quadratique(2) = 1.0; 

 for (i = 0; i < NUMPOINTS1; i++) 
 {
  H1(i) = lineaire(0) + lineaire(1)*exp(-I*(i*pi/NUMPOINTS))+ lineaire(2)*exp(-I*(2*i*pi/NUMPOINTS));
  H2(i) = quadratique(0) + quadratique(1)*exp(-I*(i*pi/NUMPOINTS))+ quadratique(2)*exp(-I*(2*i*pi/NUMPOINTS));
 }
 /* 
 Init le LM solver
 */
 
 for (s=0, i = 0; i < NUMPOINTS; i++)
  for (j = i; j <= NUMPOINTS-i; j++, s++);
 for (i = 0; i < NUMPOINTS; i++)
  for (j = i; j < NUMPOINTS-i; j++)
    for (k = j; k <= NUMPOINTS-i-j;k++,s++);

 LM_Minimiser<double> lm(4*NUMPOINTS1, 2*s);
 TVecteur<double> x0(4*NUMPOINTS1, 0);
 
 
 lm.set_f(ffx);


 for (i = 0; i < NUMPOINTS1; i++) 
 {
  x0(2*i) = 1.0;
  x0(2*i+1) = 0.0;
  x0(2*(i+NUMPOINTS1)) = 1.0;
  x0(2*(i+NUMPOINTS1)+1) = 0.0;
 }
 lm.set_x0(x0);

 for(k=0; k<NUMTEST; k++)
 {
   cout<<"k="<<k<<endl;
  /*
  Générer du signal gaussien
 */

 
 GenereVecteurInitialGaussien(x, 1.0, tampon); 
 
 /*
  Filtrage Hammerstein
 */
 
 FiltrageHammerstein(x, y, lineaire, quadratique);
 /*
  Ajouter du bruit (blanc, gaussien)
 */

 GenereVecteurInitialGaussien(x, 1.0, tampon);
 y=y+x;
 

 /*
  Estimer le bispectre et le trispectre
 */

 
 cout<< "Estimation moments ..."<<endl;
 EstimationMoment(&Moment4, &Moment3, &Moment2, &y, 0);
 mean = 0.0;
 for(i = 0; i<y.size(); i++) mean += y(i);
 mean /= y.size();
 cout <<"Estimation cumulants..."<<endl;
 EstimationCumulant(&Cumulant4, &Cumulant3, &Cumulant2, &Moment4, &Moment3, &Moment2, mean);
 cout <<"Estimation polyspectre..."<<endl;
 EstimationSpectre3(&bisp, &Cumulant3);
 EstimationSpectre4(&trisp, &Cumulant4);
 

 
 /*
 Optimisation
 */
 
 cout <<"Solving..."<<endl;
 s = lm.execute(1, 0.000000001, 10000000, 1.05, 30000, &x0);
 cout << "Return code: "<<s<<endl;

 cout <<"Writing solution to file..."<<endl;
 myfile<<k<<endl;
 tmp = complex<double>(x0(2*NUMPOINTS), x0(2*NUMPOINTS+1));
 k1 = arg(tmp);
 tmp = complex<double>(x0(4*NUMPOINTS+2), x0(4*NUMPOINTS+3));
 k2 = arg(tmp);
 for(i = 0; i<NUMPOINTS1; i++)
 {  
 /*Corriger la phase et sauvegarder les résultats*/
  tmp = complex<double>(x0(2*i), x0(2*i+1))*exp(I*(-k1*i/NUMPOINTS));
  x0(2*i) = tmp.real();
  x0(2*i+1) = tmp.imag();
  myfile<<tmp<<"        ";
  tmp = complex<double>(x0(2*(i+NUMPOINTS1)), x0(2*(i+NUMPOINTS1)+1))*exp(I*(-k2*i/NUMPOINTS));
  x0(2*(i+NUMPOINTS1)) = tmp.real();
  x0(2*(i+NUMPOINTS1)+1) = tmp.imag();
  myfile<<tmp<<endl;
 }
 if(k==0)
 {
 /*
 Pour fixer le point initial
 (Ce point devait être obtenu par un algorithme d'optimisation globale,
 mais ici nous utilisons LMA)
 */

  lm.set_x0(x0);
 }
 }
 /*
 Calculer le biais et la variance, en relisant le fichier data.txt
*/
 myfile.seekg(0, ios_base::beg);
 for(k = 0; k<NUMTEST; k++)
 {
  myfile>>s;
  for(i = 0; i<NUMPOINTS1; i++)
  { 
   myfile>>tmp;
   H1bar(i) += tmp/(double)NUMTEST;
   bias1(i) +=abs(tmp-H1(i))/NUMTEST;
   myfile>>tmp;
   H2bar(i) += tmp/(double)NUMTEST;
   bias2(i) +=abs(tmp-H2(i))/NUMTEST;
  }
 }
 myfile.seekg(0, ios_base::beg);
 for(k = 0; k<NUMTEST; k++)
 {
  myfile>>s;
  for(i = 0; i<NUMPOINTS1; i++)
  { 
   myfile>>tmp;
   var1(i) +=pow(abs(tmp-H1bar(i)),2)/NUMTEST;
   myfile>>tmp;
   var2(i) +=pow(abs(tmp-H2bar(i)),2)/NUMTEST;
  }
 }
 /*
 Afficher les resultats
 */
 cout.precision(4);
 cout.setf(ios::fixed,ios::floatfield);
 cout<<"H1, H2 estimés"<<endl;
 for(i = 0; i<NUMPOINTS1; i++)
  cout<<abs(H1bar(i))<<"<"<<arg(H1bar(i))<<"	"<<abs(H2bar(i))<<"<"<<arg(H2bar(i))<<endl;
 cout <<"Le biais:"<<endl;
 for(i = 0; i<NUMPOINTS1; i++)
  cout<<bias1(i)<<"	"<<bias2(i)<<endl;
 cout <<"La variance:"<<endl;
 for(i = 0; i<NUMPOINTS1; i++)
  cout<<var1(i)<<"	"<<var2(i)<<endl;
 myfile.close();

 return 0;
}
