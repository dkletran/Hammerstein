


#include <complex>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "Vecteur.h"
#include "Image.h"

#include "Volume.h"
#include "Moment_1D.h"
#include "Moment3_1D.h"
#include "Moment4_1D.h"

using namespace std;

#define  MOD    2147483647L     /* modulus for generator (2**31 - 1) */
#define  MAX_VALUE  (MOD-1)     /* valeurs dans [1, MAX_VALUE] */
#define  MULT        41358L     /* multiplier            */

void EstimationDensiteNoyauConjointe(TVecteur<double> *SerieChronologiqueX,
				     TVecteur<double> *SerieChronologiqueY,
				     TImage<double> *Densite,
				     TVecteur<double> *EchelleX,
				     TVecteur<double> *EchelleY,
				     double MinSerieX,
				     double MaxSerieX,
				     double MinSerieY,
				     double MaxSerieY);




void ResolutionEquation5Ordre(TVecteur <double> CoefficientPolynomes,
			      TVecteur <complex<double> > *Solution);

void ResolutionEquation4Ordre(TVecteur <double> CoefficientPolynomes,
			      TVecteur <complex<double> > *Solution);


void ResolutionEquation3Ordre(TVecteur <double> CoefficientPolynomes,
			      TVecteur <complex<double> > *Solution);



bool OrdreMA(TVecteur<double> *SerieChronologique,double alpha, int *OrdreMA);


bool OrdreMA(TVecteur<double> *SerieChronologique,double alpha,
	     int *OrdreMA, int qmax,int NombreSousSignaux);


void Rotation(TImage<double> *IR,TImage<double> *Image);



void FiltrageSOVM(TVecteur<double> *VecteurSortie,
		  TVecteur<double> * lineaire,
		  TImage<double> *quadratique,
		  double EcartTypeBruit);





void CalculCumulant4SOVM(float VarianceBruit,
			 TVecteur<double>  lineaire,
			 TImage<double>  quadratiqueentree,
			 TMoment4_1D<double> *Moment4);










void GenereVecteurInitialDoubleChi2Centre(TVecteur<double> & VecteurEntree,
					  float ParametreSkewness,struct drand48_data *tampon);

void Genere2DInitialGaussien(TImage<double> *ImageEntreeGaussien,
			     double EcartTypeBruit);


void GenereImageInitialDoubleChi2Centre(TImage<double> *ImageEntree,
					float ParametreSkewness);


void FiltrageBilineaire(TVecteur<double> & VecteurEntree,
			TVecteur<double> & VecteurSortie,
			TVecteur<double> lineaire,
			TVecteur<double> lineairefinal,
			TImage<double> quadratique,
			TImage<double> quadratiquefinal,
			float CoefficientMA);


void FiltrageBilineaire(TVecteur<double> *VecteurSortie,
		   TVecteur<double> *lineaire,
			TImage<double> *quadratique,double VarianceBruit );


void CalculParamatreSkewness(double *ParametreSkewness,
			       double Skewness);

void ResolutionMoindresCarres(TImage<double> *MatriceSysteme, 
			     TVecteur<double> *VecteurSysteme,TVecteur<double> *VecteurCoefficient);


int  Compare(const void *ele1, const void *ele2);

double CalculDeterminant(TImage<double> *Matrice);

void CalculRIF(TVecteur<double> *Correlation,
	       TVecteur<double> *RIF);

void CalculRIF(TVecteur<double> *Correlation,
	       TVecteur<double> *RIF,int ordre);

void CalculRIF(TVecteur<double> *Correlation,
	       TVecteur<double> *Initalisation,
	       TVecteur<double> *RIF);




// Noise.cc
float ran1(long *idum);

float gasdev(long *idum);
float gasdev(struct drand48_data *tampon);

TImage<float> ranGaussien(int N1,int N2);

TVecteur<float> ranGaussien(int N1);

TImage<float> ranRayleigh(int N1,int N2);

TImage<float> ranUniforme(int N1,int N2);


void  ranUniforme (TVecteur<double> *bruit);


TVecteur<double> ranUniforme(int N1);


// GenereSignal.cc
void GenereSignal(char *FichierParametre);

void GenereSignal(int NombreEnsemble,
		  int Nbrealisations,
		  long TailleSignal,
		  float CoefficientMA,
		  char *SNR,
		  float VarianceBruit,
		  char type_model[],
		  char type_entree[],
		  int PAM,
		  char FichierSimulation[],
		  TVecteur<double> lineaire,
		  TImage<double> quadratique,
		  TVecteur<double> lineairefinal,
		  TImage<double> quadratiquefinal);


void GenereSignalBis(char *FichierParametre);

void GenereSignalBis(int NombreEnsemble,
		     int Nbrealisations,
		     long TailleSignal,
		     float CoefficientMA,
		     char *SNR,
		     float VarianceBruit,
		     char type_model[],
		     char type_entree[],
		     int PAM,
		     char FichierSimulation[],
		     TVecteur<double> lineaire,
		     TImage<double> quadratique,
		     TVecteur<double> lineairefinal,
		     TImage<double> quadratiquefinal);







void GenereVecteurInitialPAM(TVecteur<double> & VecteurEntree,int PAM);

int QuantificationPAM(double a,int PAM);


void GenereVecteurInitialGaussien(TVecteur<double> & VecteurEntree,
				  float coef,
				  struct drand48_data *tampon);


void GenereVecteurInitialRayleigh(TVecteur<double> & VecteurEntree,float coef, struct drand48_data *tampon);

void GenereVecteurInitialChi2Centre(TVecteur<double> & VecteurEntree,struct drand48_data *tampon);

// FiltargeSOVM.cc
void FiltrageSOVM(TVecteur<double> *VecteurEntree,
		  TVecteur<double> *VecteurSortie,
		  TVecteur<double> lineaire,
		  TVecteur<double> lineairefinal,
		  TImage<double> quadratique,
		  TImage<double> quadratiquefinal);


void FiltrageHammerstein(TVecteur<double> &VecteurEntree,
			TVecteur<double> &VecteurSortie,
			TVecteur<double> &lineaire,
			TVecteur<double> &quadratique);
			


void FiltrageSOVM(TVecteur<double> *VecteurEntree,
		  TVecteur<double> *VecteurSortie,
		  TVecteur<double> lineaire,
		  TImage<double> quadratique);


void  Bruit(TVecteur<double> & VecteurSortie,float SNR, struct drand48_data *tampon);

void ChargeVecteur(TVolume<double> & CollVecteur,
		   TVecteur<double> VecteurSortie,
		   int NumEnsemble,
		   int insert);

int EcritureSignal(TVolume<double> & matrice,
		   TVolume<double> & matriceentree,
		   char *nom_fic, char *type_model);



//EstimationMoment.cc 
void EstimationMoment(TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset);




void EstimationMoment(TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset,long BorneSup);


void EstimationMoment(TMoment3_1D<double> *Moment3,
		      TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset);

void EstimationMoment(TMoment4_1D<double> *Moment4,
		      TMoment3_1D<double> *Moment3,
		      TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset);





//ChercheCombinaison.cc
void ChercheCombinaison(TVecteur<int> *VecteurBase,
			TImage<int> *Combinaison);

void EstimationCumulant(TMoment4_1D<double> *Cumulant4,
			TMoment3_1D<double> *Cumulant3,
			TVecteur<double> *Cumulant2,
			TMoment4_1D<double> *Moment4,
			TMoment3_1D<double> *Moment3,
			TVecteur<double> *Moment2,
			double Moyenne);


void EstimationSpectre3(TMoment3_1D<complex<double> > *bisp,
			TMoment3_1D<double> *cum
                        );

void EstimationSpectre4(TMoment4_1D<complex<double> > *trisp,
			TMoment4_1D<double> *cum
                        );





