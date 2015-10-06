
#include "general.h"





int Max(int a, int b)
{
 return (a > b)?a:b;
}
int Min(int a, int b)
{
 return (a<b)?a:b;
}

void FiltrageSOVM(TVecteur<double> *VecteurSortie,
		  TVecteur<double> * lineaire,
		  TImage<double> *quadratique,
		  double EcartTypeBruit)
{

  /////////////////////////////////////////////////////////
  // 
  // On génère un peu + de data à acuse de la longueur
  //
  /////////////////////////////////////////////////////////


  int Offset=Max(lineaire->size(), quadratique->size(0))+1;
  TVecteur<double> VecteurEntree(VecteurSortie->size()+Offset);



  ////////////////////////////////////////////////////////////////////////////
  //
  // Initisalisation des variables de bruit
  //
  ////////////////////////////////////////////////////////////////////////////


  unsigned short int seed16v[3];
  seed16v[0]=time(NULL);
  struct drand48_data *tampon=new  struct drand48_data [1];
  double resultat;
  seed48_r(seed16v,tampon);

  GenereVecteurInitialGaussien(VecteurEntree,EcartTypeBruit,tampon);

  FiltrageSOVM(&VecteurEntree,VecteurSortie,*lineaire,*lineaire,
	       *quadratique,*quadratique);

}






void FiltrageSOVM(TVecteur<double> *VecteurEntree,
		  TVecteur<double> *VecteurSortie,
		  TVecteur<double> lineaire,
		  TImage<double> quadratique)
{

  FiltrageSOVM(VecteurEntree,VecteurSortie,lineaire,lineaire,
	       quadratique,quadratique);
		   

}


void FiltrageSOVM(TVecteur<double> *VecteurEntree,
		   TVecteur<double> *VecteurSortie,
		   TVecteur<double> lineaire,
		   TVecteur<double> lineairefinal,
		   TImage<double> quadratique,
		   TImage<double> quadratiquefinal)
{

  int signal;
  
  double somme;

  long TailleSignal=VecteurSortie->size();
  long Offset=VecteurEntree->size()-TailleSignal;



  cerr <<"----------------------"<<endl;
  cerr <<"Generation de la sortie d'un modele de Volterra du second ordre"
       <<endl;
  cerr << "Partie Lineaire : " <<endl;

  for (int i=0;i<lineaire.size();i++)    
    cerr << lineaire(i) << "   ";

  cerr <<endl;

  cerr << "Partie Quadratique : " <<endl;

  for (int i=0;i<quadratique.size(0)-1;i++)
    {
      for (int j=0;j<quadratique.size(1)-1;j++)
	cerr << quadratique(i+1,j+1) << "   ";
  
      cerr <<endl;
    }


  for (int k=0;k<TailleSignal;k++)
    {
      somme=0.;

      double poids2=((double)(k-2))/((double)(TailleSignal-3));
      double poids1=1.0-poids2;

      for (int i=0;i<lineaire.size();i++)     
	somme+=(poids1*lineaire(i)
		+poids2*lineairefinal(i))*(*VecteurEntree)(k-i+Offset);

      //////////////////////////////////////////////////////////
      //
      // Le modele de noyau classique est celui donne
      // par le filtarge QAR qui impose que les coefficients du borde de la
      // matrice sont nuls. Comme ce n'est pas le cas pour les SOVM
      // on commence  à l'indice 1,1, d'ou le decalage ci-dessous
      //
      /////////////////////////////////////////////////////////


      
      for (int i=0;i<quadratique.size(0)-1;i++)
	for (int j=0;j<quadratique.size(1)-1;j++)
	  somme+=(poids1*quadratique(i+1,j+1)
		  +poids2*quadratiquefinal(i+1,j+1))
	    *(*VecteurEntree)(k-i+Offset)*(*VecteurEntree)(k-j+Offset);
      
   
      (*VecteurSortie)(k)=somme;   


      
    }
     
}



void FiltrageHammerstein(TVecteur<double> &VecteurEntree,
			TVecteur<double> &VecteurSortie,
			TVecteur<double> &lineaire,
			TVecteur<double> &quadratique)
{

  /////////////////////////////////////////////////////////////
  //
  // Procedure réécrite pour gagner du tps et éviter les double boucles
  // spécifiques aux modèles de Volterra
  //
  /////////////////////////////////////////////////////////////


  int signal;
  
  double somme;

  long TailleSignal=VecteurSortie.size();

  cerr <<"----------------------"<<endl;
  cerr <<"Generation de la sortie d'un modele de Hammestein du second ordre"
       <<endl;
  cerr << "Partie Lineaire : " <<endl;

  for (int i=0;i<lineaire.size();i++)    
    cerr << lineaire(i) << "   ";

  cerr <<endl;

  cerr << "Partie Quadratique : " <<endl;

  for (int i=0;i<quadratique.size();i++)
    cerr << quadratique(i) << "   ";
  
  cerr <<endl;
  

  int BorneInf=Min(lineaire.size(),
		   quadratique.size());

  for (int k=BorneInf;k<TailleSignal;k++)
    {
      somme=0.;

      for (int i=0;i<lineaire.size();i++)     
	somme+=lineaire(i)*VecteurEntree(k-i);


      //////////////////////////////////////////////////////////
      //
      // Le modele de noyau classique est celui donne
      // par le filtarge QARMA qui impose que les coefficients du borde de la
      // matrice sont nuls. Comme ce n'est pas le cas pour les SOVM
      // on commence  à l'indice 1,1, d'ou le decalage ci-dessous
      //
      /////////////////////////////////////////////////////////
     
      
       for (int i=0;i<quadratique.size();i++)
	somme+=quadratique(i)*VecteurEntree(k-i)*VecteurEntree(k-i);
      
   
      VecteurSortie(k)=somme;      
      
    }
     
  cerr << "Filtrage Hammerstein termine"<<endl;
}
