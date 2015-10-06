
#include "general.h"






void EstimationMoment(TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset)
{

  long FenetreMoment=Moment2->size()/2;

  EstimationMoment(Moment2,VecteurCourant,Offset,FenetreMoment);

}



void EstimationMoment(TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset,long BorneSup)

{ 


  int i,j,k,n;
  
  double cour,cour1;
   
  /////////////////////////////////////////////////////
  //
  // Cette procedure prend un vecteur vect_dep de taille
  // NbrEchantillon et en calcul les moments d ordre 3 
  // et 4 avec un algorithme par plan et ceci pour des indices
  // compris entre dep_fen et fin fen. 
  // le moment d ordre mom3(i,j) est calcule tel que j<=i et
  // mom4(i,j,k) tel que k<=j<=i ce qui constitue les
  // support non redondant.
  //
  /////////////////////////////////////////////////////


  long FenetreMoment=BorneSup;
  long NbrEchantillon=VecteurCourant->size();

  // Allocation memoire tampon pour le calcul des moments 
   
  double vec_2=0.;


  for(i=Offset;i<=FenetreMoment;i++)
    {
      // Description des indices pour lequel les
      // moments doivent etre calcules

    
      for(n=i;n<NbrEchantillon;n++)
	{  

	  // Description de tous les indices qui inetrviennent
	  // dans l'estimation des moments 
	  // calcul du produit x(n)*x(n-i) 

	  cour=(*VecteurCourant)(n)*(*VecteurCourant)(n-i);
	  vec_2+=cour;
	}

      ////////////////////////////////////////////
      //
      // lorsque tous les indices ont ete parcourus
      // le resulats est normalise 
      // puis ranges dans les structures charge de recevoir le resultat 
      // et les buffer nettoyes 
      //
      ///////////////////////////////////////////
      
      (*Moment2)(i)=vec_2/((double)(NbrEchantillon-i));

      if(i !=0)
	  (*Moment2)(Moment2->size()-i)=(*Moment2)(i);


      vec_2=0.0;

    }
   
   
  // liberation de la memoire allouee localement 
     
}  



void EstimationMoment(TMoment3_1D<double> *Moment3,
		      TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset)

{ 

 
  long i,j,k,n;
  
  double cour,cour1;
   
  /////////////////////////////////////////////////////
  //
  // Cette procedure prend un vecteur vect_dep de taille
  // NbrEchantillon et en calcul les moments d ordre 3 
  // et 4 avec un algorithme par plan et ceci pour des indices
  // compris entre dep_fen et fin fen. 
  // le moment d ordre mom3(i,j) est calcule tel que j<=i et
  // mom4(i,j,k) tel que k<=j<=i ce qui constitue les
  // support non redondant.
  //
  /////////////////////////////////////////////////////


  long FenetreMoment=(Moment3->size())/2;
  long NbrEchantillon=VecteurCourant->size();


  ////////////////////////////////////////
  //
  // Allocation memoire tampon pour le calcul des moments 
  //
  /////////////////////////////////////////


  double vec_2=0;
  double * vec_3=new double [FenetreMoment+1];

  for (i=0;i<=FenetreMoment;i++)
    vec_3[i]=0.0;

  // Fin d'allocation de memoire locale

  for(i=Offset;i<=FenetreMoment;i++)
    {
      // Description des indices pour lequel les
      // moments doivent etre calcules

      for(n=i;n<NbrEchantillon;n++)
	{  

	  // Description de tous les indices qui inetrviennent
	  // dans l'estimation des moments 
	  // calcul du produit x(n)*x(n-i) 

	  cour=(*VecteurCourant)(n)*(*VecteurCourant)(n-i);
	  for(j=0;j<=i;j++)
	    {
	      //puis calcul du produit  x(n)*x(n-i)*x(n-i-j)
	      //pour le moemnt d'ordre 3

	      cour1=cour*(*VecteurCourant)(n-i+j);
	      
	      vec_3[j]+=cour1;
	    }
	  vec_2+=cour;
	}

      ////////////////////////////////////////////
      //
      // lorsque tous les indices ont ete parcourus
      // le resultat est normalisé
      // puis rangé dans les structures chargées
      // de recevoir le resultat 
      // et les buffers sont nettoyés 
      //
      ///////////////////////////////////////////
      
      (*Moment2)(i)=vec_2/((double)(NbrEchantillon-i));

      if(i !=0)
	  (*Moment2)(2*FenetreMoment-i)=(*Moment2)(i);


      vec_2=0.0;

      for(j=0;j<=i;j++)
	{
	  (*Moment3)(i,j)=vec_3[j]/((double)(NbrEchantillon-i));
	  vec_3[j]=0.0;
	  
	}
    }
   

    
  // liberation de la memoire allouee localement 
  delete [] vec_3;
   
}  
/*
void EstimationMoment(TMoment4_1D<double> *Moment4,
		      TMoment3_1D<double> *Moment3,
		      TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset)

{ 


  long i,j,k,n;
  
  double cour,cour1;
   
  /////////////////////////////////////////////////////
  //
  // Cette procedure prend un vecteur vect_dep de taille
  // NbrEchantillon et en calcul les moments d ordre 3 
  // et 4 avec un algorithme par plan et ceci pour des indices
  // compris entre dep_fen et fin fen. 
  // le moment d ordre mom3(i,j) est calcule tel que j<=i et
  // mom4(i,j,k) tel que k<=j<=i ce qui constitue les
  // support non redondant.
  //
  /////////////////////////////////////////////////////


  long FenetreMoment=(Moment4->size())/2;
  long NbrEchantillon=VecteurCourant->size();

  // Allocation memoire tampon pour le calcul des moments 
   
  double vec_2=0;
  double * vec_3=new double [FenetreMoment+1];
  double ** vec_4=new double* [FenetreMoment+1];

  for (i=0;i<=FenetreMoment;i++)
    vec_4[i]=new double [FenetreMoment+1];

  /// Fin d'allocation de memoire locale

  for (i=0;i<=FenetreMoment;i++)
    {
      vec_3[i]=0.0;
      for(j=0;j<=FenetreMoment;j++)
	vec_4[i][j]=0.0;
    }


  for(i=Offset;i<=FenetreMoment;i++)
    {
      // Description des indices pour lequel les
      // moments doivent etre calcules

      for(n=i;n<NbrEchantillon;n++)
	{  

	  // Description de tous les indices qui inetrviennent
	  // dans l'estimation des moments 
	  // calcul du produit x(n)*x(n-i) 

	  cour=(*VecteurCourant)(n)*(*VecteurCourant)(n-i);
	  for(j=0;j<=i;j++)
	    {
	      //puis calcul du produit  x(n)*x(n-i)*x(n-i-j)
	      //pour le moemnt d'ordre 3

	      cour1=cour*(*VecteurCourant)(n-i+j);
	      for(k=0;k<=j;k++)
		{
		  // puis calcul du produit  x(n)*x(n-i)*x(n-i-j)*x(n-i-k)
		  // pour le moment d'ordre 3

		  vec_4[j][k]=vec_4[j][k]+cour1*(*VecteurCourant)(n-i+k);
		}
	      vec_3[j]=vec_3[j]+cour1;
	    }
	  vec_2+=cour;
	}

      ////////////////////////////////////////////
      //
      // lorsque tous les indices ont ete parcourus
      // le resulats est normalise 
      // puis ranges dans les structures charge de recevoir le resultat 
      // et les buffer nettoyes 
      //
      ///////////////////////////////////////////
      
      (*Moment2)(i)=vec_2/((double)(NbrEchantillon-i));

      if(i !=0)
	  (*Moment2)(2*FenetreMoment-i)=(*Moment2)(i);


      vec_2=0.0;

      for(j=0;j<=i;j++)
	{
	  (*Moment3)(i,j)=vec_3[j]/((double)(NbrEchantillon-i));
	  vec_3[j]=0.0;
	  for(k=0;k<=j;k++)
	    {
	      (*Moment4)(i,j,k)=vec_4[j][k]/((double)(NbrEchantillon-i));
	      vec_4[j][k]=0.0;
	    }
	}
    }
   
   
  // liberation de la memoire allouee localement 
  delete [] vec_3;
   
  for(i=0;i<=FenetreMoment;i++)
    delete [] vec_4[i];
   
  delete [] vec_4;
}  

*/


void EstimationMoment(TMoment4_1D<double> *Moment4,
		      TMoment3_1D<double> *Moment3,
		      TVecteur<double> *Moment2,
		      TVecteur<double> *VecteurCourant,
		      long Offset)

{ 


  long i,j,k,n;
  
  double cour,cour1;
   
  /////////////////////////////////////////////////////
  //
  // Cette procedure prend un vecteur vect_dep de taille
  // NbrEchantillon et en calcul les moments d ordre 3 
  // et 4 avec un algorithme par plan et ceci pour des indices
  // compris entre dep_fen et fin fen. 
  // le moment d ordre mom3(i,j) est calcule tel que j<=i et
  // mom4(i,j,k) tel que k<=j<=i ce qui constitue les
  // support non redondant.
  //
  /////////////////////////////////////////////////////


  long FenetreMoment=(Moment4->size())/2;
  long NbrEchantillon=VecteurCourant->size();

  // Allocation memoire tampon pour le calcul des moments 
   
  double vec_2=0;
  double * vec_3=new double [FenetreMoment+1];
  double ** vec_4=new double* [FenetreMoment+1];

  for (i=0;i<=FenetreMoment;i++)
    vec_4[i]=new double [FenetreMoment+1];

  /// Fin d'allocation de memoire locale

  for (i=0;i<=FenetreMoment;i++)
    {
      vec_3[i]=0.0;
      for(j=0;j<=FenetreMoment;j++)
	vec_4[i][j]=0.0;
    }


  for(i=Offset;i<=FenetreMoment;i++)
    {
      // Description des indices pour lequel les
      // moments doivent etre calcules

      for(n=0;n<NbrEchantillon-i;n++)
	{  

	  // Description de tous les indices qui inetrviennent
	  // dans l'estimation des moments 
	  // calcul du produit x(n)*x(n-i) 

	  cour=(*VecteurCourant)(n)*(*VecteurCourant)(n+i);
	  for(j=0;j<=i;j++)
	    {
	      //puis calcul du produit  x(n)*x(n+i)*x(n+j)
	      //pour le moemnt d'ordre 3

	      cour1=cour*(*VecteurCourant)(n+j);
	      for(k=0;k<=j;k++)
		{
		  // puis calcul du produit  x(n)*x(n+i)*x(n-i-j)*x(n-i-k)
		  // pour le moment d'ordre 3

		  vec_4[j][k]=vec_4[j][k]+cour1*(*VecteurCourant)(n+k);
		}
	      vec_3[j]=vec_3[j]+cour1;
	    }
	  vec_2+=cour;
	}

      ////////////////////////////////////////////
      //
      // lorsque tous les indices ont ete parcourus
      // le resulats est normalise 
      // puis ranges dans les structures charge de recevoir le resultat 
      // et les buffer nettoyes 
      //
      ///////////////////////////////////////////
      
      (*Moment2)(i)=vec_2/((double)(NbrEchantillon-i));

      if(i !=0)
	  (*Moment2)(2*FenetreMoment-i)=(*Moment2)(i);


      vec_2=0.0;

      for(j=0;j<=i;j++)
	{
	  (*Moment3)(i,j)=vec_3[j]/((double)(NbrEchantillon-i));
	  vec_3[j]=0.0;
	  for(k=0;k<=j;k++)
	    {
	      (*Moment4)(i,j,k)=vec_4[j][k]/((double)(NbrEchantillon-i));
	      vec_4[j][k]=0.0;
	    }
	}
    }
   
   
  // liberation de la memoire allouee localement 
  delete [] vec_3;
   
  for(i=0;i<=FenetreMoment;i++)
    delete [] vec_4[i];
   
  delete [] vec_4;
}  








