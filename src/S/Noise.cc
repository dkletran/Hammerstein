
#include "general.h"



#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran1(long *idum)
{

  
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   float temp;
  
   
   if (*idum <= 0 || !iy)
   {
      if (-(*idum)<1)
	 *idum=1;
      else
	 *idum =-(*idum);
    

      for (j=NTAB+7; j>=0; j--)
      {
	 k=(*idum)/IQ;
	 *idum=IA*(*idum-k*IQ)-IR*k;
	 if (*idum <0) *idum += IM;
	 if (j < NTAB) iv[j]= *idum;
      }
      iy=iv[0];
   }

    

   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if(*idum <0)
      *idum+=IM;
   j=iy/NDIV;
   iy=iv[j];
   

   iv[j]=*idum;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}


float gasdev(struct drand48_data *tampon)
{ 
   static int iset=0;
   static float gset;
   float fac,rsq,v1,v2;
   double resultat;


   if (iset == 0)
   {
      do {
	drand48_r(tampon, &resultat);
	v2=2.0*((float)(resultat))-1.0; 
	drand48_r(tampon, &resultat);
	v1=2.0*((float)(resultat))-1.0; 
	
	
	rsq=v1*v1+v2*v2;
      }
      while (rsq>= 1.0 || rsq==0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   }
   else
   {
      iset=0;
      return gset;
   }

    do {	
      drand48_r(tampon, &resultat);
      v2=2.0*((float)(resultat))-1.0; 
      drand48_r(tampon, &resultat);
      v1=2.0*((float)(resultat))-1.0; 
      
     rsq=v1*v1+v2*v2;
   }
   while (rsq>= 1.0 || rsq==0.0);

   return v2*sqrt(-2.0*log(rsq)/rsq);


}




TVecteur<double> ranUniforme(int N1)
{

  
  int i;
  TVecteur<double>  bruit(N1);

  unsigned short int seed16v[3];
  

  seed16v[0]=time(NULL);

  struct drand48_data *tampon=new  struct drand48_data [1];
  double resultat;
  seed48_r(seed16v,tampon);
  // coup pour rin car le debut de la seed est trs semblable
  drand48_r(tampon, &resultat);

  for(i=0;i<N1;i++) 
    {
      drand48_r(tampon, &resultat);
      bruit(i)=resultat;
    }

  return(bruit);

}



void  ranUniforme (TVecteur<double> *bruit)
{


  int i;


  unsigned short int seed16v[3];
  seed16v[0]=time(NULL);

  struct drand48_data *tampon=new  struct drand48_data [1];

  double resultat;
  seed48_r(seed16v,tampon);
 // coup pour rin car le debut de la seed est trs semblable
  drand48_r(tampon, &resultat);


  
  long N1=bruit->size();

  for(i=0;i<N1;i++) 
    {
       drand48_r(tampon, &resultat);
       (*bruit)(i)=resultat;
    }

}

TImage<float> ranUniforme(int N1,int N2)
{

  
  int i,j;

  unsigned short int seed16v[3];  
  seed16v[0]=time(NULL);

  struct drand48_data *tampon=new  struct drand48_data [1];

  double resultat;
  seed48_r(seed16v,tampon);
 // coup pour rin car le debut de la seed est trs semblable
  drand48_r(tampon, &resultat);




  TImage<float>  bruit(N1,N2);

    for(i=0;i<N1;i++) 
      for(j=0;j<N2;j++) 
	{
	  drand48_r(tampon, &resultat);
	  bruit(i,j)=resultat;
	}

    return(bruit);

}


TImage<float> ranGaussien(int N1,int N2)
{
  // Generation d'un bruit blanc gaussien;

  int i,j;
struct drand48_data *tampon=new  struct drand48_data [1];

  unsigned short int seed16v[3];
  seed16v[0]=time(NULL);
  seed48_r(seed16v,tampon);


  double resultat;
  // coup pour rin car le debut de la seed est trs semblable
  drand48_r(tampon, &resultat);

  TImage<float>  bruit(N1,N2);

    for(i=0;i<N1;i++) 
      for(j=0;j<N2;j++)
	bruit(i,j)=gasdev(tampon);

    
    return(bruit);

}


TVecteur<float> ranGaussien(int N1)
{
  ///////////////////////////////////////
  //
  // Generation d'un bruit blanc gaussien;
  //
  ////////////////////////////////////////


  int i;

  struct drand48_data *tampon=new  struct drand48_data [1];

  unsigned short int seed16v[3];
  seed16v[0]=time(NULL);
  seed48_r(seed16v,tampon);


  double resultat;
  // coup pour rin car le debut de la seed est trs semblable
  drand48_r(tampon, &resultat);



  TVecteur<float>  bruit(N1);

    for(i=0;i<N1;i++)
      bruit(i)=gasdev(tampon);
    

    return(bruit);

}









TImage<float> ranRayleigh(int N1,int N2)
{

  // Generation d'un bruit blanc Rayleigh;

  int i,j;

 
  TImage<float> bruit(N1,N2);

  TImage<float> bruitgaussien=ranGaussien(2*N1,N2);

    for(i=0;i<N1;i++) 
      for(j=0;j<N2;j++)
	bruit(i,j)=sqrt(pow(bruitgaussien(i,j),2)+
			pow(bruitgaussien(i+N1,j),2)); 

    return(bruit);

}


