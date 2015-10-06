

#include <complex>
#include <iostream>
#include "Moment4_2D.h"


template <class T>
TMoment4_2D<T>::~TMoment4_2D()
{

  int ix,iy,jx,jy,kx,ky;
 

  for (ix = 0; ix <= fenetre; ix++ )
    for (iy = 0; iy <= 2*fenetre; iy++ )
      for (jx = 0; jx <= ix; jx++ )
	for (jy = 0; jy <= 2*fenetre; jy++ )
	  for (kx = 0; kx <= jx; kx++ )
	    delete[] data[ix][iy][jx][jy][kx];


  for (ix = 0; ix <= fenetre; ix++ )
    for (iy = 0; iy <= 2*fenetre; iy++ )
      for (jx = 0; jx <= ix; jx++ )
	for (jy = 0; jy <= 2*fenetre; jy++ )
	  delete[] data[ix][iy][jx][jy];


  for (ix = 0; ix <= fenetre; ix++ )
    for (iy = 0; iy <= 2*fenetre; iy++ )
      for (jx = 0; jx <= ix; jx++ )
	delete[] data[ix][iy][jx];


  for (ix = 0; ix <= fenetre; ix++ )
    for (iy = 0; iy <= 2*fenetre; iy++ )
      delete[] data[ix][iy];


  for (ix = 0; ix <= fenetre; ix++ )
    delete[] data[ix];

  
  delete[] data;

  
  data=NULL;
  fenetre = 0;
  
}

template <class T>
TMoment4_2D<T>::TMoment4_2D() 
{
  
  fenetre = 0;
  data=NULL;
}


template <class T>
TMoment4_2D<T>::TMoment4_2D(int n) 
{

  int ix,jx,kx,iy,jy,ky;
  fenetre= n/2;
  
  
  if ((data=new T***** [fenetre+1])==NULL)
    {
      cout << "Probleme d'allocation memoire pour le cumulant sur la premiere dimension"<<endl;
      exit(1);
      
    }    
  
	
  for (ix=0;ix<=fenetre;ix++)
    {                       
      if ((data[ix]=new T**** [2*fenetre+1])==NULL)
	{
	  cout << "Probleme d'allocation memoire pour le cumulant sur la deuxieme dimension"<<endl;
	  exit(1);
	} 
    }
  
  
  for (ix=0;ix<=fenetre;ix++)
    for (iy=0;iy<=2*fenetre;iy++)
      {         
	if ((data[ix][iy]=new T*** [ix+1])==NULL)
	  { 
	    cout << "Probleme d'allocation memoire pour le cumulant sur la troisieme dimension"<<endl;
	    exit(1);
	  } 
      }      
  
  



  for (ix=0;ix<=fenetre;ix++)
    for (iy=0;iy<=2*fenetre;iy++)
      for (jx=0;jx<=ix;jx++)
	{       
	  if ((data[ix][iy][jx]=new T** [2*fenetre+1])==NULL)
	    {
	      cout << "Probleme d'allocation memoire pour le cumulant sur la quatrieme dimension"<< ix<<
		"  "<<iy<<" "<< jx<<endl;
	      exit(1);
	    }
	}




  for (ix=0;ix<=fenetre;ix++)
    for (iy=0;iy<=2*fenetre;iy++)
      for (jx=0;jx<=ix;jx++)    
	for (jy=0;jy<=2*fenetre;jy++)

	{       
	  if ((data[ix][iy][jx][jy]=new T* [jx+1])==NULL)
	    {
	      cout << "Probleme d'allocation memoire pour le cumulant sur la cinquieme dimension"<< ix<<
		"  "<<iy<<" "<< jx<<"   " << jy<< endl;
	      exit(1);
	    }
	}


  for (ix=0;ix<=fenetre;ix++)
    for (iy=0;iy<=2*fenetre;iy++)
      for (jx=0;jx<=ix;jx++)    
	for (jy=0;jy<=2*fenetre;jy++)
	  for (kx=0;kx<=jx;kx++)    
	    {       
	      if ((data[ix][iy][jx][jy][kx]=new T [2*fenetre+1])==NULL)
		{
		  cout << "Probleme d'allocation memoire pour le cumulant sur la sixieme dimension"<< ix<<
		    "  "<<iy<<" "<< jx<<endl;
		  exit(1);
		}
	    }



  for (ix=0;ix<=fenetre;ix++)
    for (iy=0;iy<=2*fenetre;iy++)
      for (jx=0;jx<=ix;jx++)    
	for (jy=0;jy<=2*fenetre;jy++)
	  for (kx=0;kx<=jx;kx++) 
	    for (ky=0;ky<=2*fenetre;ky++)
	      data[ix][iy][jx][jy][kx][ky]=0;


 

}


template class TMoment4_2D<short>;
template class TMoment4_2D<float>;
template class TMoment4_2D<double>;

