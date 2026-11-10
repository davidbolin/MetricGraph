#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

void BB(int *nEdges,double *lenEdges,int *nPointsE,double *vecNorm,double *vec)
{
     int id=0,n;
     double tj,ta,aux,T;

     for(int i=0;i<(*nEdges);i++)
     {
         T = lenEdges[i];
         n=nPointsE[i];

         vec[id]=0;
         id=id+1;
         for(int j=1;j<n;j++)
         {
		  if(j==(n-1)){vec[id]=0;}
		  else{
		         tj=j*T/(n-1);
		         ta=(j-1)*T/(n-1);
                         aux=(T-2*tj+ta)/(T-ta);
                         vec[id]=0.5*((1+aux)*vec[id-1]+vecNorm[id]*pow(T-ta-(T-2*tj+ta)*aux,0.5));
                      }
                      id=id+1;
         }
     }

}
