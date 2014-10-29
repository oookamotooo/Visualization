/*
#include <stdio.h>
#include <math.h>
//#include <lucomp.h>

#define EPS 1.0E-10

int lucomp(double fjac[3][3],
		   double fvec[3],
		   double * ds,double * dt,double * du)
{

  double dx[3],ul[3][3];
  int ierr,l[3];

  lu(fjac,ul,l,ierr);
  solve(ul,fvec,dx,l);

  *ds=dx[0]; *dt=dx[1]; *du=dx[2];

  return 0;
}

int lu(double fjac[3][3],double ul[3][3],
	   int l[3],int ierr)
{

  int i,j,k,ij,ik,dummy;
  double Max;

  for (i=0; i<=2; i++)
    for (j=0; j<=2; j++)
      ul[i][j]=fjac[i][j];

  l[0]=0; l[1]=1; l[2]=2;

  for (k=0; k<=1; k++) {
    ik=k;
    Max=ul[l[ik]][k];
    for (ij=k+1; ij<=2; ij++) {
      if (fabs(ul[l[ij]][k]) > fabs(Max)) {
	Max=ul[l[ij]][k];
	ik=ij;
      }
    }

    if (fabs(Max) <= EPS) return -1;

    if (ik != k) {
      dummy=l[k];
      l[k]=l[ik];
      l[ik]=dummy;
    }

    for (i=k+1; i<=2; i++) {
      ul[l[i]][k]=ul[l[i]][k]/Max;

      for (j=k+1; j<=2; j++)
	ul[l[i]][j]=ul[l[i]][j]-ul[l[i]][k]*ul[l[k]][j];
    }
  }

  return 0;
}

int solve(double ul[3][3],double fvec[3],
		  double dx[3],int l[3])
//double ul[3][3],fvec[3];
//double dx[3];
//int l[3];
{

  int i,j;

  for (i=0; i<=2; i++) {
    dx[i]=fvec[l[i]];
    for (j=0; j<=i-1; j++)
      dx[i]=dx[i]-ul[l[i]][j]*dx[j];
  }

  for (i=2; i>=0; i--) {
    for (j=i+1; j<=2; j++)
      dx[i]=dx[i]-ul[l[i]][j]*dx[j];
    dx[i]=dx[i]/ul[l[i]][i];
  }

  return 0;
}*/