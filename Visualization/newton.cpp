/*#include <one_basic_include.h>
#include <newton.h>
#include <simarea.h>

#if 0
#define sizeX  61
#define sizeY    61
#define sizeZ 61 
#endif*/

#include <math.h>
#include "newton.h"
#include "Vector3.h"
#include <vector>
#include "define.h"

#define NTRIAL 50000
#define EPS 1.0E-10

// kokomadayattenaiyo main to jacob dakedayo


/*
extern double x_data[sizeX][sizeY][sizeZ];
extern double y_data[sizeX][sizeY][sizeZ];
extern double z_data[sizeX][sizeY][sizeZ];
*/


double Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8;
double Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8;
double Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8;

double Vx(double s,double t,double u)
{
  return u*(1-s)*(1-t)*Vx5+u*t*(1-s)*Vx8+u*s*(1-t)*Vx6+u*s*t*Vx7
	+(1-u)*(1-s)*(1-t)*Vx1+(1-u)*(1-s)*t*Vx4+(1-u)*s*(1-t)*Vx2
	+(1-u)*s*t*Vx3;
}

double Vy(double s,double t,double u)
{
  return u*(1-s)*(1-t)*Vy5+u*t*(1-s)*Vy8+u*s*(1-t)*Vy6+u*s*t*Vy7
	+(1-u)*(1-s)*(1-t)*Vy1+(1-u)*(1-s)*t*Vy4+(1-u)*s*(1-t)*Vy2
	+(1-u)*s*t*Vy3;
}

double Vz(double s,double t,double u)
{
  return u*(1-s)*(1-t)*Vz5+u*t*(1-s)*Vz8+u*s*(1-t)*Vz6+u*s*t*Vz7
	+(1-u)*(1-s)*(1-t)*Vz1+(1-u)*(1-s)*t*Vz4+(1-u)*s*(1-t)*Vz2
	+(1-u)*s*t*Vz3;
}

double Vxs(double s,double t,double u)
{
  return -u*(1-t)*Vx5-u*t*Vx8+u*(1-t)*Vx6+u*t*Vx7-(1-u)*(1-t)*Vx1
	-(1-u)*t*Vx4+(1-u)*(1-t)*Vx2+(1-u)*t*Vx3;
}

double Vxt(double s,double t,double u)
{
  return -u*(1-s)*Vx5+u*(1-s)*Vx8-u*s*Vx6+u*s*Vx7-(1-u)*(1-s)*Vx1
	+(1-u)*(1-s)*Vx4-(1-u)*s*Vx2+(1-u)*s*Vx3;
}

double Vxu(double s,double t,double u)
{
  return (1-s)*(1-t)*Vx5+t*(1-s)*Vx8+s*(1-t)*Vx6+s*t*Vx7
	-(1-s)*(1-t)*Vx1-(1-s)*t*Vx4-s*(1-t)*Vx2-s*t*Vx3;
}

double Vys(double s,double t,double u)
{
  return -u*(1-t)*Vy5-u*t*Vy8+u*(1-t)*Vy6+u*t*Vy7-(1-u)*(1-t)*Vy1
	-(1-u)*t*Vy4+(1-u)*(1-t)*Vy2+(1-u)*t*Vy3;
}

double Vyt(double s,double t,double u)
{
  return -u*(1-s)*Vy5+u*(1-s)*Vy8-u*s*Vy6+u*s*Vy7-(1-u)*(1-s)*Vy1
	+(1-u)*(1-s)*Vy4-(1-u)*s*Vy2+(1-u)*s*Vy3;
}

double Vyu(double s,double t,double u)
{
  return (1-s)*(1-t)*Vy5+t*(1-s)*Vy8+s*(1-t)*Vy6+s*t*Vy7
	-(1-s)*(1-t)*Vy1-(1-s)*t*Vy4-s*(1-t)*Vy2-s*t*Vy3;
}

double Vzs(double s,double t,double u)
{
  return -u*(1-t)*Vz5-u*t*Vz8+u*(1-t)*Vz6+u*t*Vz7-(1-u)*(1-t)*Vz1
	-(1-u)*t*Vz4+(1-u)*(1-t)*Vz2+(1-u)*t*Vz3;
}

double Vzt(double s,double t,double u)
{
  return -u*(1-s)*Vz5+u*(1-s)*Vz8-u*s*Vz6+u*s*Vz7-(1-u)*(1-s)*Vz1
	+(1-u)*(1-s)*Vz4-(1-u)*s*Vz2+(1-u)*s*Vz3;
}

double Vzu(double s,double t,double u)
{
  return (1-s)*(1-t)*Vz5+t*(1-s)*Vz8+s*(1-t)*Vz6+s*t*Vz7
	-(1-s)*(1-t)*Vz1-(1-s)*t*Vz4-s*(1-t)*Vz2-s*t*Vz3;
}

int set_data(int i, int j, int k)
{
	int index;
	
	index = Index(i-1, j-1, k-1);
    Vx1=datas[index].x;	Vy1=datas[index].y;	Vz1=datas[index].z;
	
	index = Index(i, j-1, k-1);
    Vx2=datas[index].x;	Vy2=datas[index].y;	Vz2=datas[index].z;
	
	index = Index(i, j, k-1);
    Vx3=datas[index].x;	Vy3=datas[index].y;	Vz3=datas[index].z;
	
	index = Index(i-1, j, k-1);
    Vx4=datas[index].x;	Vy4=datas[index].y;	Vz4=datas[index].z;
	
	index = Index(i-1, j-1, k);
    Vx5=datas[index].x;	Vy5=datas[index].y;	Vz5=datas[index].z;
	
	index = Index(i, j-1, k);
    Vx6=datas[index].x;	Vy6=datas[index].y;	Vz6=datas[index].z;
	
	index = Index(i, j, k);
    Vx7=datas[index].x;	Vy7=datas[index].y;	Vz7=datas[index].z;
	
	index = Index(i-1, j, k);
    Vx8=datas[index].x;	Vy8=datas[index].y;	Vz8=datas[index].z;

    return 0;
}

int newton1(int x,int y,int z/*double *s,double *t,double *u*/)
{

	int x_c,y_c,z_c;
	double fjac[3][3],fvec[3];
	double ds,dt,du;
	double norm1,norm2;

	double bunkatsu = 500.0;

	set_data(x,y,z);

	//以下元々のソースコード
	/*
    for (x_c=0; x_c<=10; x_c++) {
		*s = (double)(x_c/10);
		for (y_c=0; y_c<=10; y_c++) {
			*t=y_c/10.0;
			for (z_c=0; z_c<=10; z_c++) {
				*u=z_c/10;
				if (fabs(Vx(*s,*t,*u)) < 0.1 && fabs(Vy(*s,*t,*u)) < 0.1 && fabs(Vz(*s,*t,*u)) < 0.1) {
					printf("find\n");
					x_c=11; y_c=11; z_c=11; break;
				}
			}
		}
    }*/
	
	for (x_c=0; x_c<=bunkatsu; x_c++) {
		s = ((double)x_c/bunkatsu);
		for (y_c=0; y_c<=bunkatsu; y_c++) {
			t = ((double)y_c/bunkatsu);
			for (z_c=0; z_c<=bunkatsu; z_c++) {
				u = ((double)z_c/bunkatsu);
				/*
				printf("%d %d %d\n", x, y, z);
				printf("%f %f %f\n", s, t, u);
				printf("%f\n", Vx(s,t,u));
				*/
				if (fabs(Vx(s,t,u)) < 0.00001 && fabs(Vy(s,t,u)) < 0.00001 && fabs(Vz(s,t,u)) < 0.00001) {
					printf("find\n");
					cp[cp_count].x = (double)x+s;	cp[cp_count].y = (double)y+t;	cp[cp_count].z = (double)z+u;
					printf("%f, %f, %f\n", cp[cp_count].x, cp[cp_count].y, cp[cp_count].z);
					cp_count++;
					x_c=bunkatsu+1; y_c=bunkatsu+1; z_c=bunkatsu+1; break;
				}
			}
		}
		//printf("%d %d %d : x = %d\n", x, y, z, x_c);
    }
	/*
	for (int count=1; count<=NTRIAL; count++) {
      fvec[0]=-Vx(s,t,u); fvec[1]=-Vy(s,t,u); fvec[2]=-Vz(s,t,u);
      fjac[0][0]=Vxs(s,t,u); fjac[0][1]=Vxt(s,t,u); fjac[0][2]=Vxu(s,t,u);
      fjac[1][0]=Vys(s,t,u); fjac[1][1]=Vyt(s,t,u); fjac[1][2]=Vyu(s,t,u);
      fjac[2][0]=Vzs(s,t,u); fjac[2][1]=Vzt(s,t,u); fjac[2][2]=Vzu(s,t,u);
      //lucomp(fjac,fvec,&ds,&dt,&du);
      norm1=sqrt(ds*ds+dt*dt+du*du);
      s=s+ds; t=t+dt; u=u+du;
      norm2=sqrt(Vx(s,t,u)*Vx(s,t,u)+Vy(s,t,u)*Vy(s,t,u)
		 +Vz(s,t,u)*Vz(s,t,u));
      if ((norm1 < 1.0E-15) || (norm2 < 1.0E-15)) {
	return 0;
      }
    }
    return -1;
	*/
	
	/*
	
    for (count=1; count<=NTRIAL; count++) {
      fvec[0]=-Vx(*s,*t,*u); fvec[1]=-Vy(*s,*t,*u); fvec[2]=-Vz(*s,*t,*u);
      fjac[0][0]=Vxs(*s,*t,*u); fjac[0][1]=Vxt(*s,*t,*u); fjac[0][2]=Vxu(*s,*t,*u);
      fjac[1][0]=Vys(*s,*t,*u); fjac[1][1]=Vyt(*s,*t,*u); fjac[1][2]=Vyu(*s,*t,*u);
      fjac[2][0]=Vzs(*s,*t,*u); fjac[2][1]=Vzt(*s,*t,*u); fjac[2][2]=Vzu(*s,*t,*u);
      lucomp(fjac,fvec,&ds,&dt,&du);
      norm1=sqrt(ds*ds+dt*dt+du*du);
      *s=*s+ds; *t=*t+dt; *u=*u+du;
      norm2=sqrt(Vx(*s,*t,*u)*Vx(*s,*t,*u)+Vy(*s,*t,*u)*Vy(*s,*t,*u)
		 +Vz(*s,*t,*u)*Vz(*s,*t,*u));
      if ((norm1 < 1.0E-15) || (norm2 < 1.0E-15)) {
	return 0;
      }
    }
    return -1;*/

	return 0;
}







