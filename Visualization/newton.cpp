#include "Vector3.h"
#define sizeX 91
#define sizeY 171
#define sizeZ 86

int newton(int x,int y,int z,double *s,double *t,double *u,
		   double bx[sizeX][sizeY][sizeZ],
           double by[sizeX][sizeY][sizeZ],
           double bz[sizeX][sizeY][sizeZ])
{
	double Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8;
	double Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8;
	double Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8;

	int count,x_c,y_c,z_c;
	double fjac[3][3],fvec[3];
	double ds,dt,du;
	double norm1,norm2;

	set_data(x,y,z,bx,by,bz,
		Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8,
		Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8,
		Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8);

    for (x_c=0; x_c<=10; x_c++) {
		*s=x_c/10;
		for (y_c=0; y_c<=10; y_c++) {
			*t=y_c/10.0;
			for (z_c=0; z_c<=10; z_c++) {
				*u=z_c/10;
				if (fabs(Vx(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8)) < 0.1 
					&& fabs(Vy(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8)) < 0.1 
				&& fabs(Vz(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8)) < 0.1) {
					x_c=11; y_c=11; z_c=11; break;
				}
			}
		}
    }
    for (count=1; count<=NTRIAL; count++) {
		fvec[0]=-Vx(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8); 
		fvec[1]=-Vy(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8); 
		fvec[2]=-Vz(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8);
		fjac[0][0]=Vxs(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8); 
		fjac[0][1]=Vxt(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8); 
		fjac[0][2]=Vxu(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8);
		fjac[1][0]=Vys(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8); 
		fjac[1][1]=Vyt(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8); 
		fjac[1][2]=Vyu(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8);
		fjac[2][0]=Vzs(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8); 
		fjac[2][1]=Vzt(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8); 
		fjac[2][2]=Vzu(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8);
		lucomp(fjac,fvec,&ds,&dt,&du);
		norm1=sqrt(ds*ds+dt*dt+du*du);
		*s=*s+ds; *t=*t+dt; *u=*u+du;
		norm2=sqrt(Vx(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8)*Vx(*s,*t,*u,Vx1,Vx2,Vx3,Vx4,Vx5,Vx6,Vx7,Vx8)
			+Vy(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8)*Vy(*s,*t,*u,Vy1,Vy2,Vy3,Vy4,Vy5,Vy6,Vy7,Vy8)
		 +Vz(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8)*Vz(*s,*t,*u,Vz1,Vz2,Vz3,Vz4,Vz5,Vz6,Vz7,Vz8));
		if ((norm1 < 1.0E-15) || (norm2 < 1.0E-15)) {
			return 0;
		}
    }
    return -1;
}

double Vx(double s,double t,double u,
		  double Vx1,double Vx2,double Vx3,double Vx4,
		  double Vx5,double Vx6,double Vx7,double Vx8)
{
  return u*(1-s)*(1-t)*Vx5+u*t*(1-s)*Vx8+u*s*(1-t)*Vx6+u*s*t*Vx7
	+(1-u)*(1-s)*(1-t)*Vx1+(1-u)*(1-s)*t*Vx4+(1-u)*s*(1-t)*Vx2
	+(1-u)*s*t*Vx3;
}

double Vy(double s,double t,double u,
		  double Vy1,double Vy2,double Vy3,double Vy4,
		  double Vy5,double Vy6,double Vy7,double Vy8)
{
  return u*(1-s)*(1-t)*Vy5+u*t*(1-s)*Vy8+u*s*(1-t)*Vy6+u*s*t*Vy7
	+(1-u)*(1-s)*(1-t)*Vy1+(1-u)*(1-s)*t*Vy4+(1-u)*s*(1-t)*Vy2
	+(1-u)*s*t*Vy3;
}

double Vz(double s,double t,double u,
		  double Vz1,double Vz2,double Vz3,double Vz4,
		  double Vz5,double Vz6,double Vz7,double Vz8)
{
  return u*(1-s)*(1-t)*Vz5+u*t*(1-s)*Vz8+u*s*(1-t)*Vz6+u*s*t*Vz7
	+(1-u)*(1-s)*(1-t)*Vz1+(1-u)*(1-s)*t*Vz4+(1-u)*s*(1-t)*Vz2
	+(1-u)*s*t*Vz3;
}

double Vxs(double s,double t,double u,
		  double Vx1,double Vx2,double Vx3,double Vx4,
		  double Vx5,double Vx6,double Vx7,double Vx8)
{
  return -u*(1-t)*Vx5-u*t*Vx8+u*(1-t)*Vx6+u*t*Vx7-(1-u)*(1-t)*Vx1
	-(1-u)*t*Vx4+(1-u)*(1-t)*Vx2+(1-u)*t*Vx3;
}

double Vxt(double s,double t,double u,
		  double Vx1,double Vx2,double Vx3,double Vx4,
		  double Vx5,double Vx6,double Vx7,double Vx8)
{
  return -u*(1-s)*Vx5+u*(1-s)*Vx8-u*s*Vx6+u*s*Vx7-(1-u)*(1-s)*Vx1
	+(1-u)*(1-s)*Vx4-(1-u)*s*Vx2+(1-u)*s*Vx3;
}

double Vxu(double s,double t,double u,
		  double Vx1,double Vx2,double Vx3,double Vx4,
		  double Vx5,double Vx6,double Vx7,double Vx8)
{
  return (1-s)*(1-t)*Vx5+t*(1-s)*Vx8+s*(1-t)*Vx6+s*t*Vx7
	-(1-s)*(1-t)*Vx1-(1-s)*t*Vx4-s*(1-t)*Vx2-s*t*Vx3;
}

double Vys(double s,double t,double u,
		  double Vy1,double Vy2,double Vy3,double Vy4,
		  double Vy5,double Vy6,double Vy7,double Vy8)
{
  return -u*(1-t)*Vy5-u*t*Vy8+u*(1-t)*Vy6+u*t*Vy7-(1-u)*(1-t)*Vy1
	-(1-u)*t*Vy4+(1-u)*(1-t)*Vy2+(1-u)*t*Vy3;
}

double Vyt(double s,double t,double u,
		  double Vy1,double Vy2,double Vy3,double Vy4,
		  double Vy5,double Vy6,double Vy7,double Vy8)
{
  return -u*(1-s)*Vy5+u*(1-s)*Vy8-u*s*Vy6+u*s*Vy7-(1-u)*(1-s)*Vy1
	+(1-u)*(1-s)*Vy4-(1-u)*s*Vy2+(1-u)*s*Vy3;
}

double Vyu(double s,double t,double u,
		  double Vy1,double Vy2,double Vy3,double Vy4,
		  double Vy5,double Vy6,double Vy7,double Vy8)
{
  return (1-s)*(1-t)*Vy5+t*(1-s)*Vy8+s*(1-t)*Vy6+s*t*Vy7
	-(1-s)*(1-t)*Vy1-(1-s)*t*Vy4-s*(1-t)*Vy2-s*t*Vy3;
}

double Vzs(double s,double t,double u,
		  double Vz1,double Vz2,double Vz3,double Vz4,
		  double Vz5,double Vz6,double Vz7,double Vz8)
{
  return -u*(1-t)*Vz5-u*t*Vz8+u*(1-t)*Vz6+u*t*Vz7-(1-u)*(1-t)*Vz1
	-(1-u)*t*Vz4+(1-u)*(1-t)*Vz2+(1-u)*t*Vz3;
}

double Vzt(double s,double t,double u,
		  double Vz1,double Vz2,double Vz3,double Vz4,
		  double Vz5,double Vz6,double Vz7,double Vz8)
{
  return -u*(1-s)*Vz5+u*(1-s)*Vz8-u*s*Vz6+u*s*Vz7-(1-u)*(1-s)*Vz1
	+(1-u)*(1-s)*Vz4-(1-u)*s*Vz2+(1-u)*s*Vz3;
}

double Vzu(double s,double t,double u,
		  double Vz1,double Vz2,double Vz3,double Vz4,
		  double Vz5,double Vz6,double Vz7,double Vz8)
{
  return (1-s)*(1-t)*Vz5+t*(1-s)*Vz8+s*(1-t)*Vz6+s*t*Vz7
	-(1-s)*(1-t)*Vz1-(1-s)*t*Vz4-s*(1-t)*Vz2-s*t*Vz3;
}

int set_data(int x,int y,int z,
			 double x_data[sizeX][sizeY][sizeZ],
             double y_data[sizeX][sizeY][sizeZ],
             double z_data[sizeX][sizeY][sizeZ],
			 double &Vx1,double &Vx2,double &Vx3,double &Vx4,
			 double &Vx5,double &Vx6,double &Vx7,double &Vx8,
			 double &Vy1,double &Vy2,double &Vy3,double &Vy4,
			 double &Vy5,double &Vy6,double &Vy7,double &Vy8,
			 double &Vz1,double &Vz2,double &Vz3,double &Vz4,
			 double &Vz5,double &Vz6,double &Vz7,double &Vz8)
{
    Vx1=x_data[z-1][y-1][x-1]; Vx2=x_data[z-1][y-1][x];
    Vx3=x_data[z-1][y][x]; Vx4=x_data[z-1][y][x-1];
    Vx5=x_data[z][y-1][x-1]; Vx6=x_data[z][y-1][x];
    Vx7=x_data[z][y][x]; Vx8=x_data[z][y][x-1];

    Vy1=y_data[z-1][y-1][x-1]; Vy2=y_data[z-1][y-1][x];
    Vy3=y_data[z-1][y][x]; Vy4=y_data[z-1][y][x-1];
    Vy5=y_data[z][y-1][x-1]; Vy6=y_data[z][y-1][x];
    Vy7=y_data[z][y][x]; Vy8=y_data[z][y][x-1];

    Vz1=z_data[z-1][y-1][x-1]; Vz2=z_data[z-1][y-1][x];
    Vz3=z_data[z-1][y][x]; Vz4=z_data[z-1][y][x-1];
    Vz5=z_data[z][y-1][x-1]; Vz6=z_data[z][y-1][x];
    Vz7=z_data[z][y][x]; Vz8=z_data[z][y][x-1];

    return 0;
}

