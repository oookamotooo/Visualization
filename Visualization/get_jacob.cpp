#include <stdio.h>
//#include <simarea.h>
#include "get_jacob.h"
#include "Vector3.h"
#include <vector>
#include "define.h"

/*
#if 0
#define LAYERS 61
#define ROWS 61
#define COLUMNS 61
#endif
*/

//global
//vector<Vector3d> datas;

//local
double Vol1,Vol2,Vol3,Vol4,Vol5,Vol6,Vol7,Vol8;

/*
extern double Vol1,Vol2,Vol3,Vol4,Vol5,Vol6,Vol7,Vol8;
extern double x_data[LAYERS][ROWS][COLUMNS];
extern double y_data[LAYERS][ROWS][COLUMNS];
extern double z_data[LAYERS][ROWS][COLUMNS];
*/

//extern int x,y,z;  //not used in this function

void volume(double x, double y, double z){
	
	double s = x - floor(x);
	double t = y - floor(y);
	double u = z - floor(z);

	Vol1=s*t*u;
	Vol2=(1.0-s)*t*u;
	Vol3=(1.0-s)*(1.0-t)*u;
	Vol4=s*(1.0-t)*u;
	Vol5=s*t*(1.0-u);
	Vol6=(1.0-s)*t*(1.0-u);
	Vol7=(1.0-s)*(1.0-t)*(1.0-u);
	Vol8=s*(1.0-t)*(1.0-u);

	return;
}

int get_jacob(double x, double y, double z, double jac[3][3])
{

	int i,j,k;
//	int index;

	double DBxDx[8],DBxDy[8],DBxDz[8];
	double DByDx[8],DByDy[8],DByDz[8];
	double DBzDx[8],DBzDy[8],DBzDz[8];

	i= floor(x); j= floor(y); k= floor(z);

	//(i,j,k)
	DBxDx[0]=0.5*(datas[Index(i+1, j, k)].x - datas[Index(i-1, j, k)].x);
	DBxDy[0]=0.5*(datas[Index(i, j+1, k)].x - datas[Index(i, j-1, k)].x);
	DBxDz[0]=0.5*(datas[Index(i, j, k+1)].x - datas[Index(i, j, k-1)].x);

	DByDx[0]=0.5*(datas[Index(i+1, j, k)].y - datas[Index(i-1, j, k)].y);
	DByDy[0]=0.5*(datas[Index(i, j+1, k)].y - datas[Index(i, j-1, k)].y);
	DByDz[0]=0.5*(datas[Index(i, j, k+1)].y - datas[Index(i, j, k-1)].y);

	DBzDx[0]=0.5*(datas[Index(i+1, j, k)].z - datas[Index(i-1, j, k)].z);
	DBzDy[0]=0.5*(datas[Index(i, j+1, k)].z - datas[Index(i, j-1, k)].z);
	DBzDz[0]=0.5*(datas[Index(i, j, k+1)].z - datas[Index(i, j, k-1)].z);

	//(i+1,j,k)
	DBxDx[1]=0.5*(datas[Index(i+2, j, k)].x - datas[Index(i, j, k)].x);
	DBxDy[1]=0.5*(datas[Index(i+1, j+1, k)].x - datas[Index(i+1, j-1, k)].x);
	DBxDz[1]=0.5*(datas[Index(i+1, j, k+1)].x - datas[Index(i+1, j, k-1)].x);

	DByDx[1]=0.5*(datas[Index(i+2, j, k)].y - datas[Index(i, j, k)].y);
	DByDy[1]=0.5*(datas[Index(i+1, j+1, k)].y - datas[Index(i+1, j-1, k)].y);
	DByDz[1]=0.5*(datas[Index(i+1, j, k+1)].y - datas[Index(i+1, j, k-1)].y);

	DBzDx[1]=0.5*(datas[Index(i+2, j, k)].z - datas[Index(i, j, k)].z);
	DBzDy[1]=0.5*(datas[Index(i+1, j+1, k)].z - datas[Index(i+1, j-1, k)].z);
	DBzDz[1]=0.5*(datas[Index(i+1, j, k+1)].z - datas[Index(i+1, j, k-1)].z);

	//(i+1, j+1, k)
	DBxDx[2]=0.5*(datas[Index(i+2, j+1, k)].x - datas[Index(i, j+1, k)].x);
	DBxDy[2]=0.5*(datas[Index(i+1, j+2, k)].x - datas[Index(i+1, j, k)].x);
	DBxDz[2]=0.5*(datas[Index(i+1, j+1, k+1)].x - datas[Index(i+1, j+1, k-1)].x);

	DByDx[2]=0.5*(datas[Index(i+2, j+1, k)].y - datas[Index(i, j+1, k)].y);
	DByDy[2]=0.5*(datas[Index(i+1, j+2, k)].y - datas[Index(i+1, j, k)].y);
	DByDz[2]=0.5*(datas[Index(i+1, j+1, k+1)].y - datas[Index(i+1, j+1, k-1)].y);

	DBzDx[2]=0.5*(datas[Index(i+2, j+1, k)].z - datas[Index(i, j+1, k)].z);
	DBzDy[2]=0.5*(datas[Index(i+1, j+2, k)].z - datas[Index(i+1, j, k)].z);
	DBzDz[2]=0.5*(datas[Index(i+1, j+1, k+1)].z - datas[Index(i+1, j+1, k-1)].z);

	//(i, j+1, k)

	DBxDx[3]=0.5*(datas[Index(i+1, j+1, k)].x - datas[Index(i-1, j+1, k)].x);
	DBxDy[3]=0.5*(datas[Index(i, j+2, k)].x - datas[Index(i, j, k)].x);
	DBxDz[3]=0.5*(datas[Index(i, j+1, k+1)].x - datas[Index(i, j+1, k-1)].x);

	DByDx[3]=0.5*(datas[Index(i+1, j+1, k)].y - datas[Index(i-1, j+1, k)].y);
	DByDy[3]=0.5*(datas[Index(i, j+2, k)].y - datas[Index(i, j, k)].y);
	DByDz[3]=0.5*(datas[Index(i, j+1, k+1)].y - datas[Index(i, j+1, k-1)].y);

	DBzDx[3]=0.5*(datas[Index(i+1, j+1, k)].z - datas[Index(i-1, j+1, k)].z);
	DBzDy[3]=0.5*(datas[Index(i, j+2, k)].z - datas[Index(i, j, k)].z);
	DBzDz[3]=0.5*(datas[Index(i, j+1, k+1)].z - datas[Index(i, j+1, k-1)].z);

	//(i, j, k+1)

	DBxDx[4]=0.5*(datas[Index(i+1, j, k+1)].x - datas[Index(i-1, j, k+1)].x);
	DBxDy[4]=0.5*(datas[Index(i, j+1, k+1)].x - datas[Index(i, j-1, k+1)].x);
	DBxDz[4]=0.5*(datas[Index(i, j, k+2)].x - datas[Index(i, j, k)].x);

	DByDx[4]=0.5*(datas[Index(i+1, j, k+1)].y - datas[Index(i-1, j, k+1)].y);
	DByDy[4]=0.5*(datas[Index(i, j+1, k+1)].y - datas[Index(i, j-1, k+1)].y);
	DByDz[4]=0.5*(datas[Index(i, j, k+2)].y - datas[Index(i, j, k)].y);

	DBzDx[4]=0.5*(datas[Index(i+1, j, k+1)].z - datas[Index(i-1, j, k+1)].z);
	DBzDy[4]=0.5*(datas[Index(i, j+1, k+1)].z - datas[Index(i, j-1, k+1)].z);
	DBzDz[4]=0.5*(datas[Index(i, j, k+2)].z - datas[Index(i, j, k)].z);

	//(i+1, j, k+1)

	DBxDx[5]=0.5*(datas[Index(i+2, j, k+1)].x - datas[Index(i, j, k+1)].x);
	DBxDy[5]=0.5*(datas[Index(i+1, j+1, k+1)].x - datas[Index(i+1, j-1, k+1)].x);
	DBxDz[5]=0.5*(datas[Index(i+1, j, k+2)].x - datas[Index(i+1, j, k)].x);

	DByDx[5]=0.5*(datas[Index(i+2, j, k+1)].y - datas[Index(i, j, k+1)].y);
	DByDy[5]=0.5*(datas[Index(i+1, j+1, k+1)].y - datas[Index(i+1, j-1, k+1)].y);
	DByDz[5]=0.5*(datas[Index(i+1, j, k+2)].y - datas[Index(i+1, j, k)].y);

	DBzDx[5]=0.5*(datas[Index(i+2, j, k+1)].z - datas[Index(i, j, k+1)].z);
	DBzDy[5]=0.5*(datas[Index(i+1, j+1, k+1)].z - datas[Index(i+1, j-1, k+1)].z);
	DBzDz[5]=0.5*(datas[Index(i+1, j, k+2)].z - datas[Index(i+1, j, k)].z);

	//(i+1, j+1, k+1)

	DBxDx[6]=0.5*(datas[Index(i+2, j+1, k+1)].x - datas[Index(i, j+1, k+1)].x);
	DBxDy[6]=0.5*(datas[Index(i+1, j+2, k+1)].x - datas[Index(i+1, j, k+1)].x);
	DBxDz[6]=0.5*(datas[Index(i+1, j+1, k+2)].x - datas[Index(i+1, j+1, k)].x);

	DByDx[6]=0.5*(datas[Index(i+2, j+1, k+1)].y - datas[Index(i, j+1, k+1)].y);
	DByDy[6]=0.5*(datas[Index(i+1, j+2, k+1)].y - datas[Index(i+1, j, k+1)].y);
	DByDz[6]=0.5*(datas[Index(i+1, j+1, k+2)].y - datas[Index(i+1, j+1, k)].y);

	DBzDx[6]=0.5*(datas[Index(i+2, j+1, k+1)].z - datas[Index(i, j+1, k+1)].z);
	DBzDy[6]=0.5*(datas[Index(i+1, j+2, k+1)].z - datas[Index(i+1, j, k+1)].z);
	DBzDz[6]=0.5*(datas[Index(i+1, j+1, k+2)].z - datas[Index(i+1, j+1, k)].z);

	//(i, j+1, k+1)

	DBxDx[7]=0.5*(datas[Index(i+1, j+1, k+1)].x - datas[Index(i-1, j+1, k+1)].x);
	DBxDy[7]=0.5*(datas[Index(i, j+2, k+1)].x - datas[Index(i, j, k+1)].x);
	DBxDz[7]=0.5*(datas[Index(i, j+1, k+2)].x - datas[Index(i, j+1, k)].x);

	DByDx[7]=0.5*(datas[Index(i+1, j+1, k+1)].y - datas[Index(i-1, j+1, k+1)].y);
	DByDy[7]=0.5*(datas[Index(i, j+2, k+1)].y - datas[Index(i, j, k+1)].y);
	DByDz[7]=0.5*(datas[Index(i, j+1, k+2)].y - datas[Index(i, j+1, k)].y);

	DBzDx[7]=0.5*(datas[Index(i+1, j+1, k+1)].z - datas[Index(i-1, j+1, k+1)].z);
	DBzDy[7]=0.5*(datas[Index(i, j+2, k+1)].z - datas[Index(i, j, k+1)].z);
	DBzDz[7]=0.5*(datas[Index(i, j+1, k+2)].z - datas[Index(i, j+1, k)].z);

	volume(x, y, z);
	
	jac[0][0]=DBxDx[0]*Vol7+DBxDx[1]*Vol8+DBxDx[2]*Vol5+DBxDx[3]*Vol6+
		DBxDx[4]*Vol3+DBxDx[5]*Vol4+DBxDx[6]*Vol1+DBxDx[7]*Vol2;

	jac[0][1]=DBxDy[0]*Vol7+DBxDy[1]*Vol8+DBxDy[2]*Vol5+DBxDy[3]*Vol6+
		DBxDy[4]*Vol3+DBxDy[5]*Vol4+DBxDy[6]*Vol1+DBxDy[7]*Vol2;

	jac[0][2]=DBxDz[0]*Vol7+DBxDz[1]*Vol8+DBxDz[2]*Vol5+DBxDz[3]*Vol6+
		DBxDz[4]*Vol3+DBxDz[5]*Vol4+DBxDz[6]*Vol1+DBxDz[7]*Vol2;

	jac[1][0]=DByDx[0]*Vol7+DByDx[1]*Vol8+DByDx[2]*Vol5+DByDx[3]*Vol6+
		DByDx[4]*Vol3+DByDx[5]*Vol4+DByDx[6]*Vol1+DByDx[7]*Vol2;

	jac[1][1]=DByDy[0]*Vol7+DByDy[1]*Vol8+DByDy[2]*Vol5+DByDy[3]*Vol6+
		DByDy[4]*Vol3+DByDy[5]*Vol4+DByDy[6]*Vol1+DByDy[7]*Vol2;

	jac[1][2]=DByDz[0]*Vol7+DByDz[1]*Vol8+DByDz[2]*Vol5+DByDz[3]*Vol6+
		DByDz[4]*Vol3+DByDz[5]*Vol4+DByDz[6]*Vol1+DByDz[7]*Vol2;

	jac[2][0]=DBzDx[0]*Vol7+DBzDx[1]*Vol8+DBzDx[2]*Vol5+DBzDx[3]*Vol6+
		DBzDx[4]*Vol3+DBzDx[5]*Vol4+DBzDx[6]*Vol1+DBzDx[7]*Vol2;

	jac[2][1]=DBzDy[0]*Vol7+DBzDy[1]*Vol8+DBzDy[2]*Vol5+DBzDy[3]*Vol6+
		DBzDy[4]*Vol3+DBzDy[5]*Vol4+DBzDy[6]*Vol1+DBzDy[7]*Vol2;

	jac[2][2]=DBzDz[0]*Vol7+DBzDz[1]*Vol8+DBzDz[2]*Vol5+DBzDz[3]*Vol6+
		DBzDz[4]*Vol3+DBzDz[5]*Vol4+DBzDz[6]*Vol1+DBzDz[7]*Vol2;

  return 0;
}