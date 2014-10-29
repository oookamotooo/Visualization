
#include<iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <gl/glut.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include "Vector3.h"
#include "newton.h"

//using namespace std;

#define x_start 0
#define y_start 0
#define z_start 0

//#define CP_NUM 20

const int sizeX = 91;
const int sizeY = 171;
const int sizeZ = 86;
const int Size = sizeX * sizeY * sizeZ;

extern vector<Vector3d> datas;
extern double s, t, u;
extern int cp_count;
extern Vector3d about_cp[CP_NUM];

Vector3i size(91, 171, 86);

int Index(int i, int j, int k)
{
	return k*sizeX*sizeY + j*sizeX + i;
	//return i*sizeY*sizeZ + j*sizeZ + k;
	//return k*sizeZ*sizeY + j*sizeZ + i;
}

bool samesign_check(double a, double b){
	if(a*b >= 0){
		return true;
	}else{
		return false;
	}
}

bool isCp_check(double test_x[8],double test_y[8],double test_z[8])
{
	for (int i = 0; i < 7; i++)
	for (int j = i+1; j < 8; j++) 
	{
		if ( samesign_check(test_x[i], test_x[j]) || 
			 samesign_check(test_y[i], test_y[j]) || 
			 samesign_check(test_z[i], test_z[j]))
			continue;
		
		
		//printf(" point1: %f %f %f point2: %f %f %f \n", test_x[i], test_y[i], test_z[i], test_x[j], test_y[j], test_z[j]);
		return true;
	}
	return false;
}

float phi = 0, theta = 0, fov = 90, radius = 20;

bool readBinary( string filename, int sizeX, int sizeY, int sizeZ, float *data)
{
	FILE *fp = fopen(filename.c_str(), "rb");
	if( fp == NULL )
	{
		cout << "can not open file " + filename << endl;
		fclose(fp);
		return false;
	}

	fseek(fp, 0, SEEK_END);

	fpos_t s;
	fgetpos(fp, &s);
	cout << s / sizeof(float) << endl;

	fseek(fp, 0, SEEK_SET);
	int num = sizeX*sizeY*sizeZ;
	size_t size = fread(data, sizeof(float), num*3, fp);
	cout << size << ", " << num*3 << endl;
	if( size < num)
	{
		cout << "指定した領域分のデータがありません" << endl;
		fclose(fp);
		return false;
	}

	fclose(fp);
	return true;
}


void normalize()
{
	for(int i=0; i<Size; i++)
		datas[i].normalize();
}

void readText(){
	int index;

	datas.resize(sizeX*sizeY*sizeZ);


	FILE *fp = fopen("bfield_near_cusp.txt", "r");

	//以下要書き換え?
	
	for(int k=0;k<sizeZ;k++){
		for(int j=0;j<sizeY;j++){
			for(int i=0;i<sizeX;i++){
				index = Index(i, j, k);

				fscanf(fp,"%lf %lf %lf", &datas[index].x, &datas[index].y, &datas[index].z);
				/*
				fscanf(fp, "%lf",&datas[index].x);
				fscanf(fp, "%lf",&datas[index].y);
				fscanf(fp, "%lf",&datas[index].z);
				*/
			}
		}
	}
	/*
	for(int i=0;i<sizeX*sizeY*sizeZ;i++){
		fscanf(fp,"%lf %lf %lf", &datas[i].x, &datas[i].y, &datas[i].z);
	}
	*/

	//normalize();

	return;
}


void inTotest(int i, int j, int k, double test_x[8], double test_y[8], double test_z[8]){
	int index;

	index = Index(i, j, k);
	test_x[0] = datas[index].x;	test_y[0] = datas[index].y;	test_z[0] = datas[index].z;

	index = Index(i, j, k-1);
	test_x[1] = datas[index].x;	test_y[1] = datas[index].y;	test_z[1] = datas[index].z;

	index = Index(i, j-1, k);
	test_x[2] = datas[index].x;	test_y[2] = datas[index].y;	test_z[2] = datas[index].z;
	
	index = Index(i-1, j, k);
	test_x[3] = datas[index].x;	test_y[3] = datas[index].y;	test_z[3] = datas[index].z;
	
	index = Index(i, j-1, k-1);
	test_x[4] = datas[index].x;	test_y[4] = datas[index].y;	test_z[4] = datas[index].z;
	
	index = Index(i-1, j-1, k);
	test_x[5] = datas[index].x;	test_y[5] = datas[index].y;	test_z[5] = datas[index].z;
	
	index = Index(i-1, j, k-1);
	test_x[6] = datas[index].x;	test_y[6] = datas[index].y;	test_z[6] = datas[index].z;

	index = Index(i-1, j-1, k-1);
	test_x[7] = datas[index].x;	test_y[7] = datas[index].y;	test_z[7] = datas[index].z;

	return;
}



void FindCP1()
{
	double test_x[8], test_y[8], test_z[8];

	cout << " Find CP start ..." << endl;
	int numOfCP = 0;
	for(int k=z_start+1; k < z_start+sizeZ; k++)
		for(int j=y_start+1; j < y_start+sizeY; j++)
			for(int i=x_start+1; i < x_start+sizeX; i++){
				/*
				for(int l = 0; l<8; l++)
				{
					int index = Index(i + (l&1), j + ((l>>1)&1), k + ((l>>2)&1) );
					test_x[l] = datas[index].x;
					test_y[l] = datas[index].y;
					test_z[l] = datas[index].z;
				}*/

				inTotest(i, j, k, test_x, test_y, test_z);
				
				if( isCp_check(test_x, test_y, test_z) ){
					cout << i << "," << j << "," << k << endl;
					numOfCP++;
				}

			}
	cout << numOfCP << endl;
	return;
}

void FindCP2(){

	double test_x[8], test_y[8], test_z[8];

	cp_count = 0; 

	for(int k=z_start+1; k<sizeZ; k++){
		for(int j=y_start+1; j<sizeY; j++){
			for(int i=x_start+1; i<sizeX; i++){

				inTotest(i, j, k, test_x, test_y, test_z);
				//cout << i << "," << j << "," << k << endl;
				if(isCp_check(test_x, test_y, test_z)){
					newton1(i, j, k);
				}

			}
		}
	}

	printf("cp count is %d\n", cp_count);

	return;

}

void display()
{
  Vector3d plot;	
  glClear(GL_COLOR_BUFFER_BIT);

  //図形の描画 
  glColor3d(0.0, 0.0, 0.0);

  cout << "diso" << endl;
  //int sx = 110, sy=55, sz=30;
  
  
  for(int k = z_start; k < z_start + sizeZ; k++)
	  for(int j = y_start; j < y_start + sizeY; j++)
		  for(int i = x_start; i < x_start + sizeX; i++){
			  if(k%10 == 0){
				int index = Index(i, j, k);
				float x = i - (double)sizeX/2;
				float y = j - (double)sizeY/2;
				float z = k - (double)sizeZ/2;
				glBegin(GL_LINES);
				glVertex3f(x, y, z);
				glVertex3f(x+datas[index].x,y+datas[index].y, z+datas[index].z);
				glEnd();
			  }
			}

	glColor3d(1.0, 0.0, 0.0);
	for(int i=0; i<CP_NUM; i++){
		printf("i = %d\n", i);
		if(about_cp[i].x == 0 && about_cp[i].y == 0 && about_cp[i].z == 0) break;

		plot.x = about_cp[i].x - (double)sizeX/2;	plot.y = about_cp[i].y - (double)sizeY/2;	plot.z = about_cp[i].z - (double)sizeZ/2;

		glBegin(GL_LINES);
		glVertex3f(plot.x-1, plot.y, plot.z);
		glVertex3f(plot.x+1, plot.y, plot.z);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(plot.x, plot.y-1, plot.z);
		glVertex3f(plot.x, plot.y+1, plot.z);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(plot.x, plot.y, plot.z-1);
		glVertex3f(plot.x, plot.y, plot.z+1);
		glEnd();
	}
  /*
  int k=0;
  for(int i = 0; i < sizeY; i++)
	  for(int j = 0; j < sizeX; j++){
		  float x = j-sizeX/2;
		  float y = i-sizeY/2;
		  //float z = 0;
		  glBegin(GL_LINES);
		  glVertex2d(y,x);
		  glVertex2d(y+datas[k].y, x+datas[k].x);
		  glEnd();
		  k++;
	  }
	*/
	glFlush();
}

int width, height;
void look()
{
	
  glViewport(0, 0, width, height);

  glLoadIdentity();
	  //glOrtho(-w / 200.0, w / 200.0, -h / 200.0, h / 200.0, -1.0, 1.0);
  cout << "look" << endl;
  gluPerspective(fov, (double)width / (double)height, 1.0, 500.0);
  
  float px = radius * cos(phi) * cos(theta);
  float py = radius * sin(phi);
  float pz = radius * cos(phi) * sin(theta);
  gluLookAt(px, py, pz,
	   0,  0, 0, 
	   0.0, 1.0, 0.0);

  glutPostRedisplay();
  
}

void resize(int w, int h)
{
	width = w;
	height = h;

  look();
}

void init(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
}


void keyboard(unsigned char key, int x, int y)
{
	cout << "keyborad" << endl;
  switch (key) {
  case 'w':
    phi += 1.0 / 180.0 * M_PI;
	break;
  case 's':
	  phi -= 1.0 / 180.0 * M_PI;
	  break;
  case 'a':
	  theta += 1.0 / 180.0 * M_PI;
	  break;
  case 'd':
	  theta -= 1.0 /180.0 *  M_PI;
	  break;
  case 'q':
	  //fov+=1;
	  radius += 1.0;
	  break;
  case 'e':
	  radius = max(1.0, radius-1.0);
	  break;
  default:
    break;
  }
  look();
}

void toVector(float *data)
{
	for(int i=0; i<Size; i++)
	{
		datas[i].x = data[i];
	}
	for(int i=0; i<Size; i++)
	{
		datas[i].y = data[i+Size];
	}
	for(int i=0; i<Size; i++)
	{
		datas[i].z = data[i+2*Size];
	}
}

void toVector2(float *data)
{
	for(int i=0; i<sizeX; i++)
	{
		for(int j=0; j<sizeY; j++)
		{
			for(int k=0; k<sizeZ; k++)
			{
				int index = Index(i,j,k);
				datas[index].x = data[3*index];
				datas[index].y = data[3*index+1];
				datas[index].z = data[3*index+2];
			}
		}
	}
}

void debugOutput()
{
	FILE *fp = fopen("test2.txt", "w");
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			for(int k=0;k<275;k++)
			{
				int index = Index(i,j,k);
				//cout << k << endl;
				fprintf(fp, "%.20lf, %.20lf, %.20lf\n",datas[index].x, datas[index].y,datas[index].z );
				//	cout << datas[index].x << "," << datas[index].y << "," << datas[index].z << endl;
			}
		}
	}
	cout << "書き出し完了"<<endl;
}

int main(int argc, char *argv[])
{
	Vector3i s;
	Vector3i a = Vector3i()+s;

	cout << size << endl;
	readText();
	

	//読み込み部分
	/*if( readBinary("fort.9005", sizeX, sizeY, sizeZ, data) )
	{
		cout << "読み込み成功" << endl;
	} 
	
	datas.resize(Size);
	toVector2(data);
	normalize();
	cout << "Vector型への変換完了" << endl;
	*/

		
	/*
	float epsilon = 0.00000000000000001;
	for(int i=0; i<sizeX; i++)
		for(int j=0; j<sizeY; j++)
			for(int k=0; k<sizeZ; k++)
			{
				index = Index(i,j,k);

				if( datas[index].x == 0 &&
					datas[index].y == 0 &&
					datas[index].z == 0 )
				{
					cout << i << "," << j << "," << k << endl;
				}
			}

			
		*/
	//about_cp配列を初期化
	for(int i=0; i<CP_NUM; i++){
		about_cp[i].x = 0;	about_cp[i].y = 0;	about_cp[i].z = 0;
 	}

	//FindCP1();

	FindCP2();
	
	system("pause");
	
	glutInitWindowSize(600, 600);
	glutInit(&argc, argv);
	//glClearColor(0.0, 0.0, 1.0, 1.0);
	//3d
	glutInitDisplayMode(GLUT_RGBA);
	//3d
	glutCreateWindow(argv[0]);
	glutDisplayFunc(display);
	//3d
	glutReshapeFunc(resize);
	//3d
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
	
	return 0;
}

/*
void DrawSphere(const Vector3d cen, const double rad){
		glPushMatrix();
		glTranslatef(cen.x, cen.y, cen.z);
		glRotated(90, 1.0, 0.0, 0.0);
		glutSolidSphere(rad, 20, 10);
		glPopMatrix();
	}
	*/