
#include<iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <gl/glut.h>
#include "FileManager.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "Vector3.h"
#include "newton.h"
#include "get_jacob.h"
#include "define.h"

//using namespace std;

#define x_start 0
#define y_start 0
#define z_start 0

//#define CP_NUM 20

/*
const int sizeX = 91;
const int sizeY = 171;
const int sizeZ = 86;
const int Size = sizeX * sizeY * sizeZ;
*/

vector<Vector3d> datas;
double s, t, u;
double jacob[3][3];
int cp_count;
Vector3d cp[CP_NUM];
Vector3d round_cp[CP_NUM][round_num];

//double runge_kutta(double x, double y, double z, int i, double step);

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

bool isCp_check(Vector3d test[8])
{
	for (int i = 0; i < 7; i++)
	for (int j = i+1; j < 8; j++) 
	{
		if ( samesign_check(test[i].x, test[j].x) || 
			 samesign_check(test[i].y, test[j].y) || 
			 samesign_check(test[i].z, test[j].z))
			continue;
		
		
		//printf(" point1: %f %f %f point2: %f %f %f \n", test_x[i], test_y[i], test_z[i], test_x[j], test_y[j], test_z[j]);
		return true;
	}
	return false;
}

float phi = 0, theta = 0, fov = 90, radius = 20;
/*
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
*/

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
	fclose(fp);
	return;
}

void read_cp(){
	FILE *fp = fopen("cp_sub.txt", "r");

	for(int i=0; i < 4;/*CP_NUM;*/ i++){
		fscanf(fp,"%lf %lf %lf", &cp[i].x, &cp[i].y, &cp[i].z);
	}

	printf("Read CP success!");

	fclose(fp);
	return;
}

void read_round(int i){
	char filepath[256];
	sprintf(filepath, "round_point%d.txt", i);
	FILE *fp = fopen(filepath, "r");
	
	for(int s=0; s<round_num; s++){
		fscanf(fp, "%lf %lf %lf", &round_cp[i][s].x, &round_cp[i][s].y, &round_cp[i][s].z);
	}

	fclose(fp);

	return;
}

void inTotest(int i, int j, int k, Vector3d test[8]){
	int index;

	index = Index(i-1, j-1, k-1);
	test[0].x = datas[index].x;	test[0].y = datas[index].y;	test[0].z = datas[index].z;

	index = Index(i, j-1, k-1);
	test[1].x = datas[index].x;	test[1].y = datas[index].y;	test[1].z = datas[index].z;

	index = Index(i, j, k-1);
	test[2].x = datas[index].x;	test[2].y = datas[index].y;	test[2].z = datas[index].z;
	
	index = Index(i-1, j, k-1);
	test[3].x = datas[index].x;	test[3].y = datas[index].y;	test[3].z = datas[index].z;
	
	index = Index(i-1, j-1, k);
	test[4].x = datas[index].x;	test[4].y = datas[index].y;	test[4].z = datas[index].z;
	
	index = Index(i, j-1, k);
	test[5].x = datas[index].x;	test[5].y = datas[index].y;	test[5].z = datas[index].z;
	
	index = Index(i, j, k);
	test[6].x = datas[index].x;	test[6].y = datas[index].y;	test[6].z = datas[index].z;

	index = Index(i-1, j, k);
	test[7].x = datas[index].x;	test[7].y = datas[index].y;	test[7].z = datas[index].z;

	return;
}


/*
void FindCP1()
{
	double test_x[8], test_y[8], test_z[8];

	cout << " Find CP start ..." << endl;
	int numOfCP = 0;
	for(int k=z_start+1; k < z_start+sizeZ; k++)
		for(int j=y_start+1; j < y_start+sizeY; j++)
			for(int i=x_start+1; i < x_start+sizeX; i++){
				
				inTotest(i, j, k, test_x, test_y, test_z);
				
				if( isCp_check(test_x, test_y, test_z) ){
					cout << i << "," << j << "," << k << endl;
					numOfCP++;
				}

			}
	cout << numOfCP << endl;
	return;
}
*/

void FindCP2(){

	Vector3d test[8];

	cp_count = 0; 

	for(int k=z_start+1; k<sizeZ; k++){
		for(int j=y_start+1; j<sizeY; j++){
			for(int i=x_start+1; i<sizeX; i++){

				inTotest(i, j, k, test);
				//cout << i << "," << j << "," << k << endl;
				if(isCp_check(test)){
					newton1(i, j, k);
					printf("Not same sign find (%d / %d)", k*sizeX*sizeY + j*sizeX + i, sizeX*sizeY*sizeZ);
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
  
  /*
  for(int k = z_start; k < z_start + sizeZ; k++)
	  for(int j = y_start; j < y_start + sizeY; j++)
		  for(int i = x_start; i < x_start + sizeX; i++){
			  //if(k%10 == 0){
				int index = Index(i, j, k);
				float x = i - (double)sizeX/2;
				float y = j - (double)sizeY/2;
				float z = k - (double)sizeZ/2;
				glBegin(GL_LINES);
				glVertex3f(x, y, z);
				glVertex3f(x+datas[index].x,y+datas[index].y, z+datas[index].z);
				glEnd();
			  //}
			}*/
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,-(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f(-(double)sizeX/2.0,-(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f(-(double)sizeX/2.0,(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f((double)sizeX/2.0,-(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,-(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f((double)sizeX/2.0,(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,-(double)sizeY/2.0,(double)sizeZ/2.0);
  glVertex3f(-(double)sizeX/2.0,(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f((double)sizeX/2.0,-(double)sizeY/2.0,(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,-(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f(-(double)sizeX/2.0,(double)sizeY/2.0,-(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f((double)sizeX/2.0,-(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,(double)sizeY/2.0,-(double)sizeZ/2.0);
  glEnd;

  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,(double)sizeY/2.0,(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,-(double)sizeY/2.0,(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,-(double)sizeY/2.0,(double)sizeZ/2.0);
  glEnd;

  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,(double)sizeY/2.0,-(double)sizeZ/2.0);
  glEnd;
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,-(double)sizeY/2.0,-(double)sizeZ/2.0);
  glVertex3f((double)sizeX/2.0,-(double)sizeY/2.0,-(double)sizeZ/2.0);
  glEnd;

  glColor3d(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(-(double)sizeX/2.0,0,0);
  glVertex3f((double)sizeX/2.0,0,0);
  glEnd;
  glColor3d(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3f(0,-(double)sizeY/2.0,0);
  glVertex3f(0,(double)sizeY/2.0,0);
  glEnd;
  glColor3d(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0,0,-(double)sizeZ/2.0);
  glVertex3f(0,0,(double)sizeZ/2.0);
  glEnd;

  Vector3d earth;
  earth = Vector3d(60.0-(double)sizeX/2.0, 82.0-(double)sizeY/2.0, -13.0-(double)sizeZ/2.0);

  glColor3d(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(earth.x-1, earth.y, earth.z);
	glVertex3f(earth.x+1, earth.y, earth.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(earth.x, earth.y-1, earth.z);
	glVertex3f(earth.x, earth.y+1, earth.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(earth.x, earth.y, earth.z-1);
	glVertex3f(earth.x, earth.y, earth.z+1);
	glEnd();

	glColor3d(1.0, 0.0, 0.0);
	for(int i=0; i<CP_NUM; i++){
		printf("i = %d\n", i);
		if(cp[i].x == 0 && cp[i].y == 0 && cp[i].z == 0) break;

		//printf("%f %f %f\n", cp[i].x, cp[i].y, cp[i].z);

		plot.x = cp[i].x - (double)sizeX/2;	plot.y = cp[i].y - (double)sizeY/2;	plot.z = cp[i].z - (double)sizeZ/2;

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


	Vector3d s_line;
	int remake;
	double step;
	double len=0.1, d_size;
	double first_len = 5.0;
	//double tmp;

	//printf("%f\n", cp[0].x);
	//printf("d-series = %f %f %f\n", dx1, dy1, dz1);
	//streamline
	for(int t=1; t<2; t++){
		if(t == 2){
			step = 0.1;
		}else{
			step = -0.1;
		}
	//int a=-40, b=40; 
	for(int s = 0; s<round_num; s++){
		s_line = (cp[t] + round_cp[t][s]*(first_len));
		glColor3d(0.0, 0.0, 1.0);
		streamline(s_line, len, step);
	}
	}

	step = -0.1;
	Vector3d e_start[4];
	e_start[0].x = earth.x+(double)sizeX/2.0;	e_start[0].y = earth.y+(double)sizeY/2.0;	e_start[0].z = earth.z +(double)sizeZ/2.0+ 20;
	e_start[1].x = earth.x+(double)sizeX/2.0 + 3;	e_start[1].y = earth.y+(double)sizeY/2.0;	e_start[1].z = earth.z +(double)sizeZ/2.0+ 20;
	e_start[2].x = earth.x+(double)sizeX/2.0 - 4;	e_start[2].y = earth.y+(double)sizeY/2.0;	e_start[2].z = earth.z +(double)sizeZ/2.0+ 20;
	e_start[3].x = earth.x+(double)sizeX/2.0 - 6;	e_start[3].y = earth.y+(double)sizeY/2.0;	e_start[3].z = earth.z +(double)sizeZ/2.0+ 20;
	glColor3d(0.0, 1.0, 0.0);
	for(int l=0; l<4; l++){
		streamline(e_start[l], len, step);
	}

	/*
	glBegin(GL_LINES);
	glVertex3f(e_start.x-1, e_start.y, e_start.z);
	glVertex3f(e_start.x+1, e_start.y, e_start.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(e_start.x, e_start.y-1, e_start.z);
	glVertex3f(e_start.x, e_start.y+1, e_start.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(e_start.x, e_start.y, e_start.z-1);
	glVertex3f(e_start.x, e_start.y, e_start.z+1);
	glEnd();
	*/
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
/*
void output_jacob(double x, double y, double z, double jacob[3][3]){
	FILE *fp = fopen("jacob.txt", "w");

	for(int a=0; a < CP_NUM; a++){
		if(cp[a].x == 0 && cp[a].y == 0 && cp[a].z == 0) break;

		get_jacob(cp[a].x, cp[a].y, cp[a].z, jacob);
		fprintf(fp)
	}
}*/

/*
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
*/

int main(int argc, char *argv[])
{
	
	cout << size << endl;
	readText();
	printf("readText done.\n");
	

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

	//cp配列を初期化
	for(int i=0; i<CP_NUM; i++){
		cp[i].x = 0;	cp[i].y = 0;	cp[i].z = 0;
 	}

	//FindCP1();

	/*
	FindCP2();
	printf("Find CP done.\n");

	for(int i=0; i<CP_NUM; i++){
		printf("cp%d = (%.10lf, %.10lf, %.10lf)\n", i, cp[i].x, cp[i].y, cp[i].z);
	}*/

	
	read_cp();
	printf("read_cp done.\n");

	for(int i=0; i<4; i++){
		read_round(i);
	}
	printf("read_round done.\n");
	

	/*
	for(int i=1; i<2; i++){
		Geodesic(cp, round_cp, i);
	}
	
	printf("Geodesic done.\n");
	*/
	
	/*
	vector<Jacobian> jacobians;
	FileManager::ReadJacobianData("p_eigen_out.txt", jacobians);
	for(auto it = jacobians.begin(); it != jacobians.end(); it++)
		cout << (*it) << endl;
	system("pause");
	*/
	
	/*
	//Jacobian導出
	FILE *fp2 = fopen("new_jacob2.txt", "w");

	for(int a=0; a < CP_NUM; a++){
		if(cp[a].x == 0 && cp[a].y == 0 && cp[a].z == 0) break;

		get_jacob(cp[a].x, cp[a].y, cp[a].z, jacob);
		
		fprintf(fp2, "%lf %lf %lf\n", cp[a].x, cp[a].y, cp[a].z);
		for(int j=0;j<3;j++){
			for(int i=0;i<3;i++){
				fprintf(fp2, "%.10lf ,",jacob[i][j]);
			}
		}
		fprintf(fp2, "\n");
	}
	fclose(fp2);
	printf("Get Jacob done.\n");
	*/

	//streamline 描画
	
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