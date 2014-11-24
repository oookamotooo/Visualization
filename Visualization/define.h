#include <iostream>
#include <vector>
#include "Vector3.h"

#define CP_NUM 20
#define PI  3.14159265358979323846

using namespace std;

const int sizeX = 91;
const int sizeY = 171;
const int sizeZ = 86;
const int Size = sizeX * sizeY * sizeZ;
//const double alpha = 30.0;
const int round_num = 20; //(int)(360.0/alpha);

const int band_size = 1000;

extern vector<Vector3d> datas;
extern double s, t, u;
extern double jacob[3][3];
extern int cp_count;
extern Vector3d cp[CP_NUM];

extern double e_val[CP_NUM][3];
extern Vector3d e_vec[CP_NUM][3];

extern Vector3d round_cp[CP_NUM][round_num];



int newton1(int x, int y, int z/*, double *s, double *t, double *u*/);

int Index(int i, int j, int k);

//int streamline(double x, double y, double z, int i);

void inTotest(int i, int j, int k, Vector3d test[8]);

//double runge_kutta(double x, double y, double z, int i, double step);
int lucomp(double fjac[3][3],double fvec[3],double ds,double dt,double du);

void Geodesic(Vector3d critical[CP_NUM], Vector3d critical_round[CP_NUM][round_num], int cpnum);

Vector3d runge_kutta(Vector3d point, double step);

bool samesign_check(double a, double b);

void streamline(Vector3d particle, double len, double step);