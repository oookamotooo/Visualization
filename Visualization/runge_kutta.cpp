#include<iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <gl/glut.h>
#include "FileManager.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "Vector3.h"
#include "define.h"

Vector3d interpolate(double is,double it,double iu, Vector3d test[8])
{
	Vector3d d;
	//printf("inter\n");
	d.x =  iu*(1-is)*(1-it)*test[4].x + iu*it*(1-is)*test[7].x + iu*is*(1-it)*test[5].x + iu*is*it*test[6].x +(1-iu)*(1-is)*(1-it)*test[0].x + (1-iu)*(1-is)*it*test[3].x + (1-iu)*is*(1-it)*test[1].x + (1-iu)*is*it*test[2].x;
	d.y =  iu*(1-is)*(1-it)*test[4].y + iu*it*(1-is)*test[7].y + iu*is*(1-it)*test[5].y + iu*is*it*test[6].y +(1-iu)*(1-is)*(1-it)*test[0].y + (1-iu)*(1-is)*it*test[3].y + (1-iu)*is*(1-it)*test[1].y + (1-iu)*is*it*test[2].y;
	d.z =  iu*(1-is)*(1-it)*test[4].z + iu*it*(1-is)*test[7].z + iu*is*(1-it)*test[5].z + iu*is*it*test[6].z +(1-iu)*(1-is)*(1-it)*test[0].z + (1-iu)*(1-is)*it*test[3].z + (1-iu)*is*(1-it)*test[1].z + (1-iu)*is*it*test[2].z;
	
	return d;
}

//単純に引数として座標(小数点含め)を与える形にしたほうが反復処理のとき楽
//与えられた座標からタイムステップ後の座標を返す関数
Vector3d streamline(Vector3d p){

	Vector3d test[8];
	int x_i, y_i, z_i;
	double x_d, y_d, z_d;

	//printf("%f %f %f\n", x, y, z);
	
	/*
	float x_i = floor(x);
	float x_d = x - x_i;
	*/

	//xyz座標からijkの値とstuの値を出す。
	
	x_i = floor(p.x);
	x_d = p.x - x_i;
	y_i = floor(p.y);
	y_d = p.y - y_i;
	z_i = floor(p.z);
	z_d = p.z - z_i;

	/*
	x_d = modf(x, &x_i);
	y_d = modf(y, &y_i);
	z_d = modf(z, &z_i);
	*/
	//それぞれに対応する線形補完値を返

	//printf("aaaaaaaaaaaaaaaaa\n");
	//printf("%d %d %d, %f %f %f\n", x_i, y_i, z_i, x_d, y_d, z_d);
	//if(x_i+1 > sizeX || y_i+1 > sizeY || z_i+1 > sizeZ) return 0;
	//for(int i=0; i < 1; i++){
	//inTotestでは-1-1-1までいれるので、+1+1+1しておく
	//printf("into %d %d %d\n", x_i+1, y_i+1, z_i+1);
	inTotest(x_i+1, y_i+1, z_i+1, test);

	return interpolate(x_d, y_d, z_d, test);
}

Vector3d runge_kutta(Vector3d point, double step){
	Vector3d k0, k1, k2, k3;
	Vector3d del;

	//printf("runge\n");
	del = streamline(point);
	k0 = step*del;
	//printf("k0:%f\n", k0);
	del = streamline(point + (k0/2.0));
	k1 = step*del;
	//printf("k1:%f\n", k1);
	del = streamline(point + (k1/2.0));
	k2 = step*del;
	//printf("k2:%f\n", k2);
	del = streamline(point + k2);
	k3 = step*del;
	//printf("k3:%f\n", k3);
	return ((k0 + 2*k1 + 2*k2 + k3)/6.0);

}