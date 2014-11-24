#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "get_jacob.h"
#include "Vector3.h"
#include <vector>
#include "define.h"

Vector3d rotate(Vector3d housen, Vector3d vec, double alpha){
	Vector3d result;
	double theta;
	
	theta = alpha * PI / 180.0;

	result.x = ((housen.x*housen.x*(1-cos(theta)) + cos(theta))*vec.x + (housen.x*housen.y*(1-cos(theta)) - housen.z*sin(theta))*vec.y + (housen.z*housen.x*(1-cos(theta)) + housen.y*sin(theta))*vec.z);
	result.y = ((housen.x*housen.y*(1-cos(theta)) + housen.z*sin(theta))*vec.x + (housen.y*housen.y*(1-cos(theta)) + cos(theta))*vec.y + (housen.y*housen.z*(1-cos(theta)) - housen.x*sin(theta))*vec.z);
	result.z = ((housen.z*housen.x*(1-cos(theta)) - housen.y*sin(theta))*vec.x + (housen.y*housen.z*(1-cos(theta)) + housen.x*sin(theta))*vec.y + (housen.z*housen.z*(1-cos(theta)) + cos(theta))*vec.z);

	return result; 
}

void make_round(double sub_val[CP_NUM][3], Vector3d sub_vec[CP_NUM][3], int n){
	int a, b;
	double h_size;
	double alpha = 360.0/round_num;
	Vector3d housen, rotated;

	FILE *fp;
	//char name = (char)n;
	char filepath[256];

	sprintf(filepath, "All_round_point%d.txt", n);
	fp = fopen(filepath, "w");
	
	/*
	FILE *fp = fopen("round_point.txt", "w");
	*/
	//同符号ベクトル探す

	//今はすべてが++-or--+なのでこれでいいが、将来的にすべて同符号パターンが生まれた場合は書き直す

	for(int i=0; i<3; i++){
		a = i;
		if(i == 2){
			b = 0;
		}else{
			b = i+1;
		}

		if(samesign_check(sub_val[n][a], sub_val[n][b])){
			printf("samesign value = (%lf, %lf)\n", sub_val[n][a], sub_val[n][b]);
			break;
		}
	}

	//法線ベクトル求める
	housen.x = (sub_vec[n][a].y*sub_vec[n][b].z - sub_vec[n][a].z*sub_vec[n][b].y);
	housen.y = (sub_vec[n][a].z*sub_vec[n][b].x - sub_vec[n][a].x*sub_vec[n][b].z);
	housen.z = (sub_vec[n][a].x*sub_vec[n][b].y - sub_vec[n][a].y*sub_vec[n][b].x);

	h_size = sqrt(housen.x*housen.x + housen.y*housen.y + housen.z*housen.z);

	housen.x = (housen.x/h_size);
	housen.y = (housen.y/h_size);
	housen.z = (housen.z/h_size);

	printf("unit vector make %lf %lf %lf -> %lf\n", housen.x, housen.y, housen.z, sqrt(housen.x*housen.x + housen.y*housen.y + housen.z*housen.z));
	//法線ベクトルを軸に回転行列かける
	rotated.x = sub_vec[n][a].x;	rotated.y = sub_vec[n][a].y;	rotated.z = sub_vec[n][a].z;
	fprintf(fp, "%lf %lf %lf\n", rotated.x, rotated.y, rotated.z);

	for(int i=0;i < (int)(360.0/alpha); i++){
		rotated = rotate(housen, rotated, alpha);
		fprintf(fp, "%lf %lf %lf\n", rotated.x, rotated.y, rotated.z);
	}

	printf("aaa\n");

	fclose(fp);
	return;
}