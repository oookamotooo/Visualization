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

void streamline(Vector3d particle, double len, double step){
	Vector3d d1,d2;
	double d_size;

	d1 = particle;

	while(true){
			if(d1.x > sizeX-5 || d1.x < 5 || d1.y > sizeY-5 || d1.y < 5 || d1.z > sizeZ-5 || d1.z < 5 ) break;
			//printf("%d %d %d\n", dx1, dy1, dz1);
			d2 = runge_kutta(d1, step);

			d_size = d2.dist(Vector3d());
				
			d2 = (d2/d_size)*len;

			glBegin(GL_LINES);
			glVertex3f(d1.x-(double)sizeX/2.0, d1.y-(double)sizeY/2.0, d1.z-(double)sizeZ/2.0);
			glVertex3f((d1.x-(double)sizeX/2.0)+d2.x, (d1.y-(double)sizeY/2.0)+d2.y, (d1.z-(double)sizeZ/2.0)+d2.z);
			glEnd();

			d1 += d2;
		}
}