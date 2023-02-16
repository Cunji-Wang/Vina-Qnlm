/*
 * common_slave.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "common_slave.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "simd.h"

double vecMultiVec(double a[3], double b[3]){
    double r;
    r = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return r;
}

double sqr(double x[3]){
    double r = vecMultiVec(x, x);
    return r;
}

double sqr_double(double x){
    return x*x;
}

double vec_distance_sqr(double *a, double *b){
	return sqr_double(a[0]-b[0]) + sqr_double(a[1]-b[1]) + sqr_double(a[2]-b[2]);
}

void vecminus(double a[3], double b[3], double r[3]){     //???

    r[0] = a[0] - b[0];
    r[1] = a[1] - b[1];
    r[2] = a[2] - b[2];
}

void vecmultiply(double a[3], double s, double r[3]){   //??

    r[0] = a[0]*s;
    r[1] = a[1]*s;
    r[2] = a[2]*s;
}

double vec_dotproduct(double a[3], double b[3]){      
    //return data[0] * v[0] + data[1] * v[1] + data[2] * v[2];
    double tmp = a[0]*b[0]+a[1]*b[1] + a[2]*b[2];
    return tmp;
}

void vecminusequal(double *a, double *b){
    a[0] -= b[0];
    a[1] -= b[1];
    a[2] -= b[2];
}

void vecplusequal(double *a, double *b){
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
}

void vecplus(double a[3], double b[3], double r[3]){

    r[0] = a[0] + b[0];
    r[1] = a[1] + b[1];
    r[2] = a[2] + b[2];
}

double *vecmultiplyequal(double *a, double b){
    a[0] *= b;
    a[1] *= b;
    a[2] *= b;
    return a;
}

bool not_max(double x){
    return (x<0.1*DBL_MAX);
}

double norm_sqr(double a[3]) {
    return sqr(a);
}
double norm(double a[3]) {
    return sqrt(norm_sqr(a));
}

void cross_product(double a[3], double b[3], double r[3]) {
 
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
}

void elementwise_product(double a[3], double b[3], double r[3]) {
    r[0] = a[0] * b[0];
    r[1] = a[1] * b[1];
    r[2] = a[2] * b[2];
}

void normalize_angle(double *x) { // subtract or add enough 2*pi's to make x be in [-pi, pi]
    if(*x > 3*pi) { // very large
        double n = ( *x - pi) / (2*pi); // how many 2*pi's do you want to subtract?
        //*x -= 2*pi*std::ceil(n); // ceil can be very slow, but this should not be called often
//       int tmp = (int)(n+0.5);    //ymx originl
//       *x -= 2*pi*tmp;
	*x-=2*pi*ceil(n);
        normalize_angle(x);
    }
    else if(*x < -3*pi) { // very small
        double n = (-*x - pi) / (2*pi); // how many 2*pi's do you want to add?
        //x += 2*pi*std::ceil(n); // ceil can be very slow, but this should not be called often
//        int tmp = (int)(n+0.5);
//        *x += 2*pi*tmp;   //ymx originl
        *x += 2*pi*ceil(n);
        normalize_angle(x);
    }
    else if(*x > pi) { // in (   pi, 3*pi]
        *x -= 2*pi;
    }
    else if(*x < -pi) { // in [-3*pi,  -pi)
        *x += 2*pi;
    }
// in [-pi, pi]
}

double normalized_angle(double x) {
    normalize_angle(&x);
    return x;
}

// mat
double *matfuncCall(int i, int j, double mat[9]) {
    //if((i<3)&&(j<3))
     return &mat[i + 3*j];
}

void matMultiply(double v[3], double mat[9], double r[3]) {
   // return vec(data[0]*v[0] + data[3]*v[1] + data[6]*v[2], data[1]*v[0] + data[4]*v[1] + data[7]*v[2],
   //             data[2]*v[0] + data[5]*v[1] + data[8]*v[2]);
    r[0] = mat[0]*v[0] + mat[3]*v[1] + mat[6]*v[2];
    r[1] = mat[1]*v[0] + mat[4]*v[1] + mat[7]*v[2];
    r[2] = mat[2]*v[0] + mat[5]*v[1] + mat[8]*v[2];
}

void matMultiplyAssign(double s, double mat[9]) {
    //for(int i=0; i<9; i++)
        //mat[i] *= s;
   //double tmp[8];
   doublev8 mat_v;
   simd_load(mat_v, &(mat[0]));
   mat_v *= s;
   simd_store(mat_v,&mat[0]);

   mat[8] *= s;
/*   mat[0] *= s;
   mat[1] *= s;
   mat[2] *= s;
   mat[3] *= s;
   mat[4] *= s;
   mat[5] *= s;
   mat[6] *= s;
   mat[7] *= s;
   mat[8] *= s;
*/
  //for(int i=0;i<8;i++){if(tmp[i]!=mat[i])printf("tmp:%lf mat:%lf\n",tmp[i],mat[i]);}
}

void curl(double *e, double *deriv, double v){

    if(*e>0 && not_max(v)){
        double tmp = (v < DBL_EPSILON)? 0 : (v/(v+*e));
        *e *= tmp;
//		if(my_id==0)
//			printf("[slave] curl e is updated, e is %lf, tmp is %lf\n", *e, tmp);
        double val = sqr_double(tmp);
        vecmultiplyequal(deriv, val);
    }
}

void curl_two(double *e, double v) {
    if(*e > 0 && not_max(v)) {
        double tmp = (v < DBL_EPSILON) ? 0 : (v / (v + *e));
        *e *= tmp;
    }
}
