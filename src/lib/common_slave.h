/*
 * common_slave.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef COMMON_SLAVE_H
#define COMMON_SLAVE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

double pi = 3.1415926535897931;
double vecMultiVec(double a[3], double b[3]);
double sqr(double x[3]);
double sqr_double(double x);
void vecminus(double a[3], double b[3], double r[3]);
void vecmultiply(double *a, double s, double *r);
double vec_dotproduct(double a[3], double b[3]);
void vecminusequal(double *a, double *b);
void vecplusequal(double *a, double *b);
void vecplus(double a[3], double b[3], double r[3]);
double *vecmultiplyequal(double *a, double b);
double vec_distance_sqr(double *a, double *b);
bool not_max(double x);
double norm_sqr(double a[3]);
double norm(double a[3]);
void cross_product(double a[3], double b[3], double r[3]);
void elementwise_product(double a[3], double b[3], double r[3]);
void normalize_angle(double *x);
double normalized_angle(double x);
double *matfuncCall(int i, int j, double mat[9]);
void matMultiply(double v[3], double mat[9], double r[3]);
void matMultiplyAssign(double s, double mat[9]);
void curl(double *e, double *deriv, double v);
void curl_two(double *e, double v);
#endif //COMMON_SLAVE_H
