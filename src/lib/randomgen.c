/*
 * randomgen.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */


#include "randomgen.h"
#include "common_slave.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "/usr/include/time.h"
#include <float.h>
#include "slave.h"

extern __thread_local int my_id;
//#define PI 3.141592654

double uniform(double min, double max, long int *seed){
    double t = 0;
    *seed = 2045 * (*seed) + 1;
    *seed = *seed - (*seed / 1048576) * 1048576;
    t = (*seed) / 1048576.0;
    t = min + (max - min) * t;
    return t;
}

//generate uniform double rand in [a, b]
double random_fl(double a, double b) { // expects a < b, returns rand in [a, b]
	
    long int s;
    double x;
	//unsigned int seed = time(NULL);
	//srand((unsigned int)time(0));
    s = rand();
   // printf("--------%d\n",s);
    x = uniform(a, b, &s);
    return x;
}
//generate uniform int rand in [a, b]
int random_integer(double a, double b) { // expects a < b, returns rand in [a, b]
    long int s;
    double x;
	s = rand();
    x = uniform(a, b, &s);
    
   // return (int)x;   ymx original
   return (int)(x+0.5);
}

void random_inside_sphere(double *var) {
    while(true) { // on average, this will have to be run about twice
        var[0] = random_fl(-1, 1);
        var[1] = random_fl(-1, 1);
        var[2] = random_fl(-1, 1);
        if(sqr(var) < 1)
        //return var;
          break;
    }
}


//generate gauss rand min: 0, sigma: 1, this algorithm is same as boost::normal_distribution
double gaussrand(){
    static double U, V;
    static int phase = 0;
    double Z;
    //generate seed first
//    unsigned int seed = time(NULL);
//	srand(seed);
    if(phase == 0){
        U = rand() / (RAND_MAX + 1.0);
        V = rand() / (RAND_MAX + 1.0);
        Z = sqrt(-2.0 * log(U))* sin(2.0 * pi * V);
    }
    else{
        Z = sqrt(-2.0 * log(U)) * cos(2.0 * pi * V);
    }
    phase = 1 - phase;
    return Z;
}

void random_init(void){
	srand(my_id * (unsigned int)time(0));
}

void random_in_box(double *corner1, double *corner2, double *out) { // expects corner1[i] < corner2[i]
	out[0] = random_fl(corner1[0], corner2[0]);
        out[1] = random_fl(corner1[1], corner2[1]);
        out[2] = random_fl(corner1[2], corner2[2]);
}

void random_orientation(double *out) {

   out[0] = gaussrand();
   out[1] = gaussrand();
   out[2] = gaussrand();
   out[3] = gaussrand();
   // fl nrm = boost::math::abs(q);
   double sqrsum = out[0]*out[0]+out[1]*out[1]+out[2]*out[2]+out[3]*out[3];
   double nrm = sqrt(sqrsum);
    if(nrm > DBL_EPSILON) {
        //q /= nrm;
        out[0] = out[0]/nrm;
        out[1] = out[1]/nrm;
        out[2] = out[2]/nrm;
        out[3] = out[3]/nrm;
//      printf("random_orientation is updated\n");
    }
}

void torsions_randomize(double *torsions, int size) {
    for(int i=0; i<size;i++)
        torsions[i] = random_fl(-pi, pi);
}
