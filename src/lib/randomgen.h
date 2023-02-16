/*
 * randomgen.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef RANDOMGEN_H
#define RANDOMGEN_H

#include <stdio.h>
#include <stdlib.h>

double gaussrand();
void random_inside_sphere(double *var);
int random_integer(double a, double b);
double random_fl(double a, double b);
void random_init(void);
void random_in_box(double *corner1, double *corner2, double *out);
void random_orientation(double *out);
void torsions_randomize(double *torsions, int size);
#endif //RANDOMGEN_H
