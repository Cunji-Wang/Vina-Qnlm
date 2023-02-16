/*
 * tree_hetero_cal.h
 * Copyright (C) 2020 yumaoxue <mxyu@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */
#ifndef TREE_HETERO_CAL_H
#define TREE_HETERO_CAL_H
#include "monte_carlo_slave.h"

void ligands_set_conf(struct output_type_s *c);
void ligands_derivative(struct change *c);
void branches_set_conf(struct flexibleHeritTree *tree, double *atoms, double *coords, double *c, int offset);
void branches_derivative(struct flexibleHeritTree *tree, double origin[3], double *coords, double *forces,
        double out[2][3], double *d, int offset, bool isrigid);
#endif //AUTODOCK_VINA_1_1_2_TREE_HETERO_CAL_H
