/*
 * model_slave.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */


#ifndef MODEL_SLAVE_H
#define MODEL_SLAVE_H

#include "monte_carlo_slave.h"
#include "tranModel.h"

int num_movable_atoms(void);
double gyration_radius(void);
double model_eval_deriv  (double v[3], struct output_type_s *c, struct change *g);
void model_set(struct output_type_s *c);
void get_heavy_atom_movable_coords(double *coords);

#endif //MODEL_SLAVE_H
