/*
 * mutate_slave.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */


#ifndef MUTATE_SLAVE_H
#define MUTATE_SLAVE_H

#include "tranModel.h"
#include "monte_carlo_slave.h"

void mutate_conf(struct output_type_s *candidate, int torsionsize, double amplitude);
void createnewchange(struct change *g_orgin, struct change *g_new, int size);
void createnewconf(struct output_type_s *x_orgin, struct output_type_s *x_new, int size);
void print_conf(struct output_type_s *x, int size);

#endif //MUTATE_SLAVE_H
