/*
 * bfgs_slave.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */
#ifndef BFGS_SLAVE_H
#define BFGS_SLAVE_H

double bfgs(struct output_type_s *candidate, struct change *g, double v[3], struct paras *para);

#endif //BFGS_SLAVE_H
