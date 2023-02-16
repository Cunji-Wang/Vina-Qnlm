/*
 * quaternion_slave.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */


#ifndef QUATERNION_SLAVE_H
#define QUATERNION_SLAVE_H


void quaternion_increment(double q[4], double rotation[3]);
void quaternion_normalize_approx(double q[4]);
void angle_to_quaternion(double q[4], double rotation[3]);
void quaternion_to_r3(double q[4], double mat[9]);
void quaternion_mul(double a[4], double b[4], double r[4]);
void angle_to_quaternion_two(double axis[3], double angle, double q[4]);
#endif //QUATERNION_SLAVE_H
