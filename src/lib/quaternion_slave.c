/*
 * quaternion_slave.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "common_slave.h"
#include "quaternion_slave.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



void angle_to_quaternion_two(double axis[3], double angle, double q[4]) { // axis is assumed to be a unit vector

    normalize_angle(&angle); // this is probably only necessary if angles can be very big
    double c = cos(angle/2);
    double s = sin(angle/2);
    //return qt(c, s*axis[0], s*axis[1], s*axis[2]);
    q[0] = c;
    q[1] = s*axis[0];
    q[2] = s*axis[1];
    q[3] = s*axis[2];
}

void angle_to_quaternion(double q[4], double rotation[3]) {
    //fl angle = tvmet::norm2(rotation);
    double angle = norm(rotation);
    if(angle > DBL_EPSILON) {
        //vec axis;
        //axis = rotation / angle;
        //vec axis = (1/angle) * rotation;
        double axis[3];
        axis[0] = (1/angle) * rotation[0];
        axis[1] = (1/angle) * rotation[1];
        axis[2] = (1/angle) * rotation[2];
        angle_to_quaternion_two(axis, angle, q);
    }
    else{
        //return qt_identity;
        q[0] = 1;
        q[1] = 0;
        q[2] = 0;
        q[3] = 0;
    }

}

double quaternion_norm_sqr(double q[4]) { // equivalent to sqr(boost::math::abs(const qt&))
    return sqr_double(q[0]) + sqr_double(q[1]) + sqr_double(q[2]) + sqr_double(q[3]);
}

void quaternion_normalize_approx(double q[4]) {

    double s = quaternion_norm_sqr(q);
    double tolerance = 1e-6;
  
    if(abs(s - 1) < tolerance){
        ; // most likely scenario}
}
    else {
        double a = sqrt(s);
        //q *= 1/a;

        q[0] = q[0] * (1/a);
        q[1] = q[1] * (1/a);
        q[2] = q[2] * (1/a);
        q[3] = q[3] * (1/a);
        // assert(quaternion_is_normalized(q));
    }
}

/*
 *             template<typename X>
            BOOST_CXX14_CONSTEXPR quaternion<T> &        operator *= (quaternion<X> const & rhs)
            {
                T    ar = static_cast<T>(rhs.R_component_1());
                T    br = static_cast<T>(rhs.R_component_2());
                T    cr = static_cast<T>(rhs.R_component_3());
                T    dr = static_cast<T>(rhs.R_component_4());

                quaternion<T> result(a*ar - b*br - c*cr - d*dr, a*br + b*ar + c*dr - d*cr,
                        a*cr - b*dr + c*ar + d*br, a*dr + b*cr - c*br + d*ar);
                swap(result);
                return(*this);
            }

 */

void quaternion_mul(double a[4], double b[4], double r[4]){
    r[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    r[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    r[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    r[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

void quaternion_increment(double q[4], double rotation[3]) {

    //q = angle_to_quaternion(rotation) * q;
    double tmp[4];
    angle_to_quaternion(tmp, rotation);

    double tmp_q[4];
    tmp_q[0] = q[0];
    tmp_q[1] = q[1];
    tmp_q[2] = q[2];
    tmp_q[3] = q[3];

    quaternion_mul(tmp, tmp_q, q);

    quaternion_normalize_approx(q); // normalization added in 1.1.2
    //quaternion_normalize(q); // normalization added in 1.1.2
}


void quaternion_to_r3(double q[4], double mat[9]) {

    double a = q[0];
    double b = q[1];
    double c = q[2];
    double d = q[3];

    double aa = a*a;
    double ab = a*b;
    double ac = a*c;
    double ad = a*d;
    double bb = b*b;
    double bc = b*c;
    double bd = b*d;
    double cc = c*c;
    double cd = c*d;
    double dd = d*d;

    // from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
    //tmp(0, 0) = (aa+bb-cc-dd);
    double *tmp = matfuncCall(0, 0, mat);
    *tmp = (aa+bb-cc-dd);
    tmp = matfuncCall(0, 1, mat);
    *tmp = 2*(-ad+bc);
    tmp = matfuncCall(0, 2, mat);
    *tmp = 2*(ac+bd);
    tmp = matfuncCall(1, 0, mat);
    *tmp = 2*(ad+bc);
    tmp = matfuncCall(1, 1, mat);
    *tmp = aa-bb+cc-dd;
    tmp = matfuncCall(1, 2, mat);
    *tmp = 2*(-ab+cd);
    tmp = matfuncCall(2, 0, mat);
    *tmp = 2*(-ac+bd);
    tmp = matfuncCall(2, 1, mat);
    *tmp = 2*(ab+cd);
    tmp = matfuncCall(2, 2, mat);
    *tmp = (aa-bb-cc+dd);
    
    
}

