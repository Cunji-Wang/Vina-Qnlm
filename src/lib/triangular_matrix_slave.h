/*
 * triangular_matrix_slave.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 *
 * Distributed under terms of the MIT license.
 */
#ifndef TRIANGULAR_MATRIX_SLAVE_H
#define TRIANGULAR_MATRIX_SLAVE_H

struct triangular_matrix{
    double *m_data;
    int m_dim;
};


inline int triangular_matrix_index(int n, int i, int j) {
    //assert(j < n);

    return i + j*(j+1)/2;
}

inline int triangular_matrix_index_permissive(int n, int i, int j) {
    return (i <= j) ? triangular_matrix_index(n, i, j)
                    : triangular_matrix_index(n, j, i);
}

inline int matrix_index(struct triangular_matrix *m, int i, int j) {
    return triangular_matrix_index(m->m_dim, i, j);
}

inline int index_permissive(struct triangular_matrix *m, int i, int j) {
    return (i < j) ? matrix_index(m, i, j) : matrix_index(m, j, i);
}


#endif //TRIANGULAR_MATRIX_SLAVE_H
