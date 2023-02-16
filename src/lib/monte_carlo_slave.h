/*
 * tranModel.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef MONTE_CARLO_SLAVE_H
#define MONTE_CARLO_SLAVE_H

#define MAX_TORSIONS_SIZE = 100;

struct change{
    double position[3];
    double orientation[3];
    double torsionsize;
    //its size depends on the input pdb file, will allocate dynamic from host passed para: torsionsize, my test is 12
    double torsions[100];

};

// change operator func: fl& operator()(sz index)
inline double* change_funcCall(struct change *g, int index) {

   // VINA_FOR_IN(i, ligands) {
        if(index < 3) return &g->position[index];
        index -= 3;

        if(index < 3) return &g->orientation[index];
        index -= 3;

        if(index < g->torsionsize) return &g->torsions[index];
        index -= g->torsionsize;


    //return &ligands[0].rigid.position[0]; // shouldn't happen, placating the compiler
	return &g->position[0];
}

#endif
