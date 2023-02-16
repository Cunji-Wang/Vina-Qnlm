/*
 * model_slave.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <math.h>
#include "tranModel.h"
#include "atom_slave.h"
#include "common_slave.h"
#include "cachegrid.h"
#include "tree_hetero_cal.h"
#include "slave.h"

extern double ***mysmooth;
extern double *smooth_slave;

extern __thread_local struct paras *para_s __attribute__ ((aligned(64)));

__thread_local extern double *ligcoords_slave;
__thread_local extern double *minus_forces_slave;
__thread_local extern struct flexibleHeritTree *heritTree_slave;

__thread_local extern double *ligandAtom_slave;
//__thread_local extern double receptorAtom_slave[100];
__thread_local extern int *ligandpair_slave; //ligand_pair 

int num_movable_atoms() {
    return para_s->m_num_movable_atoms;
}
int num_internal_pairs() {
    return para_s->pairsize;
}

double gyration_radius() {

    double acc = 0;
    unsigned counter = 0;
    for(int i=heritTree_slave->begin; i<heritTree_slave->end; i++){
        //if(globalTm->receptorAtom[i*25] != EL_TYPE_H){
        //if(receptorAtom_slave[i] != EL_TYPE_H){
        if(ligandAtom_slave[i*37] != EL_TYPE_H){
            double tmp = vec_distance_sqr(&ligcoords_slave[i*3], heritTree_slave->origin);
            acc += tmp;
            ++counter;
        }
    }
    return (counter > 0) ? sqrt(acc/counter) : 0;
}

void p_eval_deriv(int index, double r2, double *rst){

    double r2_factored = 32 * r2;
    int i1 = (int)(r2_factored);
    int i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
    double rem = r2_factored - i1;
    //const pr& p1 = smooth[i1];
    //const pr& p2 = smooth[i2];
    double p1[2], p2[2];
    //0 means first, 1 means second

    //int totalsize = para_s->matsize*para_s->pre_n*2;
    //p1[0] = smooth_slave[totalsize*_MYID + para_s->pre_n*2*index + i1*2 + 0];
    p1[0] = smooth_slave[para_s->pre_n*2*index + i1*2 + 0];
    p1[1] = smooth_slave[para_s->pre_n*2*index + i1*2 + 1];
    p2[0] = smooth_slave[para_s->pre_n*2*index + i2*2 + 0];
    p2[1] = smooth_slave[para_s->pre_n*2*index + i2*2 + 1];

/*
    mooth_slave1[0] = mysmooth[index][i1][0];
    p1[1] = mysmooth[index][i1][1];
    p2[0] = mysmooth[index][i2][0];
    p2[1] = mysmooth[index][i2][1];
*/

    //orignial code
    //fl e   = p1.first  + rem * (p2.first  - p1.first);
    //fl dor = p1.second + rem * (p2.second - p1.second);
    double e   = p1[0]  + rem * (p2[0]  - p1[0]);
    double dor = p1[1] + rem * (p2[1] - p1[1]);
    rst[0] = e;
    rst[1] = dor;

}
double eval_interacting_pairs_deriv(double v) { // adds to forces  // clean up
    double cutoff_sqr = 64;
    double e = 0;
    for(int i=0; i<para_s->pairsize; i++){
        //int *ip = globalTm->ligandpair[i];
        //int *ip = ligandpair_slave[i];
        int *ip = &ligandpair_slave[i*3];
        double *coord_b = &ligcoords_slave[ip[2]*3];
        double *coord_a = &ligcoords_slave[ip[1]*3];
        double r[3];
        vecminus(coord_b, coord_a, r);	
        double r2 = sqr(r);
        if(r2<cutoff_sqr){
            double tmp[2];
            p_eval_deriv(ip[0], r2, tmp);
            double force[3];
           // vecmultiply(r, tmp[1], force);
	    force[0]=r[0]*tmp[1];
	    force[1]=r[1]*tmp[1];
	    force[2]=r[2]*tmp[1];

            curl(&tmp[0], force, v);
            e +=tmp[0];
            vecminusequal(&minus_forces_slave[ip[1]*3], force);
            vecplusequal(&minus_forces_slave[ip[2]*3], force);
        }
    }
    return e;
}


void model_set(struct output_type_s *c) {
    ligands_set_conf(c);
    //flex   .set_conf(atoms, coords, c.flex);
}

double model_eval_deriv  (double v[3], struct output_type_s *c, struct change *g) { // clean up
    model_set(c);
/*
for(int i=0;i<135;i++){unsigned long *p;
        p=(unsigned long*)&minus_forces_slave[i];
        printf("forces:%.8lf--orientation_x:%lx\n",minus_forces_slave[i],*p);
}
printf("++++++++++++++++++++++++++\n");
*/
    double e = cache_eval_deriv(v[1]); // sets minus_forces, except inflex
    // e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
    //VINA_FOR_IN(i, ligands)
/*
        for(int i=0;i<15;i++){unsigned long *p;
        p=(unsigned long*)&minus_forces_slave[i];
        printf("forces:%.8lf--orientation_x:%lx\n",minus_forces_slave[i],*p);
}
printf("++++++++++++++++++++++++++\n");
*/
    e += eval_interacting_pairs_deriv(v[0]); // adds to minus_forces
//for(int i=0;i<135;i++)printf("%dcoords:%.8lf\n",i,ligcoords_slave[i]);
 // for(int i=0;i<4;i++)printf("coords:%.15lf\n",g->orientation[i]);
//while(1);

/*for(int i=0;i<15;i++){unsigned long *p;
        p=(unsigned long*)&minus_forces_slave[i];
        printf("forces:%.8lf--orientation_x:%lx\n",minus_forces_slave[i],*p);
}
*/  
    ligands_derivative(g);
//for(int i=0;i<3;i++)printf("coords:%.15lf\n",g->orientation[i]);
/*for(int i=0;i<4;i++){unsigned long *p;
        p=(unsigned long*)&g->orientation[i];
        printf("orientation_x:%lx\n",*p);
}

for(int i=0;i<4;i++)printf("coords:%.15lf\n",g->orientation[i]);
while(1);
*/
    //flex   .derivative(coords, minus_forces, g.flex); // inflex forces are ignored   	
	return e;
}

void get_heavy_atom_movable_coords(double *coords)  { // FIXME mv
    int num = para_s->m_num_movable_atoms;
    int torsionsize = para_s->torsionsize;
    //VINA_FOR(i, num_movable_atoms())
    for(int i=0; i<num; i++){
        //if(atoms[i].el != EL_TYPE_H)
        //saveAtoms[i*25] = (double)tmp.el;
        if(ligandAtom_slave[i*37]!= EL_TYPE_H){
            //tmp.push_back(coords[i]);
            coords[i*3 + torsionsize] = ligcoords_slave[i*3];
            coords[i*3+1 + torsionsize] = ligcoords_slave[i*3+1];
            coords[i*3+2 + torsionsize] = ligcoords_slave[i*3+2];
        }
    }

}
