/*
 * mutate_slave.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

//#include "mutate_slave.h"
#include "tranModel.h"
#include "common_slave.h"
#include "quaternion_slave.h"
#include "randomgen.h"
#include "model_slave.h"
#include <float.h>
#include "slave.h"
extern __thread_local struct flexibleHeritTree *heritTree_slave;

void mutate_conf(struct output_type_s *candidate, int torsionsize, double amplitude) {
    int mutable_entities_num = torsionsize+2;
    //int which_int = random_int(0, int(mutable_entities_num - 1), generator);
    int which_int = random_integer(0, (double)(mutable_entities_num - 1));

   // if(which_int == 0) { c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator); return; }
    double randomtmp[3];
    if(which_int == 0){
        random_inside_sphere(randomtmp);
//	for(int i=0;i<3;i++){printf("%d-mutate_randomtmp:%lf\n",my_id,randomtmp[i]);}
        randomtmp[0] = amplitude * randomtmp[0];
        randomtmp[1] = amplitude * randomtmp[1];
        randomtmp[2] = amplitude * randomtmp[2];
        vecplusequal(candidate->position, randomtmp);
        return;
    }
    --which_int;
    if(which_int == 0) {
        double gr = gyration_radius();
        if(gr > DBL_EPSILON) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
             double rotation[3];
             //rotation = amplitude / gr * random_inside_sphere(generator);
            random_inside_sphere(rotation);
            rotation[0] = amplitude / gr * rotation[0];
            rotation[1] = amplitude / gr * rotation[1];
            rotation[2] = amplitude / gr * rotation[2];

            quaternion_increment(candidate->orientation, rotation);
        }
        return;
    }
    --which_int;
    if(which_int < torsionsize) {
       // c.ligands[i].torsions[which_int] = random_fl(-pi, pi, generator);
        candidate->torsPlsCoord[which_int] = random_fl(-pi, pi);
        return;
    }

    which_int -= torsionsize;

}




void createnewchange(struct change *g_orgin, struct change *g_new, int size){
    
        g_new->position[0] = g_orgin->position[0];
        g_new->position[1] = g_orgin->position[1];
        g_new->position[2] = g_orgin->position[2];
        g_new->orientation[0] = g_orgin->orientation[0];
        g_new->orientation[1] = g_orgin->orientation[1];
        g_new->orientation[2] = g_orgin->orientation[2];

    for(int i=0; i<size; i++)
        g_new->torsions[i] = g_orgin->torsions[i];
    g_new->torsionsize = g_orgin->torsionsize;
}

void createnewconf(struct output_type_s *x_orgin, struct output_type_s *x_new, int size){
  
        x_new->position[0] = x_orgin->position[0];
        x_new->position[1] = x_orgin->position[1];
        x_new->position[2] = x_orgin->position[2];
        x_new->orientation[0] = x_orgin->orientation[0];
        x_new->orientation[1] = x_orgin->orientation[1];
        x_new->orientation[2] = x_orgin->orientation[2];
        x_new->orientation[3] = x_orgin->orientation[3];

    for(int i=0; i<size; i++)
        x_new->torsPlsCoord[i] = x_orgin->torsPlsCoord[i];
    

    x_new->e = x_orgin->e;
}

void print_conf(struct output_type_s *x, int size){
    for(int i=0; i<3; i++){
        printf("position: %lf\n", x->position[i]);
    }
    for(int i=0; i<4; i++){
        printf("orientation: %lf\n", x->orientation[i]);
    }
    for(int i=0; i<size; i++)
        printf("torsPlsCoord:%lf\n", x->torsPlsCoord[i]);
    printf("e: %lf\n", x->e);
}

