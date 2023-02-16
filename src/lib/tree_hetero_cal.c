/*
 * tree_hetero_cal.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * 
 * Distributed under terms of the MIT license.
 */
#include "monte_carlo_slave.h"
#include "tranModel.h"
#include "quaternion_slave.h"
#include "common_slave.h"
#include "tree_hetero_cal.h"
#include <stdbool.h>
#include "slave.h"

__thread_local extern double *ligcoords_slave;
__thread_local extern double *minus_forces_slave;
__thread_local extern struct segmentbranch *segment_slave;
__thread_local extern struct flexibleHeritTree *heritTree_slave;

__thread_local extern double *ligandAtom_slave;

void local_to_lab_rigidbody(struct flexibleHeritTree *tree, double local_coords[3], double rst[3]) {
    double tmp[3];
    //tmp = origin + orientation_m*local_coords;
    matMultiply(local_coords, tree->orientation_m, tmp);
    vecplus(tmp,tree->origin, rst);
}

void local_to_lab_segment(struct segmentbranch *branch, double local_coords[3], double rst[3]) {
    double tmp[3];
    //tmp = origin + orientation_m*local_coords;
    matMultiply(local_coords, branch->orientation_m, tmp);
    
    vecplus(tmp,branch->origin, rst);
}

void local_to_lab_direction(struct segmentbranch *branch, double local_direction[3], double rst[3]) {

    //tmp = orientation_m * local_direction;
    matMultiply(local_direction, branch->orientation_m, rst);
}

void local_to_lab_direction_flexible(struct flexibleHeritTree *tree, double local_direction[3], double rst[3]) {

    //tmp = orientation_m * local_direction;
    matMultiply(local_direction, tree->orientation_m, rst);
}

void set_segment_coords(struct segmentbranch *branch, double *atoms, double *coords) {

    for(int i=branch->begin; i<branch->end;++i){
        local_to_lab_segment(branch, &atoms[i*37+5], &coords[i*3]);

	}
}

void set_rigidbody_coords(struct flexibleHeritTree *tree, double *atoms, double *coords) {
    //VINA_RANGE(i, begin, end)
    //for(std::size_t VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))
    //coords[i] = local_to_lab(atoms[i].coords);
    for(int i=tree->begin; i<tree->end;++i){
        local_to_lab_rigidbody(tree, &atoms[i*37+5], &coords[i*3]);
	}
}
void set_segment_orientation(struct segmentbranch *branch, double orientation[4]) {
    // does not normalize the orientation
    //orientation_q = q;
        branch->orientation_q[0] = orientation[0];
        branch->orientation_q[1] = orientation[1];
        branch->orientation_q[2] = orientation[2];
        branch->orientation_q[3] = orientation[3];

   // quaternion_to_r3(heritTree_slave->orientation_q, branch->orientation_m); ymx oringinal
    quaternion_to_r3(branch->orientation_q, branch->orientation_m);
//    for(int i=0;i<9;i++)printf("branch->orientation_m:%lf\n",branch->orientation_m[i]);
}


void set_rigidbody_orientation(struct flexibleHeritTree *tree, double orientation[4]) {
    // does not normalize the orientation
    //orientation_q = q;
        tree->orientation_q[0] = orientation[0];
        tree->orientation_q[1] = orientation[1];
        tree->orientation_q[2] = orientation[2];
        tree->orientation_q[3] = orientation[3];

    quaternion_to_r3(tree->orientation_q, tree->orientation_m);
}

void segment_set_conf(struct segmentbranch *children, struct segmentbranch *parents, double *atoms, double *coords, double torsion) {
    double origin[3];
    //origin = parent.local_to_lab(relative_origin);
    local_to_lab_segment(parents, children->relative_origin, children->origin);
    //axis = parent.local_to_lab_direction(relative_axis);
    local_to_lab_direction(parents, children->relative_axis, children->axis);
    
    //qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();
    double q[4];
    double tmp[4];
    angle_to_quaternion_two(children->axis,torsion,q);
   // quaternion_mul(q, heritTree_slave->orientation_q, tmp);   ymx original
    
    quaternion_mul(q, parents->orientation_q, tmp);

    quaternion_normalize_approx(tmp); // normalization added in 1.1.2
    
    //quaternion_normalize(tmp); // normalization added in 1.1.2
    //set_orientation(tmp);
    set_segment_orientation(children, tmp);
 
    //set_coords(atoms, coords);
    set_segment_coords(children, atoms, coords);
}

//this is for parents is rootnode
void segment_set_conf_level1(struct segmentbranch *children, double *atoms, double *coords, double torsion) {

    double tmp_1[3];
    //origin = parent.local_to_lab(relative_origin);
    //hetritree的local_to_lab fun
    matMultiply(children->relative_origin, heritTree_slave->orientation_m , tmp_1);  
  	
	vecplus(tmp_1, heritTree_slave->origin, children->origin);   
  //  vecplus(children->origin, heritTree_slave->origin, tmp_1);

    //axis = parent.local_to_lab_direction(relative_axis);
    matMultiply(children->relative_axis, heritTree_slave->orientation_m , children->axis);
    
 
    //qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();
    double q[4];
    double tmp[4];
    angle_to_quaternion_two(children->axis, torsion, q);
    quaternion_mul(q, heritTree_slave->orientation_q, tmp);
    quaternion_normalize_approx(tmp); // normalization added in 1.1.2

    //quaternion_normalize(tmp); // normalization added in 1.1.2
    //set_orientation(tmp);
    set_segment_orientation(children, tmp);

    //set_coords(atoms, coords);
    set_segment_coords(children, atoms, coords);
}


/* ymx oringinal useless???
void tree_set_conf(struct segmentbranch *children, struct segmentbranch *parents, double *atoms, double *coords, double *c, int offset) {
    //node.set_conf(parent, atoms, coords, c);
    segment_set_conf(children, parents, atoms, coords, c);
    printf("准备调用branch_set_conf\n");
    branches_set_conf(heritTree_slave, atoms, coords, c, offset);

}
*/
/*
void branches_set_conf(struct segmentbranch *children, struct segmentbranch *parents, double *atoms, double *coords, double *c) {
    //VINA_FOR_IN(i, b)
    //b[i].set_conf(parent, atoms, coords, c);
    int size = children->size;
    printf("branches_set_conf, children->size: %d, level: %d\n", children->size, children->level);
    for(int i=0; i<size; i++)
        tree_set_conf(children, parents, atoms, coords, c);
}
*/
//offset is current segmentbranch offset
/*
void branches_set_conf(struct flexibleHeritTree *tree, double *atoms, double *coords, double *c, int offset) {
    //VINA_FOR_IN(i, b)
    //b[i].set_conf(parent, atoms, coords, c);
    int size = segment_slave[offset].size;
    int parent = 0;
    int next = offset;
    printf("branches_set_conf, current branch size: %d, offset: %d\n", size, offset);
    
    for(int i=0; i<size; i++) {
        parent = segment_slave[offset + i].parent;
        next++;
        tree_set_conf(&segment_slave[offset + i], &segment_slave[parent], atoms, coords, c, next);
    }
}
*/

void rigidbody_set_conf(struct flexibleHeritTree *tree, double *atoms, double *coords, struct output_type_s *c) {
    //origin = c.position;
        heritTree_slave->origin[0] = c->position[0];
        heritTree_slave->origin[1] = c->position[1];
        heritTree_slave->origin[2] = c->position[2];
    set_rigidbody_orientation(tree, c->orientation);
    //set_coords(atoms, coords);
    set_rigidbody_coords(tree, atoms, coords);
}

void flexible_set_conf(struct segmentbranch *children, struct flexibleHeritTree *tree, double *atoms, double *coords, double *c){
    double torsion = *c;
    ++c;
    //origin = parent.local_to_lab(relative_origin);
    local_to_lab_rigidbody(tree, children->relative_origin, children->origin);
    //axis = parent.local_to_lab_direction(relative_axis);
    local_to_lab_direction_flexible(tree, children->relative_axis, children->axis);
   // qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();
    double q[4];
    double tmp[4];
    angle_to_quaternion_two(children->axis,torsion,q);
    quaternion_mul(q, tree->orientation_q, tmp);

    quaternion_normalize_approx(tmp); // normalization added in 1.1.2
    //quaternion_normalize(tmp); // normalization added in 1.1.2
    //set_orientation(tmp);
    set_rigidbody_orientation(tree, tmp);
    set_rigidbody_coords(tree, atoms, coords);
}

extern int forbranch;
void ligands_set_conf(struct output_type_s *c) {
    //lu: if forbranch is 0 this is useless
    rigidbody_set_conf(heritTree_slave, ligandAtom_slave, ligcoords_slave, c);

    //flv::const_iterator p = c.torsions.begin();
    double *p = &c->torsPlsCoord[0];
    // first we update heritTree parents node
    //branches_set_conf(children, node, atoms, coords, p);
//printf("ligands_set_conf first p value: %lf\n", *p);
    /*for(int i=0; i<heritTree_slave->branchsize; i++) {
        flexible_set_conf(&segment_slave[i], heritTree_slave,
                          globalTm->ligandAtom, ligcoords_slave, p);
        // p should be updated according to branchsize
        printf("ligands_set_conf update P value: %lf, i: %d\n", *p, i);
    	printf("\n\nligands_set_conf 进行setconf，准备branches_set_conf\n");
        branches_set_conf(heritTree_slave,globalTm->ligandAtom, ligcoords_slave, p, i);
    }*/
    // here we update branches node
    //branches_set_conf(heritTree_slave,globalTm->ligandAtom, ligcoords_slave, p, heritTree_slave->branchsize);
    

    // luhao modify
    /*
    int *levelindex = (int *)malloc( sizeof(int) * forbranch);
    int *levelvalue = (int *)malloc( sizeof(int) * forbranch);
    //sortByleval(level);
    for(int i=0; i<forbranch; i++){
        levelvalue[i]=segment_slave[i].level;  // level
        //printf("%d\n",levelvalue[i]);
        levelindex[i]=i;// 记录下标
    }

    //bubble sort useless?
    for(int i = 0; i<forbranch-1; i++){
        for(int j=0; j<forbranch-i-1; j++){
            if(levelvalue[j] > levelvalue[j+1]){
                int tmp;
                tmp = levelvalue[j+1];
                levelvalue[j+1] = levelvalue[j];
                levelvalue[j] = tmp;

                tmp = levelindex[j+1];
                levelindex[j+1] = levelindex[j];
                levelindex[j] = tmp;  
            }
        }  
    }
*/
    // qnml_mnp_10033，have 7 branch;
    for(int i = 0; i<forbranch; i++){
        //int update = levelindex[i];
        struct segmentbranch *children = &segment_slave[i];
        int parent = segment_slave[i].parent;
	double torsion = c->torsPlsCoord[i];
        if(parent == -1){ //first level，parent is -1
            segment_set_conf_level1(children, ligandAtom_slave, ligcoords_slave, torsion);
        }else{
            struct segmentbranch *parents = &segment_slave[parent];	
 	    double torsion = c->torsPlsCoord[i];
            segment_set_conf(children, parents, ligandAtom_slave, ligcoords_slave, torsion);
        }
        
   }

}

void sum_force_and_torque_tree(struct flexibleHeritTree *tree, double *coords, double *forces, double tmp[2][3]) {

    //tmp.first.assign(0);
    //tmp.second.assign(0);
    tmp[0][0] = 0;
    tmp[0][1] = 0;
    tmp[0][2] = 0;
    tmp[1][0] = 0;
    tmp[1][1] = 0;
    tmp[1][2] = 0;
    //VINA_RANGE(i, begin, end) {
    for(int i=tree->begin; i<tree->end;++i){
        //tmp.first  += forces[i];
        vecplusequal(tmp[0], &forces[i*3]);

        //tmp.second += cross_product(coords[i] - origin, forces[i]);
        double rst[3];
	rst[0]=2;
        vecminus(&coords[i*3], tree->origin, rst);
        double cs[3];
        cross_product(rst, &forces[i*3], cs);
        vecplusequal(tmp[1], cs);	
    }
}

void sum_force_and_torque_branch(struct segmentbranch *children, double *coords, double *forces, double tmp[2][3]) {

    //tmp.first.assign(0);
    //tmp.second.assign(0);
    tmp[0][0] = 0;
    tmp[0][1] = 0;
    tmp[0][2] = 0;
    tmp[1][0] = 0;
    tmp[1][1] = 0;
    tmp[1][2] = 0;
    //VINA_RANGE(i, begin, end) {
    for(int i=children->begin; i<children->end;++i){
        //tmp.first  += forces[i];
        vecplusequal(tmp[0], &forces[i*3]);
        //tmp.second += cross_product(coords[i] - origin, forces[i]);
        double rst[3];
        vecminus(&coords[i*3], children->origin, rst);
        double cs[3];
        cross_product(rst, &forces[i*3], cs);
        vecplusequal(tmp[1], cs);
    }
}

void rigid_set_derivative(double force_torque[2][3], struct change *c) {
    //c.position     = force_torque.first;
    //c.orientation  = force_torque.second;
        c->position[0] = force_torque[0][0];
        c->position[1] = force_torque[0][1];
        c->position[2] = force_torque[0][2];
        c->orientation[0] = force_torque[1][0];
        c->orientation[1] = force_torque[1][1];
        c->orientation[2] = force_torque[1][2];
}

void frame_set_derivative(struct segmentbranch *children, double force_torque[2][3], double *c) {
    //c = force_torque.second * axis;
    *c = vec_dotproduct(force_torque[1], children->axis);
}

void tree_derivative(struct segmentbranch *children, double force_torque[2][3], double *coords, double *forces, double *p, int offset){
    //vecp force_torque = node.sum_force_and_torque(coords, forces);
    sum_force_and_torque_branch(children, coords, forces, force_torque);
    //fl& d = *p; // reference
    double d = *p;
    ++p;
    int next = offset;
    next++;
    int parent = children->parent;
    //branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
    branches_derivative(heritTree_slave, segment_slave[parent].origin, coords, 
				forces, force_torque, p, next, false);
    //node.set_derivative(force_torque, d);
    // todo: yumaoxue d will be moved by ++p? need validate!!
    frame_set_derivative(children, force_torque, &d);

}

/*
void branches_derivative(struct flexibleHeritTree *tree, double origin[3], double *coords, double *forces,
        double out[2][3], double *d, int offset, bool isrigid) {
    // adds to out
    int next = 0;
    int parent = 0;    
    int size = segment_slave[offset].size;
    if(isrigid)
        next = offset+tree->branchsize;
    else
        next = offset;
    //VINA_FOR_IN(i, b) {
    for(int i=0; i<size; i++){
        //vecp force_torque = b[i].derivative(coords, forces, d);
        next++;
        parent = segment_slave[offset + i].parent;
        printf("branches_derivative, curent chilid location: %d, next: %d\n", offset+i, next);
        double force_torque[2][3];
        tree_derivative(&segment_slave[offset+i], force_torque, coords, forces, d, next);
        //out.first  += force_torque.first;
        vecplusequal(out[0], force_torque[0]);
        //vec r; r = b[i].node.get_origin() - origin;
        double r[3];
        double parent_origin[3];
        if(isrigid){
            for(int j=0; j<3; j++)
                parent_origin[j] = tree->origin[j];
        }
        else{
            for(int j=0; j<3; j++)
                parent_origin[j] = segment_slave[parent].origin[j];
        }
        vecminus(parent_origin,segment_slave[offset+i].origin, r);
        //out.second += cross_product(r, force_torque.first) + force_torque.second;
        double rst[3], tmp[3];
        cross_product(r, force_torque[0], rst);
        vecplus(rst, force_torque[1], tmp);
        vecplusequal(out[1], tmp);
    }
}
*/

void inter_derivative(int *i, double tmp_force_torque[2][3],struct change *c){
    sum_force_and_torque_branch(&segment_slave[*i], ligcoords_slave, minus_forces_slave, tmp_force_torque);

    int offset = *i;
    int size = segment_slave[*i].childrensize;

    for(int k=0; k<size; k++){//if this node have children???    
        *i = *i+1;
        int tmp_p = *i;
        double local_force_torque[2][3];
	inter_derivative(i, local_force_torque, c);

        //out.first += force_torque.first;
        tmp_force_torque[0][0] += local_force_torque[0][0];
        tmp_force_torque[0][1] += local_force_torque[0][1];
        tmp_force_torque[0][2] += local_force_torque[0][2];
        
        // node.get_origin()
        int parent = segment_slave[tmp_p].parent;
        double origin[3];
        origin[0] = segment_slave[tmp_p].origin[0];
        origin[1] = segment_slave[tmp_p].origin[1];
        origin[2] = segment_slave[tmp_p].origin[2];
	       
        double parent_origin[3];
        if(parent == -1){
            parent_origin[0] = heritTree_slave->origin[0];
            parent_origin[1] = heritTree_slave->origin[1];
            parent_origin[2] = heritTree_slave->origin[2];
        }else{
            parent_origin[0] = segment_slave[parent].origin[0];
            parent_origin[1] = segment_slave[parent].origin[1];
            parent_origin[2] = segment_slave[parent].origin[2];
        }
        double r[3]; 
        r[0] = origin[0] - parent_origin[0];
        r[1] = origin[1] - parent_origin[1];
        r[2] = origin[2] - parent_origin[2];
 
        //out.second += cross_product(r, force_torque.first) + force_torque.second;
	double tmp[3];
        double force_torque_first[3];
        force_torque_first[0] = local_force_torque[0][0];
        force_torque_first[1] = local_force_torque[0][1];
        force_torque_first[2] = local_force_torque[0][2];

        cross_product(r, force_torque_first , tmp);

/*        this will be error
 *        for(int j=0; j<3; j++){
            tmp_force_torque[1][j] = tmp[j] + tmp_force_torque[1][j] + local_force_torque[1][j] ;
	
        }
*/
       tmp_force_torque[1][0] = tmp[0] + local_force_torque[1][0] + tmp_force_torque[1][0]; //c8
       tmp_force_torque[1][1] = tmp[1] + local_force_torque[1][1] + tmp_force_torque[1][1]; //c8
       tmp_force_torque[1][2] = tmp[2] + local_force_torque[1][2] + tmp_force_torque[1][2]; //c8
   } 
    //node.set_derivate;
    c->torsions[offset] = vecMultiVec(tmp_force_torque[1] , segment_slave[offset].axis);
    
}


void ligands_derivative(struct change *c) {
    //vecp force_torque = node.sum_force_and_torque(coords, forces);
    double force_torque[2][3];
    sum_force_and_torque_tree(heritTree_slave, ligcoords_slave, minus_forces_slave, force_torque);

    //flv::iterator p = c.torsions.begin();
    double *p = &c->torsions[0];
    int parent = 0;
    double tmp_orgin[3];
    //branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
    /*for(int i=0; i<heritTree_slave->branchsize; i++) {
        //first we need to get node.get_origin()
        //parent = segment_slave[offset + i]->parent;
        for(int j=0; j<3; j++)
            tmp_orgin[j] = heritTree_slave->origin[j];

        // here we update branches node
        branches_derivative(heritTree_slave, tmp_orgin, ligcoords_slave,
                minus_forces_slave, force_torque, p, i, true);
    }
    */

    //luhao modify
    double f_t[forbranch*3];
    double tmp_force_torque[2][3];

    for(int i=0; i<forbranch; i++){
	 tmp_force_torque[0][0]=0;
	 tmp_force_torque[0][1]=0;
	 tmp_force_torque[0][2]=0;
	 tmp_force_torque[1][0]=0;
	 tmp_force_torque[1][1]=0;
	 tmp_force_torque[1][2]=0;
        
    //sum_force_and_torque_tree(&segment_slave[i], ligcoords_slave, minus_forces_slave, tmp_force_torque);
	
	int offset = i;	
	inter_derivative(&i, tmp_force_torque, c);

        //out.first += force_torque.first; 
        force_torque[0][0] += tmp_force_torque[0][0];
        force_torque[0][1] += tmp_force_torque[0][1];
        force_torque[0][2] += tmp_force_torque[0][2];

        // node.get_origin()
        int parent = segment_slave[offset].parent;
        double origin[3];
        origin[0] = segment_slave[offset].origin[0];
        origin[1] = segment_slave[offset].origin[1];
        origin[2] = segment_slave[offset].origin[2];
        double parent_origin[3];
       
        if(parent == -1){
            parent_origin[0] = heritTree_slave->origin[0];
            parent_origin[1] = heritTree_slave->origin[1];
            parent_origin[2] = heritTree_slave->origin[2];
        }else{
            parent_origin[0] = segment_slave[parent].origin[0];
            parent_origin[1] = segment_slave[parent].origin[1];
            parent_origin[2] = segment_slave[parent].origin[2];
        }
        double r[3]; 
        r[0] = origin[0] - parent_origin[0];
        r[1] = origin[1] - parent_origin[1];
        r[2] = origin[2] - parent_origin[2];
 
        //out.second += cross_product(r, force_torque.first) + force_torque.second;
        double tmp[3];
        double force_torque_first[3];
        force_torque_first[0] = tmp_force_torque[0][0];
        force_torque_first[1] = tmp_force_torque[0][1];
        force_torque_first[2] = tmp_force_torque[0][2];
        cross_product(r, force_torque_first , tmp);
/*     this is error
 *     for(int j=0; j<3; j++){
            force_torque[1][j] =tmp[j] + force_torque[1][j] + tmp_force_torque[1][j];
        }
*/


	force_torque[1][0] = tmp[0] + tmp_force_torque[1][0] + force_torque[1][0]; //c8	
	force_torque[1][1] = tmp[1] + tmp_force_torque[1][1] + force_torque[1][1]; //c8	
	force_torque[1][2] = tmp[2] + tmp_force_torque[1][2] + force_torque[1][2]; //c8	

        //update torsion
        //*p = elementwise_product()
/*        *p = force_torque[1][0] * segment_slave[i].axis[0] + force_torque[1][1] * segment_slave[i].axis[1]
                + force_torque[1][2] * segment_slave[i].axis[2];
*/     //yxm original

    }

   // node.set_derivative(force_torque, c.rigid); 
    rigid_set_derivative(force_torque, c);
}
