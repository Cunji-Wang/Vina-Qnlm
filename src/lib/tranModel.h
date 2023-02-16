/*
 * tranModel.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef VINA_TRANMODEL_H
#define VINA_TRANMODEL_H

enum t {EL, AD, XS, SY};

struct paras{
    int coordsize;
    int pairsize;
    int torsionsize;
    int singlegridsize;
    int bfgs_steps;
    int mc_steps;
    int num_saved_mins;
    int m_num_movable_atoms;
    enum t atom_typing_used;
    double temperature;
    double min_rmsd;
    double mutation_amplitude;
    double hunt_cap[3];
    double corner1[3];
    double corner2[3];
    //luhao add
    //int grid_index[19];
    //int realgridsize;
    int ligsize;  
    int pre_n;    //precaculate
    int matsize;  //precaculate
	
   // int degrees_of_freedom;
};

struct segmentbranch{
    /* total elements are 27. vec relative_axis(3), vec relative_origin(3), vec axis（3),
     * vec origin_(3), mat orientation_m (9), qt  orientation_q (4), sz begin, sz end
     */
    //double parentNode[25];
    double relative_axis[3];
    double relative_origin[3];
    double axis[3];
    double origin[3];
    double orientation_m[9];
    double orientation_q[4];
    int begin;    
    int end;     //
    int size; //the size of this level children size      //有多少个孩子
    int childrensize; //signal the childern are in the same situation  //本节点，children size
    int parent; //signal the upper layer array num    父节点
    struct segmentbranch *segchildren;
};

// this struct flexibleHeritTree will replace the structure heterotree
struct flexibleHeritTree {
    // rigidbody have 18 element: vec origin_(3), mat orientation_m (9), qt  orientation_q (4), sz begin, sz end
    double origin[3];
    double orientation_m[9];
    double orientation_q[4];
    int begin;
    int end;
    int branchsize;     //记录有多少个孩子
    struct segmentbranch flexiblechildren[0];//f[0],f[1]//    其实应该有多个；
};

struct cachegrid {   //对应grid.h
    // gd should be gd[3][3], element is: double begin, double end, int n
    double gd[9];
    double slope;
    enum t atu;

    double m_init[3];
    double m_range[3];
    double m_factor[3];
    double m_dim_fl_minus_1[3];
    double m_factor_inv[3];

    int m_grid_data_i;
    int m_grid_data_j;
    int m_grid_data_k;
    /*vec m_init; vec m_range; vec m_factor; vec m_dim_fl_minus_1; vec m_factor_inv; total have 15 elements
     * array3d<fl> m_data: sz m_i, m_j, m_k, std::vector<double> m_data; the m_data size is dynamic allocated
     * total 18 element + m_data.size
     */
    double grids[0];
};

//this datas can be allocated to continuous memeory area
struct tranModel {    //对应model.h
    /* this receptor atom has element: sz el, ad, xs, sy, double charge, double coords[3]
     * std::vector<bond> bonds can be separated: sz bondsize, int i, bool in_grid, bool rotatable, double length
     * default, the max bond size is 4. so one atoms should have 9+4*4=25
     * the type sz, bool and int will be forced translate to double type for uniformity
     */
    double *receptorAtom; // for receptor atom;
    double *ligandAtom; // for ligand atoms;
    double *ligcoords;
    double *internal_coords;
    double *minus_forces;
    int **ligandpair;
    struct flexibleHeritTree *heritTree;
    struct paras *para;
    struct cachegrid *cachegrids;
};

struct output_type_s {     //对应out_type对应
    double e;
    double position[3];
    double orientation[4];
    //int count_step;
    //its size depends on the input pdb file, will allocate dynamic from host passed para: torsionsize, my test is 12
    //also the coords will be stored following the torsion
    // the size should be torsionsize+coordssize
    double *torsPlsCoord;
    //double torsPlsCoord[0];

};

#endif


