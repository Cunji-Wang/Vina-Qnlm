/* * cachegrid.c
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 * Checked by 2021 luhao <luhao1027@126.com>
 *
 * Distributed under terms of the MIT license.
 */

#include "atom_slave.h"
#include "model_slave.h"
#include "common_slave.h"
#include "slave.h"
#include "tranModel.h"
#include "simd.h"

extern struct tranModel *globalTm;
extern __thread_local struct paras *para_s __attribute__ ((aligned(64)));

//updata
extern __thread_local double *ligcoords_slave;
extern __thread_local double *minus_forces_slave;

//fixed
extern __thread_local struct cachegrid cachegrids_slave;
//extern __thread_local double grids_slave[19*18];
extern __thread_local double *ligandAtom_slave;
//extern double *grids_array;

extern __thread_local volatile crts_rply_t dma_rply;
extern __thread_local unsigned int D_COUNT;
extern __thread_local int my_id;


__thread_local int cache_flag=1;
int times = 0;
//__thread_local int offset;

//cache grid contants
extern __thread_local double m_init[3];
extern __thread_local double m_factor[3];
extern __thread_local double m_dim_fl_minus_1[3];
extern __thread_local double m_factor_inv[3];
extern __thread_local int dim[3];
extern __thread_local doublev8 e_x;
extern __thread_local doublev8 e_y;
extern __thread_local doublev8 e_z; 


double array3d_opr_funcCall(double *grids, int i, int j, int k) {
/* tm->cachegrids->grids[i * singlegridsize + 15] = (double) mygrids[i].m_data.dim0();
 * tm->cachegrids->grids[i * singlegridsize + 16] = (double) mygrids[i].m_data.dim1();
 * tm->cachegrids->grids[i * singlegridsize + 17] = (double) mygrids[i].m_data.dim2();
 *
*/
    int m_i = grids[15];
    int m_j = grids[16];
    return grids[18 + i + m_i*(j + m_j*k)];
    // globalTm->cachegrids->grids[offset]
    //return globalTm->cachegrids->grids[offset + 18 + i + m_i*(j + m_j*k)];
/*    printf("%lf\n",globalTm->cachegrids->grids[offset + 18 + i + m_i*(j + m_j*k)]);
    printf("----%d===%d\n",offset + 18 + i + m_i*(j + m_j*k), _MYID);
    printf("array:%lf\n",grids_array[(offset + 18 + i + m_i*(j + m_j*k))]);
    while(1);
  */  
   // printf("hello%d\n", offset + 18 + i + m_i*(j + m_j*k) + _MYID*para_s->singlegridsize*para_s->realgridsize);
    
    //return grids_array[(offset + 18 + i + m_i*(j + m_j*k)) +  _MYID*para_s->singlegridsize*para_s->realgridsize];
}

double evaluate_aux(double *grids, double location[3], double slope, double v, double deriv[3]) {
    // sets *deriv if not NULL
    //vec s  = elementwise_product(location - m_init, m_factor);
    double s[3];
    //double m_init[3];
    //double m_factor[3];
    double tmpr[3];
    double miss[3]={0};
    //double m_dim_fl_minus_1[3];
    //double m_factor_inv[3];
    //int dim[3];
    //  vec m_init; vec m_range; vec m_factor; vec m_dim_fl_minus_1; vec m_factor_inv;
/*    for(int i=0; i<3; i++){
        m_init[i] = cachegrids_slave.m_init[i];
        m_factor[i] = cachegrids_slave.m_factor[i];
        miss[i] = 0;
        m_dim_fl_minus_1[i] = cachegrids_slave.m_dim_fl_minus_1[i];
        m_factor_inv[i] = cachegrids_slave.m_factor_inv[i];
    }
*/
/*
    m_init[0] = cachegrids_slave.m_init[0];
    m_init[1] = cachegrids_slave.m_init[1];
    m_init[2] = cachegrids_slave.m_init[2];
    m_factor[0] = cachegrids_slave.m_factor[0];
    m_factor[1] = cachegrids_slave.m_factor[1];
    m_factor[2] = cachegrids_slave.m_factor[2];
    //miss[i] = 0;
    m_dim_fl_minus_1[0] = cachegrids_slave.m_dim_fl_minus_1[0];
    m_dim_fl_minus_1[1] = cachegrids_slave.m_dim_fl_minus_1[1];
    m_dim_fl_minus_1[2] = cachegrids_slave.m_dim_fl_minus_1[2];
    m_factor_inv[0] = cachegrids_slave.m_factor_inv[0];
    m_factor_inv[1] = cachegrids_slave.m_factor_inv[1];
    m_factor_inv[2] = cachegrids_slave.m_factor_inv[2];
    dim[0] = cachegrids_slave.m_grid_data_i;
    dim[1] = cachegrids_slave.m_grid_data_j;
    dim[2] = cachegrids_slave.m_grid_data_k;
*/	

    //location - m_init
    vecminus(location, m_init, tmpr);
    //vecminus(location, m_init, tmpr);
    
    elementwise_product(tmpr, m_factor, s);
    //elementwise_product(tmpr, m_factor, s);

    int region[3];
    int a[3];

    for(int i=0; i<3; i++) {
        if(s[i] < 0.0000) {
            miss[i] = -s[i];
            region[i] = -1;
            a[i] = 0;
            s[i] = 0;
        }
        else if(s[i] >= m_dim_fl_minus_1[i]) {
            miss[i] = s[i] - m_dim_fl_minus_1[i];
            region[i] = 1;
            //a[i] = m_data.dim(i) -  2;
            a[i] = dim[i] - 2;
            s[i] = 1;
        }
        else {
            region[i] = 0; // now that region is boost::array, it's not initialized
            a[i] = (int)(s[i]);
            s[i] -= (double)a[i];
        }
    }
    //double penalty = slope * (miss * m_factor_inv);
    double result = vecMultiVec(miss, m_factor_inv);
    double penalty = slope * result;
    // FIXME check that inv_factor is correctly initialized and serialized
    int x0 = a[0];
    int y0 = a[1];
    int z0 = a[2];

    //int x1 = x0+1;
    //int y1 = y0+1;
    //int z1 = z0+1;
    //get 8 grid data:f000 f100 f010 f110 f001 f101 f011 f111
/*    CRTS_dma_iget(&f8_array[0], &grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))], sizeof(double) * 8, &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
*/
    //printf("before array3d_opr_funcCall");
    double f000 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))];
    double f100 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+1];
    double f010 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+2];
    double f110 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+3];
    double f001 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+4];
    double f101 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+5];
    double f011 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+6];
    double f111 = grids[8 * (x0 + dim[0]*(y0 + dim[1]*z0))+7];
/*
    double f000 = array3d_opr_funcCall(grids, x0, y0, z0);
    double f100 = array3d_opr_funcCall(grids, x1, y0, z0);
    double f010 = array3d_opr_funcCall(grids, x0, y1, z0);
    double f110 = array3d_opr_funcCall(grids, x1, y1, z0);
    double f001 = array3d_opr_funcCall(grids, x0, y0, z1);
    double f101 = array3d_opr_funcCall(grids, x1, y0, z1);
    double f011 = array3d_opr_funcCall(grids, x0, y1, z1);
    double f111 = array3d_opr_funcCall(grids, x1, y1, z1);
*/
 //   printf("ater array3d_opr_funcCall,%lf,,,%lf,,,,%lf,,,%lf,,,,%lf,,,,%lf,,,%lf,,,%lf\n",f000,f100,f010,f110,f001,f101,f011,f111);

    double x = s[0];
    double y = s[1];
    double z = s[2];

    double mx = 1-x;
    double my = 1-y;
    double mz = 1-z;

    doublev8 result_v,x_g_v,y_g_v,z_g_v,f_v,x_v,y_v,z_v;
    f_v = simd_set_doublev8(f000,f100,f010,f110,f001,f101,f011,f111); 
    //simd_print_doublev8(f8);
    x_v = simd_set_doublev8(mx,x,mx,x,mx,x,mx,x);
    y_v = simd_set_doublev8(my,my,y,y,my,my,y,y);
    z_v = simd_set_doublev8(mz,mz,mz,mz,z,z,z,z);
    result_v = z_v*y_v*x_v*f_v;
    //simd_print_doublev8(result_v);
    double f=simd_reduc_plusd(result_v);
    if(deriv) {
 //   e_x = simd_set_doublev8(-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0);
   // e_y = simd_set_doublev8(-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0);
    //e_z = simd_set_doublev8(-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0);
    x_g_v = z_v*y_v*e_x*f_v;
    y_g_v = z_v*e_y*x_v*f_v;
    z_g_v = e_z*y_v*x_v*f_v;
    double x_g=simd_reduc_plusd(x_g_v);
    double y_g=simd_reduc_plusd(y_g_v);
    double z_g=simd_reduc_plusd(z_g_v);
/*
    double f =
        f000 *  mx * my * mz  +
        f100 *   x * my * mz  +
        f010 *  mx *  y * mz  +
        f110 *   x *  y * mz  +
        f001 *  mx * my *  z  +
        f101 *   x * my *  z  +
        f011 *  mx *  y *  z  +
        f111 *   x *  y *  z  ;
    //if(abs(f_tmp-f)>1e-6){printf("f_tmp:%.6lf f:%.6lf dif:%lf\n",f_tmp,f,f_tmp-f);}
    //printf("nonon:%.6lf =%.6lf=%.6lf=%.6lf=%.6lf=%.6lf=%.6lf=%.6lf\n",f000 *  mx * my * mz,f100 *   x * my * mz,f010 *  mx *  y * mz,f110 *   x *  y * mz,f001 *  mx * my *  z,f101 *   x * my *  z,f011 *  mx *  y *  z,f111 *   x *  y *  z);
    if(deriv) { // valid pointer
        double x_g =
            f000 * (-1)* my * mz  +
            f100 *   1 * my * mz  +
            f010 * (-1)*  y * mz  +
            f110 *   1 *  y * mz  +
            f001 * (-1)* my *  z  +
            f101 *   1 * my *  z  +
            f011 * (-1)*  y *  z  +
            f111 *   1 *  y *  z  ;


        double y_g =
            f000 *  mx *(-1)* mz  +
            f100 *   x *(-1)* mz  +
            f010 *  mx *  1 * mz  +
            f110 *   x *  1 * mz  +
            f001 *  mx *(-1)*  z  +
            f101 *   x *(-1)*  z  +
            f011 *  mx *  1 *  z  +
            f111 *   x *  1 *  z  ;


        double z_g =
            f000 *  mx * my *(-1) +
            f100 *   x * my *(-1) +
            f010 *  mx *  y *(-1) +
            f110 *   x *  y *(-1) +
            f001 *  mx * my *  1  +
            f101 *   x * my *  1  +
            f011 *  mx *  y *  1  +
            f111 *   x *  y *  1  ;
*/
        //vec gradient(x_g, y_g, z_g);
    //if(abs(x_g_tmp-x_g)>1e-6){printf("x_tmp:%.6lf x:%.6lf dif:%lf\n",x_g_tmp,x_g,x_g_tmp-x_g);}
    //if(abs(y_g_tmp-y_g)>1e-6){printf("y_tmp:%.6lf y:%.6lf dif:%lf\n",y_g_tmp,y_g,y_g_tmp-y_g);}
    //if(abs(z_g_tmp-z_g)>1e-6){printf("z_tmp:%.6lf z:%.6lf dif:%lf\n",z_g_tmp,z_g,z_g_tmp-z_g);}
        double gradient[3] = {x_g, y_g, z_g};
        curl(&f, gradient, v);
        double gradient_everywhere[3];
/*        for(int i=0; i<3; i++) {
            gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
            deriv[i] = m_factor[i] * gradient_everywhere[i] + slope * region[i];
        }
*/
	gradient_everywhere[0] = ((region[0] == 0) ? gradient[0] : 0);
	gradient_everywhere[1] = ((region[1] == 0) ? gradient[1] : 0);
	gradient_everywhere[2] = ((region[2] == 0) ? gradient[2] : 0);
	deriv[0] = m_factor[0] * gradient_everywhere[0] + slope * region[0];
	deriv[1] = m_factor[1] * gradient_everywhere[1] + slope * region[1];
	deriv[2] = m_factor[2] * gradient_everywhere[2] + slope * region[2];

	return f+penalty;
    }
    else {
        curl_two(&f, v);
        return f + penalty;
    }
}


//fl evaluate(const vec& location, fl slope, fl c)   const { return evaluate_aux(location, slope, c, NULL); }
double grid_evaluate_three(double *grid, double location[3], double slope, double c){
    return evaluate_aux(grid, location, slope, c, NULL);
}

/* fl evaluate(const vec& location, fl slope, fl c, vec& deriv)
 * const { return evaluate_aux(location, slope, c, &deriv); }
 * sets deriv
 */

double grid_evaluate_four(double *grid, double location[3], double slope, double c, double deriv[3]){
    return evaluate_aux(grid, location, slope, c, deriv);
}


double cache_eval_deriv(double v) { // needs m.coords, sets m.minus_forces
    double e = 0;
    //int nat = num_atom_types(globalTm->cachegrids->atu);
    int nat = num_atom_types(cachegrids_slave.atu);

/*    if(cache_flag==1)
	printf("cachegrid:%d\n",my_id);
    cache_flag=0;
*/
    for(int i=0; i<num_movable_atoms(); i++){
        //we transfer the atoms to all double pointer for uniformity
        double *atom = &ligandAtom_slave[i*37];
        //sz t = a.get(atu);
        //int t = get(globalTm->cachegrids->atu, atom);
        int t = get(cachegrids_slave.atu, atom);
        
        if(t >= nat) {
            //m.minus_forces[i].assign(0);
            minus_forces_slave[i*3] = 0;
            minus_forces_slave[i*3+1] = 0;
            minus_forces_slave[i*3+2] = 0;
            continue;
        }
        //const grid& g = grids[t];
       	//int index = para_s->grid_index[t];
        int gridsize = para_s->singlegridsize;


        //luhao modify
	//offset = t*gridsize*8;
	//printf("cache——index:%d\n",index);
	//double *g = &grids_slave[index*18];;
	//printf("index*18--:%d\n",index*18);
	double *g = &globalTm->cachegrids->grids[t*gridsize*8];
        double deriv[3];
        //e += g.evaluate(m.coords[i], slope, v, deriv);
        //e += grid_evaluate_four(g, &ligcoords_slave[i*3], globalTm->cachegrids->slope, v, deriv);
        e += grid_evaluate_four(g, &ligcoords_slave[i*3], cachegrids_slave.slope, v, deriv);
	
        //m.minus_forces[i] = deriv;
        minus_forces_slave[i*3] = deriv[0];
        minus_forces_slave[i*3+1] = deriv[1];
        minus_forces_slave[i*3+2] = deriv[2];
    }
    return e;
}
