/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "monte_carlo.h"
#include "coords.h"
#include "mutate.h"
#include "quasi_newton.h"
#ifdef ATHREAD_MC
extern "C"{

//#include <crts.h>
#include <athread.h>
//void slave_monteCarlo(void* para);
void slave_monteCarlo(void);
}
//extern char athread_init_called;
static int athread_init_called_flag = 0;
struct tranModel *globalTm;
//struct output_type_s *output_m;
//double *output_m_torsPlsCoord;
int *finished_sign;
int *count_step_sign;
int *result_tmp_sign;
struct output_type_s *result_tmp;
double *result_tmp_torsPlsCoord;
//double *grids_array;   //for grids, 64 grids;
#endif

#ifdef ATHREAD_MC
output_type monte_carlo::operator()(model& m, tranModel *tm, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
#else
output_type monte_carlo::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
#endif	
	output_container tmp;
#ifdef ATHREAD_MC
	this->operator()(m, tm, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator); // call the version that produces the whole container
#else
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator); // call the version that produces the whole container
#endif
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

void monte_carlo::single_run(model& m, output_type& out, const precalculate& p, const igrid& ig, rng& generator) const {
	conf_size s = m.get_size();
	change g(s);
	vec authentic_v(1000, 1000, 1000);
	out.e = max_fl;
	output_type current(out);
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	VINA_U_FOR(step, num_steps) {
		output_type candidate(current.c, max_fl);
		mutate_conf(candidate.c, m, mutation_amplitude, generator);
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		if(step == 0 || metropolis_accept(current.e, candidate.e, temperature, generator)) {
			quasi_newton_par(m, p, ig, candidate, g, authentic_v);
			current = candidate;
			if(current.e < out.e)
				out = current;
		}
	}
	quasi_newton_par(m, p, ig, out, g, authentic_v);
}

void monte_carlo::many_runs(model& m, output_container& out, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const {
	conf_size s = m.get_size();
	VINA_FOR(run, num_runs) {
		output_type tmp(s, 0);
		tmp.c.randomize(corner1, corner2, generator);
		single_run(m, tmp, p, ig, generator);
		out.push_back(new output_type(tmp));
	}
	out.sort();
}

output_type monte_carlo::many_runs(model& m, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const {
	output_container tmp;
	many_runs(m, tmp, p, ig, corner1, corner2, num_runs, generator);
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}


// out is sorted
#ifdef ATHREAD_MC
void monte_carlo::operator()(model& m, tranModel *tm, output_container& out, const precalculate& p, const igrid& ig,
        const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
#else
void monte_carlo::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
#endif

	#ifdef ATHREAD_MC
	globalTm = tm;
/*    	output_m = (struct output_type_s*)malloc((sizeof(struct output_type_s) +
            	sizeof(double)*(tm->para->coordsize*3+tm->para->torsionsize))*num_saved_mins);
*/      int slave_finshed_num = 0,tmp_index = 0;  	
	int slave_num = CRTS_athread_get_max_threads();
	int torsions_size = tm->para->torsionsize;
	int coords_size = tm->para->coordsize;
	int confsize = coords_size * 3 + torsions_size;
	//int singlegridsize = globalTm->para->singlegridsize;

	//lu: output_m size should be 20	
	//struct output_type_s output_m_arr[20 * slave_num];
	//output_m = output_m_arr;
       
	//finish sign
	int finish_flag = 0;
	finished_sign = &finish_flag;

	//count num_step 
	long sum_count_step = 0;
	int count_step_arr[slave_num]={0};
	count_step_sign = count_step_arr;
	//wcj: all slave result_tmp 
	int result_tmp_sign_arr[slave_num]={0};
	struct output_type_s result_tmp_arr[slave_num];
        double result_tmp_torsPlsCoord_arr[confsize * slave_num];	
	result_tmp_sign = result_tmp_sign_arr;
	result_tmp = result_tmp_arr;
	result_tmp_torsPlsCoord = result_tmp_torsPlsCoord_arr;	

	//creat tmp 	
	conf_size s = m.get_size();
        output_type tmp(s, 0);
	vec tmp_coord;
	//printf("coords_size:%d\n",coords_size);
        for(int k=0; k<coords_size; k++){
		tmp.coords.push_back(tmp_coord);
	}


/*	//lu: for speed,there are 64 grids
	int real_gridsize=globalTm->para->realgridsize;
	double grids_array_ent[slave_num*singlegridsize*real_gridsize];
	for(int i=0; i<slave_num; i++)
	    for(int j=0; j<singlegridsize*real_gridsize; j++)
	    	grids_array_ent[i*singlegridsize*real_gridsize + j] = globalTm->cachegrids->grids[j];
*/
	
/*	for(int i=0; i<slave_num; i++)
	    memcpy(&grids_array_ent[i*singlegridsize*5], globalTm->cachegrids->grids, sizeof(double)*singlegridsize*5);
*/	
	//printf("monetor start..\n");
	//grids_array=grids_array_ent;

        //lu: for percaculate
        globalTm->para->pre_n = p.element_size();
	globalTm->para->matsize = p.matrix_size() * (p.matrix_size()+1) /2;
//while(1);

	// lu：output_m torsPlsCoord
	//double output_m_torsPlsCoord11[20 * confsize * slave_num];
	//output_m_torsPlsCoord = output_m_torsPlsCoord11;
	//m.print_stuff();
	fl new_min_rmsd = 2.0;
	sz new_num_saved_mins = 20;
	//printf("min_rmsd, num_saved_mins:%lf=%d\n",new_min_rmsd, new_num_saved_mins);
	m.initializepara(p);
		int ret;
	//	athread_init_called = 0;
	        if(athread_init_called_flag == 0)
		{
       		    athread_init();
	            athread_init_called_flag = 1;
   		}
		//long int start=athread_time_cycle();
		
    		//ret = __real_athread_spawn((void *)slave_monteCarlo, tm->para, 0);
    		athread_spawn((void *)slave_monteCarlo, 0);
	//	ret = __real_athread_spawn((void *)slave_monteCarlo, 0, 0);
		//ret = __real_athread_spawn(0, (void *)slave_monteCarlo, 0);
		
		//ret = CRTS_athread_create(0, (void *)slave_monteCarlo,0);

		//wcj: 1.receive tmp, and add tmp to output_container
		//wcj: 2.receive slave finished signal
		while(1){
			if(result_tmp_sign[tmp_index]==1){
//master: position[3] and orientation[4]
				tmp.e = result_tmp[tmp_index].e;
				//printf("slave %d tmp.e:%lf\n",tmp_index,tmp.e);
				tmp.c.ligands[0].rigid.position[0]=result_tmp[tmp_index].position[0];
				tmp.c.ligands[0].rigid.position[1]=result_tmp[tmp_index].position[1];	
				tmp.c.ligands[0].rigid.position[2]=result_tmp[tmp_index].position[2];	
				qt qt_s(result_tmp[tmp_index].orientation[0], result_tmp[tmp_index].orientation[1], result_tmp[tmp_index].orientation[2],result_tmp[tmp_index].orientation[3]);
            			tmp.c.ligands[0].rigid.orientation = qt_s;
//matser: torsions
            			VINA_FOR(k, torsions_size){
                			tmp.c.ligands[0].torsions[k] = result_tmp_torsPlsCoord[confsize * tmp_index + k];
				}
//master: coords[n][3]
				
				VINA_FOR(k, coords_size){
            				tmp.coords[k][0] = result_tmp_torsPlsCoord[confsize * tmp_index + torsions_size + 3*k];
            				tmp.coords[k][1] = result_tmp_torsPlsCoord[confsize * tmp_index + torsions_size + 3*k+1];
            				tmp.coords[k][2] = result_tmp_torsPlsCoord[confsize * tmp_index + torsions_size + 3*k+2];

        			}
				add_to_output_container(out, tmp, new_min_rmsd, num_saved_mins); // 20 - max size                   
				//printf("out size:%d\n",out.size());
				result_tmp_sign[tmp_index] = 0;
			}
/*			if(result_tmp_sign[tmp_index]==-1){
				//printf("%d : slave %d finished\n",slave_finshed_num+1,tmp_index);
				result_tmp_sign[tmp_index] = 0;
				slave_finshed_num++;
				if(slave_finshed_num == 64){break;}
			}
*/
			if(count_step_sign[tmp_index]==1){
				sum_count_step +=100;
				if(sum_count_step > num_steps*8){
					finish_flag = 1;
					break;
				}
				count_step_sign[tmp_index]=0;
			}
			tmp_index = (tmp_index + 1) % 64;
		}

		//check slaves really finished
		bool search_flag = true;
		while(search_flag){
			int count=0;
			for(int i=0;i<slave_num;i++){
				if(result_tmp_sign[i] = -1){
					count++;
				}
                        }
			if(count==slave_num) break;
		}
		/*
 * extern int *result_tmp_sign;
 * extern struct output_type_s *result_tmp;
 * extern double *result_tmp_torsPlsCoord;    
 * */
    		athread_join();
		//long int end=athread_time_cycle();
		//printf("athread_time_cycle:%ld\n",end-start);
//		std::cout<<"在从核中出来..."<<std::endl;
		//while(1);
	#endif
	//this is original
	/*
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);
	fl best_e = max_fl;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	VINA_U_FOR(step, num_steps) {
		if(increment_me)
			++(*increment_me);
		output_type candidate = tmp;
		mutate_conf(candidate.c, m, mutation_amplitude, generator);
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		if(step == 0 || metropolis_accept(tmp.e, candidate.e, temperature, generator)) {
			tmp = candidate;

			m.set(tmp.c); // FIXME? useless?

			// FIXME only for very promising ones
			if(tmp.e < best_e || out.size() < num_saved_mins) {
				quasi_newton_par(m, p, ig, tmp, g, authentic_v);
				m.set(tmp.c); // FIXME? useless?
				tmp.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
				if(tmp.e < best_e)
					best_e = tmp.e;
			}
		}
	}
	*/
	//
	#ifdef ATHREAD_MC
	//out and out.coords szie is 0, torsions is 0, fix it
/*	vec aa;
        for(int i=0; i<20 * slave_num; i++){
            struct conf confdemo;
            out.push_back(new output_type(confdemo, i));

            for(int k=0; k<tm->para->coordsize; k++){

                out[i].coords.push_back(aa);
            }

            struct ligand_conf ligand_conf_ent;
	    out[i].c.ligands.push_back(ligand_conf_ent);
	    
	    //this is for torsionsize
	    for(int j=0; j<tm->para->torsionsize; j++){
		double tori_tmp=1.0;
		out[i].c.ligands[0].torsions.push_back(tori_tmp);
	    }
        }

//		printf("+++++++++++--------%lf\n",output_m[0].e);
//matser: e
	for(int i=0; i<num_saved_mins * slave_num; i++){
	    out[i].e = output_m[i].e;

//master: position[3] and orientation[4]
//        VINA_FOR_IN(j, out[i].c.ligands){
	
            out[i].c.ligands[0].rigid.position[0] = output_m[i].position[0];
            out[i].c.ligands[0].rigid.position[1] = output_m[i].position[1];
            out[i].c.ligands[0].rigid.position[2] = output_m[i].position[2];

            qt qt_s(output_m[i].orientation[0], output_m[i].orientation[1], output_m[i].orientation[2],
                    output_m[i].orientation[3]);
            out[i].c.ligands[0].rigid.orientation = qt_s;

//matser: torsions
// yxm original
            sz torsize = out[i].c.ligands[0].torsions.size();
            VINA_FOR(k, torsize)
                out[i].c.ligands[0].torsions[k] = output_m[i].torsPlsCoord[k];
//
            sz torsize = out[i].c.ligands[0].torsions.size();
	    VINA_FOR(k, torsize)
		out[i].c.ligands[0].torsions[k] = output_m_torsPlsCoord[k + i*confsize];

	        	
//master: coords[n][3]
        //
        VINA_FOR(k, tm->para->coordsize){
	    std::cout<<"K : "<<k<<"---i:"<<i<<std::endl;
            out[i].coords[k][0] = output_m[i].torsPlsCoord[tm->para->torsionsize+3*k];
            out[i].coords[k][1] = output_m[i].torsPlsCoord[tm->para->torsionsize+3*k+1];
            out[i].coords[k][2] = output_m[i].torsPlsCoord[tm->para->torsionsize+3*k+2];
        }//
	VINA_FOR(k, tm->para->coordsize){
            out[i].coords[k][0] = output_m_torsPlsCoord[confsize*i + torsize + 3*k];
	    out[i].coords[k][1] = output_m_torsPlsCoord[confsize*i + torsize + 3*k+1];
            out[i].coords[k][2] = output_m_torsPlsCoord[confsize*i + torsize + 3*k+2];

        }
	
//
	for(int j=0; j<tm->para->torsionszie; j++){
         
	    out[i].c.ligands[0].torsions[]	
        }
//
	}
*/
/*	for(int i=0; i<slave_num; i++){
		printf("e:%lf\n",out[i*20].e);
	}
*/		
	out.sort();
	m.freepara(p);
	//free(output_m);

       //luhao
	//free(globalTm->receptorAtom);
	free(globalTm->ligandAtom);
	free(globalTm->ligcoords);
	free(globalTm->internal_coords);
	free(globalTm->minus_forces);
//	free(tranM->ligandpair);
	for(int i=0;i<globalTm->para->pairsize;i++)
		free(globalTm->ligandpair[i]);
	free(globalTm->ligandpair);
        free(globalTm->heritTree);
        free(globalTm->para);
        free(globalTm->cachegrids);
	free(globalTm);
	#endif
/*
	for(int i=0; i<num_saved_mins * slave_num; i++){
            printf("out-%d--.e:%lf\n",i,out[i].e);
	}
*/
	VINA_CHECK(!out.empty());
//	std::cout<<"front:"<<out.front().e<<"end:"<<out.back().e<<std::endl;
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order

	
}
