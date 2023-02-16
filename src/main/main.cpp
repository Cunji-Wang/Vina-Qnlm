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
#include <iostream>
#include <string>
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include <parse_pdbqt.h>
#include "parallel_mc.h"
#include "file.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "current_weights.h"
#include "quasi_newton.h"
#include "tee.h"
#include "coords.h" // add_to_output_container
#ifdef ATHREAD_MC
#include "tranModel.h"
bool refineStruct = false;
#endif

#include <boost/optional/optional_io.hpp>
#include "LigandList.h"
#include "LigandList.cpp"
#include "mpi.h"

//#define MANAGE_NUM 1 //manage node num,must match ligandlist num
using namespace std;

LigandList lgndsList;
string scoreLigandName;//use only_score
string result;

using boost::filesystem::path;
//using boost::program_options::positional_options_description;
//using boost::program_options;

//new
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>


path make_path(const std::string& str) {
	return path(str);
}

void doing(int verbosity, const std::string& str, tee& log) {
	if(verbosity > 1) {
		log << str << std::string(" ... ");
		log.flush();
	}
}

void done(int verbosity, tee& log) {
	if(verbosity > 1) {
		log << "done.";
		log.endl();
	}
}
std::string default_output(const std::string& input_name) {
	std::string tmp = input_name;
	if(tmp.size() >= 6 && tmp.substr(tmp.size()-6, 6) == ".pdbqt")
		tmp.resize(tmp.size() - 6); // FIXME?
	return tmp + "_out.pdbqt";
}

void write_all_output(model& m, const output_container& out, sz how_many,
				  const std::string& output_name,
				  const std::vector<std::string>& remarks) {
	if(out.size() < how_many)
		how_many = out.size();
	VINA_CHECK(how_many <= remarks.size());
	//ofile f(make_path(output_name));
	VINA_FOR(i, how_many) {
		m.set(out[i].c);
		//m.write_model(f, i+1, remarks[i]); // so that model numbers start with 1
		m.write_model_new(i+1, remarks[i]);
	}
}

//Ximing resultinfo

//end resultinfo

void do_randomization(model& m,
					  const std::string& out_name,
					  const vec& corner1, const vec& corner2, int seed, int verbosity, tee& log) {
	conf init_conf = m.get_initial_conf();
	rng generator(static_cast<rng::result_type>(seed));
	if(verbosity > 1) {
		log << "Using random seed: " << seed;
		log.endl();
	}
	const sz attempts = 10000;
	conf best_conf = init_conf;
	fl best_clash_penalty = 0;
	VINA_FOR(i, attempts) {
		conf c = init_conf;
		c.randomize(corner1, corner2, generator);
		m.set(c);
		fl penalty = m.clash_penalty();
		if(i == 0 || penalty < best_clash_penalty) {
			best_conf = c;
			best_clash_penalty = penalty;
		}
	}
	m.set(best_conf);
	if(verbosity > 1) {
		log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
		log.endl();
	}
	m.write_structure(make_path(out_name));
}

void refine_structure(model& m, const precalculate& prec, non_cache& nc, output_type& out, const vec& cap, sz max_steps = 1000) {
	change g(m.get_size());
	quasi_newton quasi_newton_par;
	quasi_newton_par.max_steps = max_steps;
	const fl slope_orig = nc.slope;
	#ifdef ATHREAD_MC
//    m.initializepara(prec);
	#endif
	VINA_FOR(p, 5) {
		nc.slope = 100 * std::pow(10.0, 2.0*p);
		quasi_newton_par(m, prec, nc, out, g, cap);
		m.set(out.c); // just to be sure
		if(nc.within(m))
			break;
	}
	out.coords = m.get_heavy_atom_movable_coords();
	if(!nc.within(m))
		out.e = max_fl;
	nc.slope = slope_orig;
	#ifdef ATHREAD_MC
//    m.freepara(prec);
    
	#endif
}

std::string vina_remark(fl e, fl lb, fl ub) {
	std::ostringstream remark;
	remark.setf(std::ios::fixed, std::ios::floatfield);
	remark.setf(std::ios::showpoint);
	remark << "REMARK VINA RESULT: " 
		                << std::setw(9) << std::setprecision(1) << e
	                    << "  " << std::setw(9) << std::setprecision(3) << lb
						<< "  " << std::setw(9) << std::setprecision(3) << ub
						<< '\n';
	return remark.str();
}

output_container remove_redundant(const output_container& in, fl min_rmsd) {
	output_container tmp;
	VINA_FOR_IN(i, in)
		add_to_output_container(tmp, in[i], min_rmsd, in.size());
	return tmp;
}

#ifdef ATHREAD_MC
void do_search(model& m, const boost::optional<model>& ref, const scoring_function& sf, const precalculate& prec, const igrid& ig, const precalculate& prec_widened, const igrid& ig_widened, non_cache& nc, // nc.slope is changed
			   const std::string& out_name,
			   const vec& corner1, const vec& corner2,
			   const parallel_mc& par, fl energy_range, sz num_modes,
			   int seed, int verbosity, bool score_only, bool local_only, tee& log, const terms& t, const flv& weights, tranModel *tm) {

#else
void do_search(model& m, const boost::optional<model>& ref, const scoring_function& sf, const precalculate& prec, const igrid& ig, const precalculate& prec_widened, const igrid& ig_widened, non_cache& nc, // nc.slope is changed
			   const std::string& out_name,
			   const vec& corner1, const vec& corner2,
			   const parallel_mc& par, fl energy_range, sz num_modes,
			   int seed, int verbosity, bool score_only, bool local_only, tee& log, const terms& t, const flv& weights) {
#endif
	conf_size s = m.get_size();
	conf c = m.get_initial_conf();
	fl e = max_fl;
	const vec authentic_v(1000, 1000, 1000);
	if(score_only) {
		fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, c);
		naive_non_cache nnc(&prec); // for out of grid issues
		e = m.eval_adjusted(sf, prec, nnc, authentic_v, c, intramolecular_energy);

		log << "Intramolecular energy: " << std::fixed << std::setprecision(5)
		<< intramolecular_energy << "\n";
		log <<scoreLigandName<<" "<< "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
		log.endl();
		std::vector<flv> atom_values;
		flv term_values = t.evale_robust(m,atom_values);
		VINA_CHECK(term_values.size() == 4);
		log << "Intermolecular contributions to the terms, before weighting:\n";
		log << std::setprecision(5);
		//log << "    gauss 1     : " << term_values[0] << '\n';
		log << "    gauss       : " << term_values[0] << '\n';
		log << "    repulsion   : " << term_values[1] << '\n';
		log << "    hydrophobic : " << term_values[2] << '\n';
		log << "    Hydrogen    : " << term_values[3] << '\n';
		VINA_CHECK(weights.size() == term_values.size() + 1);
		fl e2 = 0;
		VINA_FOR_IN(i, term_values)
			e2 += term_values[i] * weights[i];
		e2 = sf.conf_independent(m, e2);
		if(e < 100 && std::abs(e2 - e) > 0.05) {
			log << "WARNING: the individual terms are inconsisent with the\n";
			log << "WARNING: affinity. Consider reporting this as a bug:\n";
			log << "WARNING: http://vina.scripps.edu/manual.html#bugs\n";
		}

		//get atom values
		vecv coords = m.get_ligand_coords();
		assert(atom_values.size() == coords.size());
		for (unsigned i = 0, n = atom_values.size(); i < n; i++) {
			//VINA_CHECK(weights.size() == term_values.size() + 1);
			log << m.ligand_atom_str(i) << " ";
			log <<"<"<< coords[i][0] <<"," <<coords[i][1] <<"," <<coords[i][2]<<">";
			for (unsigned j = 0, m= atom_values[i].size(); j < m; j++){
				log << " " << atom_values[i][j]*weights[j];
			}
			log << "\n";
		}
		log << "END\n";
	}
	else if(local_only) {
		output_type out(c, e);
		doing(verbosity, "Performing local search", log);
		refine_structure(m, prec, nc, out, authentic_v, par.mc.ssd_par.evals);
		done(verbosity, log);
		fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out.c);
		e = m.eval_adjusted(sf, prec, nc, authentic_v, out.c, intramolecular_energy);

		log <<scoreLigandName<<" "<< "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
		log.endl();
		if(!nc.within(m))
			log << "WARNING: not all movable atoms are within the search space\n";

		doing(verbosity, "Writing output", log);
		output_container out_cont;
		out_cont.push_back(new output_type(out));
		std::vector<std::string> remarks(1, vina_remark(e, 0, 0));
		write_all_output(m, out_cont, 1, out_name, remarks); // how_many == 1
		done(verbosity, log);
	}
	else {
		rng generator(static_cast<rng::result_type>(seed));
		//log << "Using random seed: " << seed;
		//log.endl();
		output_container out_cont;
		//doing(verbosity, "Performing search", log);
		//model tmp_model=m;
		double search_startTime = MPI_Wtime();
		#ifdef ATHREAD_MC
		par(m, out_cont, prec, ig, prec_widened, ig_widened, corner1, corner2, generator,tm);
		#else
		par(m, out_cont, prec, ig, prec_widened, ig_widened, corner1, corner2, generator);
		
		#endif
                double search_endTime = MPI_Wtime();
                printf("Search time : %.2lf seconds\n",search_endTime-search_startTime);

		//done(verbosity, log);

		//doing(verbosity, "Refining results", log);
		
		VINA_FOR_IN(i, out_cont){
			//printf("e:%lf\n",out_cont[i].e);
			refine_structure(m, prec, nc, out_cont[i], authentic_v, par.mc.ssd_par.evals);
		}
		if(!out_cont.empty()) {
			out_cont.sort();
			const fl best_mode_intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out_cont[0].c);
			VINA_FOR_IN(i, out_cont)
				if(not_max(out_cont[i].e))
					out_cont[i].e = m.eval_adjusted(sf, prec, nc, authentic_v, out_cont[i].c, best_mode_intramolecular_energy); 
			// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
			out_cont.sort();
		}
		const fl out_min_rmsd = 1;
		out_cont = remove_redundant(out_cont, out_min_rmsd);
		//modify:error return
		if(out_cont.empty()) {
			printf("not found model!");
			return;
		}

		//done(verbosity, log);

		//log.setf(std::ios::fixed, std::ios::floatfield);
		//log.setf(std::ios::showpoint);
		//log << '\n';
		//log << "mode |   affinity | dist from best mode\n";
		//log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
		//log << "-----+------------+----------+----------\n";

  	
		model best_mode_model = m;

		if(!out_cont.empty())
			best_mode_model.set(out_cont.front().c);
		sz how_many = 0;
		std::vector<std::string> remarks;
		VINA_FOR_IN(i, out_cont) {
			if(how_many >= num_modes || !not_max(out_cont[i].e) || out_cont[i].e > out_cont[0].e + energy_range) break; // check energy_range sanity FIXME
			++how_many;
			//log << std::setw(4) << i+1
			//	<< "    " << std::setw(9) << std::setprecision(1) << out_cont[i].e; // intermolecular_energies[i];
			m.set(out_cont[i].c);
			const model& r = ref ? ref.get() : best_mode_model;
			const fl lb = m.rmsd_lower_bound(r);
			const fl ub = m.rmsd_upper_bound(r);
			//log << "  " << std::setw(9) << std::setprecision(3) << lb
			//    << "  " << std::setw(9) << std::setprecision(3) << ub; // FIXME need user-readable error messages in case of failures

			remarks.push_back(vina_remark(out_cont[i].e, lb, ub));
			//log.endl();
		}
		//doing(verbosity, "Writing output", log);
		write_all_output(m, out_cont, how_many, out_name, remarks);
		//done(verbosity, log);

		if(how_many < 1) {
			log << "WARNING: Could not find any conformations completely within the search space.\n"
				<< "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
			log.endl();
		}
	}
}

#ifdef ATHREAD_MC
void main_procedure(model& m, const boost::optional<model>& ref, // m is non-const (FIXME?)
			     const std::string& out_name,
				 bool score_only, bool local_only, bool randomize_only, bool no_cache,
				 const grid_dims& gd, int exhaustiveness,
				 const flv& weights,
				 int cpu, int seed, int verbosity, sz num_modes, fl energy_range, tee& log, tranModel *tm) {
#else
void main_procedure(model& m, const boost::optional<model>& ref, // m is non-const (FIXME?)
			     const std::string& out_name,
				 bool score_only, bool local_only, bool randomize_only, bool no_cache,
				 const grid_dims& gd, int exhaustiveness,
				 const flv& weights,
				 int cpu, int seed, int verbosity, sz num_modes, fl energy_range, tee& log) {
#endif

	//doing(verbosity, "Setting up the scoring function", log);

	everything t;
	VINA_CHECK(weights.size() == 5);

	weighted_terms wt(&t, weights);
	precalculate prec(wt);
	const fl left  = 0.25; 
	const fl right = 0.25;
	precalculate prec_widened(prec); prec_widened.widen(left, right);

	//done(verbosity, log);

	vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
	vec corner2(gd[0].end,   gd[1].end,   gd[2].end);
	#ifdef ATHREAD_MC
    tm->para->corner1[0] = corner1[0];
    tm->para->corner1[1] = corner1[1];
    tm->para->corner1[2] = corner1[2];
    tm->para->corner2[0] = corner2[0];
    tm->para->corner2[1] = corner2[1];
    tm->para->corner2[2] = corner2[2];
/*
    printf("corner1: %lf, %lf, %lf\n", tm->para->corner1[0], tm->para->corner1[1], tm->para->corner1[2]);
    printf("corner2: %lf, %lf, %lf\n", tm->para->corner2[0], tm->para->corner2[1], tm->para->corner2[2]);
*/	
	#endif
	parallel_mc par;
	sz heuristic = m.num_movable_atoms() + 10 * m.get_size().num_degrees_of_freedom();
	par.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
	par.mc.ssd_par.evals = unsigned((25 + m.num_movable_atoms()) / 3);
	par.mc.min_rmsd = 1.0;
	par.mc.num_saved_mins = 20;
	par.mc.hunt_cap = vec(10, 10, 10);
	par.num_tasks = exhaustiveness;
	par.num_threads = cpu;
	par.display_progress = (verbosity > 1);
	#ifdef ATHREAD_MC
	//get mc parameters
    tm->para->bfgs_steps = par.mc.ssd_par.evals;
    tm->para->mc_steps = par.mc.num_steps;
    tm->para->temperature = par.mc.temperature;
    tm->para->min_rmsd = par.mc.min_rmsd;
    tm->para->mutation_amplitude = par.mc.mutation_amplitude;
    tm->para->num_saved_mins = par.mc.num_saved_mins;
    tm->para->hunt_cap[0] = 10;
    tm->para->hunt_cap[1] = 10;
    tm->para->hunt_cap[2] = 10;
	#endif
	const fl slope = 1e6; // FIXME: too large? used to be 100
	if(randomize_only) {
		do_randomization(m, out_name,
			             corner1, corner2, seed, verbosity, log);
	}
	else {
		non_cache nc        (m, gd, &prec,         slope); // if gd has 0 n's, this will not constrain anything
		non_cache nc_widened(m, gd, &prec_widened, slope); // if gd has 0 n's, this will not constrain anything
		if(no_cache) {
		#ifdef ATHREAD_MC
			do_search(m, ref, wt, prec, nc, prec_widened, nc_widened, nc,
					  out_name,
					  corner1, corner2,
					  par, energy_range, num_modes,
					  seed, verbosity, score_only, local_only, log, t, weights, tm);
		#else
			do_search(m, ref, wt, prec, nc, prec_widened, nc_widened, nc,
					  out_name,
					  corner1, corner2,
					  par, energy_range, num_modes,
					  seed, verbosity, score_only, local_only, log, t, weights);
		#endif
		}
		else {
			bool cache_needed = !(score_only || randomize_only || local_only);
			if(cache_needed) ;//doing(verbosity, "Analyzing the binding site", log);
			cache c("scoring_function_version001", gd, slope, atom_type::XS);
			if(cache_needed){ c.populate(m, prec, m.get_movable_atom_types(prec.atom_typing_used()));
			#ifdef ATHREAD_MC
			//initialize tm->cachegrids
            std::vector<grid> mygrids = c.get_grids();   //grids
            sz gridsize = mygrids.size();                //有多少个grid
            sz grid_datasize;  //有多少数据
	    //int real_gridsize=0;
	    for(int i=0;i<gridsize;i++){
		if(mygrids[i].initialized()){
		    //real_gridsize++;
		    grid_datasize = mygrids[i].m_data.get_datasize();
	        }
	    }
            tm->cachegrids = (struct cachegrid*)malloc(sizeof(struct cachegrid) + sizeof(double) *gridsize * grid_datasize * 8);
            // gd should be gd[3][3], element is: double begin, double end, int n
           // tm->cachegrids->gd = (double *)malloc(sizeof(double)*9);
            for(int i=0; i<3; i++){
                tm->cachegrids->gd[i*3] = gd[i].begin;
                tm->cachegrids->gd[i*3+1] = gd[i].end;
                tm->cachegrids->gd[i*3+2] = gd[i].n;
                //printf("gd[%d] begin: %lf, end: %lf, n: %lf\n", i, tm->cachegrids->gd[i*3],
                //        tm->cachegrids->gd[i*3+1], tm->cachegrids->gd[i*3+2]);
            }
            tm->cachegrids->slope = slope;
            tm->cachegrids->atu = XS;
            //printf("start to initialize grids\n");
            /* vec m_init; vec m_range; vec m_factor; vec m_dim_fl_minus_1; vec m_factor_inv; total have 15 elements
             * array3d<fl> m_data: sz m_i, m_j, m_k, std::vector<double> m_data; the m_data size is dynamic allocated
             * total 18 element + m_data.size
            */
            tm->para->singlegridsize = grid_datasize;
            for(int i=0;i<gridsize;i++){
                if(mygrids[i].initialized()){
	            for(int j=0; j<3; j++){
        	        tm->cachegrids->m_init[j] = mygrids[i].m_init[j];
        	        tm->cachegrids->m_range[j] = mygrids[i].m_range[j];
        	        tm->cachegrids->m_factor[j] = mygrids[i].m_factor[j];
        	        tm->cachegrids->m_dim_fl_minus_1[j] = mygrids[i].m_dim_fl_minus_1[j];
        	        tm->cachegrids->m_factor_inv[j] = mygrids[i].m_factor_inv[j];
		    }
		    tm->cachegrids->m_grid_data_i = mygrids[i].m_data.dim0();
		    tm->cachegrids->m_grid_data_j = mygrids[i].m_data.dim1();
		    tm->cachegrids->m_grid_data_k = mygrids[i].m_data.dim2();
//printf("i:%d j:%d k:%d\n",tm->cachegrids->m_grid_data_i,tm->cachegrids->m_grid_data_j,tm->cachegrids->m_grid_data_k);
		    break;
		}
	    }
	    //tm->para->realgridsize = real_gridsize;
	    //int index=0;
            for(int index=0; index<gridsize; index++){ //gridsize is 19
                if(mygrids[index].initialized()) {
                    //VINA_FOR(j, grid_datasize){
		    VINA_FOR(i, tm->cachegrids->m_grid_data_i){
			VINA_FOR(j, tm->cachegrids->m_grid_data_j){
			   VINA_FOR(k, tm->cachegrids->m_grid_data_k){
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+0] = mygrids[index].m_data(i,j,k);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+1] = mygrids[index].m_data(i+1,j,k);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+2] = mygrids[index].m_data(i,j+1,k);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+3] = mygrids[index].m_data(i+1,j+1,k);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+4] = mygrids[index].m_data(i,j,k+1);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+5] = mygrids[index].m_data(i+1,j,k+1);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+6] = mygrids[index].m_data(i,j+1,k+1);
	                    tm->cachegrids->grids[index*grid_datasize*8 + (i + tm->cachegrids->m_grid_data_i*(j + tm->cachegrids->m_grid_data_j*k))*8+7] = mygrids[index].m_data(i+1,j+1,k+1);
		           }
			}
		    }
                }
            }
		//printf("grid success----\n");
        } 
		    //
		    //tm->para->grid_index[i]=index;
                    //printf("mygrids index: %d--:%d\n",i,index);
/*
                    VINA_FOR(j, 3)
                    {tm->cachegrids->grids[index * singlegridsize + j] = mygrids[i].m_init[j];
		    }
                    VINA_FOR(j, 3)
                    tm->cachegrids->grids[index * singlegridsize + 3 + j] = mygrids[i].m_range[j];
                    VINA_FOR(j, 3)
                    tm->cachegrids->grids[index * singlegridsize + 6 + j] = mygrids[i].m_factor[j];
                    VINA_FOR(j, 3)
                    tm->cachegrids->grids[index * singlegridsize + 9 + j] = mygrids[i].m_dim_fl_minus_1[j];
                    VINA_FOR(j, 3)
                    tm->cachegrids->grids[index * singlegridsize + 12 + j] = mygrids[i].m_factor_inv[j];
                    // m_i, m_j, m_k

                    tm->cachegrids->grids[index * singlegridsize + 15] = (double) mygrids[i].m_data.dim0();
                    tm->cachegrids->grids[index * singlegridsize + 16] = (double) mygrids[i].m_data.dim1();
                    tm->cachegrids->grids[index * singlegridsize + 17] = (double) mygrids[i].m_data.dim2();
*/
		    //index++;

		#endif
            if(cache_needed) ;//done(verbosity, log);
			#ifdef ATHREAD_MC
			do_search(m, ref, wt, prec, c, prec, c, nc,
					  out_name,
					  corner1, corner2,
				          par, energy_range, num_modes,
					  seed, verbosity, score_only, local_only, log, t, weights, tm);
			#else
			do_search(m, ref, wt, prec, c, prec, c, nc,
					  out_name,
					  corner1, corner2,
					  par, energy_range, num_modes,
					  seed, verbosity, score_only, local_only, log, t, weights);
			#endif
		}
	}
}

struct usage_error : public std::runtime_error {
	usage_error(const std::string& message) : std::runtime_error(message) {}
};

struct options_occurrence {
	bool some;
	bool all;
	options_occurrence() : some(false), all(true) {} // convenience
	options_occurrence& operator+=(const options_occurrence& x) {
		some = some || x.some;
		all  = all  && x.all;
		return *this;
	}
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
	options_occurrence tmp;
	VINA_FOR_IN(i, d.options()) 
		if(vm.count((*d.options()[i]).long_name())) 
			tmp.some = true;
		else 
			tmp.all = false;
	return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
	VINA_FOR_IN(i, d.options()) {
		const std::string& str = (*d.options()[i]).long_name();
		if(!vm.count(str))
			std::cerr << "Required parameter --" << str << " is missing!\n";
	}
}

model parse_bundle(const std::string& rigid_name, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names) {
	model tmp = (flex_name_opt) ? parse_receptor_pdbqt(make_path(rigid_name), make_path(flex_name_opt.get()))
		                        : parse_receptor_pdbqt(make_path(rigid_name));
	VINA_FOR_IN(i, ligand_names)
		tmp.append(parse_ligand_pdbqt(make_path(ligand_names[i])));
	return tmp;
}

model parse_bundle(const std::vector<std::string>& ligand_names) {
	VINA_CHECK(!ligand_names.empty()); // FIXME check elsewhere
	model tmp = parse_ligand_pdbqt(make_path(ligand_names[0]));
	VINA_RANGE(i, 1, ligand_names.size())
		tmp.append(parse_ligand_pdbqt(make_path(ligand_names[i])));
	return tmp;
}

model parse_bundle(const boost::optional<std::string>& rigid_name_opt, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names) {
	
	if(rigid_name_opt){
		return parse_bundle(rigid_name_opt.get(), flex_name_opt, ligand_names);
	}
	else
		return parse_bundle(ligand_names);
}

model parse_bundle_content(const std::string& rigid_name, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names, std::vector<std::string>& receptor_content,std::vector<std::string>& ligand_content) {
        model tmp = (flex_name_opt) ? parse_receptor_pdbqt(make_path(rigid_name), make_path(flex_name_opt.get()))
                                        : parse_receptor_pdbqt_content(make_path(rigid_name),receptor_content);
	VINA_FOR_IN(i, ligand_names)
                tmp.append(parse_ligand_pdbqt_content(make_path(ligand_names[i]),ligand_content));
        return tmp;
}

model parse_bundle_content(const std::vector<std::string>& ligand_names,std::vector<std::string>& ligand_content) {
        VINA_CHECK(!ligand_names.empty()); // FIXME check elsewhere
        model tmp = parse_ligand_pdbqt(make_path(ligand_names[0]));
        VINA_RANGE(i, 1, ligand_names.size())
                tmp.append(parse_ligand_pdbqt_content(make_path(ligand_names[i]),ligand_content));
        return tmp;
}
model parse_bundle_content(const boost::optional<std::string>& rigid_name_opt, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names, std::vector<std::string>& receptor_content,std::vector<std::string>& ligand_content) {

        if(rigid_name_opt){
                return parse_bundle_content(rigid_name_opt.get(), flex_name_opt, ligand_names,receptor_content,ligand_content);
        }
        else
                return parse_bundle_content(ligand_names,ligand_content);
}



int main_ori(int argc, char*argv[], std::string proten, std::string ligand, int rank, int num_caculted){
			std::string save_name;//log save_name
			try{
			//std::cout<<ligand<<std::endl;
			std::string rigid_name,ligand_name, flex_name, config_name, out_name, log_name;
			fl center_x, center_y, center_z, size_x, size_y, size_z;
			int cpu = 1, seed, exhaustiveness, verbosity = 2, num_modes = 9, manageNum;
			fl energy_range = 2.0;
			// -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 0.05846
			//fl weight_gauss1      = -0.035579;
			fl weight_gauss       = -0.045000;
			fl weight_repulsion   =  0.800000;
			fl weight_hydrophobic = -0.035000;
			fl weight_hydrogen    = -0.600000;
			fl weight_rot         =  0.000000;
			bool score_only = false, local_only = false, randomize_only = false, help = false, help_advanced = false, version = false; // FIXME

		//	positional_options_description positional; // remains empty
			boost::program_options::positional_options_description positional; // remains empty

			boost::program_options::options_description inputs("Input");
			inputs.add_options()
				("receptor",boost::program_options::value<std::string>(&rigid_name), "rigid part of the receptor (PDBQT)")
				("flex",boost::program_options::value<std::string>(&flex_name), "flexible side chains, if any (PDBQT)")
				("ligand", boost::program_options::value<std::string>(&ligand_name), "ligand (PDBQT)")
		                ("manageNum", boost::program_options::value<int>(&manageNum)->default_value(1), "the number of manage core ,to read ligandlist and manage work process")
			;
			
			
			//options_description search_area("Search area (required, except with --score_only)");
			boost::program_options::options_description search_area("Search space (required)");
			search_area.add_options()
				("center_x",boost::program_options::value<fl>(&center_x), "X coordinate of the center")
				("center_y", boost::program_options::value<fl>(&center_y), "Y coordinate of the center")
				("center_z",boost::program_options::value<fl>(&center_z), "Z coordinate of the center")
				("size_x", boost::program_options::value<fl>(&size_x), "size in the X dimension (Angstroms)")
				("size_y", boost::program_options::value<fl>(&size_y), "size in the Y dimension (Angstroms)")
				("size_z", boost::program_options::value<fl>(&size_z), "size in the Z dimension (Angstroms)")
			;
			
			//options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
			boost::program_options::options_description outputs("Output (optional)");
			outputs.add_options()
				("out", boost::program_options::value<std::string>(&out_name), "output models (PDBQT), the default is chosen based on the ligand file name")
				("log", boost::program_options::value<std::string>(&log_name), "optionally, write log file")
			;
			boost::program_options::options_description advanced("Advanced options (see the manual)");
			advanced.add_options()
				("score_only",     boost::program_options::bool_switch(&score_only),     "score only - search space can be omitted")
				("local_only",     boost::program_options::bool_switch(&local_only),     "do local search only")
				("randomize_only", boost::program_options::bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")
				//("weight_gauss1", value<fl>(&weight_gauss1)->default_value(weight_gauss1),                "gauss_1 weight")
				("weight_gauss", boost::program_options::value<fl>(&weight_gauss)->default_value(weight_gauss),                "gauss weight")
				("weight_repulsion", boost::program_options::value<fl>(&weight_repulsion)->default_value(weight_repulsion),       "repulsion weight")
				("weight_hydrophobic", boost::program_options::value<fl>(&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
				("weight_hydrogen", boost::program_options::value<fl>(&weight_hydrogen)->default_value(weight_hydrogen),          "Hydrogen bond weight")
				("weight_rot", boost::program_options::value<fl>(&weight_rot)->default_value(weight_rot),                         "N_rot weight")
			;
			
			boost::program_options::options_description misc("Misc (optional)");
			misc.add_options()
				("cpu", boost::program_options::value<int>(&cpu), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
				("seed", boost::program_options::value<int>(&seed), "explicit random seed")
				("exhaustiveness", boost::program_options::value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
				("num_modes", boost::program_options::value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
				("energy_range", boost::program_options::value<fl>(&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
			;
			boost::program_options::options_description config("Configuration file (optional)");
			config.add_options()
				("config", boost::program_options::value<std::string>(&config_name), "the above options can be put here")
			;
			boost::program_options::options_description info("Information (optional)");
			info.add_options()
				("help",          boost::program_options::bool_switch(&help), "display usage summary")
				("help_advanced", boost::program_options::bool_switch(&help_advanced), "display usage summary with advanced options")
				("version",       boost::program_options::bool_switch(&version), "display program version")
			;
			boost::program_options::options_description desc, desc_config, desc_simple;
			desc       .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
			desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
			desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);
			
			
			boost::program_options::variables_map vm;
			
			try {
				//std::cout<<"ligand = "<<argv[4]<<std::endl;
				//store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
				store(boost::program_options::command_line_parser(argc, argv)
					.options(desc)
					.style(boost::program_options::command_line_style::default_style ^ boost::program_options::command_line_style::allow_guessing)
					.positional(positional)
					.run(), 
					vm);
				notify(vm); 
				
			}
			catch(boost::program_options::error& e) {
				std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
				return 1;
			}
			
			if(vm.count("config")) {
				
				try {
					path name = make_path(config_name);
					ifile config_stream(name);
					//printf("debug2++++=\n");
					store(parse_config_file(config_stream, desc_config), vm);
					notify(vm);
				}
				catch(boost::program_options::error& e) {
					std::cerr << "Configuration file parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
					return 1;
				}
			}
			
			if(help) {
				std::cout << desc_simple << '\n';
				return 0;
			}
			if(help_advanced) {
				std::cout << desc << '\n';
				return 0;
			}
			if(version) {
				std::cout << "version 1." << '\n';
				return 0;
			}

			bool search_box_needed = !score_only; // randomize_only and local_only still need the search space
			bool output_produced   = !score_only; 
			bool receptor_needed   = !randomize_only;
			
/*	
			if(receptor_needed) {
				if(vm.count("receptor") <= 0) {
					std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
					return 1;
				}
			}
			if(vm.count("ligand") <= 0) {
				std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
				return 1;
			}
*/		    
			if(cpu < 1) 
				cpu = 1;
			if(vm.count("seed") == 0) 
				seed = auto_seed();
			if(exhaustiveness < 1)
				throw usage_error("exhaustiveness must be 1 or greater");
			if(num_modes < 1)
				throw usage_error("num_modes must be 1 or greater");
			sz max_modes_sz = static_cast<sz>(num_modes);
			
			boost::optional<std::string> rigid_name_opt;
			if(vm.count("receptor"))
				rigid_name_opt = rigid_name;

			boost::optional<std::string> flex_name_opt;
			if(vm.count("flex"))
				flex_name_opt = flex_name;

			if(vm.count("flex") && !vm.count("receptor"))
				throw usage_error("Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

			tee log;
			if(vm.count("log") > 0)
/*				log.init(log_name);

			if(search_box_needed) { 
				options_occurrence oo = get_occurrence(vm, search_area);
				if(!oo.all) {
					check_occurrence(vm, search_area);
					std::cerr << "\nCorrect usage:\n" << desc_simple << std::endl;
					return 1;
				}
				if(size_x <= 0 || size_y <= 0 || size_z <= 0)
					throw usage_error("Search space dimensions should be positive");
			}
*/
			//log << cite_message << '\n';

			if(search_box_needed && size_x * size_y * size_z > 27e3) {
				log << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
			}

			// if(score_only|local_only){
			// 	scoreLigandName.clear();
			// 	int length=strlen(lignadName);
   //                      	for(int i=0;i<length;i++){
   //                              	scoreLigandName+=lignadName[i];
   //                              	scoreLigandName[length]='\0';

   //                      }
			// }
			scoreLigandName = ligand_name;
			//if(rank == 444) std::cout<<proten<<std::endl;
			//printf("-------------------\n");
			//this is modify
			std::vector<std::string> proten_content;
			std::vector<std::string> ligand_content;
			boost::split(proten_content, proten, boost::is_any_of("\n"), boost::token_compress_on);
			boost::split(ligand_content, ligand, boost::is_any_of("\n"), boost::token_compress_on);
			
			while(proten_content[0]==""){
                                proten_content.erase(proten_content.begin());
                        }

			while(ligand_content[0]==""){
				ligand_content.erase(ligand_content.begin());
			}
			while(ligand_content[0]=="\r"){
                                ligand_content.erase(ligand_content.begin());
                        }
                        while(ligand_content[0]=="\n"){
                                ligand_content.erase(ligand_content.begin());
                        }	
			//ligand_name = ligand_content[0];
			ligand_name = "hello.txt";

			std::vector<std::string> configs;
			for(int i=0;i<7;i++){
				boost::split(configs, proten_content[i], boost::is_any_of("="), boost::token_compress_on);
				proten_content[i] = configs[1];
			}
			//rigid_name_opt = proten_content[0];
			rigid_name_opt = "vina.in";
			out_name = "useless";
                        //out_name =  std::to_string(rank/2048) + "/" + std::to_string(rank%2048);   //eg. 0/8.data
                        save_name = proten_content[0] + "&&&&" + ligand_content[0];   //eg. 10gs_proten.pdbqt_ZINC1001.pdbqt
			result = result + save_name + '\n';
			
			std::cout<<"save name : "<<save_name<<std::endl;
			center_x = atof(proten_content[1].c_str());
			center_y = atof(proten_content[2].c_str());
			center_z = atof(proten_content[3].c_str());
			size_x = atof(proten_content[4].c_str());
			size_y = atof(proten_content[5].c_str());
			size_z = atof(proten_content[6].c_str());
/*
			printf("ligand size:%d\n",ligand_content.size());
			for(int i=0;i<ligand_content.size(); i++){
                                std::cout<<ligand_content[i]<<std::endl;
                        }
                        printf("ligand==========================\n");
*/
			ligand_content.erase(ligand_content.begin());
			ligand_content.pop_back();
			for(int i=0;i<7;i++)
				proten_content.erase(proten_content.begin());
			proten_content.pop_back();
			//end modify
			//test code
			//std::cout<<"x:"<<center_x<<std::endl;
			//std::cout<<"name:"<<save_name<<std::endl;
			//std::cout<<"rank:"<<rank<<"|num:"<<num_caculted<<"|name:"<<save_name<<std::endl;
/*			
			for(int i=0;i<ligand_content.size(); i++){
				std::cout<<ligand_content[i]<<std::endl;
			}
			printf("==========================\n");
*/
/*
			for(int i=0;i<proten_content.size();i++){
				std::cout<<proten_content[i]<<std::endl;
			}
*/
/*
			for(int =7;i<ligand_content.size(); i++){
                                std::cout<<ligand_conten[i]<<std::endl;
                        }
*/
			//while(1);

			grid_dims gd; // n's = 0 via default c'tor

			flv weights;
			//weights.push_back(weight_gauss1);
			weights.push_back(weight_gauss);
			weights.push_back(weight_repulsion);
			weights.push_back(weight_hydrophobic);
			weights.push_back(weight_hydrogen);
			weights.push_back(weight_rot);
			//weights.push_back(5 * weight_rot / 0.1 - 1); // linearly maps onto a different range, internally. see everything.cpp

			if(search_box_needed) { 
				const fl granularity = 0.375;
				//vec span(size_x,   size_y,   size_z);
				//vec center(center_x, center_y, center_z);
				//modify 
				vec span(size_x,   size_y,   size_z);
				vec center(center_x, center_y, center_z);
				VINA_FOR_IN(i, gd) {
					gd[i].n = sz(std::ceil(span[i] / granularity));
					fl real_span = granularity * gd[i].n;
					gd[i].begin = center[i] - real_span/2;
					gd[i].end = gd[i].begin + real_span;
				}
			}
			

			//doing(verbosity, "Reading input", log);
			
			//model m       = parse_bundle(rigid_name_opt, flex_name_opt, std::vector<std::string>(1, ligand_name));
			double parse_startTime = MPI_Wtime();
			model m       = parse_bundle_content(rigid_name_opt, flex_name_opt, std::vector<std::string>(1, ligand_name),proten_content ,ligand_content);	
                        double parse_endTime = MPI_Wtime();
                        printf("\nParse bundle time : %.2lf seconds\n",parse_endTime-parse_startTime);

			#ifdef ATHREAD_MC
        	/* this receptor atom has element: sz el, ad, xs, sy, double charge, double coords[3]
        	* sz bondsize, int i, bool in_grid, bool rotatable, double length
        	* default, the max bond size is 4. so one atoms should have 9+4*4=25
        	*/
        	tranModel tranM = reform_transfer_model(m);

		int ligand_atom_num = m.get_ligand_coords().size();
		if(212 < (tranM.para->torsionsize + 2 * ligand_atom_num)){
		//if(175 < (tranM.para->torsionsize + 2 * ligand_atom_num)){
		//if(10 < tranM.para->torsionsize||63 < ligand_atom_num){
			result = result + "ligand size over big,not calculation!\n";
			//create directory and write out
                        string dir_ligandNotCalcu ="./ligandNotCalcu";
                        if(!boost::filesystem::exists(dir_ligandNotCalcu)){
                                boost::filesystem::path dir(dir_ligandNotCalcu);
                                if(boost::filesystem::create_directory(dir)){
                                }
                        }
                        string ligandNotCalcu_path = dir_ligandNotCalcu + "/" + std::to_string(rank/1024) + ".data";
			ofstream fout(ligandNotCalcu_path,ios::app);
			fout << save_name << "&&&&torsion:" << tranM.para->torsionsize <<endl;
			fout.close();
			//ofstream fout("ligandNotCalcu.txt",ios::app);
			//fout << save_name << ",torsion:" << tranM.para->torsionsize <<endl;
			//fout.close();
			return 0;
		}
			#endif
			boost::optional<model> ref;
			//done(verbosity, log);
			#ifdef ATHREAD_MC
			main_procedure(m, ref, 
					out_name,
					score_only, local_only, randomize_only, false, // no_cache == false
					gd, exhaustiveness,
					weights,
					cpu, seed, verbosity, max_modes_sz, energy_range, log, &tranM);
			#else
			main_procedure(m, ref, 
						out_name,
						score_only, local_only, randomize_only, false, // no_cache == false
						gd, exhaustiveness,
						weights,
						cpu, seed, verbosity, max_modes_sz, energy_range, log);
			#endif
			}
			catch(file_error& e) {
			std::cerr << "\n\nError: could not open \"" << e.name.filename().native() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
			return 1;
			}
			catch(boost::filesystem::filesystem_error& e) {
			std::cerr << "\n\nFile system error: " << e.what() << '\n';
			return 1;
			}
			catch(usage_error& e) {
			std::cerr << "\n\nUsage error: " << e.what() << ".\n";
			return 1;
			}
			catch(parse_error& e) {
			//std::cerr << "\n\nParse error on line " << e.line << " in file \"" << e.file.filename().native()  << "\": " << e.reason << '\n';
			return 1;
			}
			catch(std::bad_alloc&) {
			std::cerr << "\n\nError: insufficient memory!\n";
			return 1;
			}

			// Errors that shouldn't happen:

			catch(std::exception& e) { 
			std::cerr << "\n\nAn error occurred: " << e.what() << ". " << std::endl;
			return 1; 
			}
			catch(internal_error& e) {
			std::cerr << "\n\nAn internal error occurred in " << e.file << "(" << e.line << "). " << save_name << std::endl;
			return 1;
			}
			catch(...) {
			std::cerr << "\n\nAn unknown error occurred. " << save_name << std::endl;
			return 1;
			}
}


int main(int argc,char*argv[]){

	std::string tmp_ligand_name,tmp_config_name;	
	int manageNum;
	boost::program_options::options_description genm("generate manage core");
	genm.add_options()
		("ligand",boost::program_options::value<std::string>(&tmp_ligand_name), "part of the ligand (PDBQT)")
		("config", boost::program_options::value<std::string>(&tmp_config_name), "the above options can be put here")
		("manageNum", boost::program_options::value<int>(&manageNum)->default_value(1), "the number of manage core ,to read ligandlist and manage work process")
	;
        boost::program_options::variables_map vm;
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argc,argv).options(genm).run(),vm);
		notify(vm);
		//store(parse_command_line(argc, argv, manage_core, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
                //store(boost::program_options::command_line_parser(argc, argv).options(manage_core).style(boost::program_options::command_line_style::default_style ^ boost::program_options::command_line_style::allow_guessing).positional(positional).run(),vm);
	}catch(...){
		cerr << "There are undefined options in the parameters entered." << endl;
		return 1;
	}
	if(manageNum  < 1)
		throw usage_error("manage_num must be 1 or greater");

	int rank;
	int numProces;
        int manage_num = manageNum;

        int ligand_c_size;
        char *ligand_c;
	double startTime = 0.0, endTime = 0.0;

	MPI_Init (&argc, &argv );
   	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  //获得进程ID
   	MPI_Comm_size(MPI_COMM_WORLD, &numProces);  //获得进程个数
	
	//1. read ligand.pqbqt
	if(rank==0){

		startTime = MPI_Wtime();

		//create result directory
                string result_path ="./result";
                if(!boost::filesystem::exists(result_path)){
                        boost::filesystem::path dir(result_path);
                        if(boost::filesystem::create_directory(dir)){
                        }
                }

		double readLigand_startTime = MPI_Wtime();//read ligand start
		//1.1 read file
		ifstream fin("./ligand/ligand", std::ios::binary);
		if(!fin.is_open()){
                        std::cout<<"error: can not open ligand file."<<std::endl;
                }
		std::vector<char> buf(static_cast<unsigned int>(fin.seekg(0, std::ios::end).tellg()));        
                fin.seekg(0, std::ios::beg).read(&buf[0], static_cast<std::streamsize>(buf.size()));
                fin.close();		
		double readLigand_endTime = MPI_Wtime();//read ligand end
		
		double sendLigand_startTime = MPI_Wtime();//send ligand start
		char ligand_array[buf.size()];
                std::copy(buf.begin(), buf.end(), ligand_array);
                ligand_c = ligand_array;
                ligand_c_size = buf.size();

                MPI_Bcast(&ligand_c_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(ligand_c, buf.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
		double sendLigand_endTime = MPI_Wtime();//send ligand end
		printf("readLigand time:%.2lfs | sendLigand time:%.2lfs\n",readLigand_endTime - readLigand_startTime,sendLigand_endTime - sendLigand_startTime);
                                                                       
		// receive calculation finish signal
		MPI_Status mStatus;
		double countPlstime[3] = {0};
        	double tmp_cPlst[3];
		int front_process = -1, tmp_process;
		ofile f(make_path("process_status.txt"));
		for(int i=0;i<manage_num;i++){
			MPI_Recv(tmp_cPlst, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mStatus);
			//view nodes process
			tmp_process = (i+1)*100/manage_num;
			if(tmp_process!=front_process){
				printf("finished process:%d%\n",tmp_process);
				front_process = tmp_process;
			}
        		countPlstime[0] = countPlstime[0] + tmp_cPlst[0];
        		if(countPlstime[1]<tmp_cPlst[1]){
                        	countPlstime[1] = tmp_cPlst[1];
                	}
			countPlstime[2] = countPlstime[2] + tmp_cPlst[2];
                        f<<mStatus.MPI_SOURCE;
			f<<"\n"; 
		}	
		f.close();
		printf("done\n");

		endTime = MPI_Wtime(); // end timer.
		fflush(stdout);
	//sleep(000);
        printf("\n\n..........................................\n");
        //printf("   Number of workers       = %d \n", numProces - 1);
        printf("   总共对接次数（Number of docking）      = %u \n", (int)countPlstime[0]);
        printf("   总共运行时间（Total time required）    = %.2lf seconds.\n", endTime - startTime);
        printf("..........................................\n\n");
	printf("   总共计算核数（Total number caculate cores) = %u.\n", (int)(numProces - manage_num - 1));
        printf("   总的核心计算时间（Total caculate time) = %.2lf seconds.\n",  countPlstime[2]);
	printf("   平均单核计算时间（Average caculate time) = %.2lf seconds.\n",  countPlstime[2]/(numProces - manage_num - 1));
        printf("   单主核最长计算时间（Max time core）    = %.2lf seconds.\n",  countPlstime[1]);
        printf("..........................................\n\n");
		printf("...........FINISHED BY VINA@QNLM..........\n");
        fflush(stdout);

	
	}
	else if(rank<manage_num+1){
	//	std::cout<<"Hello world from process "<<rank<<" of "<<numProces<<std::endl;
		//1. receive ligand
		MPI_Bcast(&ligand_c_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		char ligand_c_m[ligand_c_size];
		ligand_c = ligand_c_m; 
                MPI_Bcast(ligand_c, ligand_c_size, MPI_CHAR, 0, MPI_COMM_WORLD);	
		//std::cout<<"Send ligand finished..."<<std::endl;

		//2.1 read receptor file
		double readReceptor_startTime = MPI_Wtime();//read receptor start
		MPI_Status mStatus;
		vector<std::string> receptors;
		std::string ReceptorName = "./protein/proteinlist" + std::to_string(rank);

		ifstream receptor_handler(ReceptorName, std::ios::binary);
                if(!receptor_handler.is_open()){
                        std::cout<<"error: can not open receptor file."<<std::endl;
                }
                std::vector<char> buf(static_cast<unsigned int>(receptor_handler.seekg(0, std::ios::end).tellg()));
                receptor_handler.seekg(0, std::ios::beg).read(&buf[0], static_cast<std::streamsize>(buf.size()));
                receptor_handler.close();
		double readReceptor_endTime = MPI_Wtime();//read receptor end
		//printf("Rank %d read ligandfile time = %.2lf seconds.\n", rank,readLigand_endTime - readLigand_startTime);

		//2.2 char to string
		double sendReceptor_startTime = MPI_Wtime();//send receptor start
		char receptor_array[buf.size()];
                std::copy(buf.begin(), buf.end(), receptor_array);
                std::string receptor_str(receptor_array, buf.size());
		boost::split(receptors, receptor_str, boost::is_any_of("$$$$"), boost::token_compress_on);
		receptors.pop_back();

		int receptor_num = receptors.size();	
		int total_receptor = manage_num * receptor_num;    // eg. 100*180000
		//bool falg;	
		int send_num;
		int control_num;  //eg. should send to all worker ,receptor number
		string tmp="";
		control_num=(numProces-1)/manage_num -1;
		
		//2.3 send to worker
		//this, we assum numProcess must be manage_num*n + 1, n is 1,2,3...
		if(receptor_num % control_num == 0){
            		send_num = receptor_num / control_num;
			int loop=0;
			int control_id=0;
			for(int i=0; i<receptors.size(); i++){
				loop++;
				//if((443<i)&&(i<445)){std::cout<<receptors[i]<<std::endl;}
				tmp = tmp+receptors[i]+"$$$$";
				if(loop==send_num){
					int size=tmp.size();
					int offset = manage_num + 1 + (rank-1)*control_num + control_id;
					MPI_Send(&size, 1, MPI_INT, offset, 0, MPI_COMM_WORLD);
					MPI_Send(tmp.c_str(), size, MPI_CHAR, offset, 0, MPI_COMM_WORLD);
					control_id++;
					tmp = "";
					loop=0;
				}
			}
        }
		else{
            send_num = receptor_num / control_num +1;
			int loop=0;
                        int control_id=0;
                        for(int i=0; i<receptors.size(); i++){
                                loop++;
				//if(2702<i&&i<2711){std::cout<<"receptor-"<<i<<receptors[i]<<std::endl;}
                                tmp = tmp+receptors[i]+"$$$$";
                                if(loop==send_num){
                                        int size=tmp.size();
                                        int offset = manage_num + 1 + (rank-1)*control_num + control_id;
                                        MPI_Send(&size, 1, MPI_INT, offset, 0, MPI_COMM_WORLD);
                                        MPI_Send(tmp.c_str(), size, MPI_CHAR, offset, 0, MPI_COMM_WORLD);
                                        control_id++;
					tmp = "";
					loop=0;
                                }
                        }
			//tmp have remaining receptors	
			for(int i=control_id; i<control_num; i++){
				if(tmp == ""){
					int nu=0; 	
					int offset = manage_num + 1 + (rank-1)*control_num + i;
					MPI_Send(&nu, 1, MPI_INT, offset , 0, MPI_COMM_WORLD);
				}
				else{
					int tmp_size =tmp.size(); 
					int offset = manage_num + 1 + (rank-1)*control_num + i;
					MPI_Send(&tmp_size, 1, MPI_INT, offset, 0, MPI_COMM_WORLD);
                                        MPI_Send(tmp.c_str(), tmp_size, MPI_CHAR, offset, 0, MPI_COMM_WORLD);
					tmp = "";
				}
			}
        }
	double sendReceptor_endTime = MPI_Wtime();//send receptor end
	printf("Rank %d readReceptor time:%.2lfs | sendReceptor time:%.2lfs\n", rank,readReceptor_endTime - readReceptor_startTime,sendReceptor_endTime - sendReceptor_startTime);
	printf("working , please waiting...\n");

        //receive count and time
        MPI_Status mStatus_1;
        double countPlstime[3] = {0};
        double tmp_cPlst[2];
	
        //std::cout<<"control_num:"<<control_num<<std::endl;
        for(int i=0;i<control_num;i++){
        	//MPI_Recv(&result_count, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mStatus_1);            
        	MPI_Recv(tmp_cPlst, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mStatus_1);
        	countPlstime[0] = countPlstime[0] + tmp_cPlst[0];
		if(countPlstime[1]<tmp_cPlst[1]){
			countPlstime[1] = tmp_cPlst[1];
		}
        	countPlstime[2] = countPlstime[2] + tmp_cPlst[1];      
        }	
		
		//send to 0
		MPI_Send(&countPlstime, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		
	}
	else{
		int control_num,offset;
		//1. receive ligand
		MPI_Bcast(&ligand_c_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
                char ligand_c_m[ligand_c_size];
                ligand_c = ligand_c_m;
                MPI_Bcast(ligand_c, ligand_c_size, MPI_CHAR, 0, MPI_COMM_WORLD);
		string single_ligands(ligand_c_m, ligand_c_size);
		std::vector<std::string> ligands;
		
                boost::split(ligands, single_ligands, boost::is_any_of("$$$$"), boost::token_compress_on);
		
		ligands.pop_back();
                //if(rank==2) std::cout<<"ligand:"<<ligands[0]<<std::endl;
		
		//2. receive receptor
		MPI_Status mStatus;	
		int all_char_num;	
		MPI_Recv(&all_char_num, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mStatus);
		if(all_char_num==0){
			//std::cout<<"Rank "<<rank<<" is sleeping..."<<std::endl;
                        //control_num=(numProces-1)/manage_num -1;
                        //offset = (rank-manage_num-1)/control_num+1;
                        //MPI_Send(&all_char_num, 1, MPI_INT, offset, 0, MPI_COMM_WORLD);
			//MPI_Send(&all_char_num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			double countPlstime[2]={0};
                        MPI_Send(countPlstime, 2, MPI_DOUBLE, mStatus.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}
		else{
		//	std::cout<<"worker:"<<all_char_num<<std::endl;
			char receptor_c[all_char_num];
			MPI_Recv(receptor_c, all_char_num, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mStatus);
			//std::cout<<"Send receptor finished..."<<std::endl;
			string single_receptors(receptor_c, all_char_num);
		
			std::vector<std::string> receptors;
			
			
        	boost::split(receptors, single_receptors, boost::is_any_of("$$$$"), boost::token_compress_on);
			receptors.pop_back();
			//std::cout<<"rank:"<<rank<<"receptor size():"<<receptors.size()<<"end:"<<std::endl;
			//create directory,and outpath
                        string dir_path ="./result/" + std::to_string(rank/2048);
                        if(!boost::filesystem::exists(dir_path)){
                                boost::filesystem::path dir(dir_path);
                                if(boost::filesystem::create_directory(dir)){
                                }
                        }
	
			//working
			double startTime = MPI_Wtime(); 
			int result_count = 0;
			for(int i=0; i<receptors.size(); i++){
				for(int j=0; j<ligands.size(); j++){
					//calculation
					main_ori(argc, argv, receptors[i], ligands[j], rank,i*ligands.size()+j+1);
					result = result + "####" + '\n';
					result_count++;
					//if(rank == 1000){printf("正在计算中，请等待。。。\n");}
					//std::cout<<result<<std::endl;
				}	
			}
			double endTime = MPI_Wtime();

			string out_path = dir_path + "/" + std::to_string(rank%2048) + ".data";
			try{
			ofile f(make_path(out_path));
			f<<result;
			f.close();                       
			}catch(file_error& e) {
                        std::cerr << "\n\nError: could not open \"" << out_path << "\" for " << (e.in ? "reading" : "writing") << ".\n";
                        std::cout<<"file__error"<<std::endl;
			return 1;
                        };
			//std::cout<<"Rank "<<rank<<" finished write "<<result_count<<" times docking result!"<<std::endl;
			
			//send to rank0 work finish signal 
			double countPlstime[2];
			countPlstime[0] = (double)result_count;
			countPlstime[1] = endTime-startTime;
			//MPI_Send(&result_count, 1, MPI_INT, offset, 0, MPI_COMM_WORLD);
			//MPI_Send(countPlstime, 2, MPI_DOUBLE, offset, 0, MPI_COMM_WORLD);
			MPI_Send(countPlstime, 2, MPI_DOUBLE, mStatus.MPI_SOURCE, 0, MPI_COMM_WORLD);
			
		}
	}

	MPI_Finalize();
	return 0;
}
