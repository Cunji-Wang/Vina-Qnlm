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

#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <boost/utility.hpp> // for noncopyable 
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include "parse_pdbqt.h"
#include "atom_constants.h"
#include "file.h"
#include "convert_substring.h"
#include "parse_error.h"
#ifdef ATHREAD_MC
#include "tranModel.h"
#endif

sz forbranch;
sz total;
//double receptorAtom[400]; 

struct stream_parse_error {
	unsigned line;
	std::string reason;
	stream_parse_error(unsigned line_, const std::string& reason_) : line(line_), reason(reason_) {}
	parse_error to_parse_error(const path& name) const {
		return parse_error(name, line, reason);
	}
};

struct parsed_atom : public atom {
	unsigned number; 
	parsed_atom(sz ad_, fl charge_, const vec& coords_, unsigned number_) : number(number_) {
		ad = ad_;
		charge = charge_;
		coords = coords_;
	}
};

void add_context(context& c, std::string& str) {
	c.push_back(parsed_line(str, boost::optional<sz>()));
}

std::string omit_whitespace(const std::string& str, sz i, sz j) {
	if(i < 1) i = 1;
	if(j < i-1) j = i-1; // i >= 1
	if(j < str.size()) j = str.size();

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	// omit trailing whitespace
	while(i <= j && std::isspace(str[j-1]))
		--j;

	VINA_CHECK(i-1 < str.size());
	VINA_CHECK(j-i+1 < str.size());

	return str.substr(i-1, j-i+1);
}

struct atom_syntax_error {
	std::string nature;
	atom_syntax_error(const std::string& nature_) : nature(nature_) {}
};

template<typename T>
T checked_convert_substring(const std::string& str, sz i, sz j, const std::string& dest_nature) {
	VINA_CHECK(i >= 1);
	VINA_CHECK(i <= j+1);
	if(j > str.size()) throw atom_syntax_error("The line is too short");

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	const std::string substr = str.substr(i-1, j-i+1);
	try {
		return boost::lexical_cast<T>(substr);
	}
	catch(...) {
		throw atom_syntax_error(std::string("\"") + substr + "\" is not a valid " + dest_nature);
	}
}

parsed_atom parse_pdbqt_atom_string(const std::string& str) {
	unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "atom number");
	vec coords(checked_convert_substring<fl>(str, 31, 38, "coordinate"),
			   checked_convert_substring<fl>(str, 39, 46, "coordinate"),
			   checked_convert_substring<fl>(str, 47, 54, "coordinate"));
	fl charge = 0;
	if(!substring_is_blank(str, 69, 76))
		charge = checked_convert_substring<fl>(str, 69, 76, "charge");
	std::string name = omit_whitespace(str, 78, 79);
	sz ad = string_to_ad_type(name);
	parsed_atom tmp(ad, charge, coords, number);
	if(is_non_ad_metal_name(name))
		tmp.xs = XS_TYPE_Met_D;
	if(tmp.acceptable_type()) 
		return tmp;
	else 
		throw atom_syntax_error(std::string("\"") + name + "\" is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.");
}

struct atom_reference {
	sz index;
	bool inflex;
	atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

struct movable_atom : public atom {
	vec relative_coords;
	movable_atom(const atom& a, const vec& relative_coords_) : atom(a) {
		relative_coords = relative_coords_;
	}
};

struct rigid {
	atomv atoms;
};

typedef std::vector<movable_atom> mav;

struct non_rigid_parsed {
	vector_mutable<ligand> ligands;
	vector_mutable<residue> flex;

	mav atoms;
	atomv inflex;

	distance_type_matrix atoms_atoms_bonds;
	matrix<distance_type> atoms_inflex_bonds;
	distance_type_matrix inflex_inflex_bonds;

	distance_type_matrix mobility_matrix() const {
		distance_type_matrix tmp(atoms_atoms_bonds);
		tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
		return tmp;
	}
};

struct parsing_struct {
	// start reading after this class
	template<typename T> // T == parsing_struct
	struct node_t {
		sz context_index;
		parsed_atom a;
		std::vector<T> ps;
		node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

		// inflex atom insertion
		void insert_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
			nr.inflex.push_back(a);
		}
		void insert_immobiles_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile_inflex(nr);
		}

		// insertion into non_rigid_parsed
		void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
			vec relative_coords; relative_coords = a.coords - frame_origin;
			c[context_index].second = nr.atoms.size();
			nr.atoms.push_back(movable_atom(a, relative_coords));
		}
		void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile(nr, c, frame_origin);
		}
	};

	typedef node_t<parsing_struct> node;
	boost::optional<sz> immobile_atom; // which of `atoms' is immobile, if any
	boost::optional<atom_reference> axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
	boost::optional<atom_reference> axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
	std::vector<node> atoms;

	void add(const parsed_atom& a, const context& c) { 
		VINA_CHECK(c.size() > 0);
		atoms.push_back(node(a, c.size()-1)); 
	}
	const vec& immobile_atom_coords() const {
		VINA_CHECK(immobile_atom);
		VINA_CHECK(immobile_atom.get() < atoms.size());
		return atoms[immobile_atom.get()].a.coords;
	}
	// inflex insertion
	void insert_immobile_inflex(non_rigid_parsed& nr) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.inflex.size(), true);
			atoms[immobile_atom.get()].insert_inflex(nr);
		}
	}

	// insertion into non_rigid_parsed
	void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.atoms.size(), false);
			atoms[immobile_atom.get()].insert(nr, c, frame_origin);
		}
	}

	bool essentially_empty() const { // no sub-branches besides immobile atom, including sub-sub-branches, etc
		VINA_FOR_IN(i, atoms) {
			if(immobile_atom && immobile_atom.get() != i)
				return false;
			const node& nd = atoms[i];
			if(!nd.ps.empty())
				return false; // FIXME : iffy
		}
		return true;
	}
};

unsigned parse_one_unsigned(const std::string& str, const std::string& start, unsigned count) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp;
	in_str >> tmp;
	if(!in_str || tmp < 0) 
		throw stream_parse_error(count, "Syntax error");
	return unsigned(tmp);
}

void parse_two_unsigneds(const std::string& str, const std::string& start, unsigned count, unsigned& first, unsigned& second) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp1, tmp2;
	in_str >> tmp1;
	in_str >> tmp2;
	if(!in_str || tmp1 < 0 || tmp2 < 0) 
		throw stream_parse_error(count, "Syntax error");
	first = unsigned(tmp1);
	second = unsigned(tmp2);
}

void parse_pdbqt_rigid(const path& name, rigid& r) {
	ifile in(name);
	printf("hello2\n");
	unsigned count = 0;
	std::string str;
	while(std::getline(in, str)) {
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "TER")) {} // ignore 
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				r.atoms.push_back(parse_pdbqt_atom_string(str));
			}
			catch(atom_syntax_error& e) {
				throw parse_error(name, count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw parse_error(name, count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_rigid_content(const path& name, std::vector<std::string>& receptor_content, rigid& r) {	
	std::string str;
	//while(std::getline(in, str)) {
	for(int i=0;i<receptor_content.size();i++) {
		
		std::string str=receptor_content[i];
		
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "TER")) {} // ignore 
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				r.atoms.push_back(parse_pdbqt_atom_string(str));
			}
			catch(atom_syntax_error& e) {
				throw parse_error(name, i, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw parse_error(name, i, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(i, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, i, "Unknown or inappropriate tag");
	}
	//printf("r.atom.size::%d\n",r.atoms.size());
}

void parse_pdbqt_root_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				p.add(parse_pdbqt_atom_string(str), c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "ENDROOT")) return;
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_root_aux_content(std::vector<std::string>& ligand_content, unsigned& count, parsing_struct& p, context& c) {

        for(;count<ligand_content.size();) {
                std::string str=ligand_content[count];
                //std::cout<<str<<std::endl;
                add_context(c, str);
                ++count;
                //std::cout<<count<<std::endl;
                if(str.empty()) {} // ignore ""
                else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
                else if(starts_with(str, "REMARK")) {} // ignore
                else if(starts_with(str, "ATOM     ") || starts_with(str, "HETATM")) {
                        try {
                                p.add(parse_pdbqt_atom_string(str), c);
                        }
                        catch(atom_syntax_error& e) {
                                throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
                        }
                        catch(...) {
                                throw stream_parse_error(count, "ATOM syntax incorrect");
                        }
                }
                else if(starts_with(str, "ENDROOT")) return;
                else if(starts_with(str, "MODEL"))
                        throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
                else throw stream_parse_error(count, "Unknown or inappropriate tag");
        }
}

void parse_pdbqt_root(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ROOT")) {
			parse_pdbqt_root_aux(in, count, p, c);
			break;
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_root_content(std::vector<std::string>& ligand_content, unsigned& count, parsing_struct& p, context& c) {

        for(;count<ligand_content.size();) {
                std::string str=ligand_content[count];
                //std::cout<<str<<std::endl;
                add_context(c, str);
                ++count;
                //std::cout<<count<<std::endl;
                if(str.empty()) {} // ignore
                else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
                else if(starts_with(str, "REMARK")) {} // ignore
                else if(starts_with(str, "ROOT")) {
                        parse_pdbqt_root_aux_content(ligand_content, count, p, c);
                        break;
                }
                else if(starts_with(str, "MODEL"))
                        throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
                else throw stream_parse_error(count, "Unknown or inappropriate tag");
        }
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to); // forward declaration

void parse_pdbqt_branch_aux(std::istream& in, unsigned& count, const std::string& str, parsing_struct& p, context& c) {
	unsigned first, second;
	parse_two_unsigneds(str, "BRANCH", count, first, second); 
	sz i = 0;
	for(; i < p.atoms.size(); ++i)
		if(p.atoms[i].a.number == first) {
			p.atoms[i].ps.push_back(parsing_struct());
			parse_pdbqt_branch(in, count, p.atoms[i].ps.back(), c, first, second);
			break;
		}
	if(i == p.atoms.size())
		throw stream_parse_error(count, "No atom number " + boost::lexical_cast<std::string>(first) + " in this branch");
}

void parse_pdbqt_branch_content(std::vector<std::string>& ligand_content, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to); // forward declaration

void parse_pdbqt_branch_aux_content(std::vector<std::string>& ligand_content, unsigned& count, const std::string& str, parsing_struct& p, context& c) {
        unsigned first, second;
        parse_two_unsigneds(str, "BRANCH", count, first, second);
        sz i = 0;
        for(; i < p.atoms.size(); ++i)
                if(p.atoms[i].a.number == first) {
                        p.atoms[i].ps.push_back(parsing_struct());
                        parse_pdbqt_branch_content(ligand_content, count, p.atoms[i].ps.back(), c, first, second);
                        break;
                }
        if(i == p.atoms.size())
                throw stream_parse_error(count, "No atom number " + boost::lexical_cast<std::string>(first) + " in this branch");
}

void parse_pdbqt_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue) {
	parse_pdbqt_root(in, count, p, c);

	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(!residue && starts_with(str, "TORSDOF")) {
			if(torsdof) throw stream_parse_error(count, "TORSDOF can occur only once");
			torsdof = parse_one_unsigned(str, "TORSDOF", count);
		}
		else if(residue && starts_with(str, "END_RES")) return; 
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_aux_content(std::vector<std::string>& ligand_content, unsigned& count, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue) {
        parse_pdbqt_root_content(ligand_content, count, p, c);
        for(;count<ligand_content.size();) {
                std::string str=ligand_content[count];
                //std::cout<<str<<std::endl;
                add_context(c, str);
                ++count;
                //std::cout<<count<<std::endl;
                if(str.empty()) {} // ignore ""
                else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
                else if(starts_with(str, "REMARK")) {} // ignore
                else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux_content(ligand_content, count, str, p, c);
                else if(!residue && starts_with(str, "TORSDOF")) {
                        if(torsdof) throw stream_parse_error(count, "TORSDOF can occur only once");
                        torsdof = parse_one_unsigned(str, "TORSDOF", count);
                }
                else if(residue && starts_with(str, "END_RES")) return;
                else if(starts_with(str, "MODEL"))
                        throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
                else throw stream_parse_error(count, "Unknown or inappropriate tag");
        }
}

void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
	if(atm)
		VINA_RANGE(i, r.begin, r.end) {
			atom_reference& ar = atm.get();
			if(ar.inflex) 
				nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED; //(max_unsigned); // first index - atoms, second index - inflex
			else
				nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED; // (max_unsigned);
		}
}

void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
	if(axis_begin && axis_end) {
		atom_reference& r1 = axis_begin.get();
		atom_reference& r2 = axis_end  .get();
		if(r2.inflex) {
			VINA_CHECK(r1.inflex); // no atom-inflex rotors
			nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
		}
		else
			if(r1.inflex)
				nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (atoms, inflex)
			else
				nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
	}
}

typedef std::pair<sz, sz> axis_numbers;
typedef boost::optional<axis_numbers> axis_numbers_option;

void nr_update_matrixes(non_rigid_parsed& nr) {
	// atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)

	nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);  
	nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // first index - inflex, second index - atoms
	nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
}

template<typename B> // B == branch or main_branch or flexible_body 
void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
	b.node.begin = nr.atoms.size();
	VINA_FOR_IN(i, p.atoms) {  // postprocess atoms into 'b.node'
		parsing_struct::node& p_node = p.atoms[i];
		if(p.immobile_atom && i == p.immobile_atom.get()) {} // skip immobile_atom - it's already inserted in "THERE"
		else p_node.insert(nr, c, b.node.get_origin());
		p_node.insert_immobiles(nr, c, b.node.get_origin());
	}
	b.node.end = nr.atoms.size();

	nr_update_matrixes(nr);
	add_bonds(nr, p.axis_begin, b.node); // b.node is used as atom_range
	add_bonds(nr, p.axis_end  , b.node); // b.node is used as atom_range
	set_rotor(nr, p.axis_begin, p.axis_end);

	VINA_RANGE(i, b.node.begin, b.node.end)
		VINA_RANGE(j, i+1, b.node.end)
			nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME


	VINA_FOR_IN(i, p.atoms) { 	// postprocess children
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile already inserted // FIXME ?!
				b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch will assign begin and end
				postprocess_branch(nr, ps, c, b.children.back());
			}
		}
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
	VINA_CHECK(!p.atoms.empty());
	nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch will assign begin and end
	postprocess_branch(nr, p, c, nr.ligands.back());
	nr_update_matrixes(nr); // FIXME ?
}

void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		p_node.insert_inflex(nr);
		p_node.insert_immobiles_inflex(nr);
	}
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile atom already inserted // FIXME ?!
				nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords))); // postprocess_branch will assign begin and end
				postprocess_branch(nr, ps, c, nr.flex.back());
			}
		}
	}
	nr_update_matrixes(nr); // FIXME ?
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	printf("hello1\n");
	unsigned count = 0;
	parsing_struct p;
	boost::optional<unsigned> torsdof;
	try {
		parse_pdbqt_aux(in, count, p, c, torsdof, false);
		if(p.atoms.empty()) 
			throw parse_error(name, count, "No atoms in the ligand");
		if(!torsdof)
			throw parse_error(name, count, "Missing TORSDOF");
		postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
	}
	catch(stream_parse_error& e) {
		throw e.to_parse_error(name);
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_ligand_content(const path& name, non_rigid_parsed& nr, context& c,std::vector<std::string>& ligand_content) {
        //ifile in(name);
        unsigned count = 0;
        parsing_struct p;
        boost::optional<unsigned> torsdof;
        try {
                parse_pdbqt_aux_content(ligand_content, count, p, c, torsdof, false);
                if(p.atoms.empty())
                        throw parse_error(name, count, "No atoms in the ligand");
                if(!torsdof)
                        throw parse_error(name, count, "Missing TORSDOF");
                postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
        }
        catch(stream_parse_error& e) {
                throw e.to_parse_error(name);
        }
        VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_residue(std::istream& in, unsigned& count, parsing_struct& p, context& c) { 
	boost::optional<unsigned> dummy;
	parse_pdbqt_aux(in, count, p, c, dummy, true);
}

void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	unsigned count = 0;
	std::string str;

	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BEGIN_RES")) {
			try {
				parsing_struct p;
				parse_pdbqt_residue(in, count, p, c);
				postprocess_residue(nr, p, c);
			}
			catch(stream_parse_error& e) {
				throw e.to_parse_error(name);
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} //ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(starts_with(str, "ENDBRANCH")) {
			unsigned first, second;
			parse_two_unsigneds(str, "ENDBRANCH", count, first, second);
			if(first != from || second != to) 
				throw stream_parse_error(count, "Inconsistent branch numbers");
			if(!p.immobile_atom) 
				throw stream_parse_error(count, "Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch");
			return;
		}
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				parsed_atom a = parse_pdbqt_atom_string(str);
				if(a.number == to)
					p.immobile_atom = p.atoms.size();
				p.add(a, c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_branch_content(std::vector<std::string>& ligand_content, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to) {

        for(;count<ligand_content.size();) {
                std::string str=ligand_content[count];
                //std::cout<<str<<std::endl;
                add_context(c, str);
                ++count;
                if(str.empty()) {} //ignore ""
                else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
                else if(starts_with(str, "REMARK")) {} // ignore
                else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux_content(ligand_content, count, str, p, c);
                else if(starts_with(str, "ENDBRANCH")) {
                        unsigned first, second;
                        parse_two_unsigneds(str, "ENDBRANCH", count, first, second);
                        if(first != from || second != to)
                                throw stream_parse_error(count, "Inconsistent branch numbers");
                        if(!p.immobile_atom)
                                throw stream_parse_error(count, "Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch");                                return;                 }
                else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
                        try {
                                parsed_atom a = parse_pdbqt_atom_string(str);
                                if(a.number == to)
                                        p.immobile_atom = p.atoms.size();
                                p.add(a, c);
                        }
                        catch(atom_syntax_error& e) {
                                throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
                        }
                        catch(...) {
                                throw stream_parse_error(count, "ATOM syntax incorrect");
                        }
                }
                else if(starts_with(str, "MODEL"))
                        throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
                else throw stream_parse_error(count, "Unknown or inappropriate tag");
        }
}



//////////// new stuff //////////////////


struct pdbqt_initializer {
	model m;
	void initialize_from_rigid(const rigid& r) { // static really
		VINA_CHECK(m.grid_atoms.empty());
		m.grid_atoms = r.atoms;
	}
	void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
		VINA_CHECK(m.ligands.empty());
		VINA_CHECK(m.flex   .empty());

		m.ligands = nrp.ligands;
		m.flex    = nrp.flex;

		VINA_CHECK(m.atoms.empty());

		sz n = nrp.atoms.size() + nrp.inflex.size();
		m.atoms.reserve(n);
		m.coords.reserve(n);

		VINA_FOR_IN(i, nrp.atoms) {
			const movable_atom& a = nrp.atoms[i];
			atom b = static_cast<atom>(a);
			b.coords = a.relative_coords;
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_FOR_IN(i, nrp.inflex) {
			const atom& a = nrp.inflex[i];
			atom b = a;
			b.coords = zero_vec; // to avoid any confusion; presumably these will never be looked at
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_CHECK(m.coords.size() == n);

		m.internal_coords.resize(m.coords.size(), zero_vec); // FIXME

		m.minus_forces = m.coords;
		m.m_num_movable_atoms = nrp.atoms.size();

		if(is_ligand) {
			VINA_CHECK(m.ligands.size() == 1);
			m.ligands.front().cont = c;
		}
		else
			m.flex_context = c;

	}
	void initialize(const distance_type_matrix& mobility) {
		m.initialize(mobility);
	}
};

model parse_ligand_pdbqt  (const path& name) { // can throw parse_error
	non_rigid_parsed nrp;
	context c;
	parse_pdbqt_ligand(name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_nrp(nrp, c, true);
	tmp.initialize(nrp.mobility_matrix());
	return tmp.m;
}
model parse_ligand_pdbqt_content(const path& name,std::vector<std::string>& ligand_content) { // can throw parse_error
        non_rigid_parsed nrp;
        context c;
        //std::cout<<"1-------------------------"<<std::endl;
        parse_pdbqt_ligand_content(name, nrp, c,ligand_content);

        pdbqt_initializer tmp;
        tmp.initialize_from_nrp(nrp, c, true);
        tmp.initialize(nrp.mobility_matrix());
        return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name, const path& flex_name) { // can throw parse_error
	rigid r;
	non_rigid_parsed nrp;
	context c;
	parse_pdbqt_rigid(rigid_name, r);
	parse_pdbqt_flex(flex_name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	tmp.initialize_from_nrp(nrp, c, false);
	tmp.initialize(nrp.mobility_matrix());
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name) { // can throw parse_error
	rigid r;
	parse_pdbqt_rigid(rigid_name, r);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	distance_type_matrix mobility_matrix;
	tmp.initialize(mobility_matrix);
	return tmp.m;
}

model parse_receptor_pdbqt_content(const path& name,std::vector<std::string>& receptor_content) {
	rigid r;
	parse_pdbqt_rigid_content(name,receptor_content,r);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	distance_type_matrix mobility_matrix;
	tmp.initialize(mobility_matrix);
	return tmp.m;
}


#ifdef ATHREAD_MC
void copyatoms(atomv& modelatoms, double *saveAtoms){
    VINA_FOR_IN(i, modelatoms){
        atom tmp = modelatoms[i];
        int j = 0;
        saveAtoms[i*37] = (double)tmp.el;
        saveAtoms[i*37+1] = (double)tmp.ad;
        saveAtoms[i*37+2] = (double)tmp.xs;
        saveAtoms[i*37+3] = (double)tmp.sy;
        saveAtoms[i*37+4] = tmp.charge;
        saveAtoms[i*37+5] = tmp.coords[0];
        saveAtoms[i*37+6] = tmp.coords[1];
        saveAtoms[i*37+7] = tmp.coords[2];
        sz bondsize = tmp.bonds.size();
		if(bondsize>4){
			std::cout<<"copyatoms fun: bondsize is: "<<bondsize<<std::endl;
			//assert(bondsize<=4);
		}
        saveAtoms[i*37+8] = bondsize;
        while(bondsize!=0){
            saveAtoms[i*37+8+j*4+1] = (double)tmp.bonds[j].connected_atom_index.i;
            saveAtoms[i*37+8+j*4+2] = (double)tmp.bonds[j].connected_atom_index.in_grid;
            saveAtoms[i*37+8+j*4+3] = (double)tmp.bonds[j].rotatable;
            saveAtoms[i*37+8+j*4+4] = tmp.bonds[j].length;
           // printf("receptorAtom in_grid: %ld\n", (std::size_t)m->receptorAtom[i*25+8+j*4+2]);
           // printf("original in_grid: %ld\n", (std::size_t)tmp.bonds[j].connected_atom_index.in_grid);
            bondsize--;
            j++;
        }
    }
}
/*
void copyatoms(atomv& modelatoms, double *saveAtoms){
    VINA_FOR_IN(i, modelatoms){
        atom tmp = modelatoms[i];
        int j = 0;
        saveAtoms[i*25] = (double)tmp.el;
        saveAtoms[i*25+1] = (double)tmp.ad;
        saveAtoms[i*25+2] = (double)tmp.xs;
        saveAtoms[i*25+3] = (double)tmp.sy;
        saveAtoms[i*25+4] = tmp.charge;
        saveAtoms[i*25+5] = tmp.coords[0];
        saveAtoms[i*25+6] = tmp.coords[1];
        saveAtoms[i*25+7] = tmp.coords[2];
        sz bondsize = tmp.bonds.size();
                if(bondsize>4){
                        std::cout<<"copyatoms fun: bondszie is: "<<bondsize<<std::endl;
                        assert(bondsize<=4);
                }
        saveAtoms[i*25+8] = bondsize;
        while(bondsize!=0){
            saveAtoms[i*25+8+j*4+1] = (double)tmp.bonds[j].connected_atom_index.i;
            saveAtoms[i*25+8+j*4+2] = (double)tmp.bonds[j].connected_atom_index.in_grid;
            saveAtoms[i*25+8+j*4+3] = (double)tmp.bonds[j].rotatable;
            saveAtoms[i*25+8+j*4+4] = tmp.bonds[j].length;
           // printf("receptorAtom in_grid: %ld\n", (std::size_t)m->receptorAtom[i*25+8+j*4+2]);
           //            // printf("original in_grid: %ld\n", (std::size_t)tmp.bonds[j].connected_atom_index.in_grid);
            bondsize--;
            j++;
        }
    }
}
*/
          
//entry is equal to segchild
//1: 父节点，2：一个segment指针，3：root; 4：大小; 5：层数；6：父节点的坐标
sz iter_save_children(branch& entry, struct segmentbranch *segchild, struct flexibleHeritTree *heritTree, sz size, int level, int parent) {
    sz offset = size;
    sz nextoffset = offset+1;
	//sz nextoffset = offset;
    static sz finaloffset = 0;
    // printf("\niter_save_children current offset %d, children size: %ld\n",offset, entry.children.size());
    VINA_FOR_IN(k, entry.children){
		total=total+1;
		forbranch = forbranch + 1;
        branch child = entry.children[k];
        sz childsize = entry.children.size();
        sz curoffset = offset+k;
        //nextoffset = nextoffset+k+1;
		nextoffset = nextoffset+k;
        // printf("current children k :%ld, children size: %ld\n", k, entry.children.size());
        vec relative_axis = child.node.get_relative_axis();
        vec relative_origin = child.node.get_relative_origin();
        vec axis = child.node.get_axis();
        vec orgin = child.node.get_origin();
        VINA_FOR_IN(j, relative_axis)
        heritTree->flexiblechildren[total].relative_axis[j] = relative_axis[j];
        VINA_FOR_IN(j, relative_origin)
        heritTree->flexiblechildren[total].relative_origin[j] = relative_origin[j];
        VINA_FOR_IN(j, axis)
        heritTree->flexiblechildren[total].axis[j] = axis[j];
        VINA_FOR_IN(j, orgin)
        heritTree->flexiblechildren[total].origin[j] = orgin[j];
        mat orientation_m = child.node.get_orientation_m();
        VINA_FOR(j, 9)
        heritTree->flexiblechildren[total].orientation_m[j] = orientation_m.data[j];

        qt mqt = child.node.orientation();
        heritTree->flexiblechildren[total].orientation_q[0] = mqt.R_component_1();
        heritTree->flexiblechildren[total].orientation_q[1] = mqt.R_component_2();
        heritTree->flexiblechildren[total].orientation_q[2] = mqt.R_component_3();
        heritTree->flexiblechildren[total].orientation_q[3] = mqt.R_component_4();
        heritTree->flexiblechildren[total].begin = child.node.begin;
        heritTree->flexiblechildren[total].end = child.node.end;
        heritTree->flexiblechildren[total].size = childsize;
        //heritTree->flexiblechildren[total].level = level+1;
		heritTree->flexiblechildren[total].childrensize = child.children.size();
        heritTree->flexiblechildren[total].parent = parent;  //说明父节点是哪一个
        //make the link
        //todo: if there are more than one children which have same parenets, the previous will be override
        //segchild->segchildren = &heritTree->flexiblechildren[offset];
        // do recursion
        sz curlevel = level+1;
        /*
        printf("current offset is %ld, it's parent: %d, next offset: %ld, level is %d\n", curoffset,parent,nextoffset, curlevel);
        printf("offset is %d, &heritTree->flexiblechildren addr: 0x%x\n",
               curoffset, &heritTree->flexiblechildren[curoffset]);
        */
        if(finaloffset < nextoffset)
            finaloffset = nextoffset;
        //iter_save_children(child, &heritTree->flexiblechildren[curoffset], heritTree, nextoffset, curlevel, curoffset);
		iter_save_children(child, &heritTree->flexiblechildren[curoffset], heritTree, nextoffset, curlevel, total);

    }
	finaloffset = offset;
    //printf("final offset: %ld, forbranch:%d\n", finaloffset, forbranch);
	return forbranch;
    //return finaloffset;
}


tranModel reform_transfer_model( model& m){

    /* one atom has element: sz el, ad, xs, sy, double charge, double coords[3]
     * sz bondsize, int i, bool in_grid, bool rotatable, double length
     * default, the max bond size is 4. so one atoms should have 9+4*4=25,
     * but the max bond size sometime is 6,so we modify max bond is 7,thus one atom 9+4*7=37
     */
    total = 0;
    forbranch = 0;
    tranModel *tranM;
    atomv mygrid_atoms = m.get_grid_atoms();
    sz recsize = mygrid_atoms.size();
    tranM = (struct tranModel*)malloc(sizeof(struct tranModel));
    //printf("recsize:%d\n",recsize);
    //tranM->receptorAtom = (double *)malloc(sizeof(double)*recsize*37);
    //copyatoms(mygrid_atoms, tranM->receptorAtom);
    // save the ligand atoms
    atomv ligand_atoms = m.get_ligand_atoms();
    sz ligsize = ligand_atoms.size();
    tranM->ligandAtom = (double *)malloc(sizeof(double)*ligsize*37);
    copyatoms(ligand_atoms, tranM->ligandAtom);


    // save the ligand coordinate by using
    vecv lig_coords = m.get_ligand_coords();
    sz coordsize = lig_coords.size();
    //printf("reform_transfer_model coordsize is %ld\n", coordsize);
    tranM->ligcoords = (double *)malloc(sizeof(double)*coordsize*3);
    tranM->internal_coords = (double *)malloc(sizeof(double)*coordsize*3);
    tranM->minus_forces = (double *)malloc(sizeof(double)*coordsize*3);
    VINA_FOR_IN(i, lig_coords){
        vec tmp = lig_coords[i];
        tranM->ligcoords[i*3] = tmp[0];
        tranM->minus_forces[i*3] = tmp[0];
        tranM->internal_coords[i*3] = 0;

        tranM->ligcoords[i*3+1] = tmp[1];
        tranM->minus_forces[i*3+1] = tmp[1];
        tranM->internal_coords[i*3+1] = 0;

        tranM->ligcoords[i*3+2] = tmp[2];
        tranM->minus_forces[i*3+2] = tmp[2];
        tranM->internal_coords[i*3+2] = 0;
    }

    sz pairsize = m.num_internal_pairs();
    //printf("reform_transfer_model pairsize is %ld\n", pairsize);
    tranM->ligandpair = (int **)malloc(sizeof(int*)*pairsize);

    interacting_pairs mpairs = m.get_ligand_pairs();
    VINA_FOR_IN(j, mpairs){
        tranM->ligandpair[j] = (int *)malloc(sizeof(int)*3);
        const interacting_pair ip = mpairs[j];
        tranM->ligandpair[j][0] = ip.type_pair_index;
        tranM->ligandpair[j][1] = ip.a;
        tranM->ligandpair[j][2] = ip.b;
    }

    //for hereotree reform
    conf_size cs = m.get_size();
    sz ligand_torsion_size;
    VINA_FOR_IN(i, cs.ligands)
		ligand_torsion_size = cs.ligands[i];     
	
        
    // get the branches  //Hertirtree
    branches mbranch = m.get_branches();
    sz branchessize = mbranch.size();

    tranM->heritTree = (struct flexibleHeritTree*)malloc(sizeof(struct flexibleHeritTree) +
                            sizeof(struct segmentbranch)*ligand_torsion_size);
/*   
 printf("reform_transfer_model ligand_torsion_size is %ld, branches size: %ld, segment size: %d\n",
            ligand_torsion_size, branchessize, sizeof(struct segmentbranch));
*/
    struct rigid_body rootnode = m.get_rootbody();
    vec origin = rootnode.get_origin();
    //first initialize rootnode
    // rigidbody have 18 element: vec origin_(3), mat orientation_m (9), qt  orientation_q (4), sz begin, sz end
    tranM->heritTree->origin[0] = origin[0];
    tranM->heritTree->origin[1] = origin[1];
    tranM->heritTree->origin[2] = origin[2];
    mat orientation_m = rootnode.get_orientation_m();
    VINA_FOR(i, 9)
        tranM->heritTree->orientation_m[i] = orientation_m.data[i];

    qt mqt = rootnode.orientation();
    tranM->heritTree->orientation_q[0] = mqt.R_component_1();
    tranM->heritTree->orientation_q[1] = mqt.R_component_2();
    tranM->heritTree->orientation_q[2] = mqt.R_component_3();
    tranM->heritTree->orientation_q[3] = mqt.R_component_4();
    tranM->heritTree->begin = rootnode.begin;
    tranM->heritTree->end = rootnode.end;
    tranM->heritTree->branchsize = branchessize;     //记录第一层的大小

    /* save the first level children
     * total elements are 27. vec relative_axis(3), vec relative_origin(3), vec axis（3),
     * vec origin_(3), mat orientation_m (9), qt  orientation_q (4), sz begin, sz end
     */
    int k = 0;
    sz offset = branchessize;
//	printf("mbranch szie：%d\n,offset %d\n",mbranch.size(), offset);
	forbranch = forbranch +offset;
    VINA_FOR_IN(k, mbranch){
        branch child = mbranch[k];
        vec relative_axis = child.node.get_relative_axis();
        vec relative_origin = child.node.get_relative_origin();
        vec axis = child.node.get_axis();
        vec orgin = child.node.get_origin();
        VINA_FOR_IN(j, relative_axis)
            tranM->heritTree->flexiblechildren[total].relative_axis[j] = relative_axis[j];
        VINA_FOR_IN(j, relative_origin)
            tranM->heritTree->flexiblechildren[total].relative_origin[j] = relative_origin[j];
        VINA_FOR_IN(j, axis)
            tranM->heritTree->flexiblechildren[total].axis[j] = axis[j];
        VINA_FOR_IN(j, orgin){
            tranM->heritTree->flexiblechildren[total].origin[j] = orgin[j];
	    
	}	
        mat orientation_m = child.node.get_orientation_m();
        VINA_FOR(j, 9)
            tranM->heritTree->flexiblechildren[total].orientation_m[j] = orientation_m.data[j];

        qt mqt = child.node.orientation();
        tranM->heritTree->flexiblechildren[total].orientation_q[0] = mqt.R_component_1();
        tranM->heritTree->flexiblechildren[total].orientation_q[1] = mqt.R_component_2();
        tranM->heritTree->flexiblechildren[total].orientation_q[2] = mqt.R_component_3();
        tranM->heritTree->flexiblechildren[total].orientation_q[3] = mqt.R_component_4();
        tranM->heritTree->flexiblechildren[total].begin = child.node.begin;
        tranM->heritTree->flexiblechildren[total].end = child.node.end;
        tranM->heritTree->flexiblechildren[total].size = branchessize;    //又记录了本层的大小
//		std::cout<<"这是第一层的，第"<<k+1<<"个branch"<< "   tatol:"<<total<<std::endl;
        //tranM->heritTree->flexiblechildren[total].level = 1;     //记录它是在第几层
		tranM->heritTree->flexiblechildren[total].childrensize = child.children.size();
        tranM->heritTree->flexiblechildren[total].parent = -1;
//        printf("offset is %d, tranM->heritTree->flexiblechildren addr: 0x%x\n\n", offset, &tranM->heritTree->flexiblechildren[k]);
        //save the branch and children info to structure, update offset of the children
        offset = iter_save_children(child, &tranM->heritTree->flexiblechildren[k], tranM->heritTree, forbranch, 1, total);
		total=total+1;
//		std::cout<<"forbranch:"<<offset<<std::endl;
		
    }

    tranM->para = (struct paras*)malloc(sizeof(struct paras));
    tranM->para->coordsize = coordsize;
    tranM->para->pairsize = pairsize;
    tranM->para->torsionsize = ligand_torsion_size;
    tranM->para->m_num_movable_atoms = m.num_movable_atoms();
    tranM->para->atom_typing_used = XS;
    //luhao add
    tranM->para->ligsize = ligsize;

    //luhao modify: globalTm->receptorAtom[i*25] 
/*
    for(int i=rootnode.begin; i<rootnode.end; i++){
	printf("i----:%d\n",i);
	receptorAtom[i] = tranM->receptorAtom[i*37];
    }
*/

/*
	std::cout<<"\n\nprint tree parent index:";
	for(int i = 0;i<forbranch; i++){
		std::cout<<tranM->heritTree->flexiblechildren[i].parent<<" ";
	}
*/

//	std::cout<<"get " <<forbranch<<" tree data"<<std::endl;
	//while(1);
/*	
	std::cout<<tranM->heritTree->flexiblechildren[0].size<<std::endl;
	std::cout<<tranM->heritTree->flexiblechildren[1].size<<std::endl;
	std::cout<<tranM->heritTree->flexiblechildren[1].segchildren->size<<std::endl;
	if(tranM->heritTree->flexiblechildren[100].size == NULL)
		std::cout<<"没有100\n"<<std::endl;
	else{
		std::cout<<"100的size："<<tranM->heritTree->flexiblechildren[100].size<<std::endl;
		std::cout<<"200的size："<<tranM->heritTree->flexiblechildren[200].size<<std::endl;
	}
	*/	


    return *tranM;
}

#endif
