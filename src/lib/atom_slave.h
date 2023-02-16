/*
 * atom_slave.h
 * Copyright (C) 2020 yumaoxue <yumaoxue@qnlm.ac>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef ATOM_SLAVE_H
#define ATOM_SLAVE_H
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include "triangular_matrix_slave.h"
#include "tranModel.h"
// based on SY_TYPE_* but includes H
const int EL_TYPE_H    =  0;
const int EL_TYPE_C    =  1;
const int EL_TYPE_N    =  2;
const int EL_TYPE_O    =  3;
const int EL_TYPE_S    =  4;
const int EL_TYPE_P    =  5;
const int EL_TYPE_F    =  6;
const int EL_TYPE_Cl   =  7;
const int EL_TYPE_Br   =  8;
const int EL_TYPE_I    =  9;
const int EL_TYPE_Met  = 10;
const int EL_TYPE_SIZE = 11;

// AutoDock4
#define AD_TYPE_C   	  0
#define AD_TYPE_A   	  1
#define AD_TYPE_N   	  2
#define AD_TYPE_O   	  3
#define AD_TYPE_P	  4
#define AD_TYPE_S         5
#define AD_TYPE_H         6 // non-polar hydrogen
#define AD_TYPE_F         7
#define AD_TYPE_I         8
#define AD_TYPE_NA        9
#define AD_TYPE_OA       10
#define AD_TYPE_SA       11
#define AD_TYPE_HD       12
#define AD_TYPE_Mg       13
#define AD_TYPE_Mn       14
#define AD_TYPE_Zn       15
#define AD_TYPE_Ca       16
#define AD_TYPE_Fe       17
#define AD_TYPE_Cl       18
#define AD_TYPE_Br       19
#define AD_TYPE_SIZE   20

// X-Score
const int XS_TYPE_C_H   =  0;
const int XS_TYPE_C_P   =  1;
const int XS_TYPE_C_A   =  2;
const int XS_TYPE_N_P   =  3;
const int XS_TYPE_N_D   =  4;
const int XS_TYPE_N_A   =  5;
const int XS_TYPE_N_DA  =  6;
const int XS_TYPE_O_P   =  7;
const int XS_TYPE_O_D   =  8;
const int XS_TYPE_O_A   =  9;
const int XS_TYPE_O_DA  = 10;
const int XS_TYPE_S_P   = 11;
const int XS_TYPE_S_H   = 12;
const int XS_TYPE_P_P   = 13;
const int XS_TYPE_F_H   = 14;
const int XS_TYPE_Cl_H  = 15;
const int XS_TYPE_Br_H  = 16;
const int XS_TYPE_I_H   = 17;
const int XS_TYPE_Met_D = 18;
const int XS_TYPE_SIZE  = 19;

// DrugScore-CSD
const int SY_TYPE_C_3   =  0;
const int SY_TYPE_C_2   =  1;
const int SY_TYPE_C_ar  =  2;
const int SY_TYPE_C_cat =  3;
const int SY_TYPE_N_3   =  4;
const int SY_TYPE_N_ar  =  5;
const int SY_TYPE_N_am  =  6;
const int SY_TYPE_N_pl3 =  7;
const int SY_TYPE_O_3   =  8;
const int SY_TYPE_O_2   =  9;
const int SY_TYPE_O_co2 = 10;
const int SY_TYPE_S     = 11;
const int SY_TYPE_P     = 12;
const int SY_TYPE_F     = 13;
const int SY_TYPE_Cl    = 14;
const int SY_TYPE_Br    = 15;
const int SY_TYPE_I     = 16;
const int SY_TYPE_Met   = 17;
const int SY_TYPE_SIZE  = 18;

struct atom_kind {
    char name[2];
    double radius;
    double depth;
    double solvation;
    double volume;
    double covalent_radius; // from http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
};

// generated from edited AD4_parameters.data using a script,
// then covalent radius added from en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
struct atom_kind atom_kind_data[] = { // name, radius, depth, solvation parameter, volume, covalent radius
        { "C",    2.00000,    0.15000,   -0.00143,   33.51030,   0.77}, //  0
        { "A",    2.00000,    0.15000,   -0.00052,   33.51030,   0.77}, //  1
        { "N",    1.75000,    0.16000,   -0.00162,   22.44930,   0.75}, //  2
        { "O",    1.60000,    0.20000,   -0.00251,   17.15730,   0.73}, //  3
        { "P",    2.10000,    0.20000,   -0.00110,   38.79240,   1.06}, //  4
        { "S",    2.00000,    0.20000,   -0.00214,   33.51030,   1.02}, //  5
        { "H",    1.00000,    0.02000,    0.00051,    0.00000,   0.37}, //  6
        { "F",    1.54500,    0.08000,   -0.00110,   15.44800,   0.71}, //  7
        { "I",    2.36000,    0.55000,   -0.00110,   55.05850,   1.33}, //  8
        {"NA",    1.75000,    0.16000,   -0.00162,   22.44930,   0.75}, //  9
        {"OA",    1.60000,    0.20000,   -0.00251,   17.15730,   0.73}, // 10
        {"SA",    2.00000,    0.20000,   -0.00214,   33.51030,   1.02}, // 11
        {"HD",    1.00000,    0.02000,    0.00051,    0.00000,   0.37}, // 12
        {"Mg",    0.65000,    0.87500,   -0.00110,    1.56000,   1.30}, // 13
        {"Mn",    0.65000,    0.87500,   -0.00110,    2.14000,   1.39}, // 14
        {"Zn",    0.74000,    0.55000,   -0.00110,    1.70000,   1.31}, // 15
        {"Ca",    0.99000,    0.55000,   -0.00110,    2.77000,   1.74}, // 16
        {"Fe",    0.65000,    0.01000,   -0.00110,    1.84000,   1.25}, // 17
        {"Cl",    2.04500,    0.27600,   -0.00110,   35.82350,   0.99}, // 18
        {"Br",    2.16500,    0.38900,   -0.00110,   42.56610,   1.14}  // 19
};

double metal_solvation_parameter = -0.00110;

double metal_covalent_radius = 1.75; // for metals not on the list // FIXME this info should be moved to non_ad_metals

int atom_kinds_size =  sizeof(atom_kind_data) / sizeof(struct atom_kind);

struct atom_equivalence {
    char name[2];
    char to[1];
};

struct atom_equivalence atom_equivalence_data[] = {
        {"Se",  "S"}
};

int atom_equivalences_size = sizeof(atom_equivalence_data) / sizeof(struct atom_equivalence);

struct acceptor_kind {
    int ad_type;
    double radius;
    double depth;
};

struct acceptor_kind acceptor_kind_data[] = { // ad_type, optimal length, depth
       // {AD_TYPE_NA, 1.9, 5.0}, this type can't be complied cause const int still is variable for compiler
        {9, 1.9, 5.0},
        //{AD_TYPE_OA, 1.9, 5.0},
        {10, 1.9, 5.0},
       //{AD_TYPE_SA, 2.5, 1.0}
        {11, 2.5, 1.0}
};

int acceptor_kinds_size = sizeof(acceptor_kind_data) / sizeof(struct acceptor_kind);

inline bool ad_is_hydrogen(int ad) {
    return ad == AD_TYPE_H || ad == AD_TYPE_HD;
}

inline bool ad_is_heteroatom(int ad) { // returns false for ad >= AD_TYPE_SIZE
    return ad != AD_TYPE_A && ad != AD_TYPE_C  &&
           ad != AD_TYPE_H && ad != AD_TYPE_HD &&
           ad < AD_TYPE_SIZE;
}

inline int ad_type_to_el_type(int t) {
    switch(t) {
        case AD_TYPE_C    : return EL_TYPE_C;
        case AD_TYPE_A    : return EL_TYPE_C;
        case AD_TYPE_N    : return EL_TYPE_N;
        case AD_TYPE_O    : return EL_TYPE_O;
        case AD_TYPE_P    : return EL_TYPE_P;
        case AD_TYPE_S    : return EL_TYPE_S;
        case AD_TYPE_H    : return EL_TYPE_H;
        case AD_TYPE_F    : return EL_TYPE_F;
        case AD_TYPE_I    : return EL_TYPE_I;
        case AD_TYPE_NA   : return EL_TYPE_N;
        case AD_TYPE_OA   : return EL_TYPE_O;
        case AD_TYPE_SA   : return EL_TYPE_S;
        case AD_TYPE_HD   : return EL_TYPE_H;
        case AD_TYPE_Mg   : return EL_TYPE_Met;
        case AD_TYPE_Mn   : return EL_TYPE_Met;
        case AD_TYPE_Zn   : return EL_TYPE_Met;
        case AD_TYPE_Ca   : return EL_TYPE_Met;
        case AD_TYPE_Fe   : return EL_TYPE_Met;
        case AD_TYPE_Cl   : return EL_TYPE_Cl;
        case AD_TYPE_Br   : return EL_TYPE_Br;
        case AD_TYPE_SIZE : return EL_TYPE_SIZE;
    }
    return EL_TYPE_SIZE; // to placate the compiler in case of warnings - it should never get here though
}

double xs_vdw_radii[] = {
    2.000000, // C_H
	2.000000, // C_P
	1.900000, // C_A
	1.700000, // N_P
	1.700000, // N_Ds
	1.700000, // N_A
	1.700000, // N_DA
	1.600000, // O_P
	1.600000, // O_D
	1.600000, // O_A
	1.600000, // O_DA
	2.000000, // S_P
	2.000000, // S_A
	2.100000, // P_P
	1.500000, // F_H
	1.800000, // Cl_H
	2.000000, // Br_H
	2.200000, // I_H
	1.200000  // Met_D
};

inline double xs_radius(int t) {
    int n = sizeof(xs_vdw_radii) / sizeof(double);
    return xs_vdw_radii[t];
}

char non_ad_metal_names[][2] = { // expand as necessary
        "Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};


inline bool xs_is_hydrophobic(int xs) {
    return xs == XS_TYPE_C_H  ||
		   xs == XS_TYPE_C_A  ||
		   xs == XS_TYPE_S_H  ||
		   xs == XS_TYPE_F_H  ||
		   xs == XS_TYPE_Cl_H ||
		   xs == XS_TYPE_Br_H || 
		   xs == XS_TYPE_I_H;
}

inline bool xs_is_acceptor(int xs) {
    return xs == XS_TYPE_N_A ||
           xs == XS_TYPE_N_DA ||
           xs == XS_TYPE_O_A ||
           xs == XS_TYPE_O_DA;
}

inline bool xs_is_donor(int xs) {
    return xs == XS_TYPE_N_D ||
           xs == XS_TYPE_N_DA ||
           xs == XS_TYPE_O_D ||
           xs == XS_TYPE_O_DA ||
           xs == XS_TYPE_Met_D;
}

inline bool xs_donor_acceptor(int t1, int t2) {
    return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(int t1, int t2) {
    return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

inline struct atom_kind ad_type_property(int i) {
    return atom_kind_data[i];
}

inline double max_covalent_radius() {
    double tmp = 0;
    for(int i=0; i<atom_kinds_size; i++)
        if(atom_kind_data[i].covalent_radius > tmp)
            tmp = atom_kind_data[i].covalent_radius;
    return tmp;
}
// transfer atom_type.h here
/*we store the atoms in double:
 *saveAtoms[i*25] = (double)tmp.el;
 *saveAtoms[i*25+1] = (double)tmp.ad;
 *saveAtoms[i*25+2] = (double)tmp.xs;
 *saveAtoms[i*25+3] = (double)tmp.sy;

struct atom_type {
    enum t {
        EL, AD, XS, SY
    };
    int el, ad, xs, sy;
};
*/
/*enum t {
    EL, AD, XS, SY
};
*/
inline int get(enum t atom_typing_used, double *atoms) {
    switch(atom_typing_used) {
        case EL: return (int)atoms[0]; //el
        case AD: return (int)atoms[1]; //ad
        case XS: return (int)atoms[2]; //xs
        case SY: return (int)atoms[3]; //sy
        default: return DBL_EPSILON;
    }
}

inline bool is_hydrogen(double *atoms) {
    return ad_is_hydrogen((int)atoms[1]); //ad
}

inline bool is_heteroatom(double *atoms) {
    return ad_is_heteroatom((int)atoms[1]) || (int)atoms[2] == XS_TYPE_Met_D;
}

inline bool acceptable_type(double *atoms) {
    return (int)atoms[1] < AD_TYPE_SIZE || (int)atoms[2] == XS_TYPE_Met_D;
}

inline void assign_el(double *atoms) {
    //type->el = ad_type_to_el_type(type->ad);
    atoms[0] = (double)ad_type_to_el_type((int)atoms[1]);
    if((int)atoms[1] == AD_TYPE_SIZE && (int)atoms[2] == XS_TYPE_Met_D)
        atoms[0] = (double)EL_TYPE_Met;
}

inline bool same_element(double *atomsb, double *atomsa) { // does not distinguish metals or unassigned types
    //return type->el == a->el;
    return (atomsb[0] == atomsa[0]);
}

inline double covalent_radius(double *atoms) {
    if((int)atoms[1] < AD_TYPE_SIZE)
        return ad_type_property((int)atoms[1]).covalent_radius;
    else if((int)atoms[2] == XS_TYPE_Met_D)
        return metal_covalent_radius;
    return 0; // never happens - placating the compiler
}

inline double optimal_covalent_bond_length(double *atoms) {
    return covalent_radius(atoms) + covalent_radius(atoms);
}

inline int num_atom_types(enum t atom_typing_used) {
    switch(atom_typing_used) {
        case EL: return EL_TYPE_SIZE;
        case AD: return AD_TYPE_SIZE;
        case XS: return XS_TYPE_SIZE;
        case SY: return SY_TYPE_SIZE;
        default: return INT_MAX;
    }
}

inline int get_type_pair_index(enum t atom_typing_used, double *a, double *b) {
    // throws error if any arg is unassigned in the given typing scheme
    int n = num_atom_types(atom_typing_used);

    int i = get(atom_typing_used, a);
    int j = get(atom_typing_used, b);

    if(i <= j)
        return triangular_matrix_index(n, i, j);
    else
        return triangular_matrix_index(n, j, i);
}

#endif //ATOM_SLAVE_H
