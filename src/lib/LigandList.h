#ifndef LIGANDLIST_H
#define LIGANDLIST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LIGAND_NAME_LENGTH 25

struct ligands {
    char ligandName[MAX_LIGAND_NAME_LENGTH];
    struct ligands *next;
};

typedef struct ligands Ligand;

struct ligandList {
    Ligand *head;
};

typedef struct ligandList LigandList;

/* initialize an empty list */
void InitLigandList (LigandList *lgndList);

/* insert val at front */
void InsertFront(LigandList *lgndList, char lName[MAX_LIGAND_NAME_LENGTH]);

/* insert val at back */
void InsertBack(LigandList *lgndList, char lName[MAX_LIGAND_NAME_LENGTH]);

/* returns list length */
unsigned int length (LigandList *lgndList);

/* deletes list and clear the memory occupied by list item */
void ClearLigandList (LigandList *lgndList);

/* print the list items */
void PrintLigands(LigandList *lgndList);

#endif // LIGANDLIST_H
