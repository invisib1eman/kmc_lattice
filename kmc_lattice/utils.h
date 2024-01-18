//NANOROD: utils.h Utilities Heade (Revision Date: Oct 27, 2023)
#ifndef _UTILS_H
#define _UTILS_H
#include "header.h"
#include "molecule.h"
#include "lattice.h"
#include "aggregate.h"
#include "hbond.h"
#include "xyz.h"
vector<int> Lattice_index2_4index(int LatticeID,int Ng,int Nbasis);
int Lattice_4index2_index(vector<int> Lattice4ID,int Ng,int Nbasis);
XYZ Lattice4ID2XYZ(vector<int> Lattice4ID,vector<XYZ> BasisPoints,vector<XYZ> LatticeVectors);
vector<int> Translate_Lattice4ID(vector<int> Lattice4ID,vector<int> Lattice4ID_translate,int Ng);
int bondarm(int arm);
int neighborarm(int arm);
void DFSUtil(int index_ag,int v,bool visited[],vector<Molecule>& M,Aggregate& new_aggregate);
int newvertype(int oldvertype, int changevertype);
bool bondstate(int arm1,int arm2);


#endif
