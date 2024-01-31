//NANOROD: molecule.h Molecule Class (Revision Date:Oct 27, 2023)
#ifndef _MOLECULE_H
#define _MOLECULE_H
#include "header.h"
#include "hbond.h"
class Molecule
{
public:
    int MOL_ID;//molecule ID
    int N_VER; //Number of vertices
    vector<int> Lattice4ID;//lattice id in form of (nx,ny,nz,nbasis)
    int LatticeID;//lattice id(n=Nbasis*(Ng*Ng*nx+Ng*ny+nz)+nbasis)
    int AID;//in which aggregate
    int AsubID;//id in the aggregate
    int nbonds;//number of bonds
    vector<hbond> hbond_list;//list of hbonded neighbors and information
    vector<bool> bondstate;
    Molecule() //HARD CODE BASED ON VERTEX INFORMATION,unit length=1nm,arm length=1.1nm
    {
        N_VER=6;
        Lattice4ID.resize(4);
        LatticeID=0;
        nbonds=0;
        Lattice4ID={0,0,0,0};
        bondstate.push_back(false);
        bondstate.push_back(false);
        bondstate.push_back(false);
        bondstate.push_back(false);
        bondstate.push_back(false);
        bondstate.push_back(false);
        
    }
};
#endif

