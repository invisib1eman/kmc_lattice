//Define lattice
#ifndef _LATTICE_H
#define _LATTICE_H
#include "header.h"
class Lattice
{
public:
    vector<int> molid_list;//molid of the molecule that occupy the lattice point,-1 if empty
    int Nbasis;//number of basis points in one unit lattice
    int Ng;//number of lattices in one direction
    Lattice()
    {
        Nbasis=1;
        Ng=1;
        molid_list.resize(Nbasis*Ng*Ng*Ng,-1);
    }
    void Resize_lattice()
    {
        molid_list.resize(Nbasis*Ng*Ng*Ng,-1);
    }
    /*XYZ index_to_xyz(Latticeindex m)
    {
        return XYZ(LatticeVectors[0].x*m.n1+LatticeVectors[1].x*m.n2+LatticeVectors[2].x*m.n3+BasisPoints[m.n4].x,LatticeVectors[0].y*m.n1+LatticeVectors[1].y*m.n2+LatticeVectors[2].y*m.n3+BasisPoints[m.n4].y,LatticeVectors[0].z*m.n1+LatticeVectors[1].z*m.n2+LatticeVectors[2].z*m.n3+BasisPoints[m.n4].z);
    }*/
    
    
};
#endif
