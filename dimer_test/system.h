//NANOROD: system.h System Class (Revision Date: Oct 27,2023)
#ifndef _SYSTEM_H
#define _SYSTEM_H
#include "header.h"
#include "molecule.h"
#include "utils.h"
#include "lattice.h"
#include "aggregate.h"
#include "hbond.h"
#include "xyz.h"
class System
{
  public:
    //Define lattice_type
    string Lattice_Type;
    //define lists/vectors of the systems
    vector<Molecule> M; //List of molecules
    Lattice lattice;//create lattice
    list<hbond> H;//hbond list
    vector<Aggregate> Ag;//list of aggregates
    //define constants of systems
    const gsl_rng_type * gsl_T;
    gsl_rng * gsl_r;
    string Description;
    int NMOL; //Number of molecules
    int NAg;//Number of aggregates
    int Ng; //Length of box
    int Nbasis=2; //Number of basis points in one unit lattice
    int Nlattice;//number of lattice points
    vector<XYZ> BasisPoints;
    vector<XYZ> LatticeVectors;
    vector<vector<int>> NeighborLista;//list of neighbor lattice points on a type sublattice
    vector<vector<int>> NeighborListb;//list of neighbor lattice points on b type sublattice
    unsigned int GSL_SEED=static_cast<unsigned>(time(nullptr)); //Seed of random number generator
    int nsweep=0; //Number of MC sweeps
    int MCstep=1; //MC step size
    double deltat=0; //Timestep
    double E_1=10;//hbond dis enthalpy
    double free_bond_freeenergy=-1;//free bond entropy
    double D=1;//D=kT/6pietaR is constant
    double K_Dp=4*D*2/pow((2.26*sqrt(3)),2);//Diffusion rate inside the plane
    double K_Dz=4*D*1/pow(2.12,2);//Diffusion rate in the z-direction
    double K_Dtotal=K_Dp+K_Dz;//total diffusion rate
    double K_swap=0.01;//arm rotation rate (for initialization of arm rotation type bond breaking)
    double K_Break=0.01;//Bond Break rate (translational type)
    double K_total=K_Dtotal+K_Break+K_swap;//total rate
    double P_swap=K_swap/K_total;//probability of arm rotation
    double P_break=K_Break/K_total;//probability of bond break
    double P_Dz=K_Dz/K_Dtotal;//probability of diffusion in z-direction
    //double P_form=0.25;//probability of bond formation
    double total_time=0;
    void ReadInput(int argc, char *argv[])
    {
        
        
        options_description desc("Usage:\nNANOROD <options>");
    
        desc.add_options()
        ("help,h", "print usage message")
        ("Lattice_type,l", value<string>(&Lattice_Type)->default_value("tilted_honeycomb"), "Lattice type (default tilted_honeycomb)")
        ("NMOL,N", value<int>(&NMOL)->default_value(400), "#molecules (default 400)")
        ("box_length,L", value<int>(&Ng)->default_value(40.0), "length of box (default 40.0)")
        ("time,t", value<double>(&total_time)->default_value(500.0), "total time in tau_0 units (default 500.0)")
        ("E_1,a", value<double>(&E_1)->default_value(10.0), "E_association (default 10)")
        ("K_break,b",value<double>(&K_Break)->default_value(0.01),"Break frequency (default 0.01)")
        ("K_swap,s",value<double>(&K_swap)->default_value(0.01),"Arm rotation frequency (default 0.01)")
        //("GSL_SEED,g", value<int>(&GSL_SEED)->default_value(10), "seed for the RNG (default 10)")
        ("Description,D", value<string>(&Description)->default_value("nanorod"), "Description (default nanorod)");
        
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
        
        if (vm.count("help"))
        {
            cout << desc << "\n";
            exit(1);
        }
        
        gsl_rng_env_setup();
          
        gsl_T = gsl_rng_default;
        gsl_r = gsl_rng_alloc(gsl_T);
        gsl_rng_set(gsl_r, GSL_SEED);
        deltat=1.0/K_total;
        cout<<deltat<<endl;
        nsweep=int(ceil(total_time/deltat));
        K_total=K_Dtotal+K_Break+K_swap;
        P_swap=K_swap/K_total;//probability of arm rotation
        P_break=K_Break/K_total;//probability of bond break
        Nlattice=Ng*Ng*Ng*Nbasis;
        lattice.Ng=Ng;
        lattice.Nbasis=Nbasis;
        lattice.Resize_lattice();
        BasisPoints.push_back(XYZ(0,0,0));
        BasisPoints.push_back(XYZ(2.26*sqrt(3)/2,-2.26/2,1.06));
        LatticeVectors.push_back(XYZ(2.26*sqrt(3),0,0));
        LatticeVectors.push_back(XYZ(2.26*sqrt(3)/2,2.26*3/2,0));
        LatticeVectors.push_back(XYZ(0,0,2.12));
        vector<int> Na1{0,0,0,1};
        NeighborLista.push_back(Na1);
        vector<int> Na2{0,0,-1,1};
        NeighborLista.push_back(Na2);
        vector<int> Na3{-1,0,0,1};
        NeighborLista.push_back(Na3);
        vector<int> Na4{-1,0,-1,1};
        NeighborLista.push_back(Na4);
        vector<int> Na5{-1,1,0,1};
        NeighborLista.push_back(Na5);
        vector<int> Na6{-1,1,-1,1};
        NeighborLista.push_back(Na6);
        for(vector<int> Nai:NeighborLista)
        {
            vector<int> Nbj(4);
            Nbj[0]=-Nai[0];
            Nbj[1]=-Nai[1];
            Nbj[2]=-Nai[2];
            Nbj[3]=-Nai[3];
            NeighborListb.push_back(Nbj);
        }
    }
    void Create();
    //void WriteMol2(int timestep);
    void WriteDump(int timestep);
    void WriteBond(int timestep);
    void WriteAggregate(int timestep);
    void WriteLattice(int timestep);
    
};
#endif
