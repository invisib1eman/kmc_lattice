//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    Mnew=S.M;
    S.WriteDump(0);
    double accept=0.0;
    
    int nsample=S.nsweep/100;
    if(nsample<1)
      nsample=1;
    
    //energy=TotalEnergy();
   
    
    WriteTemplate();
	LogProfile(0,accept);
    cout<<S.nsweep<<endl;
    
    for(int i=1; i<=S.nsweep; i++)
    {
        time+=S.deltat;
        accept+=MoveMolecule();
        
        if(i%nsample==0)
        {
            
            LogProfile(i,accept);
            //WriteEnergy(i);
            S.WriteLattice(i);
            S.WriteDump(i);
            S.WriteBond(i);
            S.WriteAggregate(i);
            accept=0.0;
            bondform=0;
            bondbreak=0;
        }
    }
   S.WriteDump(S.nsweep);
}

void MC::WriteTemplate()
{
    ofstream out;
    char FileName[200];
    sprintf(FileName,"N%d_L%d_E_dis%.1f_S%.4f_T%.3f_B%.3fMC.log",S.NMOL,S.Ng,S.E_1,S.K_swap,S.total_time,S.K_Break);
    out.open(FileName);
    out<<setw(12)<<"E_a"<<setw(12)<<S.E_1<<"K_break"<<setw(12)<<S.K_Break<<"Box_length"<<setw(12)<<S.Ng<<"N_particle"<<setw(12)<<S.NMOL<<"Simulation time"<<setw(12)<<S.total_time<<endl;
    out<<setw(12)<<"sweep"<<"\t"<<setw(12)<<"time"<<"\t"<<setw(12)<<"energy"<<"\t"<<setw(8)<<"accept"<<setw(8)<<"N_aggregate"<<setw(8)<<"bond_forming"<<setw(8)<<"bond_breaking"<<endl;
    out.close();
}

void MC::LogProfile(int i, double accept)
{
    ofstream out;
    char FileName[200];
    sprintf(FileName,"N%d_L%d_E_dis%.1f_S%.4f_T%.3f_B%.3fMC.log",S.NMOL,S.Ng,S.E_1,S.K_swap,S.total_time,S.K_Break);
    out.open(FileName, ios::app);
    out<<setw(12)<<i<<"\t"<<setw(12)<<time<<"\t"<<setw(12)<<energy<<"\t"<<setw(8)<<accept<<setw(8)<<S.NAg<<setw(8)<<bondform<<setw(8)<<bondbreak<<endl;
    
    out.close();
}



//Returns acceptance fraction
double MC::MoveMolecule()
{
    double accept=0.0;//accept events
    int N_break=0;
    ofstream out;
    out.open("events.txt",ios::app);
    ofstream out2;
    out2.open("Error.txt",ios::app);
    
    //break bonds
    list<hbond>::iterator it;
    double rand=gsl_rng_uniform(S.gsl_r);
    int nbond_mol=0;
    int nbond=S.H.size();
    for(int i=0;i<S.NMOL;i++)
    {
        nbond_mol+=S.M[i].hbond_list.size();
    }
    nbond_mol=nbond_mol/2;
    if(nbond_mol!=nbond)
    {
        cout<<"Hbondlist no match"<<endl;
        exit(3);
    }
    int nmol_aggregate=0;
    for(int i=0;i<S.NAg;i++)
    {
        nmol_aggregate+=S.Ag[i].M_A.size();
    }
    if(nmol_aggregate!=S.NMOL)
    {
        cout<<"Aggregate no match"<<endl;
        exit(3);
    }
    if(rand<S.P_swap)
    {
        for(int i=0;i<S.NMOL*3;i++)
        {
            int molindex=gsl_rng_uniform_int(S.gsl_r,S.NMOL);
            int arm_index=gsl_rng_uniform_int(S.gsl_r,3);
            int arm_index1=arm_index*2;
            int arm_index2=arm_index1+1;
            

            
            if(S.M[molindex].bondstate[arm_index1]==false && S.M[molindex].bondstate[arm_index2]==false)
            {   
                int rotate_1=2*gsl_rng_uniform_int(S.gsl_r,2)-1;
                int rotate_2=2*gsl_rng_uniform_int(S.gsl_r,2)-1;
                int newvertype1=newvertype(S.M[molindex].vertype[arm_index1],rotate_1);
                int newvertype2=newvertype(S.M[molindex].vertype[arm_index2],rotate_2);
                S.M[molindex].vertype[arm_index1]=newvertype1;
                S.M[molindex].vertype[arm_index2]=newvertype2;
            } 
            

            
        }
    }
    if(rand>S.P_swap&&rand<(S.P_break+S.P_swap))//the relative frequency of breaking bond to frequency of diffusion step
    {
        for(it=S.H.begin();it!=S.H.end();)
        {
            hbond old_hbond=*it;
            //calculate bond_dissociation energy
            double E_dis=0;
            
            int free_bonds=0;
            E_dis+=S.E_1;//the basic enthalpy change of one bond
            //count # of freed bonds
            //find neighbor arms,first the one of M1, then the one of M2
            int neighborarm1=neighborarm(old_hbond.arm1);
            free_bonds+=2;
            int bonded_index1;
            vector<hbond> old_hbondlist=S.M[old_hbond.M1].hbond_list;
            for(int p=0;p<old_hbondlist.size();p++)
            {
                if(old_hbondlist[p].arm1==neighborarm1)
                {
                    free_bonds-=2;
                }
                if(old_hbondlist[p].arm1==old_hbond.arm1)
                {
                    bonded_index1=p;
                }
            } 
            int neighborarm2=neighborarm(old_hbond.arm2);
            free_bonds+=2;
            vector<hbond> bonded_neighbor_hbondlist=S.M[old_hbond.M2].hbond_list;
            int bonded_index2;
            for(int p=0;p<bonded_neighbor_hbondlist.size();p++)
            {
                if(bonded_neighbor_hbondlist[p].arm1==neighborarm2)
                {
                    free_bonds-=2;
                }
                if(bonded_neighbor_hbondlist[p].arm1==old_hbond.arm2)
                    bonded_index2=p;
                
            }    
                         
            E_dis+=free_bonds*S.free_bond_freeenergy;
                
            if(Arrhenius(1,E_dis,gsl_rng_uniform(S.gsl_r)))
            {
                //break bond
                out<<"Break bond"<<setw(12)<<old_hbond.M1<<setw(12)<<old_hbond.M2<<setw(12)<<old_hbond.arm1<<setw(12)<<old_hbond.arm2<<endl;
                bondbreak+=1;
                S.M[old_hbond.M1].bondstate[old_hbond.arm1]=false;
                S.M[old_hbond.M1].hbond_list[bonded_index1]=S.M[old_hbond.M1].hbond_list.back();
                S.M[old_hbond.M1].nbonds-=1;
                S.M[old_hbond.M1].hbond_list.pop_back();
                S.M[old_hbond.M2].bondstate[old_hbond.arm2]=false;
                S.M[old_hbond.M2].hbond_list[bonded_index2]=S.M[old_hbond.M2].hbond_list.back();
                S.M[old_hbond.M2].nbonds-=1;
                S.M[old_hbond.M2].hbond_list.pop_back();
                it=S.H.erase(it);        
            }
            else
            {
                ++it;
            }
        }
        //recalculate clusters
        bool* visited=new bool[S.NMOL];
        for(int i=0;i<S.NMOL;i++)
        {
            visited[i]=false;
        }
        S.Ag.clear();
        int index_ag=0;
        for (int v = 0; v < S.NMOL; v++) 
        {
            if (visited[v] == false) 
            {
                //Create new aggregate
                
                Aggregate new_aggregate;
                
                // Visited v and its neighbors
                
                DFSUtil(index_ag,v, visited,S.M,new_aggregate);
                S.Ag.push_back(new_aggregate);
                index_ag++;
            }
            
        }
        

        S.NAg=S.Ag.size();
        
    }
    nbond_mol=0;
    nbond=S.H.size();
    for(int i=0;i<S.NMOL;i++)
    {
        nbond_mol+=S.M[i].hbond_list.size();
    }
    nbond_mol=nbond_mol/2;
    if(nbond_mol!=nbond)
    {
        cout<<"Hbondlist no match"<<endl;
        exit(3);
    }
    nmol_aggregate=0;
    for(int i=0;i<S.NAg;i++)
    {
        nmol_aggregate+=S.Ag[i].M_A.size();
    }
    if(nmol_aggregate!=S.NMOL)
    {
        cout<<"Aggregate no match"<<endl;
        exit(3);
    }
    //diffusion limited bond formation
    //Aggregate diffusion
    //vector<Aggregate>::iterator ita;
    if(rand>(S.P_break+S.P_swap))
    {
        for(int a=0;a<S.NAg;a++)
        {
            
            //check total number of mol in lattice
            /*
            int Nmol_lattice=0;
            for(int i=0;i<S.Nlattice;i++)
            {
                if(S.lattice.molid_list[i]!=-1)
                    Nmol_lattice+=1;
            }
            if(Nmol_lattice!=S.NMOL)
            {
                cout<<Nmol_lattice<<"Wrong total molecules"<<endl;
                exit(3);
            }
            for(int i=0;i<S.NMOL;i++)
            {
                if(S.lattice.molid_list[S.M[i].LatticeID]!=i)
                {
                    cout<<S.lattice.molid_list[S.M[i].LatticeID]<<"\t"<<i<<"Wrong lattice"<<endl;
                    exit(3);
                }
            }
            */
            //Randomly choose one aggregate by equal chance, the diffusion coefficient will be rescaled in the acception rate according to the aggregate size
            int index=gsl_rng_uniform_int(S.gsl_r,S.NAg);
            Aggregate old_aggregate=S.Ag[index];
            //Rescale diffusion coefficient
            rand=gsl_rng_uniform(S.gsl_r);
            if(rand>(1/pow(old_aggregate.n,1/3)))
            {
                continue;
            }
            //Choose Diffusion direction
            vector<int> Movement(4,0);
            //Choose Diffusion in the plane or in the z-direction
            rand=gsl_rng_uniform(S.gsl_r);
            if(rand<S.P_Dz)//Diffusion in z-direction
            {
                //Choose Diffusion direction: Up or down
                rand=gsl_rng_uniform(S.gsl_r);
                if(rand<0.5)//Up
                {
                    
                    Movement[2]=1;
                }
                else//Down
                {
                    vector<int> Movement(4,0);
                    Movement[2]=-1;
                }
            }
            else
            {
                int nrand=gsl_rng_uniform_int(S.gsl_r,6);
                if(nrand==0)
                    Movement[0]=1;
                else if(nrand==1)
                    Movement[0]=-1;
                else if(nrand==2)
                    Movement[1]=1;
                else if(nrand==3)
                    Movement[1]=-1;
                else if(nrand==4)
                {
                    Movement[0]=1;
                    Movement[1]=-1;
                }    
                else if(nrand==5)
                {
                    Movement[0]=-1;
                    Movement[1]=1;
                }
            }
            //check if the new lattice points are occupied or will have unmatched bonds
            //ghost all particles
            for(int i=0;i<old_aggregate.n;i++)
            {
                S.lattice.molid_list[S.M[old_aggregate.M_A[i]].LatticeID]=-1;
            }
            bool occupied=false;
            bool unmatched=false;
            for(int i=0;i<old_aggregate.n;i++)
            {
                vector<int> old_Lattice4ID=S.M[old_aggregate.M_A[i]].Lattice4ID;
                vector<int> new_Lattice4ID=Translate_Lattice4ID(old_Lattice4ID,Movement,S.Ng);
                //check if the new lattice point is occupied
                int new_LatticeID=Lattice_4index2_index(new_Lattice4ID,S.Ng,S.Nbasis);
                if(S.lattice.molid_list[new_LatticeID]!=-1)
                {
                    occupied=true;
                    break;
                }   
                //check if will have unmatched bonds
                vector<vector<int>> Neighborlist;
                if(new_Lattice4ID[3]==0)
                {
                    Neighborlist=S.NeighborLista;
                }
                else
                {
                    Neighborlist=S.NeighborListb;
                }
                for(int k=0;k<6;k++)
                {
                    vector<int> neighbor_Lattice4ID=Translate_Lattice4ID(new_Lattice4ID,Neighborlist[k],S.Ng);
                    int neighbor_LatticeID=Lattice_4index2_index(neighbor_Lattice4ID,S.Ng,S.Nbasis);
                    if(S.lattice.molid_list[neighbor_LatticeID]!=-1)
                    {
                        if(!bondstate(S.M[old_aggregate.M_A[i]].vertype[k],S.M[S.lattice.molid_list[neighbor_LatticeID]].vertype[k]))
                            unmatched=true;
                    }
                    
                }
            }
            if(occupied==false && unmatched==false)
            {
                //update lattice
                for(int i=0;i<old_aggregate.n;i++)
                {
                    int mol_id=old_aggregate.M_A[i];
                    vector<int> old_Lattice4ID=S.M[old_aggregate.M_A[i]].Lattice4ID;
                    vector<int> new_Lattice4ID=Translate_Lattice4ID(old_Lattice4ID,Movement,S.Ng);
                    int new_LatticeID=Lattice_4index2_index(new_Lattice4ID,S.Ng,S.Nbasis);
                    S.lattice.molid_list[new_LatticeID]=old_aggregate.M_A[i];
                    S.M[mol_id].LatticeID=new_LatticeID;
                    S.M[mol_id].Lattice4ID=new_Lattice4ID;
                }
                //Check if form new bonds
                //iterate about new neighbor vertices
                for(int i=0;i<old_aggregate.n;i++)
                {
                    
                    int mol_id=old_aggregate.M_A[i];
                    if(S.M[mol_id].bonded==false)
                    {
                        vector<int> new_Lattice4ID=S.M[mol_id].Lattice4ID;
                        vector<vector<int>> Neighborlist;
                        
                        //Determine which sublattice the molecule is on
                        if(new_Lattice4ID[3]==0)
                        {
                            Neighborlist=S.NeighborLista;
                        }
                        else
                        {
                            Neighborlist=S.NeighborListb;
                        }
                        for(int k=0;k<6;k++)
                        {
                            
                            
                            if(S.M[mol_id].bondstate[k]==true)
                                continue;
                            else
                            {
                                vector<int> neighbor_Lattice4ID=Translate_Lattice4ID(new_Lattice4ID,Neighborlist[k],S.Ng);
                                int neighbor_LatticeID=Lattice_4index2_index(neighbor_Lattice4ID,S.Ng,S.Nbasis);
                                if(neighbor_Lattice4ID[3]==new_Lattice4ID[3])
                                {
                                    cout<<"Wrong bond"<<endl;
                                    exit(3);
                                }
                                if(neighbor_LatticeID>S.Nlattice)
                                {
                                    cout<<"wrong neighbor lattice id"<<endl;
                                    exit(3);
                                }
                                if(S.lattice.molid_list[neighbor_LatticeID]!=-1)
                                {
                                    if(S.M[S.lattice.molid_list[neighbor_LatticeID]].bonded==false)
                                    {
                                    rand=gsl_rng_uniform(S.gsl_r);
                                    if(bondstate(S.M[old_aggregate.M_A[i]].vertype[k],S.M[S.lattice.molid_list[neighbor_LatticeID]].vertype[k]))
                                    {
                                        //form bond
                                        bondform+=1;
                                        //Update molecules
                                        int neighbormol_id=S.lattice.molid_list[neighbor_LatticeID];
                                        S.M[mol_id].bondstate[k]=true;
                                        S.M[mol_id].nbonds+=1;
                                        S.M[mol_id].bonded=true;
                                        S.M[mol_id].hbond_list.push_back(hbond(mol_id,neighbormol_id,k,bondarm(k)));
                                        S.M[neighbormol_id].bondstate[bondarm(k)]=true;
                                        S.M[neighbormol_id].nbonds+=1;
                                        
                                        S.M[neighbormol_id].bonded=true;
                                        
                                        S.M[neighbormol_id].hbond_list.push_back(hbond(neighbormol_id,mol_id,bondarm(k),k));
                                        //Update hbondlist
                                        S.H.push_back(hbond(mol_id,neighbormol_id,k,bondarm(k)));
                                        out<<"form a hbond"<<setw(12)<<mol_id<<setw(12)<<neighbormol_id<<setw(12)<<k<<setw(12)<<bondarm(k)<<endl;
                                        //Update aggregate  
                                        if(S.M[mol_id].AID!=S.M[neighbormol_id].AID)//multiple bonds form between two aggregates
                                        {
                                            int AID1=S.M[mol_id].AID;
                                        
                                            
                                            int AID2=S.M[neighbormol_id].AID;
                                            if(AID1>AID2)
                                            {
                                            int tmp=AID1;
                                                AID1=AID2;
                                                AID2=tmp; 
                                            }
                                            Aggregate old_aggregate1=S.Ag[AID1];
                                            Aggregate old_aggregate2=S.Ag[AID2];
                                            //merge two aggregate
                                            Aggregate new_aggregate=old_aggregate1;
                                            new_aggregate.n+=old_aggregate2.n;
                                            for(int j=0;j<old_aggregate2.M_A.size();j++)
                                            {
                                                S.M[old_aggregate2.M_A[j]].AID=AID1;
                                                S.M[old_aggregate2.M_A[j]].AsubID+=old_aggregate1.n;
                                            }
                                            new_aggregate.M_A.insert(new_aggregate.M_A.end(),old_aggregate2.M_A.begin(),old_aggregate2.M_A.end());
                                            S.Ag[AID1]=new_aggregate;
                                            if(AID2<S.NAg-1)
                                            {
                                                S.Ag[AID2]=S.Ag.back();
                                                for(int j=0;j<S.Ag[AID2].M_A.size();j++)
                                                {
                                                    S.M[S.Ag[AID2].M_A[j]].AID=AID2;    
                                                }
                                            }
                                            S.Ag.pop_back();
                                            S.NAg-=1;
                                        }
                                        }

                                    }
                                }
                        }
                    }

                    }
                   
                }
                    
                    
                
                accept+=1.0;
            }
            else
            {
                //Deghost all particles
                for(int i=0;i<old_aggregate.n;i++)
                {
                    S.lattice.molid_list[S.M[old_aggregate.M_A[i]].LatticeID]=old_aggregate.M_A[i];
                }
            }

        }
    }   
            
    return accept/double(S.NMOL);                    
}
              
            
    
    
    
    
    



bool MC::Arrhenius(double A,double delta, double rand)
{
    
    if(A*(exp(-delta))>rand)
        {
        return true;}
    else
        return false;
    
    
}
/*int MC::bond_energy()
{
    double ebond=0;
    int nbond=S.H.size();
    //ebond=-nbond*S.E_1;
    return nbond;
}
int MC::bond_freeze_freenerngy()
{
    double efreeze=0;
    int nfreeze=0;
    for(int i=0;i<S.NMOL;i++)
    {
        Molecule m=S.M[i];
        for(int l=0;l<m.N_VER;l++)
        {
            if(m.vertype[l]=='I'||m.vertype[neighborarm(l)]=='I')
            nfreeze+=1;

        }
    }
    efreeze=nfreeze*S.free_bond_freeenergy;
    return nfreeze;
}
double MC::TotalEnergy()
{
    double totalenergy=0.0;
    /*for(int i;i<S.NMOL;i++)
    {
        for(int j=0;j<S.M[i].nbonds;j++)
        {
            totalenergy+=E.hbonde(S.M[i],S.M[S.M[i].hbond_list[j].M2],j);
            //totalenergy+=E.LJ(S.M[i],S.M[j]);
        }
        
    }*//*
    return totalenergy;
}
double MC::WriteEnergy(int timestep)
{
    ofstream out;
    out.open("energy.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"WCA"<<'\t'<<"FENE"<<'\t'<<"Angle"<<'\t'<<"Dihedral"<<'\t'<<"Bond"<<'\t'<<"Bond_freeze"<<endl;
    out<<'\t'<<bond_energy()<<'\t'<<bond_freeze_freenerngy()<<endl;
    out.close();
}
*/
