//NANOROD: system.cpp System Class Function Definitions (Revision Date: Oct 27, 2023)
#include "system.h"
void System::Create()
{
    
    cout<<"Creating System"<<endl;
    int i,j,k,l,n;
   
    //Fill particles Use random permutation
    bool flag;
    H.clear();//init hbond list
    //int NPair=NMOL/2;
    
    int count=0;
    /*
    for(int i=0;i<NPair;i++)
    {
        Molecule m1;//initial test for 2NMOL bonded particles
        m1.MOL_ID=2*i;
        m1.AID=i;
        m1.AsubID=0;
        m1.centre=XYZ(gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L);
        m1.gID=GridIndex_xyz(m1.centre,NGRID,GRIDL,L);
        G[m1.gID].n+=1;
        G[m1.gID].plist.push_back(m1.MOL_ID);
        m1.UpdateVertices();
        m1.nbonds=1;
        m1.hbond_list.push_back(hbond(2*i,2*i+1,0,3));
        m1.vertype[0]='I';
        M.push_back(m1);
        Molecule m2;
        m2.MOL_ID=2*i+1;
        m2.AID=i;
        m2.AsubID=1;
        m2.centre=image((XYZ(2.26,0,1.06)+m1.centre),L);
        m2.gID=GridIndex_xyz(m2.centre,NGRID,GRIDL,L);
        G[m2.gID].n+=1;
        G[m2.gID].plist.push_back(m2.MOL_ID);
        m2.orientation=angle_to_quarternion(M_PI/3,0,0);
        m2.UpdateVertices();
        m2.nbonds=1;
        m2.hbond_list.push_back(hbond(2*i+1,2*i,3,0));
        m2.vertype[3]='I';
        M.push_back(m2);
        if(min_d2(m1.ver[0],m2.ver[3],L)>0.58*0.58)
        {
            cout<<"wrong initialization"<<endl;
            exit(3);
        }
        Aggregate A;
        A.n=2;
        A.rg2=2.2;
        A.M_A.push_back(m1.MOL_ID);
        A.M_A.push_back(m2.MOL_ID);
        A.L=L;
        Ag.push_back(A);
        H.push_back(hbond(2*i,2*i+1,0,3));

    }*/
    //initialize a ring
    int N_ring=1;
    vector<Molecule> unit;
    vector<Molecule> ring;
    Molecule m1;
    m1.MOL_ID=0;
    m1.AID=0;//at first each molecule is one aggregate
    m1.AsubID=0;//every molecule is the only element of the aggregate
    int trial_lattice;
    while (-1)
    {
        trial_lattice=gsl_rng_uniform_int(gsl_r,Nlattice);
        if(lattice.molid_list[trial_lattice]==-1)
        {
            break;
        }
    }
    
    m1.LatticeID=trial_lattice;
    m1.Lattice4ID=Lattice_index2_4index(m1.LatticeID,Ng,Nbasis);
    m1.hbond_list.clear();
    m1.nbonds=4;
    //initialize vertypes of arms
    m1.bondstate[0]=true;
    m1.bondstate[1]=true;
    m1.bondstate[4]=true;
    m1.bondstate[5]=true;
    m1.vertype[0]=0;
    m1.vertype[neighborarm(0)]=0;
    m1.vertype[2]=1;
    m1.vertype[3]=3;
    m1.vertype[4]=0;
    m1.vertype[5]=0;
    m1.hbond_list.push_back(hbond(0,1,0,0));
    m1.hbond_list.push_back(hbond(0,5,4,4));
    m1.hbond_list.push_back(hbond(0,-5,1,1));
    m1.hbond_list.push_back(hbond(0,-1,5,5));

    
    Molecule m2;
    m2.MOL_ID=1;
    m2.AID=0;//at first each molecule is one aggregate
    m2.AsubID=1;//every molecule is the only element of the aggregate
    
    
    
    
    m2.Lattice4ID=Translate_Lattice4ID(m1.Lattice4ID,NeighborLista[0],Ng);
    m2.LatticeID=Lattice_4index2_index(m2.Lattice4ID,Ng,Nbasis);
    m2.hbond_list.clear();
    m2.nbonds=4;
    //initialize vertypes of arms
    m2.bondstate[0]=true;
    m2.bondstate[2]=true;
    m2.bondstate[1]=true;
    m2.bondstate[3]=true;
    m2.vertype[0]=3;
    m2.vertype[neighborarm(0)]=3;
    m2.vertype[2]=3;
    m2.vertype[3]=3;
    m2.vertype[4]=0;
    m2.vertype[5]=0;
    m2.hbond_list.push_back(hbond(1,0,0,0));
    m2.hbond_list.push_back(hbond(1,2,2,2));
    m2.hbond_list.push_back(hbond(1,6,1,1));
    m2.hbond_list.push_back(hbond(1,7,3,3 ));
    Molecule m3;
    m3.MOL_ID=2;
    m3.AID=0;//at first each molecule is one aggregate
    m3.AsubID=2;//every molecule is the only element of the aggregate
    
    
    
    vector<int> Movement(4,0);
    Movement[0]=1;
    m3.Lattice4ID=Translate_Lattice4ID(m1.Lattice4ID,Movement,Ng);
    m3.LatticeID=Lattice_4index2_index(m3.Lattice4ID,Ng,Nbasis);
    m3.hbond_list.clear();
    m3.nbonds=4;
    //initialize vertypes of arms
    m3.bondstate[2]=true;
    m3.bondstate[3]=true;
    m3.bondstate[4]=true;
    m3.bondstate[5]=true;
    m3.vertype[0]=3;
    m3.vertype[neighborarm(0)]=3;
    m3.vertype[2]=0;
    m3.vertype[3]=0;
    m3.vertype[4]=0;
    m3.vertype[5]=0;
    m3.hbond_list.push_back(hbond(2,1,2,2));
    m3.hbond_list.push_back(hbond(2,3,4,4));
    m3.hbond_list.push_back(hbond(2,-5,3,3));
    m3.hbond_list.push_back(hbond(2,-3,5,5));
    

    Molecule m4;
    m4.MOL_ID=3;
    m4.AID=0;//at first each molecule is one aggregate
    m4.AsubID=3;//every molecule is the only element of the aggregate
    
    
    
    vector<int> Movement2(4,0);
    Movement2[1]=1;
    m4.Lattice4ID=Translate_Lattice4ID(m2.Lattice4ID,Movement2,Ng);
    m4.LatticeID=Lattice_4index2_index(m4.Lattice4ID,Ng,Nbasis);
    m4.hbond_list.clear();
    m4.nbonds=4;
    //initialize vertypes of arms
    m4.bondstate[4]=true;
    m4.bondstate[5]=true;
    m4.bondstate[0]=true;
    m4.bondstate[1]=true;
    m4.vertype[0]=3;
    m4.vertype[neighborarm(0)]=3;
    m4.vertype[2]=0;
    m4.vertype[3]=0;
    m4.vertype[4]=3;
    m4.vertype[5]=3;
    m4.hbond_list.push_back(hbond(3,2,4,4));
    m4.hbond_list.push_back(hbond(3,4,0,0));
    m4.hbond_list.push_back(hbond(3,8,5,5));
    m4.hbond_list.push_back(hbond(3,10,1,1));
    Molecule m5;
    m5.MOL_ID=4;
    m5.AID=0;//at first each molecule is one aggregate
    m5.AsubID=4;//every molecule is the only element of the aggregate
    
    
    
    
    m5.Lattice4ID=Translate_Lattice4ID(m1.Lattice4ID,Movement2,Ng);
    m5.LatticeID=Lattice_4index2_index(m5.Lattice4ID,Ng,Nbasis);
    m5.hbond_list.clear();
    m5.nbonds=4;
    //initialize vertypes of arms
    m5.bondstate[0]=true;
    m5.bondstate[1]=true;
    m5.bondstate[2]=true;
    m5.bondstate[3]=true;
    m5.vertype[0]=0;
    m5.vertype[neighborarm(0)]=0;
    m5.vertype[2]=0;
    m5.vertype[3]=0;
    m5.vertype[4]=0;
    m5.vertype[5]=0;
    m5.hbond_list.push_back(hbond(4,3,0,0));
    m5.hbond_list.push_back(hbond(4,5,2,2));
    m5.hbond_list.push_back(hbond(4,-3,1,1));
    m5.hbond_list.push_back(hbond(4,-1,3,3));

    Molecule m6;
    m6.MOL_ID=5;
    m6.AID=0;//at first each molecule is one aggregate
    m6.AsubID=5;//every molecule is the only element of the aggregate
    
    
    
    vector<int> Movement3(4,0);
    Movement3[0]=-1;
    m6.Lattice4ID=Translate_Lattice4ID(m4.Lattice4ID,Movement3,Ng);
    m6.LatticeID=Lattice_4index2_index(m6.Lattice4ID,Ng,Nbasis);
    m6.hbond_list.clear();
    m6.nbonds=4;
    //initialize vertypes of arms
    m6.bondstate[2]=true;
    m6.bondstate[4]=true;
    m6.bondstate[3]=true;
    m6.bondstate[5]=true;
    m6.vertype[0]=3;
    m6.vertype[neighborarm(0)]=3;
    m6.vertype[2]=3;
    m6.vertype[3]=3;
    m6.vertype[4]=3;
    m6.vertype[5]=3;
    m6.hbond_list.push_back(hbond(5,4,2,2));
    m6.hbond_list.push_back(hbond(5,0,4,4));
    m6.hbond_list.push_back(hbond(5,10,3,3));
    m6.hbond_list.push_back(hbond(5,6,5,5));
    

    unit.push_back(m1);
    unit.push_back(m2);
    unit.push_back(m3);
    unit.push_back(m4);
    unit.push_back(m5);
    unit.push_back(m6);
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            Molecule newmolecule=unit[j];
            newmolecule.MOL_ID=i*6+j;
            newmolecule.AsubID=i*6+j;
            newmolecule.Lattice4ID[2]+=i;
            newmolecule.LatticeID=Lattice_4index2_index(newmolecule.Lattice4ID,Ng,Nbasis);
            for(int k=0;k<4;k++)
            {
            
                newmolecule.hbond_list[k].M1+=6*i;
                newmolecule.hbond_list[k].M2+=6*i;
            }
            ring.push_back(newmolecule);
            
        }
            
    }
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            int index=6*i+j;
            
            for(int k=0;k<4;k++)
            {
            
                
                if(ring[index].hbond_list[k].M2<0)
                {
                    ring[index].nbonds-=1;
                    ring[index].hbond_list[k]=ring[index].hbond_list.back();
                    ring[index].hbond_list.pop_back();
                }
                
                if(ring[index].hbond_list[k].M2>35)
                {
                    ring[index].nbonds-=1;
                    ring[index].hbond_list[k]=ring[index].hbond_list.back();
                    ring[index].hbond_list.pop_back();
                }
            }
            
            
        }
            
    }
    //Push molecules to M
    for(int i=0;i<36;i++)
    {
        M.push_back(ring[i]);
    }
    Aggregate A;
    A.n=6;
    A.rg2=1.1;
    for(int i=0;i<36;i++)
    {   
        A.M_A.push_back(i);
        cout<<M[i].LatticeID<<endl;
    }
    Ag.push_back(A);
    
    
    //Update lattice
    for(int i=0;i<36;i++)
    {
        lattice.molid_list[ring[i].LatticeID]=i;
    }
    //Update hbondlist
    for(int i=0;i<36;i++)
    {
        Molecule molecule=ring[i];
        for(int j=0;j<molecule.nbonds;j++)
        {
            hbond newbond=molecule.hbond_list[j];
            if(newbond.M1<newbond.M2)
                H.push_back(newbond);
        }
    }
    cout<<H.size()<<endl;
    for(int i=0; i<NMOL-36; i++)
    {
        Molecule m;
        m.MOL_ID=i+36;
        
        m.AID=i+1;//at first each molecule is one aggregate
        m.AsubID=0;//every molecule is the only element of the aggregate
        int trial_lattice;
        while (-1)
        {
            trial_lattice=gsl_rng_uniform_int(gsl_r,Nlattice);
            if(lattice.molid_list[trial_lattice]==-1)
            {
                break;
            }
        }
        
        m.LatticeID=trial_lattice;
        m.Lattice4ID=Lattice_index2_4index(m.LatticeID,Ng,Nbasis);
        m.hbond_list.clear();
        m.nbonds=0;
        //initialize vertypes of arms
        for(int i=0;i<m.N_VER/2;i++)
        {
            m.vertype[i]=gsl_rng_uniform_int(gsl_r,4);//random assign one of 4 types
            m.vertype[neighborarm(i)]=(2*gsl_rng_uniform_int(gsl_r,2)+m.vertype[i])%4;//neighbor arm must be the same or differ by pi
        }
        M.push_back(m);
        Aggregate A;
        A.n=1;
        A.rg2=1.1;
        A.M_A.push_back(m.MOL_ID);
        
        Ag.push_back(A);
        //Update lattice
        lattice.molid_list[m.LatticeID]=m.MOL_ID;
    }
    NAg=Ag.size();
    cout<<NAg<<endl;
    
    
}
/*
void System::WriteMol2(int timestep)
{
    ofstream out;
    out.open("conf.mol2",ios::app);
    out<<timestep<<endl;
    int n_ver=M[0].N_VER;
    //define filename from variable
    //char FileName[100];

    //sprintf(FileName,"_NG_%d_l0_%lf_lr_%d_t_%lf_f_%lf_e_%lf_C_%lf.mol2",NGRID,lambda0,lambdar,theta, tau_off, eta_off, cofactor);
    out<<"@<TRIPOS>MOLECULE"<<endl;
    out<<"Nanorod"<<endl;
    out<<NMOL*(n_ver+1)<<"\t"<<NMOL*n_ver+H.size()<<endl;
    out<<"SMALL"<<endl;
    out<<"NO_CHARGES"<<endl;

    out<<"@<TRIPOS>ATOM"<<endl;

    string name,type;
    
    int count=0;
    
    XYZ im,shift;
    
    for(int i=0; i<NMOL; i++)
    {
        
        out<<setw(6)<<++count<<"\t"<<"1"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"N.3"<<endl;
        
        for(int j=0; j<n_ver; j++)
        {
        
            out<<setw(6)<<++count<<"\t"<<"2"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"C.1"<<endl;
        }
    }
          
    out<<"@<TRIPOS>BOND"<<endl;
    
    count=0;
    for(int i=0; i<NMOL; i++)
    {
        for(int j=0; j<n_ver; j++)
        {
            out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*i+1<<"\t"<<setw(8)<<(n_ver+1)*i+j+2<<"\t"<<setw(2)<<"1"<<endl;
        }
    }
    list<hbond>::iterator it;
    for(it=H.begin();it!=H.end();it++)
    {
        hbond writebond=*it;
        out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*writebond.M1+writebond.arm1+2<<"\t"<<setw(8)<<(n_ver+1)*writebond.M2+writebond.arm2+2<<"\t"<<setw(2)<<"1"<<endl;
            
        
        
    }
    out<<"@<TRIPOS>FF_PBC"<<endl;
    out<<setw(8)<<"v1.0"<<"\t"<<setw(8)<<"1"<<"\t"<<setw(8)<<-L/2<<"\t"<<setw(8)<<-L/2<<"\t"<<setw(8)<<-L/2<<"\t"<<setw(8)<<L/2<<"\t"<<setw(8)<<L/2<<"\t"<<setw(8)<<L/2<<endl;   
    

    out.close();
    //   cout<<"Mol2 input written in\t"<<FileName<<endl;
    return;
}
*/
void System::WriteDump(int timestep)
{
    char FileName[200];
    sprintf(FileName,"N%d_L%d_E_dis%.1f_S_%.4f_T_%.3f_B%.3f_Dump.lammpstrj",NMOL,Ng,E_1,K_swap,total_time,K_Break);
    ofstream out;
    out.open(FileName,ios::app);
    
    out<<"ITEM: TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"ITEM: NUMBER OF ATOMS"<<endl;
    int NT=NMOL;
    out<<NT<<endl;
    //Calculate parameters for the hexagonal box

    out<<"ITEM: BOX BOUNDS xy xz yz"<<endl;
    out<<-0.5*Ng*LatticeVectors[0].x<<"\t"<<0.5*Ng*LatticeVectors[0].x+Ng*LatticeVectors[1].x<<"\t"<<Ng*LatticeVectors[1].x<<endl;
    out<<-0.5*Ng*LatticeVectors[1].y<<"\t"<<0.5*Ng*LatticeVectors[1].y<<"\t"<<0<<endl;
    out<<-0.5*Ng*LatticeVectors[2].z<<"\t"<<0.5*Ng*LatticeVectors[2].z<<"\t"<<0<<endl;
    out<<"ITEM: ATOMS index type x y z"<<endl;
    /*int NT=NMOL*7;
    out<<NT<<endl;
    out<<"Time="<<timestep<<endl;*/

    for(int i=0;i<NMOL;i++)
    {   
        XYZ position=Lattice4ID2XYZ(M[i].Lattice4ID,BasisPoints,LatticeVectors);

        out<<setw(6)<<i<<"\t"<<1<<"\t"<<setw(8)<<position.x<<"\t"<<setw(8)<<position.y<<"\t"<<setw(8)<<position.z<<endl;

    }
    out.close();
}
void System::WriteBond(int timestep)
{
    ofstream out;
    char FileName[200];
    sprintf(FileName,"N%d_L%d_E_dis%.1f_S_%.4f_T_%.3f_B%.3f_Bonlist.txt",NMOL,Ng,E_1,K_swap,total_time,K_Break);
    out.open(FileName,ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<setw(12)<<"molecule1"<<"\t"<<setw(12)<<"molecule2"<<"\t"<<setw(12)<<"arm1"<<"\t"<<setw(8)<<"arm2"<<endl;
    list<hbond>::iterator it;
    for(it=H.begin();it!=H.end();it++)
    {
        hbond writebond=*it;
        out<<setw(12)<<writebond.M1<<setw(12)<<writebond.M2<<setw(12)<<writebond.arm1<<setw(12)<<writebond.arm2<<endl;
            
        
        
    }
    out.close();
}


void System::WriteAggregate(int timestep)
{
    ofstream out;
    char FileName[200];
    sprintf(FileName,"N%d_L%d_E_dis%.1f_S_%.4f_T_%.3f_B%.3f_Aggregate.txt",NMOL,Ng,E_1,K_swap,total_time,K_Break);
    out.open(FileName,ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<setw(12)<<"AID"<<setw(12)<<"moleculeid"<<endl;
    for(int i=0;i<Ag.size();i++)
    {
        Aggregate ag=Ag[i];
        for(int j=0;j<ag.n;j++)
        out<<setw(12)<<i<<setw(12)<<ag.M_A[j]<<endl;
            
        
        
    }
    out.close();
}
void System::WriteLattice(int timestep)
{
    ofstream out;
    char FileName[200];
    sprintf(FileName,"N%d_L%d_E_dis%.1f_S_%.4f_T_%.3f_B%.3f_Lattice.txt",NMOL,Ng,E_1,K_swap,total_time,K_Break);
    out.open(FileName,ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<setw(12)<<"LatticeID"<<setw(12)<<"moleculeid"<<endl;
    for(int i=0;i<Nlattice;i++)
    {
        
        if(lattice.molid_list[i]!=-1)
        {
            out<<setw(12)<<i<<setw(12)<<lattice.molid_list[i]<<endl;
        }

            
        
        
    }
    out.close();
}


