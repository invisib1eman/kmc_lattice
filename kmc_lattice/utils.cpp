//NANOROD: utils.cpp Utilities Functions (Revision Date: Oct 27, 2023)
#include "utils.h"
//do not assume where the particle is (particle may have crossed the box twice or more)
vector<int> Lattice_index2_4index(int LatticeID,int Ng,int Nbasis)
{
  vector<int> Lattice4ID(4);
  Lattice4ID[3]=LatticeID%Nbasis;
  Lattice4ID[0]=(LatticeID/Nbasis)%Ng;
  Lattice4ID[1]=(LatticeID/Nbasis/Ng)%Ng;
  Lattice4ID[2]=(LatticeID/Nbasis/Ng/Ng)%Ng;
  return Lattice4ID;  
}
//ID=Nbasis*(IDx+Ng*IDy+Ng*Ng*IDz)+IDbasis
int Lattice_4index2_index(vector<int> Lattice4ID,int Ng,int Nbasis)
{
  return Nbasis*(Lattice4ID[0]+Ng*Lattice4ID[1]+Ng*Ng*Lattice4ID[2])+Lattice4ID[3];
}
XYZ Lattice4ID2XYZ(vector<int> Lattice4ID,vector<XYZ> BasisPoints,vector<XYZ> LatticeVectors)
{
    return XYZ(LatticeVectors[0].x*Lattice4ID[0]+LatticeVectors[1].x*Lattice4ID[1]+LatticeVectors[2].x*Lattice4ID[2]+BasisPoints[Lattice4ID[3]].x,LatticeVectors[0].y*Lattice4ID[0]+LatticeVectors[1].y*Lattice4ID[1]+LatticeVectors[2].y*Lattice4ID[2]+BasisPoints[Lattice4ID[3]].y,LatticeVectors[0].z*Lattice4ID[0]+LatticeVectors[1].z*Lattice4ID[1]+LatticeVectors[2].z*Lattice4ID[2]+BasisPoints[Lattice4ID[3]].z);
}
vector<int> Translate_Lattice4ID(vector<int> Lattice4ID,vector<int> Lattice4ID_translate,int Ng)
{
  vector<int> new_Lattice4ID(4,0);
  std::transform (Lattice4ID.begin(), Lattice4ID.end(), Lattice4ID_translate.begin(), new_Lattice4ID.begin(), std::plus<int>());
  for(int i=0;i<3;i++)
  {
    if(new_Lattice4ID[i]<0)
      new_Lattice4ID[i]+=Ng;
    if(new_Lattice4ID[i]>=Ng)
      new_Lattice4ID[i]-=Ng;
  }
  return new_Lattice4ID;
}
int bondarm(int arm)
{
  return arm;//match in our definition
}
//return armindex of neighbor arm of arm n
int neighborarm(int n)
{
  return n+1-2*(n%2);
}
void DFSUtil(int index_ag,int v,bool visited[],vector<Molecule>& M,Aggregate& pnew_aggregate)
{
  visited[v]=true;
  M[v].AID=index_ag;
  (pnew_aggregate).M_A.push_back(M[v].MOL_ID);
  (pnew_aggregate).n+=1;
  vector<hbond>::iterator i;
  for(i=M[v].hbond_list.begin();i!=M[v].hbond_list.end();++i)
  {
    hbond j=*i;
    int next=j.M2;
    if(!visited[next])
    {
      DFSUtil(index_ag,next,visited,M,pnew_aggregate);
    }
  }
  
}
int newvertype(int oldvertype,int changevertype)
{
  int newvertype=oldvertype+changevertype;
  if(newvertype>=4)
    newvertype-=4;
  if(newvertype<0)
   newvertype+=4;
  return newvertype;
}
bool bondstate(int armstate1,int armstate2)
{
  if(armstate1+armstate2==3)
    return true;
  else
    return false;
}
