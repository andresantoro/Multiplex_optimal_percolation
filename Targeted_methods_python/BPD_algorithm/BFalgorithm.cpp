  /*
********************************************************************************
TAbyFVSbpdV03.cpp

STATEMENT:
This code is aimed at approximately solving an important problem of network
science, namely to divide a network into small components by deleting vertices
in a most economic way.
The author (Hai-Jun Zhou) of this code is against any form of terrorism. This
code should only be used for pure academic purposes and for practical purposes
that improve the welfare of humanity. The author strongly condemn any intention
of adapting this code for destructive purposes.

DESCRIPTION:
Belief propagation guided decimation (BPD) as a solver for the network targeted
attack problem. The goal is to construct a minum set S such that, if all the
vertices in this set and all the attached edges are deleted from the graph,
the remaining subgraph will be formed by disconnected components of size not
exceeding a given value Sthreshold.

This program is applicable on a single graph instance. The imput graph file has
the following form:

N   M                % first row specifies N (vertex number) and M (edge number)

i_1  j_1                                            % undirected edge (i_1, j_1)
i_2  j_2                                            % undirected edge (i_2, j_2)
.    .
.    .
.    .
i_M  j_M                                            % undirected edge (i_M, j_M)

The program reads only the first EdgeNumber edges from the imput graph, with
EdgeNumber being explicitly specified (EdgeNumber <= M, of course).

Each vertex has two states (unoccupied, 0; occupied, 1). An unoccupied vertex
belongs to the constructed target set S, while an occupied vertex remains in the
network.

The target set S is contructed by three steps:
1) a feedback vertex set S0 is constructed by BPD, all the vertices in S0 are
deleted from the input graph G.
2) if a tree component of the remaining subgraph contains more than Sthreshold
vertices, one of its vertices is deleted.
3) some deleted vertices are added back to the remaining network as long as
the netowrk component sizes do not exceed Sthreshold.


For more details about the algorithm and the replica-symmetric mean field
theory and the spin glass model on which this algorithm was based, please
consult the following references:
[1] Salomon Mugisha, Hai-Jun Zhou, "Identifying optimal targets of network
attack by belief propagation", arXiv:1603.05781 (2016).
[2] Hai-Jun Zhou, "Spin glass approach to the feedback vertex set problem",
European Physical Journal B 86: 455 (2013).
[3] Hai-Jun Zhou, "Spin Glasses and Message Passing" (Science Press, Beijing,
2015), chapter 7, pages 218--238.

To generate the executive file, you can simply compile as
* c++ -O3 -o tabpd.exe TAbyFVSbpdV03.cpp
!!! please make sure some of the key parameters, such as input graph name,
number of edges, number of vertices, output file names, are appropriately
specified in the TAbyFVSbpdV03.cpp file.

After you successfully compiled the code, you can then run the algorithm as
* tabpd.exe

Good luck!

LOG:
28.03.2016: TAbyFVSbpdV03.cpp (publicly accessible version).
28.03.2016: copied TAbyFVSbpdV02.cpp to TAbyFVSbpdV03.cpp.

PROGRAMMER:
Hai-Jun Zhou
Institute of Theoretical Physics, Chinese Academy of Sciences
Zhong-Guan-Cun East Road 55, Beijing 100190
email: zhouhj@itp.ac.cn
webpage: power.itp.ac.cn/~zhouhj/
********************************************************************************
*/

#include <exception>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <valarray>
#include <time.h>
#include <sstream>


#include "zhjrandom.h"                             //a random number generator

using namespace std;

/*---               random real uniformly distributed in [0,1)            ---*/
double u01prn(ZHJRANDOMv3 *rd)
{
  return rd->rdflt();
}


struct IntPair
{
  int first;
  int second;
  IntPair(void);
  IntPair(int, int);
};

IntPair::IntPair(void)
{
  first=0;
  second=0;
  return;
}

IntPair::IntPair(int a, int b)
{
  first  = a;
  second = b;
  return ;
}

bool operator<(IntPair a, IntPair b)
{
  if(a.first < b.first)
    return true;
  else if(a.first > b.first)
    return false;
  else                                                      //a.first = b.first
    {
      if(a.second < b.second)
	return true;
      else
	return false;
    }
}

bool operator==(IntPair a, IntPair b)
{
  return (a.first == b.first) && (a.second == b.second);
}


struct message                                   //cavity message from a vertex
{
  struct vstruct *v_ptr;                        //pointer to the sending vertex
  double q_0;                                 //probability of taking state A=0
  double q_root;                                //probability of being the root
  message(void);
  message(double, double);
};

message::message(void)
{
  v_ptr  = 0;
  q_0    = 0;
  q_root = 0;
  return;
}

message::message(double a, double b)
{
  v_ptr  = 0;
  q_0    = a;
  q_root = b;
  return;
}


struct mpointer                                            //pointer to message
{
  struct message *m_ptr;
  mpointer(void);
};

mpointer::mpointer(void)
{
  m_ptr=0;
  return;
}


struct vstruct                                                  //vertex struct
{
  int index;                                //index of vertex, positive integer
  int degree;                                  //number of neighboring vertices
  int active_degree;                    //number of active neighboring vertices
  int c_index;                                             //index of component
  bool occupied;                        //=true (not deleted); =false (deleted)
  bool active;                  //=true (need to be considered in BP iteration)
  double q_0;                                               //empty probability
  struct message *im_ptr;                     //start position of input message
  struct mpointer *omptr_ptr;   //start position of output message address list
  vstruct(void);
};

vstruct::vstruct(void)
{
  index=0;
  degree=0;
  active_degree=0;
  c_index=0;
  occupied=true;                          //initially all vertices are occupied
  active=true;                              //initially all vertices are active
  q_0=0;
  im_ptr=0;
  omptr_ptr=0;
  return;
}


class FVS                                                 //feedback vertex set
{
public:
  FVS(ZHJRANDOMv3* );                                             //constructor
  void SetX(double);                             //set re-weighting parameter X
  void SetDampingFactor(double);                           //set damping factor
  bool Graph(string&, int,int);                     //read graph connection pattern

  void DegreeRank(string&);                       //rank plot of active degrees
  void DeleteMaxVertex(int);            //delete most highly connected vertices
  void ComponentRefinement(int);                     //refinement of components
  bool CheckComponents(string&, int);        //check the size of each component
  void AttackEffect(string&);         //the accumulated effect of node deletion

  bool ReadFVS(string&);                          //read a FVS into the program
  void MaxComponentSize(string&);      //max component size versus attack scale
  void Initialization(void);                            //initialize population
  bool CheckFVS(string&);           //check whether occupation pattern is a FVS
  bool BeliefPropagation(double, int);                    //population updating
  void EmptyProbability(void);       //empty probability for each active vertex
  void Thermodynamics(string& );                     //thermodynamic quantities
  int  Fix0(void);          //fix some variables to be un-occupied and simplify
  void Simplify(struct vstruct* );        //simplify the graph by removing leaf
  int return_max_degree(void);

private:
  int VertexNumber;                                  //total number of vertices
  int ActiveVertexNumber;                     //total number of active vertices
  int ActiveVertexNumber0;   //total # of active vertices before this BPD round

  int EdgeNumber;                                       //total number of edges
  int MaxDegree;                           //maximal vertex degree in the graph

  int MaxFixNumber;             //maximal number of fixed vertices in one round
  int MinRange;                 //minimal range of empty_prob of fixed vertices
  float FixFraction;           /* fraction of active vertices to be fixed to be
				  un-occupied in each round of decimation */
  double X;                                            //re-weighting parameter
  double Weight0;                                                    //=exp(-X)
  double DampingFactor;                    //damping factor of message updating

  ZHJRANDOMv3 *PRNG;                                  //random number generator

  stack<int> Targets;                                     //target attack nodes

  valarray<vstruct> Vertex;
  valarray<message> InMessage;
  valarray<mpointer> OutMessageAddress;

  valarray<int> CandidateVertices;     //list of candidate vertices to be fixed
  valarray<int> CandidateSize;     //num of candidates in each empty_prob range
  valarray<int> Permutation;         //array used in random sequential updating
  valarray<double> WeightA;       //auxiliary array needed for message updating
  valarray<double> WeightB;       //auxiliary array needed for message updating
  valarray<double> WeightC;       //auxiliary array needed for message updating
  valarray<bool>   IsActive;      //auxiliary array needed for message updating

  void UpdateMessage(struct vstruct *, double&);
};

int main(int argc, char ** argv)
{
  //                                     random number generator initialization
  int rdseed=91276791;                 //you can set this seed to another value
  ZHJRANDOMv3 rdgenerator(rdseed);
  int prerun=12000000;                        //you can set it to another value
  for(int i=0; i<prerun; ++i)
    rdgenerator.rdflt();

  FVS system( &rdgenerator);

  // read graph, please check these parameters according to your graph instance
  string graphname= argv[1];                     //graph file's name
  int VertexNumber = atoi(argv[2]);   //number of vertices in the graph
  int EdgeNumber   = atoi(argv[3]);   //number of edges you want to read from the file
  int job_number = atoi(argv[4]);
  int iteration = atoi(argv[5]);
  char str1[100];
  //cout << "name of the graph: " << graphname << " -- number of nodes: " << VertexNumber << " -- Edge number: " << EdgeNumber << endl;
  bool succeed=system.Graph(graphname, VertexNumber,EdgeNumber);
  if(succeed==false)
    return -1;

  /* If you already have a FVS solution, you don't need to run BPD. Please
     use the following codes instead. */
  /*
    string FVSfile0="ERn100kM1mg100.FVS";
    succeed=system.ReadFVS(FVSfile0);
    if(succeed==false)
    return -1;
    int Csize0;
    do
    {
    cout<<"Maximal allowed size of component =";
    cout.flush();
    cin>>Csize0;
    } while(Csize0 <= 2);
    system.ComponentRefinement(Csize0);
    string Afile0="ERn100kM1mg100.TAset";
    string Bfile0="ERn100kM1mg100.TAeffect";
    if( system.CheckComponents(Afile0, Csize0) == true)
    {
    cout<<"An set of targeted attack nodes constructed.\n";
    system.AttackEffect(Bfile0);
    return 1;
    }
    else
    {
    cerr<<"Some components still too big. Problem in program.\n";
    return -1;
    }
  */

  /* In case you don't have a FVS solution at hand */
  string degreefile="ERn100kM1mg100.degree";                      //degree rank
  sprintf(str1, "_job%d_iteration%d", job_number,iteration);
  degreefile.append(str1);
  system.DegreeRank(degreefile);
  int Degree0;   /* Delete vertices whose active degree is larger than Degree0.
		    You should not set Degree0 to be too small.
		    If you don't want to delete vertices at this stage (for
		    example, if the maximal vertex degree is less than 1000),
		    you can just set Degree0 to be very large.
		    If the maximal degree of the network is very large, you can
		    set Degree0 to be a value that approximately 0.5 percent of
		    the vertices have degree beyond Degree0. */
  // do {
  //   cout<<"Please check "<<degreefile<<" file and choose Degree threshold=";
  //   cout.flush();
  //   cin>>Degree0;
  // } while(Degree0 <= 0);
  Degree0 = system.return_max_degree();
  system.DeleteMaxVertex(Degree0);
  const time_t t1=time(NULL);
  //                            initialize cavity messages in a most random way
  system.Initialization();
  //                                                         set BPD parameters
  double X=12.0e0;
  system.SetX(X);
  double DampingFactor=0.9e0;             //you don't need to change this value
  system.SetDampingFactor(DampingFactor);
  //                                                      initial BP iterations
  float epsilon = 1.0e-10;
  int eq_iterations = 100;      //please don't make it too small (i.e., >= 100)
  system.BeliefPropagation(epsilon, eq_iterations);
  string thermfile="ERn100kM1mg100.RS";
  thermfile.append(str1);
  system.Thermodynamics(thermfile);                       //compute free energy
  //                                                                        BPD
  epsilon = 1.0e-7;
  int iterations = 20;  /* the final result not sensitive to this parameter, as
			   long as it is not too small, e.g., iterations>=10 */
  int ActiveVertexNumber=0;
  do {
    system.BeliefPropagation(epsilon, iterations);
    system.EmptyProbability();
    ActiveVertexNumber=system.Fix0();
  } while(ActiveVertexNumber>0);
  //           report feedback vertex set, please change the name as you prefer
  string FVSfile="ERn100kM1mg100.FVS";
  FVSfile.append(str1);
  if( system.CheckFVS(FVSfile) == false)
    {
      cerr<<"Not a feedback vertex set.\n";
      return -1;
    }
  //          You can set Csize to be 1 percent of VertexNumber of even smaller
  int Csize = int(sqrt(VertexNumber));                       //assumed to be 1 percent
  system.ComponentRefinement(Csize);
  string Afile="ERn100kM1mg100.TAset";
  Afile.append(str1);
  if( system.CheckComponents(Afile, Csize) == true)
    {
      /* The following two lines report size evolution of the largest component
	 string Bfile="ERn100kM1mg100.TAeffect";
	 system.AttackEffect(Bfile);
      */
      const time_t t2=time(NULL);
      //cout<<"Total time used = "<<t2-t1<< "s"<<endl;
      //cout.flush();
      return 1;
    }
  else
    {
      cerr<<"Some components still too big. Problem in program.\n";
      return -1;
    }
}


/*---                       constructor of FVS cluster                    ---*/
FVS::FVS(ZHJRANDOMv3* rd)
{
  PRNG = rd;                                          //random number generator
  FixFraction = 0.001;                    //fix 1 percent of the active vertices
  return;
}

/* -                           Read graph from file                           -
   gname: the input graph name.
   enumber: read the first enumber edges only.
*/
bool FVS::Graph(string& gname, int vertex_number,int enumber)
{
  ifstream graphf(gname.c_str());
  if( !graphf.good() )
    {
      cerr<<"Graph probably non-existant.\n";
      return false;
    }
  while( !Targets.empty() )
    {
      Targets.pop();
    }
  //first read of input graph
  // graphf >>VertexNumber
	//  >>EdgeNumber;
  VertexNumber = vertex_number;
  EdgeNumber = enumber;
  // cout << "HEREEE" << endl;
  // cerr<< "Vertex Number: " << VertexNumber <<  " -- enumber: " << enumber << endl;
  if(EdgeNumber<enumber)
    {
      cerr<<"No so many edges in the graph.\n";
      graphf.close();
      return false;
    }
  EdgeNumber=enumber;                 //only the first enum edges are read into
  try { Vertex.resize(VertexNumber+1); } catch(bad_alloc)
    {
      cerr<<"Vertex construction failed.\n";
      graphf.close();
      return false;
    }
  try { Permutation.resize(VertexNumber+1); } catch(bad_alloc)
    {
      cerr<<"Permutation construction failed.\n";
      return false;
    }
  bool succeed=true;
  set<IntPair> EdgeSet;
  // cout << "HEREEE" << endl;
  for(int eindex=0; eindex<EdgeNumber && succeed; ++eindex)
    {
      int v1, v2;
      graphf >>v1
	     >>v2;
      // cout << v1 <<" "  << v2 << endl;
      if(v1>v2)
	{
	  int v3=v1;
	  v1=v2;
	  v2=v3;
	}
      if(v1==v2 || v1==0 || v1 >VertexNumber || v2==0 || v2 >VertexNumber)
	{
    if(v1==v2)
      cerr<<"Graph incorrect at line "<<eindex+1 << "because" << v1 << " == " << v2 <<endl;  
    if(v1==0)
      cerr<<"Graph incorrect at line "<<eindex+1 << "because" << v1 << " == 0" <<endl;  
    if(v2==0)
      cerr<<"Graph incorrect at line "<<eindex+1 << "because" << v2 << " == 0" <<endl;  
    if(v1 >VertexNumber)
      cerr<<"Graph incorrect at line "<<eindex+1 << "because" << v1 << " > " << VertexNumber <<endl;  
    if(v2 >VertexNumber)
      cerr<<"Graph incorrect at line "<<eindex+1 << "because" << v2 << " > " << VertexNumber <<endl;  
	  cerr<<"Graph incorrect at line "<<eindex+1<<endl;
	  succeed=false;
	}
      else if(EdgeSet.find(IntPair(v1,v2) ) != EdgeSet.end() )
	{
	  cerr<<"Multiple edges.\n";
	  succeed=false;
	}
      else
	{
	  EdgeSet.insert(IntPair(v1,v2));
	  ++(Vertex[v1].degree);
	  ++(Vertex[v2].degree);
	}
    }
  graphf.close();
  if(succeed==false)
    return false;
  EdgeSet.clear();
  try { InMessage.resize( 2*EdgeNumber ); } catch(bad_alloc)
    {
      cerr<<"InMessage construction failed.\n";
      return false;
    }
  try { OutMessageAddress.resize( 2*EdgeNumber ); } catch(bad_alloc)
    {
      cerr<<"OutMessageAddress construction failed.\n";
      return false;
    }
  int position=0;
  MaxDegree=0;
  struct vstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    {
      v_ptr->index =v;
      v_ptr->occupied=true;
      v_ptr->active =true;
      v_ptr->active_degree = v_ptr->degree;
      v_ptr->im_ptr = &InMessage[position];
      v_ptr->omptr_ptr = &OutMessageAddress[position];
      position += v_ptr->degree;
      if(v_ptr->degree>MaxDegree)
	MaxDegree = v_ptr->degree;
      v_ptr->degree=0;
    }
  //second read of input graph
  graphf.open(gname.c_str());
  // graphf >>VertexNumber
	//  >>EdgeNumber;
  // EdgeNumber=enumber;
  struct mpointer *omptr_ptr=&OutMessageAddress[0];
  for(int eindex=0; eindex<EdgeNumber; ++eindex)
    {
      int v1,v2;
      graphf >>v1
	     >>v2;
      if(v1>v2)
	{
	  int v3=v1;
	  v1=v2;
	  v2=v3;
	}
      omptr_ptr=Vertex[v1].omptr_ptr + Vertex[v1].degree;
      omptr_ptr->m_ptr = Vertex[v2].im_ptr + Vertex[v2].degree;
      omptr_ptr->m_ptr->v_ptr = &Vertex[v1];
      omptr_ptr=Vertex[v2].omptr_ptr + Vertex[v2].degree;
      omptr_ptr->m_ptr = Vertex[v1].im_ptr + Vertex[v1].degree;
      omptr_ptr->m_ptr->v_ptr = &Vertex[v2];
      ++(Vertex[v1].degree);
      ++(Vertex[v2].degree);
    }
  graphf.close();
  // cout << "Graph: N= "<<VertexNumber
  //      <<",  M= "<< EdgeNumber
  //      <<",  Mean degree= "<<(2.0*EdgeNumber)/(1.0*VertexNumber)<<endl;
  //cout.flush();
  try { WeightA.resize(MaxDegree+1); } catch(bad_alloc)
    {
      cerr<<"WeightA construction failed.\n";
      return false;
    }
  try { WeightB.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"WeightB construction failed.\n";
      return false;
    }
  try { WeightC.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"WeightC construction failed.\n";
      return false;
    }
  try { IsActive.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"IsActive construction failed.\n";
      return false;
    }
  ActiveVertexNumber=VertexNumber;
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->active && v_ptr->active_degree <=1)
      {
	v_ptr->active = false;
	v_ptr->occupied = true;                                //being occupied
	--ActiveVertexNumber;
	Simplify(v_ptr);
      }
  if(ActiveVertexNumber==0)
    {
      /*cout<<"The graph has no loops. Done! \n";
      cout.flush();*/
    }
  v_ptr=&Vertex[1];
  ActiveVertexNumber=0;
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->active)
      {
	Permutation[ActiveVertexNumber]=v_ptr->index;
	++ActiveVertexNumber;
      }
  ActiveVertexNumber0=ActiveVertexNumber;
  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0)
    MaxFixNumber = 1;
  //the empty probability q_0 is divided into bins of width 0.01
  try { CandidateVertices.resize( MaxFixNumber * 101 ); } catch( bad_alloc )
    {
      cerr<<"CandidateVertices construction failed.\n";
      return false;
    }
  CandidateSize.resize(101);
  MinRange=0;
  return true;
}

/*                simplify after fixing variable                             */
void FVS::Simplify(struct vstruct *v_ptr)
{
  struct message *im_ptr=v_ptr->im_ptr;
  for(int d=0; d < v_ptr->degree; ++d, ++im_ptr)
    {
      if( --(im_ptr->v_ptr->active_degree) <= 1)
	{
	  if(im_ptr->v_ptr->active)
	    {
	      im_ptr->v_ptr->active   = false;
	      im_ptr->v_ptr->occupied = true;
	      --ActiveVertexNumber;
	      Simplify(im_ptr->v_ptr);
	    }
	}
    }
  return;
}

int FVS::return_max_degree(void)
{
  return MaxDegree;
}


/* Rank plot about the active degrees of the active vertices */
void FVS::DegreeRank(string& rankfile)
{
  //    Notice!  In this subroutine WeightA is not used for its default purpose
  WeightA=0;
  struct vstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->active)
      WeightA[v_ptr->active_degree] += 1.0;
  ofstream output(rankfile.c_str());
  double myrank=1;
  int size0p01=VertexNumber/100;
  int degree0p01;      //lowest degree of the 1 percent highest degree vertices
  for(int d=MaxDegree; d>=0; --d)
    if(WeightA[d]>0)
      {
	output<<myrank<<'\t'<<d<<endl;                         //rank of degree
	myrank += WeightA[d];
	if(myrank<=size0p01)
	  degree0p01=d;
	if(WeightA[d]>1)
	  output<<(myrank-1)<<'\t'<<d<<endl;
      }
  output.close();
  // cout<<"Maximal degree ="<<MaxDegree<<"   "
  //     <<"Degree of top one-percent="<<degree0p01<<endl;
  // cout.flush();
  return;
}

/* Delete the most highly connected active vertices */
void FVS::DeleteMaxVertex(int Dthreshold)
{
  //                                       Dthreshold is the degree upper bound
  set<IntPair> candidates;
  struct vstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->active && v_ptr->active_degree > Dthreshold)
      candidates.insert(IntPair(- v_ptr->active_degree, v_ptr->index) );

  int Dnumber=0;                        //number of externally deleted vertices
  int vtx;
  while(!candidates.empty() )
    {
      bool succeed=false;
      do
	{
	  set<IntPair>::const_iterator sci = candidates.begin();
	  vtx  = sci->second;
	  int Dmax = - sci->first;
	  candidates.erase(IntPair(-Dmax, vtx) );
	  if(Vertex[vtx].active)
	    {
	      succeed=true;
	      if(Dmax != Vertex[vtx].active_degree)
		{
		  succeed=false;
		  if(Vertex[vtx].active_degree>Dthreshold)
		    candidates.insert(IntPair(-Vertex[vtx].active_degree,vtx));
		}
	    }
	}
      while(succeed==false && !candidates.empty() );
      if(succeed)
	{
	  v_ptr = &Vertex[vtx];
	  v_ptr->active=false;
	  v_ptr->occupied=false;                                      //deleted
    cout << "Removed node: " << v_ptr->index <<endl;
	  Targets.push(v_ptr->index);
	  ++Dnumber;
	  --ActiveVertexNumber;
	  Simplify(v_ptr);
	}
    }
  v_ptr=&Vertex[1];
  ActiveVertexNumber=0;
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->active)
      {
	Permutation[ActiveVertexNumber]=v_ptr->index;
	++ActiveVertexNumber;
      }
  ActiveVertexNumber0=ActiveVertexNumber;
  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0)
    MaxFixNumber = 1;
  //cout<<Dnumber<<" highestly connected vertices are externally deleted.\n";
  //cout.flush();
  return;
}

/*---                set the value of X                                   ---*/
void FVS::SetX(double xval)
{
  X = xval;                                            //re-weighting parameter
  Weight0=exp(-X);
  return;
}

/*---                set the value of DampingFactor                       ---*/
void FVS::SetDampingFactor(double dp)
{
  DampingFactor = dp;                                          //damping factor
  if(DampingFactor >= 1.0e0)
    DampingFactor=1.0e0;
  else if(DampingFactor<=0.0001e0)
    DampingFactor=0.0001e0;
  return;
}

/*---                      initialize population                          ---*/
void FVS::Initialization(void)
{
  struct mpointer *omptr_ptr = &OutMessageAddress[0];
  for(int v=1; v<=VertexNumber; ++v)
    {
      int degree=Vertex[v].degree;
      for(int d=0; d<degree; ++d)
	{
	  omptr_ptr->m_ptr->q_0 = 1.0e0/(1.0e0*degree+2.0e0);
	  omptr_ptr->m_ptr->q_root = 1.0e0/(1.0e0*degree+2.0e0);
	  ++omptr_ptr;
	}                                                //uniform distribution
    }
  return;
}

/* -                Belief propagation                                     - */
bool FVS::BeliefPropagation(double error, int count)
{
  if(ActiveVertexNumber==0)                                //no active vertices
    return true;
  int iter=0;
  double max_error=0;
  for(int quant=ActiveVertexNumber0-1; quant>=0; --quant)
    if(Vertex[ Permutation[quant] ].active == false)
      {
	--ActiveVertexNumber0;
	Permutation[quant]=Permutation[ActiveVertexNumber0];
      }
  if(ActiveVertexNumber0 != ActiveVertexNumber)
    {
      cerr<<"Unexpected mistake.\n";
      return false;
    }
  do
    {
      max_error=0;
      for(int quant=ActiveVertexNumber; quant>0; --quant)
	{
	  int iii = static_cast<int>(quant * u01prn(PRNG) );
	  int v = Permutation[iii];
	  Permutation[iii] = Permutation[quant-1];
	  Permutation[quant-1] = v;
	  UpdateMessage(&Vertex[v], max_error);
	}
      // if((iter % 10) == 0)
	// {
	//   cerr<<' '<<max_error<<' ';
	//   cerr.flush();
	// }
      ++iter;
    }
  while(max_error>error && iter<count);
  if(max_error<=error)
    {
      // cerr<<' '<<max_error<<"  :-)\n";
      return true;
    }
  else
    {
      // cerr<<' '<<max_error<<"  :-(\n";
      return false;
    }
}

/* -         Update the output messages from a vertex                         -
   It is assumed that central vertex is active and has at least two active
   neighbors, namely v_ptr->active = true,  v_ptr->active_degree >= 2        */
void FVS::UpdateMessage(struct vstruct *v_ptr, double& maxdiff)
{
  int degree=v_ptr->degree;
  struct message *im_ptr = v_ptr->im_ptr;
  for(int j=0; j<degree; ++j, ++im_ptr)
    {
      if(im_ptr->v_ptr->active)
	{
	  IsActive[j]=true;
	  WeightA[j]=1.0e0; /* probability all neighbors of central vertex i
			       (except vertex j) are empty or are root in
			       cavity graph */
	  WeightB[j]=0.0e0; /* probability that one neighbor (k not equal to j)
			       of central vertex i are not empty in the cavity
			       graph, all the other neighbors (except j) of
			       central vertex i are either empty or root in
			       the cavity graph */
	  WeightC[j]=Weight0;
	}
      else
	IsActive[j]=false;
    }
  im_ptr = v_ptr->im_ptr;
  for(int j=0; j<degree; ++j, ++im_ptr)
    if(IsActive[j])
      {
	double q_0     = im_ptr->q_0;
	double q_root  = im_ptr->q_root;
	double q_0root = q_0+q_root;
	for(int k=0; k<degree; ++k)
	  if( (k != j) && IsActive[k])
	    {
	      WeightB[k]  = WeightA[k]*(1.0e0-q_0) + WeightB[k]*q_0root;
	      WeightA[k] *= q_0root;
	      double maxval=max(WeightA[k], WeightB[k]);
	      if(maxval<WeightC[k])
		maxval=WeightC[k];
	      WeightA[k] /= maxval;                          //avoid underflow
	      WeightB[k] /= maxval;                          //avoid underflow
	      WeightC[k] /= maxval;                          //avoid underflow
	    }
      }
  struct mpointer *omptr_ptr = v_ptr->omptr_ptr;
  for(int j=0; j<degree; ++j, ++omptr_ptr)
    if(IsActive[j])
      {
	double norm=WeightA[j]+WeightB[j]+WeightC[j];
	double diff_q_0    = WeightC[j]/norm - omptr_ptr->m_ptr->q_0 ;
	double diff_q_root = WeightA[j]/norm - omptr_ptr->m_ptr->q_root;
	double diff = abs(diff_q_0)+abs(diff_q_root);
	if(diff>maxdiff)
	  maxdiff=diff;
	omptr_ptr->m_ptr->q_0    += DampingFactor * diff_q_0;
	omptr_ptr->m_ptr->q_root += DampingFactor * diff_q_root;
      }
  return;
}

/*---                       Empty probability of each vertex              ---*/
void FVS::EmptyProbability(void)
{
  if(ActiveVertexNumber == 0)
    return ;
  double q_0, q_root, q_0root, norm, weight_a, weight_b, weight_c;
  struct vstruct *v_ptr=&Vertex[1];
  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0)
    MaxFixNumber = 1;
  for(int r=0; r<=100; ++r) CandidateSize[r] = 0;   /* number of candidate
						       vertices in each q_0 bin
						       of width 0.01 */
  MinRange=0;   /* the maximal value of r0 such that
		   sum_{r>=r0} CandidateSize[r] >= MaxFixNumber is satisfied */
  /*         thermodynamic quantities of edges                               */
  struct message *im_ptr = &InMessage[0];
  for(int v_index=0; v_index < ActiveVertexNumber; ++v_index)
    {
      v_ptr = &Vertex[ Permutation[v_index] ];
      int degree=v_ptr->degree;
      /*  the central vertex is assumed to have at least two or more active
	  neighbors. If it has zero or only one active neighbor, it should have
	  been removed from the active subgraph */
      weight_a=1.0e0; /*          probability all neighbors of central vertex i
				  are empty or are root in cavity graph */
      weight_b=0.0e0; /*     probability that one neighbor (j)of central vertex
			     i are not empty in the cavity graph, all the other
			     neighbors of central vertex i are either empty or
			     root in the cavity graph */
      weight_c=Weight0;
      im_ptr = v_ptr->im_ptr;
      for(int j=0; j<degree; ++j, ++im_ptr)
	if(im_ptr->v_ptr->active)         // the neighboring vertex is active
	  {
	    q_0     = im_ptr->q_0;
	    q_root  = im_ptr->q_root;
	    q_0root = q_0+q_root;
	    weight_b  = weight_a * (1.0e0-q_0)+weight_b * q_0root;
	    weight_a *= q_0root;
	    double maxval=max(weight_a, weight_b);
	    if(maxval<weight_c)
	      maxval=weight_c;
	    weight_a /= maxval;
	    weight_b /= maxval;
	    weight_c /= maxval;
	  }
      norm = weight_a+weight_b+weight_c;
      q_0  = weight_c/norm;
      v_ptr->q_0 = q_0;
      int rrr = static_cast<int>(q_0*100);
      if(rrr>=MinRange)
	{
	  if(CandidateSize[rrr]<MaxFixNumber)
	    {
	      CandidateVertices[rrr*MaxFixNumber + CandidateSize[rrr] ]
		= v_ptr->index;
    cout <<  v_ptr->index << " " <<  v_ptr->q_0  << endl;
	      ++CandidateSize[rrr];
	    }
	  else
	    MinRange=rrr;
	}
    }
  return;
}

/*---                       thermodynamic quantities                      ---*/
void FVS::Thermodynamics(string& ofile)
{
  if(ActiveVertexNumber == 0)
    return;
  int degree;
  double q_0, q_root, q_0out, q_rootout, q_0root, norm, weight_a, weight_b,
    weight_c, logval;
  double phi_vtx=0;                       //vertex contribution to free entropy
  double phi_edge=0;                        //edge contribution to free entropy
  double rho_vtx=0;                                 //vertex occupation density
  /*         thermodynamic quantities of edges                               */
  struct message *im_ptr = &InMessage[0];
  struct mpointer *omptr_ptr = &OutMessageAddress[0];
  for(int v_index=0; v_index < ActiveVertexNumber; ++v_index)
    {
      struct vstruct *v_ptr= &Vertex[ Permutation[v_index] ];
      int degree=v_ptr->degree;
      /*  the central vertex is assumed to have at least two or more active
	  neighbors. If it has zero or only one active neighbor, it should have
	  been removed from the active subgraph */
      weight_a=1.0e0; /*          probability all neighbors of central vertex i
				  are empty or are root in cavity graph */
      weight_b=0.0e0; /*     probability that one neighbor (j)of central vertex
			     i are not empty in the cavity graph, all the other
			     neighbors of central vertex i are either empty or
			     root in the cavity graph */
      weight_c=Weight0;
      logval=0;
      im_ptr = v_ptr->im_ptr;
      omptr_ptr = v_ptr->omptr_ptr;
      for(int j=0; j<degree; ++j)
	{
	  q_0out = omptr_ptr->m_ptr->q_0;
	  q_rootout=omptr_ptr->m_ptr->q_root;
	  if(im_ptr->v_ptr->active)          //the neighboring vertex is active
	    {
	      q_0     = im_ptr->q_0;
	      q_root  = im_ptr->q_root;
	      q_0root = q_0+q_root;
	      weight_b  = weight_a * (1.0e0-q_0)+weight_b * q_0root;
	      weight_a *= q_0root;
	      // the following lines added on 07.09.2015 for avoiding underflow
	      double maxval = max(weight_a, weight_b);
	      if(maxval<weight_c)
		maxval=weight_c;
	      weight_a /= maxval;
	      weight_b /= maxval;
	      weight_c /= maxval;
	      logval += log(maxval);
	      phi_edge += log(q_0 + q_0out * (1.0e0-q_0) + q_root *
			      (1.0e0 - q_0out) + q_rootout * (1.0e0 - q_0) );
	    }
	  ++im_ptr;
	  ++omptr_ptr;
	}
      norm = weight_a+weight_b+weight_c;
      q_0  = weight_c/norm;
      v_ptr->q_0 = q_0;
      rho_vtx += 1.0e0-q_0;
      phi_vtx += X+log(norm)+logval;
    }
  phi_vtx  /= ActiveVertexNumber;
  phi_edge /= 2.0e0*ActiveVertexNumber;
  rho_vtx  /= ActiveVertexNumber;
  double phi = (phi_vtx - phi_edge)/X;
  ofstream output(ofile.c_str() );
  output<< X << '\t'
	<< ActiveVertexNumber << '\t'
	<< rho_vtx <<'\t'
	<< phi <<'\t'
	<< X * (phi-rho_vtx) << endl;
  //cout<<"Estimated FVS size at X="<<X<<" is "<<(1-rho_vtx)*VertexNumber<<endl;
  //cout.flush();
  return ;
}

/*-- externally fixing some variables to be empty and simplify the system --*/
int FVS::Fix0(void)
{
  int rank = 100; /* q_0 is distributed to 101 bins [0,0.01), [0.01,0.02), ...,
		     [0.99,1), [1,1] */
  int num_fixed_empty = 0;  //number of vertices externally fixed to unoccupied
  double mean_emptyprob = 0;        // ... and mean q_0 value of these vertices
  int num_examined_vertices = 0;        //number of examined candidate vertices
  struct vstruct *v_ptr;
  while(num_examined_vertices < MaxFixNumber && ActiveVertexNumber > 0)
    {
      int size=CandidateSize[rank];
      if(size>0)
	{
	  int *i_ptr = &CandidateVertices[rank * MaxFixNumber];
	  if(num_examined_vertices + size <= MaxFixNumber)
	    {
	      num_examined_vertices += size;
	      for(int s=0; s<size; ++s)
		{
		  v_ptr = &Vertex[ *i_ptr ];
		  if(v_ptr->active)
		    {
		      ++num_fixed_empty;
		      v_ptr->active = false;
		      v_ptr->occupied = false;
          exit(-1);
          cout << "Removed node: " << v_ptr->index <<endl;
		      Targets.push(v_ptr->index);
		      mean_emptyprob += v_ptr->q_0;
		      --ActiveVertexNumber;
		      Simplify(v_ptr) ;
		    }
		  ++i_ptr;
		}
	    }
	  else
	    {
	      while(num_examined_vertices < MaxFixNumber )
		{
		  ++num_examined_vertices;
		  v_ptr = &Vertex[ *i_ptr ];
		  if(v_ptr->active)
		    {
		      v_ptr->active = false;
		      v_ptr->occupied = false;
          exit(-1);
          cout << "Removed node: " << v_ptr->index <<endl;
		      Targets.push(v_ptr->index);
		      ++num_fixed_empty;
		      mean_emptyprob += v_ptr->q_0;
		      --ActiveVertexNumber;
		      Simplify(v_ptr) ;
		    }
		  ++i_ptr;
		}
	    }
	}
      --rank;
    }
  //cout<<" - number of active vertices: "<<ActiveVertexNumber<<endl;
  //cout.flush();
  return ActiveVertexNumber;
}

/* --- check whether the final occupation pattern corresponds to a FVS   --- */
bool FVS::CheckFVS(string& filename)
{
  int num_occup =0;
  int num_empty =0;                                                   //FVS size
  int num_active_edge =0;
  /* Consider the subgraph induced by all the occupied vertices and the edges
     between pairs of such vertices.
     -For each vertex in this subgraph, first determine its number of neighbors
     in this subgraph.
     - Then check whether all the edges of this subnetwork can be completely
     deleted by removing vertices of degree one iteratively.  */
  queue<int> LeafVertices;
  struct vstruct *v_ptr = &Vertex[1];
  struct message *im_ptr = &InMessage[0];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    {
      if(v_ptr->occupied == false)
	++num_empty;                                           //vertex deleted
      else
	{
	  ++num_occup;
	  im_ptr = v_ptr->im_ptr;
	  v_ptr->active_degree=0;                  //number of active neighbors
	  for(int d=0; d < v_ptr->degree; ++d, ++im_ptr)
	    if(im_ptr->v_ptr->occupied)
	      {
		++(v_ptr->active_degree);
		++num_active_edge;
	      }
	  if( v_ptr->active_degree == 1)
	    LeafVertices.push(v_ptr->index);
	}
    }
  num_active_edge /= 2;
  while( !LeafVertices.empty() )
    {
      v_ptr = &Vertex[ LeafVertices.front() ];
      LeafVertices.pop();
      im_ptr = v_ptr->im_ptr;
      for(int d=0; d < v_ptr->degree; ++d, ++im_ptr)
	if(im_ptr->v_ptr->occupied)
	  {
	    if( --(im_ptr->v_ptr->active_degree) == 1)
	      LeafVertices.push(im_ptr->v_ptr->index);
	  }
    }
  v_ptr = &Vertex[1];
  int TwoCoreSize=0;
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->occupied && v_ptr->active_degree != 0)
      ++TwoCoreSize;
  if(TwoCoreSize==0)
    {
      ofstream pfile(filename.c_str() );
      pfile <<"FVS "<<num_empty<<endl<<endl;
      //cout<< "FVS size = "<<num_empty<<endl;
      //cout.flush();
      v_ptr=&Vertex[1];
      for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
	if(v_ptr->occupied==false)
	  pfile << v_ptr->index<<endl;
      pfile.close();
      return true;
    }
  else
    {
      cerr<<"Not a proper FVS. The final two-core size is "<<TwoCoreSize<<endl;
      return false;
    }
}

/*
  Refinement of the components.
  Sthreshold is the maximal allowed component size.
*/
void FVS::ComponentRefinement(int Sthreshold)
{
  struct vstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    v_ptr->c_index=0;
  int max_comp_size=0;

  //           determine the size of each tree component, and break giant trees
  int c_index=0;                                           //index of component
  int c_size=0;                                             //size of component
  int NumberDelete=0;         //number of deleted vertices to break giant trees
  // Notice! In this subroutine, Permutation stores the size of each component.
  Permutation=0;
  v_ptr=&Vertex[1];
  struct message *im_ptr;
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->occupied && v_ptr->c_index==0)
      {
	++c_index;
	c_size=1;
	Vertex[v_ptr->index].c_index=c_index;
	queue<int> members;
	members.push(v_ptr->index);
	while( !members.empty() )
	  {
	    int vtx=members.front();
	    members.pop();
	    im_ptr = Vertex[vtx].im_ptr;
	    for(int d=0; d<Vertex[vtx].degree; ++d, ++im_ptr)
	      if(im_ptr->v_ptr->occupied)
		{
		  int vtx2=im_ptr->v_ptr->index;
		  if(Vertex[vtx2].c_index==0)
		    {
		      Vertex[vtx2].c_index=c_index;
		      members.push(vtx2);
		      c_size += 1;
		    }
		}
	  }
	if(c_size<=Sthreshold)
	  {
	    Permutation[c_index]=c_size;
	    if(c_size>max_comp_size)
	      max_comp_size=c_size;
	  }
	else
	  {
	    /* break giant tree component. Choose a vertex 'cvtx' such that its
	       deletion decreases the tree component size maximally. */
	    int cvtx;                                   // vertex to be deleted
	    double min_max_csize=c_size;//minimal possible size of max-sub-tree
	    queue<int> leafnodes;
	    set<int> gtree;
	    int vtx=v_ptr->index;
	    gtree.insert(vtx);
	    members.push(vtx);
	    while( !members.empty() )
	      {
		vtx=members.front();
		members.pop();
		Vertex[vtx].active=true;
		Vertex[vtx].active_degree=0;
		im_ptr = Vertex[vtx].im_ptr;
		for(int d=0; d<Vertex[vtx].degree; ++d, ++im_ptr)
		  if(im_ptr->v_ptr->c_index==c_index)
		    {
		      Vertex[vtx].active_degree += 1;
		      int vtx2=im_ptr->v_ptr->index;
		      if(gtree.find(vtx2)==gtree.end() )
			{
			  members.push(vtx2);
			  gtree.insert(vtx2);
			}
		    }
		if(Vertex[vtx].active_degree==1)
		  leafnodes.push(vtx);
	      }
	    while( !leafnodes.empty() )
	      {
		vtx=leafnodes.front();
		leafnodes.pop();
		Vertex[vtx].active=false;
		double maxbsize=0;                        //maximal branch size
		double psum=0;
		int d0=-1;                        //position of active neighbor
		im_ptr = Vertex[vtx].im_ptr;
		for(int d=0; d<Vertex[vtx].degree; ++d, ++im_ptr)
		  if(im_ptr->v_ptr->occupied)
		    {
		      im_ptr->v_ptr->active_degree -= 1;
		      if(im_ptr->v_ptr->active == false)
			{
			  double bsize=im_ptr->q_root;
			  if(maxbsize<bsize)
			    maxbsize=bsize;
			  psum += bsize;
			}
		      else
			{
			  d0=d;
			  if(im_ptr->v_ptr->active_degree ==1)
			    leafnodes.push(im_ptr->v_ptr->index);
			}
		    }
		if(d0>=0)
		  {
		    double bsize=c_size-1-psum;
		    if(maxbsize<bsize)
		      maxbsize=bsize;
		    struct mpointer *omptr_ptr=Vertex[vtx].omptr_ptr+d0;
		    omptr_ptr->m_ptr->q_root=psum+1;
		  }
		if(maxbsize<min_max_csize)
		  {
		    cvtx=vtx;
		    min_max_csize=maxbsize;
		  }
	      }
	    Vertex[cvtx].occupied=false;                       //vertex deleted
      cout << cvtx  << " 100000" << endl;
	    Targets.push(cvtx);
	    ++NumberDelete;
	    for(set<int>::const_iterator sci=gtree.begin(); sci != gtree.end();
		++sci)
	      Vertex[ *sci ].c_index=0;
	    --c_index;
	    --v;
	    --v_ptr;
	  }
      }
  // cout<<"Giant trees cut by deleting "<<NumberDelete
  //     <<" vertices.  Max-tree size: "<<max_comp_size<<endl;
  // cout.flush();
  //        add some vertices to small components and/or merge small components
  max_comp_size=0;
  int NumberAddition=0;
  set<IntPair> candidates;
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->c_index==0)
      {
	set<int> cgroups;
	int vtx=v_ptr->index;
	im_ptr = v_ptr->im_ptr;
	for(int d=0; d<v_ptr->degree; ++d, ++im_ptr)
	  {
	    int vtx2=im_ptr->v_ptr->index;
	    if(Vertex[vtx2].c_index !=0)
	      cgroups.insert( Vertex[vtx2].c_index);
	  }
	c_size=1;
	for(set<int>::const_iterator sci=cgroups.begin(); sci != cgroups.end();
	    ++sci)
	  c_size += Permutation[ *sci ];
	if(c_size<=Sthreshold)
	  candidates.insert(IntPair(c_size, vtx) );
      }
  while( !candidates.empty() )
    {
      c_index=0;
      set<IntPair>::const_iterator sci=candidates.begin();
      int csize0=sci->first;
      int vtx=sci->second;
      candidates.erase(IntPair(csize0, vtx) );
      set<int> cgroups;
      im_ptr = Vertex[vtx].im_ptr;
      for(int d=0; d<Vertex[vtx].degree; ++d, ++im_ptr)
	{
	  int vtx2=im_ptr->v_ptr->index;
	  if(Vertex[vtx2].c_index !=0)
	    {
	      cgroups.insert( Vertex[vtx2].c_index);
	      if(c_index<Vertex[vtx2].c_index)
		c_index=Vertex[vtx2].c_index;
	    }
	}
      c_size=1;
      for(set<int>::const_iterator sci2=cgroups.begin(); sci2 != cgroups.end();
	  ++sci2)
	c_size += Permutation[ *sci2 ];
      if(c_size<=Sthreshold)
	{
	  if(c_size>csize0)
	    candidates.insert( IntPair(c_size, vtx) );
	  else
	    {
	      ++NumberAddition;
	      Vertex[vtx].occupied=true;
	      Permutation[c_index]=c_size;
	      if(max_comp_size<c_size)
		max_comp_size=c_size;
	      queue<int> members;
	      members.push(vtx);
	      Vertex[vtx].c_index=c_index;
	      while( !members.empty() )
		{
		  int vtx2=members.front();
		  members.pop();
		  im_ptr = Vertex[vtx2].im_ptr;
		  for(int d=0; d<Vertex[vtx2].degree; ++d, ++im_ptr)
		    {
		      int vtx3=im_ptr->v_ptr->index;
		      if(Vertex[vtx3].c_index !=c_index &&
			 Vertex[vtx3].occupied)
			{
			  Permutation[Vertex[vtx3].c_index]=0;
			  Vertex[vtx3].c_index=c_index;
			  members.push(vtx3);
			}
		    }
		}
	    }
	}
    }
  // cout<<"Enlarge or merge small components by adding "
  //     <<NumberAddition      <<" vertices. Max-component size: "<<max_comp_size
  //     <<endl;
  // cout.flush();
  return;
}

/*  --- check the size of the components    --- */
bool FVS::CheckComponents( string& filename, int Sthreshold)
{
  int num_occup =0;
  int num_empty =0;                                            //number deleted
  //                                       determine the size of each component
  //  Notice! In this subroutine, Permutation stores the size of each component
  Permutation=0;
  struct vstruct *v_ptr = &Vertex[1];
  struct message *im_ptr = &InMessage[0];
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->occupied == false)
	++num_empty;                                           //vertex deleted
      else
	{
	  ++num_occup;
	  v_ptr->active=true;
	}
      ++v_ptr;
    }
  int c_index=0;
  int c_size=0;
  int max_c_size=0;
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->occupied && v_ptr->active)
      {
	++c_index;
	c_size=1;
	v_ptr->active=false;
	queue<int> members;
	members.push(v_ptr->index);
	while( !members.empty() )
	  {
	    int vtx=members.front();
	    members.pop();
	    im_ptr = Vertex[vtx].im_ptr;
	    for(int d=0; d < Vertex[vtx].degree; ++d, ++im_ptr)
	      if(im_ptr->v_ptr->occupied && im_ptr->v_ptr->active)
		{
		  im_ptr->v_ptr->active=false;
		  members.push(im_ptr->v_ptr->index);
		  ++c_size;
		}
	  }
	Permutation[c_index]=c_size;
	if(max_c_size<c_size)
	  max_c_size=c_size;
      }
  ///cout<<"Maximal component size "<<max_c_size
  //    <<";  Maximal allowed component size "<<Sthreshold<<endl;
  // cout.flush();
  ofstream pfile(filename.c_str() );
  // while (!Targets.empty()) {
  //       cout << ' ' << Targets.top();
  //       Targets.pop();
  //   }
  // cout << endl;
  // pfile <<"Targets  "<<num_empty<<endl<<endl;
  // cout <<"Targeted attack set  "<<num_empty<<endl;
  // cout.flush();
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->occupied==false)
      pfile << v_ptr->index<<endl;
  pfile.close();
  if(max_c_size<=Sthreshold)
    return true;
  else
    return false;
}

/*  --- read into a FVS pattern and check whether it is really a FVS --- */
bool FVS::ReadFVS(string& filename)
{
  while( !Targets.empty() )
    {
      Targets.pop();
    }
  int num_occup =VertexNumber;
  struct vstruct *v_ptr = &Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    v_ptr->occupied=true;
  int num_empty;                                                     //FVS size
  ifstream pfile(filename.c_str() );
  string tmpstring;
  pfile>>tmpstring;
  pfile>>num_empty;
  for(int l=0; l<num_empty; ++l)
    {
      int vtx;
      pfile>>vtx;
      Vertex[vtx].occupied=false;
      Targets.push(vtx);
    }
  pfile.close();
  num_occup -= num_empty;
  int num_active_edge =0;
  /* Consider the subgraph induced by all the occupied vertices and the edges
     between pairs of such vertices.
     -For each vertex in this subgraph, first determine its number of neighbors
     in this subgraph.
     - Then check whether all the edges of this subnetwork can be completely
     deleted by removing vertices of degree one iteratively.  */
  queue<int> LeafVertices;
  v_ptr = &Vertex[1];
  struct message *im_ptr = &InMessage[0];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->occupied)
      {
	im_ptr = v_ptr->im_ptr;
	v_ptr->active_degree=0;                    //number of active neighbors
	for(int d=0; d < v_ptr->degree; ++d, ++im_ptr)
	  if(im_ptr->v_ptr->occupied)
	    {
	      ++(v_ptr->active_degree);
	      ++num_active_edge;
	    }
	if( v_ptr->active_degree == 1)
	  LeafVertices.push(v_ptr->index);
      }
  num_active_edge /= 2;
  while( !LeafVertices.empty() )
    {
      v_ptr = &Vertex[ LeafVertices.front() ];
      LeafVertices.pop();
      im_ptr = v_ptr->im_ptr;
      for(int d=0; d < v_ptr->degree; ++d, ++im_ptr)
	if(im_ptr->v_ptr->occupied)
	  {
	    if( --(im_ptr->v_ptr->active_degree) == 1)
	      LeafVertices.push(im_ptr->v_ptr->index);
	  }
    }
  v_ptr = &Vertex[1];
  int TwoCoreSize=0;
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    if(v_ptr->occupied && v_ptr->active_degree != 0)
      ++TwoCoreSize;
  if(TwoCoreSize==0)
    {
      cerr<<"Input is indeed a FVS.\n";
      return true;
    }
  else
    {
      cerr<<"Not a proper FVS. The final two-core size is "<<TwoCoreSize<<endl;
      return false;
    }
}

/*
  Report the attack order and its effect.
*/
void FVS::AttackEffect(string& attackfile)
{
  struct vstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    v_ptr->c_index=0;
  set<int> finaltargets;
  int max_comp_size=0;
  int num_empty=0;
  //                                       determine the size of each component
  int c_index=0;                                           //index of component
  int c_size=0;                                             //size of component
  //  Notice! In this subroutine, Permutation stores the size of each component
  Permutation=0;
  v_ptr=&Vertex[1];
  struct message *im_ptr;
  for(int v=1; v<=VertexNumber; ++v, ++v_ptr)
    {
      if(v_ptr->occupied==false)
	{
	  ++num_empty;
	  finaltargets.insert(v_ptr->index);
	}
      else if(v_ptr->c_index==0)
	{
	  ++c_index;
	  c_size=1;
	  Vertex[v_ptr->index].c_index=c_index;
	  queue<int> members;
	  members.push(v_ptr->index);
	  while( !members.empty() )
	    {
	      int vtx=members.front();
	      members.pop();
	      im_ptr = Vertex[vtx].im_ptr;
	      for(int d=0; d<Vertex[vtx].degree; ++d, ++im_ptr)
		if(im_ptr->v_ptr->occupied)
		  {
		    int vtx2=im_ptr->v_ptr->index;
		    if(Vertex[vtx2].c_index==0)
		      {
			Vertex[vtx2].c_index=c_index;
			members.push(vtx2);
			c_size += 1;
		      }
		  }
	    }
	  Permutation[c_index]=c_size;
	  if(c_size>max_comp_size)
	    max_comp_size=c_size;
	}
    }
  ofstream output(attackfile.c_str() );
  output<< (1.0e0*num_empty)/(1.0e0*VertexNumber) << '\t'
	<< num_empty << '\t'
	<< (1.0e0*max_comp_size)/(1.0e0*VertexNumber) << '\t'
	<< max_comp_size << '\t'
	<< 0 <<endl;
  while( !Targets.empty() )
    {
      int vtx=Targets.top();
      Targets.pop();
      if(Vertex[vtx].occupied==false)
	{
	  int c_index=0;
	  set<int> cgroups;
	  im_ptr = Vertex[vtx].im_ptr;
	  for(int d=0; d<Vertex[vtx].degree; ++d, ++im_ptr)
	    {
	      int vtx2=im_ptr->v_ptr->index;
	      if(Vertex[vtx2].c_index !=0)
		{
		  cgroups.insert( Vertex[vtx2].c_index);
		  if(c_index<Vertex[vtx2].c_index)
		    c_index=Vertex[vtx2].c_index;
		}
	    }
	  int c_size=1;
	  for(set<int>::const_iterator sci2=cgroups.begin();
	      sci2 != cgroups.end();  ++sci2)
	    c_size += Permutation[ *sci2 ];
	  Vertex[vtx].occupied=true;
	  --num_empty;
	  Permutation[c_index]=c_size;
	  if(max_comp_size<c_size)
	    max_comp_size=c_size;
	  queue<int> members;
	  members.push(vtx);
	  Vertex[vtx].c_index=c_index;
	  while( !members.empty() )
	    {
	      int vtx2=members.front();
	      members.pop();
	      im_ptr = Vertex[vtx2].im_ptr;
	      for(int d=0; d<Vertex[vtx2].degree; ++d, ++im_ptr)
		{
		  int vtx3=im_ptr->v_ptr->index;
		  if(Vertex[vtx3].c_index !=c_index &&
		     Vertex[vtx3].occupied)
		    {
		      Permutation[Vertex[vtx3].c_index]=0;
		      Vertex[vtx3].c_index=c_index;
		      members.push(vtx3);
		    }
		}
	    }
	  output<< (1.0e0*num_empty)/(1.0e0*VertexNumber) << '\t'
		<< num_empty << '\t'
		<< (1.0e0*max_comp_size)/(1.0e0*VertexNumber) << '\t'
		<< max_comp_size <<'\t'
		<< vtx << endl;
	}
    }
  output.close();
  for(set<int>::const_iterator sci=finaltargets.begin();
      sci!=finaltargets.end(); ++sci)
    Vertex[ *sci ].occupied=false;
  return;
}
