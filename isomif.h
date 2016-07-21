// IsoMIF is a program to identify molecular interaction field similarities between proteins
// Copyright (C) 2015 - Matthieu Chartier (under the supervision or Rafael Najmanovich)

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iterator>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <map>

// #ifdef _OPENMP
// # include <omp.h>
// // multi-processors
// # define CHUNKSIZE 1
// # define SCHEDULE dynamic
// // #  define SCHEDULE static
// # define NUM_THREADS 12
// #endif

using namespace std;

typedef struct clefts cleft;
struct clefts{
  float min_x,max_x,min_y,max_y,min_z,max_z;
  int width,height,depth;
  int vrtx_count;
};

struct nodeI{
  int id;
  int neibrs;
};

struct pwRun{
  string mif1;
  string mif2;
  string rnc1;
  string rnc2;
  int getrmsd;
};

struct vertex {
  float coor[3];
  float ncoor[3];
  int grid[4];
  vector<int> cg;
  vector<int> pb;
  vector<float> nrg;
  vector<float> ang;
  vector<int> m;
  vector<int> ol;
  int id;
};

//Atom struc
typedef struct atoms atom;
struct atoms{
  float coor[3];
  float ncoor[3];
  string resn;
  string atomn;
  string chain;
  string alt;
  int resnb; //residue ID in PDB
  int atomnb; //atomnb in PDB
  int id;
  int bs;
  int mif;
};

//Atom struc
struct pseudoC{
  float coor[3];
  string type;
  int id;
};

struct mif {
  vector<vertex> mif;
  vector<atom> prot;
  vector<int> ss;
  vector<int> ssm;
  int caSize;
  vector<pseudoC> pseudoL;
  string rnc;
  vector<atom> lig;
};

typedef struct nodes node;
typedef node *pNode;
struct nodes{
  vertex* a;
  vertex* b;
  atom* ca;
  atom* cb;
  pseudoC* pa;
  pseudoC* pb;
  float cosim;
  int nbi;
  int nbiw;
  float nrg;
};

typedef struct CliqueStruct Clique;
struct CliqueStruct{
  int cg;
  vector<node> nodes;
  vector<vertex> va;
  vector<vertex> vb;
  int nbNodes;
  int nbNodesM;
  float nbNodesMW;
  float normNodes;
  float normNodesRMSD;
  float tani;
  float taniM;
  vector<int> pbweight;
  vector<float> nrgsum;
  vector<int> angCount;
  vector<float> angSum;
  float taniMW;
  float taniNorm;
  float rmsd;
  float ligRMSD;
  float nrg;
  gsl_matrix *mat_r;
  float cen_a[3];
  float cen_b[3];
  double det;
  double detOri;
};

int *compsub = NULL;
int c,numNodes;
int numEdges=0;
int Clique_threshold=3;
int stopBk = 0;
float topT=-1.0;
float topN=-1;
int bkAll = 0;
int nCliques=0;
int nCliquesExplored=0;
int maxCliques=200;
vector<Clique> cliques;

int emptOut=0;
char outH[4000];
int printDetails=0;
string pairwiseF;
string nrg_file1;
string nrg_file2;
char exePath[150];
char cmdLine[550];
char cmdArgs[550];
int flagpp=0; //To store if cmdline is -pp or -p1 -p2
char outbase[200];
int commonInt=1;

int skipDet=1;
int pc=0;
int wc=0;
int ol=-1;
float olDist=1.0;
float olDistsq=olDist*olDist;
int jttt=5;
int jtt[20][20];
int nb_of_probes=6;
int cg2=1;
char out_file[150]; //Output filename, constructed using outbase in get_info()
char outPDB[150];
char outPDB2[150];
int maxNodes=100000;
float dDist=3.0;
float ca_dDist=3.0;
float ps_dDist=1.0;
float neibr_dDist=3.0;
int cg_start=-1;
float cosdT=0.90;
vector<vertex> mif1;
vector<vertex> mif2;
vector<atom> prot1;
vector<atom> prot2;
vector<pseudoC> pseudoL1;
vector<pseudoC> pseudoL2;
string tag1;
string tag2;
int caSize1=0;
int caSize2=0;
vector<int> list1;
vector<int> list2;
vector<int> steps;
int wrfn=0; //write sim in filename
vector<int> ss1;
vector<int> ss2;
vector<int> ss1m;
vector<int> ss2m;
int totalVrtx=0;
map<int,int> topCliques;
map<string,mif> mifs;
vector<atom> lig1;
vector<atom> lig2;
string rnc1;
string rnc2;
int getrmsd=0;
vector<pwRun> pw;
FILE* fpout;
stringstream matchFileOut;

void getPairwise();
void sortArray(int *&, int nn, bool*&);
bool myfunction (nodeI,nodeI);
void AddNewClique(int, int*, int, vector<node>&);
// bool compareSim(node ,node);
void bk(int, vector<node>&, bool*&);
void adjMat(vector<node>&, bool*&, int);
void Extend(int* old,int ne,int ce, int cg, vector<node>&, bool*&, int);
int open_file_ptr(FILE**, char*, int);
float dist3d(float[], float[]);
float dist3dnosqrt(float[], float[]);
void clearStep(int);
void printNodes();
int read_commandline(int, char*[]);
bool fexists(const char*);
int createVrtxVec(string, vector<vertex>&, vector<atom>&, vector<int>&, vector<int>&, int&, vector<pseudoC>&, string, vector<atom>&);
void createNodes(int, vector<node>&, int);
void rem_spaces(char*);
int get_info(string, string);
double calcRot(vector<float>, vector<float>, float*, float*, gsl_matrix *, double&);
double SupSVD(gsl_matrix *, double&);
double gsl_matrix_Det3D(gsl_matrix *);
long long ConnID(int, int);
void setJTT(int);
int nam2n(string);
int samePseudo(pseudoC&,pseudoC&);