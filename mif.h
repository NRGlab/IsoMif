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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <math.h> 
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iterator>

using namespace std;

#define NEW(p,type)     if ((p=(type *) malloc (sizeof(type))) == NULL) { \
        printf ("Could not allocate ptr!\n");                                    \
        exit(0);                                                        \
}

#define PI 3.14159265

//Atom struc
typedef struct atoms atom;
struct atoms{
  float x,y,z,xr,yr,zr;
  string resn;
  string atomn;
  string chain;
  string alt;
  string pseudo;
  string at;
  int resnb; //residue ID in PDB
  int atomnb; //atomnb in PDB
  int dir;
  int rDir;
  int h; //if its a hydrogen 0 or 1
  int bs; //if its a binding site atom
  int mif; //if it should be used to calculated MIFS
};

//Grid point struc
typedef struct vertexes vertex;
struct vertexes{
  float x,y,z;
  int p;
  int bu;
  vector<int> ints;
  vector<float> nrgs;
  vector<float> angles;
  int grid[4];
  map<string,float> env;
  int id;
  int modulo;
};

struct angRef{
    string r1;
    string r2;
    int rDir;
    int ring;
};

struct pseudovrtx{
    float x,y,z;
    float dist;
    string type;
};

struct pbVrtx{
  float dist;
  float nrg;
};

struct pwRun{
  string pdbF;
  string cleftF;
  string gridF;
  string rnc;
  string ligF;
};
vector<pwRun> pw;

string pairwiseF="";
string cmdLine="";
string cleftFile="";
string gridFile="";
string proteinFile="";
string outBase="";
string outGridBase="";
string tag="";
string chain="";
string basePath="";
string resnumc="";
string resnumcShort="";
string matrixF="";
string probesF="";
string ligFile="";
string ff="original";
string statsF="original";
float distpbV=2.0;
float stepsize=0.5;
float maxGridDist=4.0;
float minGridDist=2.5;
float atmPbMaxDist=8.0;
float gridLigDist=2.0;
float caT=5.0;
int uID=0;
int smoothDist=0;
int printDetails=0;
int printAtoms=1;
int* probesList;
float* pbDistTmax;
float* pbDistTmin;
int nbOfAts=0;
int nbOfAtoms=0;
int nbOfProbes=0;
int ss[4];
int ssm[4];
int zip=-1;
int bul=0;
int buD=40;
float* epsilons;
float angThresh=60.0;
map<string,string> atomTypes;
map<string,angRef> atomRef;
map<int,string> idAt;
map<string,string> pseudoC;
map<string,string> ligAt;
map<string,int> eps;
map<string,int> hyd;
map<string,int> arm;
map<string,int> don;
map<string,int> acc;
map<string,int> chr;
map<string,float> minD;
map<string,float> maxD;
map<string,float> nrgT;
vector<string> probes;
vector<string> aa;
vector<pseudovrtx> pseudoList;
float min_x, min_y, min_z, max_x, max_y, max_z;
int width, height, depth;

class  Protein{
  public:

    Protein(string);
    ~Protein(void);
    void readPDB(string);
    void getAtomDir();
    int getRefAtom(float&, float&, float&, string, int, string, string, int, float, float, float, string);
    vector<atom> PROTEIN;
    vector<float> LIGAND;
    vector<atom> LIGATOMS;

  private:

};

class  Grid{
  public:

    Grid(string, Protein&);
    ~Grid(void);
    int readGetCleft(string, vector<atom>&, vector<float>&);
    int generateID(int, int, int, int, int);
    int buildGrid(vector<atom>&);
    void getProtVrtx(vector<atom>&);
    void createProtVrtx(vector<atom>&);
    void getMinMax(string);
    void getBuriedness();
    void writeGrid();
    int readGrid();
    int getDiag(int, int, int, int&);
    int inGridRes(vertex&, float);
    void smooth();
    void writeMif(vector<atom>&,vector<atom>&);
    map<int,vertex> GRID;
    vector<int> vrtxIdList;

  private:
    int vrtx050,vrtx100,vrtx150,vrtx200;
};

int readCmdLine(int, char**);
float dist_3d(float, float, float, float, float, float);
void getMif(map<int,vertex>&, vector<atom>&,vector<int>&);
void getEnv(map<int,vertex>&, vector<atom>&,vector<int>&);
void getPseudo(map<int,vertex>&, vector<atom>&,vector<int>&);
int is_coord_in_cube(float,float,float,float,float,float,float);
void stripSpace(string &);
void getStats(map<int,vertex>&, vector<atom>&, vector<int>&);
float roundCoord(float, int);
double calcNrg(vertex&, atom&, int, int&, float&);
void getAtomRef();
void getPseudoC();
void getEpsilons();
void getAtomTypes();
void readLigFile();
void getProbes();
void getaa();
bool compByNrg(pbVrtx,pbVrtx);
void getPairwise();