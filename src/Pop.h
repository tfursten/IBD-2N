#ifndef POP_H_INCLUDED
#define POP_H_INCLUDED

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <map>
#include <boost/foreach.hpp>
#include <algorithm>
#include "xorshift64.h"
#include "rexp.h"
#include "disperse.h"
#include "disk.h"
#include "aliastable.h"


typedef pair<int,int> xyCoord;
typedef vector<int> gam;
struct individual
{
    unsigned int nWeight_1;
    unsigned int nWeight_2;
    gam vgamete_1;
    gam vgamete_2;
    int nParent_1;
    int nParent_2;
    int ngameteID_1;
    int ngameteID_2;
    xyCoord xy;
    int raw_count;
    vector<int> parent_count;
    int last_parent;

    individual(xyCoord x, unsigned int weight,gam & gamete_1,gam &gamete_2, int parent_1, int parent_2) {
		nWeight_1 = weight;
		nWeight_2 = weight;
		vgamete_1 = gamete_1;
		vgamete_2 = gamete_2;
		nParent_1 = parent_1;
		nParent_2 = parent_2;
		ngameteID_1 = 0;
		ngameteID_2 = 0;
        xy = x;
        raw_count = 0;
        last_parent = -1;
    }
};


class Population
{
private:
    int m_nMaxX;
    int m_nMaxY;
    std::string m_sBound;
    int m_nOffspring;
    double m_dSigma;
    double m_dTotMut;
    int m_nMutCount;
	int m_nIndividuals;
	int m_nMarkers;
	int m_nSample;
	int m_nDistClass;
	int m_nPairs;
	int m_nTotGametes;
	int m_nAlleles;
	xorshift64 m_myrand;
	std::ofstream & pout;
	std::ofstream & nbout;
	bool verbose;
	Dispersal disp;
	alias_table alias_mut;
	map<int,int> m_mMutations;
	std::vector<double> m_vdMut;
	std::vector<int> m_vSample;
	std::vector<individual> m_vPop1;
	std::vector<individual> m_vPop2;
	int m_nAlleleID;
	int setMutCount();
	void disperse_step(int parent);
	void mutation_step();
	void reproduction_step(int offspring);
	void mutate_iam(gam & gamete);
	void mutate_smm(gam & gamete);
	void mutate_kam(gam & gamete);
	void sampleIBD(int gen);
	void samplePop();
	void sampleDist();
	void sampleNb(int gen);
	float nbSize(vector<int> & counts);
	gam make_gamete(gam &parent1, gam &parent2);

protected:
    int(Population::*disperse)(int,int);
    void(Population::*mutate)(gam &gamete);

public:
    
    Population(std::ofstream &p,std::ofstream &nb, bool v): pout(p), nbout(nb), verbose(v) {};
    void initialize(int nMaxX, int nMaxY, int nOffspring, int nMarkers, double dSigma,
    	vector<double> vdMut, unsigned int seed, int nSample, 
    	string dist_name, float param, 
    	bool fast, int nclass, int npairs, string mut_type, int nAlleles);
	void evolve(int m_nGenerations, int m_nBurnIn);
};

#endif // POP_H_INCLUDED
