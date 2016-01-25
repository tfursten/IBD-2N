#include "Pop.h"

//urandom dev/urandom if it exists use it else use create random seed
using namespace std;
// random seed generator
inline unsigned int create_random_seed() {
	unsigned int v;
	//ifstream urandom("/dev/urandom", ios::in|ios::binary);
	//if(urandom.good()) {
		//urandom >> v;
	//} else {
	v = static_cast<unsigned int>(getpid());
	v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
	v^=(v<<17);
	v^=(v>>13);
	v^=(v<<5);
	v += static_cast<unsigned int>(time(NULL));
	v^=(v<<17);
	v^=(v>>13);
	v^=(v<<5);
	//}
    return (v == 0) ? 0x6a27d958 : (v & 0x7FFFFFFF); // return at most a 31-bit seed
}



inline xyCoord i2xy(int i, int mx, int my)
{
    assert(0 <= i && i < mx*my);
    return make_pair((i/my),(i%my));
}

inline int xy2i(int x, int y, int mx, int my) {
	assert(0 <= x && x < mx);
	assert(0 <= y && y < my);
	return x*my+y;
}

inline int xy2i(xyCoord xy, int mx, int my) {
	return xy2i(xy.first,xy.second,mx,my);
}

void Population::initialize(int nMaxX, int nMaxY, int nOffspring, int nMarkers, double dSigma, vector<double> vdMut,
                            unsigned int seed, int nSample, int nPopSample,
                            string dist_name, float param, bool fast, int nclass, int npairs)
{
    ostringstream out;

    //set Random seed
    if (seed)
        out << "User set PRNG seed to: " << seed << ".\n";
    else {
        seed = create_random_seed();
        out << "Using Generated PRNG Seed: "<< seed << endl;
    }
    m_myrand.seed(seed);
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_nOffspring = nOffspring;
    m_nDistClass = nclass;
    m_nPairs = npairs;
    m_vdMut = vdMut;
    m_dTotMut = 1.0;
	for(vector<double>::iterator it = m_vdMut.begin();it!=m_vdMut.end();++it){
		m_dTotMut *= (1-*it);
	}
	m_dTotMut = -log(m_dTotMut);
	for(vector<double>::iterator it = m_vdMut.begin();it!=m_vdMut.end();++it){
		*it = *it/m_dTotMut;
	}
	alias_mut.create(m_vdMut.begin(),m_vdMut.end());
    m_nIndividuals = nMaxX * nMaxY;
    m_nTotGametes = m_nIndividuals * m_nOffspring;
    m_nAlleleID = 0;
    m_nMarkers = nMarkers;
    m_nMutCount = setMutCount();
    m_nSample = nSample;
    m_nPopSample = nPopSample;
    m_dSigma = dSigma;
    disp.initialize(dist_name, m_nMaxX, m_nMaxY, fast, "torus", dSigma, param);
    out << "Dispersal distribution set to " << disp.getName() << ".\n" ;
    out << "Extra parameter set to " << param << ".\n";


    // Initialize Population: each individual has unique allele
    for(int iii=0; iii<m_nIndividuals; iii++) {
        gam h1;
        gam h2;
        gam h3 (m_nMarkers,0);
        for(int iii=0; iii<m_nMarkers; iii++){
        	h1.push_back(m_nAlleleID++);
        	h2.push_back(m_nAlleleID++);
        }
        xyCoord xy = i2xy(iii, m_nMaxX, m_nMaxY);
        m_vPop1.emplace_back(xy,1,h1,h2,iii,iii);
        m_vPop2.emplace_back(xy,0,h3,h3,0,0);
    }

    //write settings to screen and file
    cout << out.str();
    pout << out.str();

}

int Population::setMutCount() {
    return floor(rand_exp(m_myrand, m_dTotMut));
}


//Each mutational event creates a new allele unlike any other allele currently in the population
//so that identity in state for two or more alleles is always an indication of identity by descent.
void Population::mutate(gam & gamete)
{
    gamete[alias_mut(m_myrand.get_uint64())] = m_nAlleleID ++;
}

void Population::evolve(int m_nBurnIn, int m_nGenerations)
{
    //run burn-in period
    for(int ggg=0;ggg<m_nBurnIn;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
        {
            disperse_step(parent);
        }
        
        mutation_step();
        
        for(int offspring=0; offspring<m_nIndividuals; offspring++)
        {
        	reproduction_step(offspring);
        }
        swap(m_vPop1,m_vPop2);
    }
    //outfile headers
    nbout << "Gen\tX\tY\tNb\tNb_raw\tp1_d\tp2_d\t";
    for(int m=1; m <= m_nMarkers; m++){
    	nbout << "M" << m << (m==m_nMarkers ? "\n" : "\t");
    }


    for(int ggg=0;ggg<m_nGenerations;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
        {
            disperse_step(parent);
        }
        
        mutation_step();

        for(int offspring=0; offspring<m_nIndividuals;offspring++)
        {
        	reproduction_step(offspring);
        }

        if (ggg % m_nSample == 0)
            sampleNb(ggg);
        //}
        //if (ggg % m_nPopSample == 0)
        //    samplePop();
        swap(m_vPop1,m_vPop2);
    }
    //sampleDist();
}


void Population::disperse_step(int parent)
{

    individual &parentHere = m_vPop1[parent];
    //check if there is a parent here
    if(parentHere.nWeight_1 == 0 || parentHere.nWeight_2 == 0){
    	parentHere.nWeight_1 = 0;
    	parentHere.nWeight_2 = 0;
    	parentHere.parent_count.clear();
    	parentHere.raw_count = 0;
    	parentHere.last_parent = -1;
    	return;
    }
    //clear out weight/counts for next generation
    parentHere.nWeight_1 = 0;
    parentHere.nWeight_2 = 0;
    parentHere.parent_count.clear();
    parentHere.last_parent = -1;
    parentHere.raw_count = 0;
    for (int off=0; off<m_nOffspring; off++)
    {
        //disperse
        int nNewCell = disp(m_myrand,parentHere.xy.first, parentHere.xy.second);
        //cout << "NEW CELL: " << nNewCell << endl;
        //for absorbing boundary: check if individual dispersed off the grid
        if (nNewCell == -1)
        {
            //mutate(parentHere.nAllele);
            continue;
        }
        individual &parentThere = m_vPop2[nNewCell];
        unsigned int nSeedWeight = m_myrand.get_uint32();

        //cout << "Seed Weight Outter" << nSeedWeight << endl;
        //Count how many unique parents are competing for cell
        if(parent > parentThere.last_parent){
        	parentThere.raw_count ++;
        	parentThere.last_parent = parent;
        	parentThere.parent_count.push_back(1);
        }
        else{
        	parentThere.parent_count.back()++;
        }

        //competition for cell
        unsigned int min_weight = min(parentThere.nWeight_1, parentThere.nWeight_2);
        //cout << parentThere.nWeight_1 << " " << parentThere.nWeight_2 << endl;
        bool replace = (parentThere.nWeight_1 <= parentThere.nWeight_2 ? 0 : 1);
        if(nSeedWeight > min_weight)
        {
        //	cout << "MIN WEIGHT: " << min_weight << endl;
        //	cout << "SEED WEIGHT: " << nSeedWeight << endl;
        //	cout << "IN" << replace << endl;
        	if(replace){
        //		cout << "R2" << endl;
        		parentThere.nWeight_2 = nSeedWeight;
        		parentThere.nParent_2 = parent;
        		parentThere.ngameteID_2 = parent*m_nOffspring + off;
        	}
        	else{
        //		cout << "R1" << endl;
        		parentThere.nWeight_1 = nSeedWeight;
        		parentThere.nParent_1 = parent;
        		parentThere.ngameteID_1 = parent*m_nOffspring + off;
        //		cout << parentThere.nWeight_1 << endl;
        	}
        }
    }
}

void Population::mutation_step(){
	//Make a map of gamete ID's that need a mutation
	//it is possible to draw the same gamete twice
	m_mMutations.clear();
	int count = m_nMutCount;
	int n = 0;
	while(count<m_nTotGametes){
		n = m_nTotGametes - count;
		std::map<int,int>::iterator it;
		it = m_mMutations.find(count);
		if (it !=m_mMutations.end())
			it->second += 1;
		else
			m_mMutations[count] = count;
		//m_vMutation.push_back(count);
		count+= setMutCount();
	m_nMutCount = count - n;
	}
}

gam Population::make_gamete(gam &parent1, gam &parent2){
	vector< gam > genotype;
	gam newGamete;
	genotype.push_back(parent1);
	genotype.push_back(parent2);
	//use top 32 bits
	uint64_t r = m_myrand.get_uint64() >> 32;
	if(m_nMarkers < 32){
		for(int m=0; m<m_nMarkers; m++){
			newGamete.push_back(genotype[r&1][m]);
			r >>= 1;
		}
	}
	else{
		for(int m=0; m<m_nMarkers; m++){
			if(!m%32)
				r = m_myrand.get_uint64() >> 32;
			newGamete.push_back(genotype[r&1][m]);
			r >>= 1;
		}
	}
	return newGamete;
}

void Population::reproduction_step(int offspring)
{
	individual & offspringHere = m_vPop2[offspring];
	//make sure 2 gametes landed in cell
	//cout << offspringHere.nWeight_1 << " " << offspringHere.nWeight_2 << endl;
	if(offspringHere.nWeight_1 == 0 || offspringHere.nWeight_2 == 0){
		cout << "WARNING: INCOMPLETELY FILLED CELL" << endl;
		return;
	}
	//generate a gamete from both parents
	individual & parent1 = m_vPop1[offspringHere.nParent_1];
	individual & parent2 = m_vPop1[offspringHere.nParent_2];
	offspringHere.vgamete_1 = make_gamete(parent1.vgamete_1,parent1.vgamete_2);
	offspringHere.vgamete_2 = make_gamete(parent2.vgamete_1,parent2.vgamete_2);
	//check for mutation on both gametes
	std::map<int,int>::iterator it;
	it = m_mMutations.find(offspringHere.ngameteID_1);
	if(it != m_mMutations.end()){
		for(int iii=0; iii<it->second; iii++){
			mutate(offspringHere.vgamete_1);
		}
	}
	it = m_mMutations.find(offspringHere.ngameteID_2);
	if(it != m_mMutations.end()){
		for(int iii=0; iii<it->second; iii++){
			mutate(offspringHere.vgamete_2);
		}
	}
}

double minEuclideanDist(xyCoord xy1, xyCoord xy2, int mx, int my){
    double dx = abs(1.0*(xy1.first - xy2.first));
    double dy = abs(1.0*(xy1.second - xy2.second));
    dx = (dx < mx*0.5) ? dx : mx-dx;
    dy = (dy < my*0.5) ? dy : my-dy;
    return (dx*dx+dy*dy);
}


float Population::nbSize(vector<int> & counts){
	assert(!counts.empty);
	int tot = 0;
	int psum = 0;
	for(std::vector<int>::iterator it = counts.begin(); it != counts.end(); ++it){
    	tot += *it;
    	psum += *it * *it;
	}
	return 1/(psum/(float)(tot*tot));
}

void Population::sampleNb(int gen){
	map<int,int> alleles;
    vector<individual> pop(m_vPop2);
    random_shuffle(pop.begin(),pop.end());
    vector<int> vN(m_nDistClass,0);
    while(pop.size()){
        individual ind1 = pop[0];
        if(ind1.nWeight_1==0 || ind1.nWeight_2==0){
        	pop.erase(pop.begin());
        	continue;
        }
        int x1 = ind1.xy.first;
        int y1 = ind1.xy.second;
        bool found = false;
        for(unsigned int i = 1; i < pop.size(); i++){
            individual ind2 = pop[i];
            if(ind2.nWeight_1==0 || ind2.nWeight_2==0){
            	continue;
            }
            int x2 = ind2.xy.first;
            int y2 = ind2.xy.second;
            int d;
            if(x1 == x2)
                d = abs(y1-y2);
            else if(y1 == y2)
                d = abs(x1-x2);
            else continue;
            if(d > m_nDistClass || vN[d-1] >= m_nPairs)
                continue;

            double sig1 = minEuclideanDist(ind1.xy,m_vPop1[ind1.nParent_1].xy, m_nMaxX, m_nMaxY);
            double sig2 = minEuclideanDist(ind1.xy,m_vPop1[ind1.nParent_2].xy, m_nMaxX, m_nMaxY);
            double sig3 = minEuclideanDist(ind2.xy,m_vPop1[ind2.nParent_1].xy, m_nMaxX, m_nMaxY);
            double sig4 = minEuclideanDist(ind2.xy,m_vPop1[ind2.nParent_2].xy, m_nMaxX, m_nMaxY);

            nbout << gen << "\t"<< x1 <<"\t"<< y1
             << "\t" << nbSize(ind1.parent_count) << "\t" << ind1.raw_count << "\t"
             << sig1 << "\t" << sig2 << "\t";
            for(int m=0; m<m_nMarkers; m++){
            	nbout << ind1.vgamete_1[m] << "/" << ind1.vgamete_2[m]
            	<< ((m<m_nMarkers-1) ? "\t": "\n");
            	if(verbose){
            		alleles[ind1.vgamete_1[m]] ++;
            		alleles[ind1.vgamete_2[m]] ++;
            	}
            }
            nbout << gen << "\t" << x2 << "\t" << y2
            << "\t" << nbSize(ind2.parent_count) << "\t" << ind2.raw_count << "\t"
            << sig3 << "\t" << sig4 << "\t";
            for(int m=0; m<m_nMarkers; m++){
            	nbout << ind2.vgamete_1[m] << "/" << ind2.vgamete_2[m]
            	<< ((m<m_nMarkers-1) ? "\t": "\n");
            	if(verbose){
            		alleles[ind2.vgamete_1[m]] ++;
            		alleles[ind2.vgamete_2[m]] ++;
            	}
            }
            vN[d-1] += 1;
            pop.erase(pop.begin()+i);
            pop.erase(pop.begin());
            found = true;
            break;
        }
        if(!found)
            pop.erase(pop.begin());
    }
    if(verbose){
    	cout << "Gen: " << gen << " Ko: " << alleles.size() << endl;
    	pout << "Gen: " << gen << " Ko: " << alleles.size() << endl;
    }

}