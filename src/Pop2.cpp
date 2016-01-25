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

void Population::initialize(int nMaxX, int nMaxY, int nOffspring, double dSigma,
                            unsigned int seed,
                            string dist_name, string bound, float param, bool fast)
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
    m_sBound = bound;
    m_nOffspring = nOffspring;
    m_nIndividuals = nMaxX * nMaxY;

    m_dSigma = dSigma;

    disp.initialize(dist_name, m_nMaxX, m_nMaxY, fast, m_sBound, dSigma, param);
    out << "Dispersal distribution set to " << disp.getName() << ".\n" ;
    out << "Extra parameter set to " << param << ".\n";


    // Initialize Population: each individual has unique allele
    for(int iii=0; iii<m_nIndividuals; iii++) {
        xyCoord xy = i2xy(iii, m_nMaxX, m_nMaxY);
        m_vPop1.emplace_back(xy, 1,iii+1,iii);
        m_vPop2.emplace_back(xy);
    }

    //write settings to screen and file
    cout << out.str();
    pout << out.str();

}

void Population::setMutCount() {
    m_nMutCount = floor(rand_exp(m_myrand, m_dMut));
}


//Each mutational event creates a new allele unlike any other allele currently in the population
//so that identity in state for two or more alleles is always an indication of identity by descent.
int Population::mutate(int allele)
{
    if (--m_nMutCount > 0)
        return allele;
    setMutCount();
    return m_nAlleleID++;
}

void Population::evolve(int m_nGenerations)
{
    for(int ggg=0;ggg<m_nGenerations;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++){
            step(parent);
        }
        swap(m_vPop1,m_vPop2);
    }
}


void Population::step(int parent)
{
    individual &parentHere = m_vPop1[parent];
    for (int off=0; off<m_nOffspring; off++)
    {
        //disperse
        int nNewCell = disp(m_myrand,parentHere.xy.first, parentHere.xy.second);
    }
}

xyCoord minEuclideanDist(int i, int j, int mx, int my){
    //calculate minimum distance between two positions for periodic boundary
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dx = abs(1.0*(xy1.first - xy2.first));
    double dy = abs(1.0*(xy1.second - xy2.second));
    dx = (dx < mx*0.5) ? dx : mx-dx;
    dy = (dy < my*0.5) ? dy : my-dy;
    return make_pair(dx,dy);
}

double minEuclideanDist2(int i, int j, int mx, int my) {
    xyCoord xy = minEuclideanDist(i,j,mx,my);
	return (xy.first*xy.first+xy.second*xy.second);
}

double minEuclideanDist3(int i, int j, int mx, int my){
    xyCoord xy = minEuclideanDist(i,j,mx,my);
    return (xy.first*xy.first*xy.first+xy.second*xy.second*xy.second);
}

xyCoord euclideanDist(int i, int j, int mx, int my){
    //calculate distance between two positions for absorbing boundary
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dx = abs(1.0*(xy1.first - xy2.first));
    double dy = abs(1.0*(xy1.second - xy2.second));
    return make_pair(dx,dy);
}

double euclideanDist2(int i, int j, int mx, int my) {
    xyCoord xy = euclideanDist(i,j,mx,my);
    return (xy.first*xy.first+xy.second*xy.second);
}

double euclideanDist3(int i, int j, int mx, int my){
    xyCoord xy = euclideanDist(i,j,mx,my);
    return (xy.first*xy.first*xy.first+xy.second*xy.second*xy.second);
}


void Population::sampleIBD(int gen)
{
    vector<int> vIBD(m_nLenTrans,0);
    vector<int> vgIBD(m_nLenTrans,0);
    vector<int> vpIBD(m_nLenTrans,0);
    vector<int> vN(m_nLenTrans,0);
    typedef map<int,int> mapType;
    mapType alleleMap;
    int szSample = 0;
    double dSigma2 = 0.0;
    double dSigma3 = 0.0;
    double ko = 0.0;
    double ke = 0.0;
    //loop through inner 50% of transect
    for(int i = m_nTransIdx; i < m_nTransIdx+m_nLenTrans; ++i) {
    	individual & ind = m_vPop2[i];
        if(ind.nWeight == 0)
    		continue;
    	szSample += 1;
    	alleleMap[ind.nAllele] += 1;
    	int & p = ind.nParent_id;
        if (m_sBound == "torus")
        {   
            dSigma2 += minEuclideanDist2(i,p,m_nMaxX,m_nMaxY);
            dSigma3 += minEuclideanDist3(i,p,m_nMaxX,m_nMaxY);
        }
        else
        {
            dSigma2 += euclideanDist2(i,p,m_nMaxX,m_nMaxY);
            dSigma3 += euclideanDist3(i,p,m_nMaxX,m_nMaxY);        
        }
 
        for(int j=i; j < m_nTransIdx+m_nLenTrans; ++j) {
            individual & ind2 = m_vPop2[j];
            if(ind2.nWeight == 0)
                continue;
            int k = j-i;
            //if(m_sBound == "torus")
                //***NOT NECESSARY WHEN USING 50% TRANSECT
                //k = (k <= m_nMaxX/2) ? k : m_nMaxX-k;
            //within individual
            if(k==0){
                //check for presence of competition allele
                if(ind.nWeight == 0)
                    continue;
                // count IIS
                if(ind.nAllele ==ind.nAllele2)
                    vIBD[k] += 1;
                // count grandparental IBD
                if(m_vPop1[p].nParent_id == m_vPop1[ind.nParent_id2].nParent_id)
                    vgIBD[k] += 1;
                // count parental IBD
                if(p == ind.nParent_id2)
                    vpIBD[k] += 1;
                vN[k] += 1;
            }
            //between individuals
            else{
                // count IIS
                if(ind.nAllele == ind2.nAllele)
                    vIBD[k] += 1;
                // count grandparental IBD
                if(m_vPop1[p].nParent_id == m_vPop1[ind2.nParent_id].nParent_id) 
                    vgIBD[k] += 1;
                // count parental IBD
                if(p == ind2.nParent_id)
                    vpIBD[k] += 1;
                vN[k] += 1;
            }
            
        }


    }
    int df = 0;
    BOOST_FOREACH(mapType::value_type & vv, alleleMap){
            df += vv.second*vv.second;
        }
    ko = (double)alleleMap.size();
    double f = df/(double)(szSample*szSample);
    ke = 1.0/f;
    //verbose option for checking burn-in
    if(verbose)
        cout << "Gen: " << gen << " Ko: " << ko << " Ke: " << ke << endl;
    demout << gen << "\t" << dSigma2/(2.0*szSample) <<"\t" << dSigma3/(2.0*szSample)<<"\t"<<ko<<"\t"<<ke<<"\t"<<f<<"\n";
    for(unsigned int k=0; k<vIBD.size();++k){
        dout << vIBD[k] << "/" << vN[k] << ((k < vIBD.size()-1) ? "\t" : "\n");
        gout << vgIBD[k] << "/" << vN[k] << ((k < vgIBD.size()-1) ? "\t" : "\n");
        iout << vpIBD[k] << "/" << vN[k] << ((k < vpIBD.size()-1) ? "\t" : "\n");
    }

}

void Population::samplePop(){
    for(int i = 0; i < m_nIndividuals; ++i) {
        if(m_vPop2[i].nWeight <= 0)
            popout << -1;
        else
            popout << m_vPop2[i].nAllele;
        popout << ((i<m_nIndividuals-1) ? "," : "\n");
    }
}

void Population::sampleDist(){
    disp.set_test(true);
    for(int i=m_nBoxIdxMin; i<m_nBoxIdxMax; i++){
        for(int j=m_nBoxIdxMin; j<m_nBoxIdxMax; j++){
            for(int k=0; k<m_nOffspring; k++){
                disp(m_myrand,i,j);
            }
        }
    }
    disp.set_test(false);
}

void Population::sampleNb(){
    vector<individual> pop(m_vPop2);
    random_shuffle(pop.begin(),pop.end());
    vector<int> vIBD(m_nDistClass,0);
    vector<int> vN(m_nDistClass,0);
    while(pop.size()){
        individual ind1 = pop[0];
        int x1 = ind1.xy.first;
        int y1 = ind1.xy.second;
        unsigned int i = 1;
        bool found = false;
        for(; i < pop.size(); i++){
            individual ind2 = pop[i];
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
            if(ind1.nAllele == ind2.nAllele){
                vIBD[d-1] += 1;
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
    for(unsigned int k=0; k<vIBD.size();++k)
        nbout << vIBD[k] << ((k< vIBD.size()-1) ? "\t" : "\n");
}