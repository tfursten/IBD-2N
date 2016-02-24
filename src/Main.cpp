#include "Main.h"


int main(int ac, char** av)
{
	namespace po = boost::program_options;
    using namespace std;

    static int nGenerations, nMaxX, nMaxY, nOffspring, nBurnIn, nSample, ndClass, nPairs, nMarkers, nAlleles;
    unsigned int seed;
    //static double dMut;
    static vector<double> vdMut;
    static double dSigma;
    static float param;
    bool f;
    string dist_name, infile, outfileName, mut_type;
    bool verbose;

    ostringstream out;
    try
    {
        po::options_description generic("General Options");
        generic.add_options()
            ("help", "Produce help message")
            ;

        po::options_description config("Configuration");
        config.add_options()
            ("maxX,x", po::value<int>(&nMaxX)->default_value(100),"Set X dimension")
            ("maxY,y", po::value<int>(&nMaxY)->default_value(100),"Set Y dimension")
            ("generations,g", po::value<int>(&nGenerations)->default_value(10), "Set number of Generations to run after burn-in")
            ("offspring,o", po::value<int>(&nOffspring)->default_value(10), "Set number of offspring per individual")
            ("markers,m", po::value<int>(&nMarkers)->default_value(1), "Set number of markers")
            ("mut", po::value<vector<double>>(&vdMut)->multitoken()->default_value(vector<double>(1,0.0001),"0.0001"), "Set mutation rates")
            ("distribution,d", po::value<string>(&dist_name)->default_value("triangular"), "Set Dispersal Distribution")
            ("sigma,s", po::value<double>(&dSigma)->default_value(2.0), "Set dispersal parameter")
            ("burn,b", po::value<int>(&nBurnIn)->default_value(0),"Set Burn-in Period")
            ("sample,t", po::value<int>(&nSample)->default_value(1),"Sample every n generations after burn-in")
            ("output_file,f", po::value<string>(&outfileName)->default_value(string("data")),"Output File Name")
            ("seed", po::value<unsigned int>(&seed)->default_value(0), "Set PRNG seed, 0 to create random seed")
            ("verbose", po::value<bool>(&verbose)->default_value(false),"Print data to screen")
            ("sparam", po::value<float>(&param)->default_value(0),"Extra Parameter for dispersal")
            ("fast", po::value<bool>(&f)->default_value(true),"Use fast dispersal when available")
            ("ndistClass", po::value<int>(&ndClass)->default_value(20),"Number of distance classes for Nb estimate")
            ("nPairs", po::value<int>(&nPairs)->default_value(20),"Number of pairs for Nb estimate")
            ("mut-type", po::value<string>(&mut_type)->default_value(string("IAM")),"Mutation Model (IAM, KAM or SMM)")
            ("nAllele", po::value<int>(&nAlleles)->default_value(20),"Number of alleles under SMM and KAM")
            ;

        po::options_description hidden("Hidden Options");
        hidden.add_options()
            ("input-file", po::value<string>(&infile), "input file")
            ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config);

        po::options_description visible("Allowed Options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac,av).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);


        if (vm.count("help"))
        {
            cout << visible << "\n";
            return 0;
        }


        if (!infile.empty())
        {
            ifstream ifs(infile.c_str());
            if (!ifs)
            {
                cout << "can not open config file: "<< infile << "\n";
                return 0;
            }
            else
            {
                po::store(parse_config_file(ifs, config_file_options), vm);
                po::notify(vm);
            }
        }
        if((int)(vdMut.size())!=nMarkers){
            if(int(vdMut.size()) == 1){
                vdMut.resize(nMarkers);
                fill(vdMut.begin(),vdMut.end(),vdMut[0]);
            }
            else{
                throw;
            }
        if(2*nPairs*ndClass > nMaxX*nMaxY){
            cout << "Sample size is larger than population size" << endl;
            throw;
        }
        if(mut_type != "IAM" && mut_type != "SMM" && mut_type != "KAM"){
            cout << "Not a valid mutation model" << endl;
            throw;
        }
            
        }
        out << "X dimension set to " << nMaxX << ".\n"
        << "Y dimension set to " << nMaxX << ".\n"
        << "Run for " << nGenerations << " generations.\n"
        << "Burn " << nBurnIn << " generation(s).\n"
        << "Collect data every " << nSample << " Generation(s).\n"
        << "Number of Offspring set to " << nOffspring << ".\n"
        << "Number of markers set to " << nMarkers << ".\n"
        << "Dispersal parameter set to " << dSigma << ".\n"
        << "Number of distances classes for Nb estimate set to " << ndClass << ".\n"
        << "Number of pairs collected for Nb estimate set to " << nPairs << ".\n"
        << "Mutation model set to " << mut_type << ".\n";
        if(mut_type == "SMM" || mut_type == "KAM")
            out << "Number of alleles set to " << nAlleles << ".\n";
        out << "Mutation rate(s) set to ";
        for(auto i=vdMut.begin();i!=vdMut.end();++i){
            out << *i << " ";
        }
        out<<".\n";


    }

    catch(exception& e)
    {
        cout<< e.what() << "\n";
        return 1;
    }
    nMaxY = nMaxX; //override maxY, landscape needs to be square at the moment
    string param_file = outfileName+"_settings.txt";
    string nb_data = outfileName+"_nb.txt";
    cout << "Parameters saved to: " << param_file << endl;
    cout << "Nb size estimation data saved to: " << nb_data << endl;
    ofstream pout;
    ofstream nbout;
    pout.open(param_file);
    nbout.open(nb_data);
    pout << out.str();
    cout << out.str();
	//Initialize Population
    clock_t start = clock();
    Population pop(pout, nbout, verbose);
	pop.initialize(nMaxX,nMaxY,nOffspring,nMarkers,dSigma,vdMut,seed, nSample,\
    dist_name, param, f, ndClass, nPairs, mut_type, nAlleles);
	//Run Simulation
	pop.evolve(nBurnIn, nGenerations);
	clock_t end = clock();
	float seconds = (float)(end-start)/ CLOCKS_PER_SEC;
    cout << "TIME: " << seconds << endl;
    pout << "TIME: " << seconds << endl;
    pout.close();
    nbout.close();



	return 0;

}

