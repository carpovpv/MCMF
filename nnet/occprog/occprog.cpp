
/*
   The program for virtual screening based on one-class classification
   approach. A diabolo networks are used to predict the probability of
   a structure under consideration to be active.

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/fingerprint.h>

#include "../Net.h"

using namespace OpenBabel;
using namespace std;
using namespace NeuralNetwork;

double tanimoto(const double *a, const double  *b, unsigned int N = 1024 )
{
    double p = 0.0, q=0, r =0;
    unsigned int i;
    for (i=0; i<N; i++)
    {
        if (a[i] >=0.5 && b[i]>= 0.5) p++;
        if (a[i] <0.5 && b[i] >= 0.5) q++;
        if (a[i] >= 0.5 && b[i] < 0.5) r++;
    }
    return p/ (p+q+r);
}

int main(int argc , char **argv)
{
    if (argc == 1 )
    {
        cerr << "Usage: ./occ-diab model1 model2 model3." << endl;
        cerr << "Structures to be predicted are taken from the input stdin. " << endl;
        return EXIT_FAILURE;
    }

    OBConversion conv(&cin);
    conv.SetInFormat("SDF");

    OBMol mol;
    OBFingerprint *ob  = OBFingerprint::FindFingerprint("FP2");

    const int ndescr = 1024;
    unsigned int one =  1;
    one = one  << (8 * sizeof(unsigned int) -1) ;

    double descrs[ndescr];
    double res[ndescr];

    vector<unsigned int> fp;
    vector<Net *> nets;

    for (int i=1; i< argc; ++i)
    {
        ifstream in (argv[i], std::ios::binary);
        if (!in)
        {
            cerr << "Can't open file " << argv[i] << " for reading. " << endl;
            return EXIT_FAILURE;
        }

        Net *net = new Net (in);
        if (!net)
        {
            cerr << "Error in read model " << argv[i] << "." << endl;
            return EXIT_FAILURE;
        }

	nets.push_back(net);
	in.close();
    }

    const int nnsize = nets.size();
    unsigned int molnum = 0;

    while (conv.Read(&mol))
    {

        fp.clear();
        ob->GetFingerprint(&mol,fp, ndescr);
        int N = 0;

        for (int i=0; i< fp.size(); ++i)
        {
            unsigned int temp = fp[i];
            for (int j=0; j< (8 * sizeof(unsigned int)); ++j)
            {
                descrs[N++] = (temp & one) ? 1.0 : 0.0;
                temp = temp << 1  ;
            }
        }

	double av = 0.0;
	for(int n = 0; n< nnsize; ++n)
	{
        	nets[n]->run(descrs, res);
	        double tan = tanimoto(descrs, res);
                av += tan;       
	}

	printf("%5d %f\n", molnum++, av/nnsize);

    }
  
    return EXIT_SUCCESS;
}

