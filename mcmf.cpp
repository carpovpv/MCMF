
/*
 * mcmf.cpp
 * Copyright (C) Carpov Pavel   2010 <carpovpv@qsar.chem.msu.ru>
                 Baskin Igor I. 2010 <igbaskin@gmail.com>
 *
 * MCMF is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MCMF is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "seal/seal.h"
#include "seal/ks_seal.h"

#include "cmfa.h"
#include "kernels/gauss.h"
#include "kernels/tanimoto.h"
#include "descfact.h"
#include "kernels/electro.h"
#include "descrs/fp2.h"
#include "kernels/hydropho.h"
#include "kernels/steric.h"

#include "machine.h"
#include "machines/oneclasssvm.h"

#include <errno.h>
#include <boost/shared_ptr.hpp>

const int EX_USAGE = 127;
const int MAX_PARAMS = 20;

void help()
{
    printf("The progrma for building models for virtual screening \n"
           "based on continuous molecular field analysis (MCMF). \n");
}

int main(int argc, char ** argv)
{

    SEAL * train = NULL, * test = NULL;
    FILE *fres = NULL;

    char * sdf_test = NULL;
    char * sdf_train = NULL; //decoys and ligands

    int do_help = 0;
    char * save_model = NULL;
    char * file_res = NULL;

    struct option longopts[] = {
        {"sdf-train", required_argument, NULL, 't'},
        {"sdf-test", optional_argument, NULL, 'v'},
        {"help", no_argument, & do_help, 1},
        {"model", required_argument, NULL , 'm'},
        {"results", required_argument, NULL, 'r'},
        {"h", optional_argument, NULL, 'h'},
        {"kernel", required_argument, NULL, 'k'},
        {"machine", required_argument, NULL, 'a'},
        {0, 0, 0, 0}
    };

    //for init global variables of OpenBabel
    OBConversion obconversion;
    Machine * machine = NULL;

    //default descriptor block
    boost::shared_ptr<FingerPrints2> fp2( new FingerPrints2());

    boost::shared_ptr<CMFA> cmfa (new CMFA());

    boost::shared_ptr<ElectroStaticKernel> electro (new ElectroStaticKernel());
    boost::shared_ptr<HydrophobicKernel> hydrophobic (new HydrophobicKernel());
    boost::shared_ptr<StericKernel> steric (new StericKernel());

    boost::shared_ptr<GaussKernel> gauss( new GaussKernel(fp2.get()));
    boost::shared_ptr<TanimotoKernel > tanimoto (new TanimotoKernel(fp2.get()));

    double usep [MAX_PARAMS];
    double * userp = usep;

    usep[0] = 0.04;
    for(int i =1; i< MAX_PARAMS; ++i)
        usep[i] = (i%2) ? 0.001 : 0.2 ;

    int c;
    while( (c = getopt_long(argc, argv, "a:t:v:m:r:k:", longopts, NULL)) != -1)
    {
        double k;
        char *p;
        switch(c)
        {
        case 'k':
            if(!strcmp(optarg,"electrostatic"))
                cmfa->addKernel(electro.get());
            else if(!strcmp(optarg, "hydrophobic"))
                cmfa->addKernel(hydrophobic.get());
            else if(!strcmp(optarg, "steric"))
                cmfa->addKernel(steric.get());
            else if(!strcmp(optarg, "mcmf"))
            {
                cmfa->addKernel(electro.get());
                cmfa->addKernel(hydrophobic.get());
                cmfa->addKernel(steric.get());
            }
            else if(!strcmp(optarg, "gauss"))
                cmfa->addKernel(gauss.get());
            else if(!strcmp(optarg, "tanimoto"))
                cmfa->addKernel(tanimoto.get());
            else
            {
                fprintf(stderr,"Unknown kernel.\n");
                return EX_USAGE;
            }
            break;
        case 'a':
            if(!strcmp(optarg, "1-svm"))
            {
                if(machine != NULL)
                {
                    fprintf(stderr,"Several machines are not supported.\n");
                    return EX_USAGE;
                }

                OneClassSVM * svm_1 = new OneClassSVM();
                svm_1->setCMFA(cmfa.get());
                machine = svm_1;
            }
            else
            {
                fprintf(stderr,"Unknown machine!\n");
                return EX_USAGE;
            }
            break;
        case 't':
            sdf_train = optarg;
            break;
        case 'v':
            sdf_test = optarg;
            break;
        case 'r':
            file_res = optarg;
            fres = fopen(file_res, "w");
            if(fres == NULL)
            {
                fprintf(stderr, "Unable to create the file %s\n", file_res);
                return EX_USAGE;
            }
            break;
        case 'm':
            save_model = optarg;
            break;
        case 'h':
            p = optarg;
            while (true)
            {
                sscanf (p, "%lf", &k);
                printf("h: %g\n", k);

                if(userp - usep > MAX_PARAMS)
                {
                    fprintf(stderr, "Max params exceeded!\n");
                    return EX_USAGE;
                }

                *userp++ = k;

                while (isdigit (*p))
                    p++;
                if (*p == ',')
                    p++;
                else if(*p == '.')
                {
                    p++;
                    while (isdigit(*p)) p++;
                    if(*p == '\0')
                        break;
                    p++;
                }
                else if (*p == '\0')
                    break;
            }

            break;
        case 0:
            break;
        default:
            return EX_USAGE;
            break;
        }
    }

    if(do_help)
    {
        help();
        return 0;
    }

    if(cmfa->count() == 0)
    {
        fprintf(stderr,"Select a kernel please.\n");
        return EX_USAGE;
    }

    if(file_res == NULL)
    {
        fprintf(stderr,"Select a result file please.\n");
        return EX_USAGE;
    }

    if(machine == NULL)
    {
        fprintf(stderr,"Please select a machine.\n");
        return EX_USAGE;
    }

    train = NULL;
    test = NULL;

    try
    {
        train = new SEAL(sdf_train);
    }
    catch(std::exception &exc)
    {
        fprintf(stderr, "%s\n", exc.what());
        return EXIT_FAILURE;
    }

    if(sdf_test)
    {
        try
        {
            test = new SEAL(sdf_test);
        }
        catch(std::exception &exc)
        {
            fprintf(stderr, "%s\n", exc.what());
            return EXIT_FAILURE;
        }
    }

    train->go();
    if(test) test->go();

    const int N = train->getNumberOfMolecules();
    if(!test)
        printf("Number of structures in the loaded database: %d.\n", N);
    else
        printf("Number of structures in the loaded database: %d (ligands) %d (decoys).\n", N, test->getNumberOfMolecules());

    if(!machine->setData(train, test))
        return EX_USAGE;

    machine->setOutput(fres);
    machine->setParameters(usep);
    machine->create();


    if(fres != NULL) fclose(fres);

    delete machine;

    return 0;
}

