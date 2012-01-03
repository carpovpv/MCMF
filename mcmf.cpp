
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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "seal/seal.h"
#include "fields.h"

#include "cmfa.h"
#include "kernels/gauss.h"
#include "kernels/tanimoto.h"
#include "kernels/electro.h"
#include "kernels/hydropho.h"
#include "kernels/steric.h"
#include "kernels/linear.h"
#include "kernels/hydrophov.h"
#include "kernels/sterick.h"
#include "kernels/abraham.h"

#include "descfact.h"
#include "descrs/fp2s.h"
#include "descrs/spectrophores.h"
#include "descrs/mnadescr.h"

#include "machine.h"
#include "machines/oneclasssvm.h"
#include "machines/svr.h"

#include "parser.h"

#include <errno.h>
#include <boost/shared_ptr.hpp>

const int EX_USAGE = 127;
const int MAX_PARAMS = 20;

void help()
{
    printf("The program for building models for virtual screening \n"
           "based on continuous molecular field analysis (MCMF). \n");
}

int main(int argc, char ** argv)
{

    parse_command_line("");

    std::cout << cond.machine  << std::endl;

    srand(time(NULL));

    SEAL * train = NULL;
    SEAL * test = NULL;

    FILE * fres = NULL;

    char * sdf_test = NULL;  //decoys
    char * sdf_train = NULL; //ligands

    int do_help = 0;
    int do_prognosis = 0;
    int max_iter = 3;
    int cv = 10;

    char * save_model = NULL;
    char * file_res = NULL;

    std::vector< std::string> props;

    struct option longopts[] = {
        {"sdf-train", required_argument, NULL, 't'},
        {"sdf-test", optional_argument, NULL, 'v'},
        {"help", no_argument, & do_help, 1},
        {"model", required_argument, NULL , 'm'},
        {"results", required_argument, NULL, 'r'},
        {"h", optional_argument, NULL, 'h'},
        {"kernel", required_argument, NULL, 'k'},
        {"machine", required_argument, NULL, 'a'},
        {"property", optional_argument, NULL, 'p'},
        {"prognosis", optional_argument, & do_prognosis, 1},
        {"max-iter", optional_argument, NULL, 'i'},
        {"cv", optional_argument, NULL, 'c'},
        {0, 0, 0, 0}
    };

    //for init global variables of OpenBabel
    OBConversion obconversion;
    Machine * machine = NULL;

    //default descriptor block
    boost::shared_ptr<FingerPrints2s> fp2s( new FingerPrints2s());
    boost::shared_ptr<Spectrophores> Spectr( new Spectrophores());
    boost::shared_ptr<MnaDescr> mna( new MnaDescr());

    boost::shared_ptr<CMFA> cmfa (new CMFA());

    boost::shared_ptr<ElectroStaticKernel> electro (new ElectroStaticKernel());
    boost::shared_ptr<HydrophobicKernel> hydrophobic (new HydrophobicKernel());
    boost::shared_ptr<StericKernel> steric (new StericKernel());
    boost::shared_ptr<LinearKernel> linear (new LinearKernel(fp2s.get()));

    boost::shared_ptr<HydrophobicKernelV> hydrophobicv (new HydrophobicKernelV());
    boost::shared_ptr<StericKernelK> sterick (new StericKernelK());

    boost::shared_ptr<AbrahamKernelA> abrahama (new AbrahamKernelA());
    boost::shared_ptr<AbrahamKernelB> abrahamb (new AbrahamKernelB());
    boost::shared_ptr<AbrahamKernelE> abrahame (new AbrahamKernelE());
    boost::shared_ptr<AbrahamKernelS> abrahams (new AbrahamKernelS());

    boost::shared_ptr<GaussKernel> gauss( new GaussKernel(mna.get()));
    boost::shared_ptr<GaussKernel> gaussSpect( new GaussKernel(Spectr.get()));
    boost::shared_ptr<TanimotoKernel > tanimoto (new TanimotoKernel(fp2s.get()));

    double usep [MAX_PARAMS];
    double * userp = usep;

    for(int i =0; i< MAX_PARAMS; ++i)
        usep[i] = UNKNOWN_VALUE;


    ////start old parsing
    int c;
    while( (c = getopt_long(argc, argv, "", longopts, NULL)) != -1)
    {
        double k;
        char *p, *s;
        char prop[1024];
        int i =0,l;

        switch(c)
        {
        case 'i':
            max_iter=atoi(optarg);
            break;
        case 'c':
            cv=atoi(optarg);
            break;
        case 'p':
            p = optarg;
            do
            {
                i = 0;
                while(*p!='\0' && *p != ',')
                    prop[i++] = *p++;
                prop[i]='\0';
                props.push_back(prop);

            } while(*p++!='\0');

            printf("Properties:\n");
            for(i=0; i< props.size(); ++i)
                printf("\t%s\n", props[i].c_str());

            break;
        case 'k':
            sprintf(prop,"%s", optarg);
            s = p = prop;
            l=0;
            while(true)
            {
                if(*p == ',' || *p == '\0')
                {
                    if(*p=='\0') l = 1;

                    *p='\0';
                    p++;
                    if(!strcmp(s,"electrostatic"))
                        cmfa->addKernel(electro.get());
                    else if(!strcmp(s, "hydrophobic"))
                        cmfa->addKernel(hydrophobic.get());
                    else if(!strcmp(s, "steric"))
                        cmfa->addKernel(steric.get());
                    else if(!strcmp(s, "linear"))
                        cmfa->addKernel(linear.get());
                    else if(!strcmp(s, "gauss"))
                        cmfa->addKernel(gauss.get());
                    else if(!strcmp(s, "tanimoto"))
                        cmfa->addKernel(tanimoto.get());
                    else if(!strcmp(s, "hydrophobicv"))
                        cmfa->addKernel(hydrophobicv.get());
                    else if(!strcmp(s, "sterick"))
                        cmfa->addKernel(sterick.get());
                    else if(!strcmp(s,"abrahama"))
                        cmfa->addKernel(abrahama.get());
                    else if(!strcmp(s,"abrahamb"))
                        cmfa->addKernel(abrahamb.get());
                    else if(!strcmp(s,"abrahame"))
                        cmfa->addKernel(abrahame.get());
                    else if(!strcmp(s,"abrahams"))
                        cmfa->addKernel(abrahams.get());
                    else if(!strcmp(s,"gaussspectr"))
                        cmfa->addKernel(gaussSpect.get());
                    else
                    {
                        fprintf(stderr,"Unknown kernel %s.\n", s);
                        return EX_USAGE;
                    }
                    s = p;
                    if(l)
                        break;
                }
                else
                    p++;
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
            else if(!strcmp(optarg, "svr"))
            {
                if(machine != NULL)
                {
                    fprintf(stderr,"Several machines are not supported.\n");
                    return EX_USAGE;
                }

                Svr * svr = new Svr();
                svr->setCMFA(cmfa.get());
                svr->setProps(&props);
                machine = svr;
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


    ///// end parsing command line

    if(do_help)
    {
        help();
        return 0;
    }

    if(do_prognosis)
    {
        printf("Prognosis mode...\n");
        FILE * fp = fopen(save_model, "rb");
        if(fp == NULL)
        {
            fprintf(stderr, "Can't open file %s.\n", save_model);
            return 1;
        }

        char str[255];
        do
        {
            fgets(str, 255, fp);
        } while(str[0] == '#');

        std::stringstream sin;

        int N = 0;
        sscanf(str, "Structures: %d\n", &N);
        if(N == 0)
        {
            fprintf(stderr, "Number of structures in the model must be greater than zero.\n");
            return 1;
        }


        std::string buf;
        for(int i =0; i< N && !feof(fp); i++)
        {

            do {
                fgets(str, 255, fp);
                buf += str;
            } while(strncmp("$$$$", str, 4));
        }

        sin.str(buf);
        obconversion.SetInStream(&sin);
        obconversion.SetInFormat("sdf");

        SEAL * mols = new SEAL(&sin, &props);

        mols->go();
        printf("Loaded: %d of model's molecules.\n", mols->getNumberOfMolecules());

        N = 0;

        fscanf(fp, "Kernels: %d\n", &N);
        if(N == 0)
        {
            fprintf(stderr, "The number of kernels must be greater than zero.\n");
            return 1;
        }

        cmfa->clear();

        printf("Kernels: %d\n", N);
        for ( int i = 0; i < N; ++i)
        {
            char kernel[100], descr[100];
            fgets(str, 255, fp);

            descr[0] = '\0';
            char * p = strrchr(str,':');
            if(p)
            {
                sprintf(descr,"%s", p+1);
                descr[strlen(descr) - 1] = '\0';
                *p='\0';

            }
            else
                str[strlen(str) - 1 ] = '\0';

            sprintf(kernel,"%s", str);
            printf("Descriptors: %s. Kernel: %s.\n", descr, kernel);

            DescriptorFactory * descrfact;

            if(!strcmp("Spectrophores", descr))
                descrfact = Spectr.get();
            else if(!strcmp("FP2", descr))
                descrfact = fp2s.get();

            if(!strcmp(kernel, "Gaussian"))
            {
                printf("Loading gauss %s\n", descrfact->getName().c_str());
                cmfa->addKernel(new GaussKernel(descrfact));
            }
            else if(!strcmp(kernel, "Tanimoto"))
            {
                printf("Loading tanimoto %s\n", descrfact->getName().c_str());
                cmfa->addKernel(new TanimotoKernel(descrfact));
            }
            else if(!strcmp(kernel, "Electro-Static"))
            {
                printf("loading Electro-Static kernel.\n");
                cmfa->addKernel(new ElectroStaticKernel());
            }
            else if(!strcmp(kernel, "Steric"))
            {
                printf("loading Steric kernel.\n");
                cmfa->addKernel(new StericKernel());
            }
            else if(!strcmp(kernel, "Hydrophobic"))
            {
                printf("loading Hydrophobic kernel.\n");
                cmfa->addKernel(new HydrophobicKernel());
            }
            else if(!strcmp(kernel, "Linear"))
            {
                printf("loading Linear kernel.\n");
                cmfa->addKernel(new LinearKernel(descrfact));
            }
        }

        N= 0;
        fscanf(fp, "Parameters: %d\n", &N);
        if(N == 0)
        {
            fprintf(stderr, "The number of parameters must be greater than zero.\n");
            return 1;
        }

        for(int i=0; i< N; ++i)
        {
            fscanf(fp, "%lf", &usep[i]);
        }

        cmfa->setParameters(usep);

        for(int i =0; i< N; ++i)
            printf("%g ", usep[i]);
        printf("\n");

        fscanf(fp, "\nMachine: %s\n", str);
        printf("Machine: %s\n", str);

        Machine * machine;

        if(!strcmp("1-SVM", str))
        {
            OneClassSVM * svm = new OneClassSVM();

            if(svm->load(fp))
            {
                printf("1-SVM loaded.\n");
                svm->setCMFA(cmfa.get());
                svm->init();
                svm->setData(mols,NULL);
                svm->setParameters(usep);
                cmfa->setParameters(usep);
                machine = svm;
            }
        }
        else if(!strcmp("SVR", str))
        {
            Svr * svm = new Svr();

            if(svm->load(fp))
            {
                printf("SVR loaded.\n");
                svm->setCMFA(cmfa.get());
                svm->init();
                svm->setData(NULL, mols);
                svm->setParameters(usep);
                cmfa->setParameters(usep);
                machine = svm;
            }
        }
        fclose(fp);

        //model read

        std::ifstream ifs(sdf_test);

        obconversion.SetInFormat("SDF");
        obconversion.SetInStream(&ifs);

        OBMol mol;

        while(obconversion.Read(&mol))
        {

            OBAtom *atom;
            FOR_ATOMS_OF_MOL(atom, mol)
            {
                Fields * ff = new Fields();
                ff->calcValues(&*atom);
                atom->SetData(ff);
            }

            machine->predict(&mol);            
        }

        return 0;
    }

    //building a model

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
        train = new SEAL(sdf_train, &props);
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
    machine->set_CV(cv);
    machine->init();

    machine->setParameters(usep);

    //machine->create();
    machine->create_random(max_iter);

    machine->save(save_model);

    if(fres != NULL) fclose(fres);

    delete machine;

    return 0;
}

