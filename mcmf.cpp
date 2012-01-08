
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

#ifdef WIN32
#include <windows.h>
#endif

void help()
{
    std::cout << "The program for building models for virtual screening \n"
              "based on continuous molecular field analysis (MCMF). \n";
}

DescriptorFactory * createDescriptor(int code)
{
    switch(code)
    {
    case D_MNA:
        return new MnaDescr();
    case D_FP2:
        return new FingerPrints2s();
    case D_SPECTROPHORES:
        return new Spectrophores();
    case D_UNKNOWNDESCR:
        return NULL;
    }
    std::cerr << "Unknown descriptor block!" << std::endl;
    return NULL;

}

CKernel * createKernel(int kernel, DescriptorFactory * descr)
{
    switch(kernel)
    {
    case D_GAUSS:
        return new GaussKernel(descr);
    case D_LINEAR:
        return new LinearKernel(descr);
    case D_TANIMOTO:
        return new TanimotoKernel(descr);
    case D_ELECTROSTATIC:
        return new ElectroStaticKernel();
    case D_STERIC:
        return new StericKernel();
    case D_HYDROPHOBIC:
        return new HydrophobicKernel();
    default:
        std::cerr << "Unknown kernel block!\n";
    }
    return NULL;
}

int do_prognosis(const std::string &)
{
    //init global variables of OpenBabel
    OBConversion obconversion;

    double usep [MAX_PARAMS];

    for(int i =0; i< MAX_PARAMS; ++i)
        usep[i] = UNKNOWN_VALUE;

    if(cond.sdf_test.empty())
    {
        fprintf(stderr, "Please, supply sdf-test file!\n");
        return EX_USAGE;
    }

    if(cond.model.empty())
    {
        fprintf(stderr, "Please, supply a model file!\n");
        return EX_USAGE;
    }

    FILE * fp = fopen(cond.model.c_str(), "rb");
    if(fp == NULL)
    {
        fprintf(stderr, "Can't open model's file %s.\n", cond.model.c_str());
        return EX_USAGE;
    }

    //комментарии в файле модели начинаются с символа #
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
        fclose(fp);

        return EX_USAGE;
    }

    //в модельном файле структуры должны быть всегда в SDF
    //потому что нужны координаты атомов.

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

    boost::shared_ptr<SEAL> mols ( new SEAL(&sin, NULL)) ;

    mols->go();
    fprintf(stderr, "Loaded: %d of model's molecules.\n", mols->getNumberOfMolecules());

    N = 0;

    fscanf(fp, "Kernels: %d\n", &N);
    if(N == 0)
    {
        fprintf(stderr, "The number of kernels must be greater than zero.\n");
        fclose(fp);

        return EX_USAGE;
    }

    boost::shared_ptr<CMFA> cmfa (new CMFA());
    fprintf(stderr, "Kernels: %d\n", N);

    for ( int i = 0; i < N; ++i)
    {
        fgets(str, 255, fp);
        int descrcode;

        sscanf(str,"Descriptors: %d", &descrcode);

        DescriptorFactory * descr = createDescriptor(descrcode);
        if(descr)
        {
            descr->load(fp);
            descr->setPrognosis(true);
        }

        fgets(str, 255, fp);
        int kerncode;

        sscanf(str,"Kernel: %d", &kerncode);
        CKernel * kernel = createKernel(kerncode, descr);

        cmfa->addKernel(kernel);
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

    fprintf(stderr, "Parameters: ");
    for(int i =0; i< N; ++i)
        fprintf(stderr, "%g ", usep[i]);
    fprintf(stderr,"\n");

    fscanf(fp, "\nMachine: %s\n", str);
    cond.machine = str;

    boost::shared_ptr<Machine> machine;
    if(cond.machine == "1-SVM")
        machine = boost::shared_ptr<Machine> (new OneClassSVM());
    else if(cond.machine == "SVR")
        machine = boost::shared_ptr<Machine> (new Svr());
    else
    {
        std::cerr << "Unknown machine! Available are: 1-svm, svr.\n" << std::endl;
        fclose(fp);
        return EX_USAGE;
    }

    machine->load(fp);

    machine->setCMFA(cmfa.get());
    machine->init();
    machine->setData(mols.get(), NULL);
    machine->setParameters(usep);

    fclose(fp);

    //do prognosis

    std::ifstream ifs(cond.sdf_test.c_str());

    obconversion.SetInFormat("SDF");
    obconversion.SetInStream(&ifs);

    OBMol mol;

    long struc = 0;

    while(obconversion.Read(&mol))
    {

        OBAtom *atom;
        FOR_ATOMS_OF_MOL(atom, mol)
        {
            Fields * ff = new Fields();
            ff->calcValues(&*atom);
            atom->SetData(ff);
        }

        //ставим метку, какую по счёту мы обрабатываем молекулу.
        //по ней дескрипторные вектора определяют надо вычислять
        //дескрипторы или можно взять уже вычисленные.

        char st[100];
        sprintf(st, "%ld", struc++);

        OBPairData  * data = new OBPairData();

        data->SetAttribute("prognosis");
        data->SetValue(st);
        data->SetOrigin(userInput);

        mol.SetData(data);

        machine->predict(&mol);
    }

    machine->after_prognosis();

    return 0;
}

int do_model(const std::string & params)
{
    //init global variables of OpenBabel
    OBConversion obconversion;

    boost::shared_ptr<Machine> machine;
    if(cond.machine == "1-svm")
        machine = boost::shared_ptr<Machine> (new OneClassSVM());
    else if(cond.machine == "svr")
        machine = boost::shared_ptr<Machine> (new Svr());
    else
    {
        std::cerr << "Unknown machine! Available are: 1-svm, svr.\n" << std::endl;
        return EX_USAGE;
    }

    if(cond.kernels.size() == 0)
    {
        std::cerr << "Please, select at least one kernel!\n" << std::endl;
        return EX_USAGE;
    }

    boost::shared_ptr<CMFA> cmfa (new CMFA());

    for(int i=0; i<  cond.kernels.size(); ++i)
    {
        DescriptorFactory * descr = NULL;
        try
        {
            descr = createDescriptor(cond.kernels[i].descr);
        } catch( DescrFailed &f)
        {
            std::cerr << f.what() << std::endl;
            return EX_USAGE;
        }

        CKernel * kernel = NULL;
        try
        {
            kernel = createKernel(cond.kernels[i].kernel, descr);
        }
        catch(KernelFailed & k)
        {
            std::cerr << k.what() << std::endl;
            if(descr)
                delete descr;

            return EX_USAGE;
        }

        cmfa->addKernel(kernel);
    }

    if(cmfa->count() == 0)
    {
        std::cerr <<  "Please, select at least one kernel!\n" << std::endl;
        return EX_USAGE;
    }

    machine->setCMFA(cmfa.get());

    if(cond.sdf_train.empty())
    {
        std::cerr <<  "Sdf-train filename is empty!\n" << std::endl;
        return EX_USAGE;
    }

    std::cerr <<  "Sdf-test: " << cond.sdf_test << std::endl;
    std::cerr <<  "Sdf-train: "<< cond.sdf_train << std::endl;

    srand(time(NULL));

    if(cond.results.empty())
    {
        std::cerr << "Please, select a result file!" << std::endl;
        return EX_USAGE;
    }

    boost::shared_ptr<SEAL> train;
    boost::shared_ptr<SEAL> test;

    int max_iter = cond.max_iter;
    int cv = cond.cv > 0 ? cond.cv : 10;

    std::vector< std::string> props;

    double usep [MAX_PARAMS];

    for(int i =0; i< MAX_PARAMS; ++i)
        usep[i] = UNKNOWN_VALUE;

    for(int i=0; i< cond.params.size() && i< MAX_PARAMS; ++i)
        usep[i] = cond.params[i];

    if(cond.results.empty())
    {
        fprintf(stderr,"Select a result file please.\n");
        return EX_USAGE;
    }

    if(machine == NULL)
    {
        fprintf(stderr,"Please select a machine.\n");
        return EX_USAGE;
    }

    try
    {
        train = boost::shared_ptr<SEAL> (new SEAL(cond.sdf_train.c_str(), &props) );
    }
    catch(std::exception &exc)
    {
        fprintf(stderr, "%s\n", exc.what());
        return EXIT_FAILURE;
    }

    if(!cond.sdf_test.empty())
    {
        try
        {
            test = boost::shared_ptr<SEAL> (new SEAL(cond.sdf_test.c_str()));
        }
        catch(std::exception &exc)
        {
            fprintf(stderr, "%s\n", exc.what());
            return EXIT_FAILURE;
        }
    }

    train->go();
    if(test.get()) test->go();

    const int N = train->getNumberOfMolecules();
    if(!test)
        printf("Number of structures in the loaded database: %d.\n", N);
    else
        printf("Number of structures in the loaded database: %d (ligands) %d (decoys).\n", N, test->getNumberOfMolecules());

    if(!machine->setData(train.get(), test.get()))
        return EX_USAGE;

    FILE * fres = fopen(cond.results.c_str(), "w");
    fprintf(fres, "MCMF RESULT FILE.\n");
    fprintf(fres, "Command: %s\n", params.c_str());

    char outstr[200];
    time_t t;
    struct tm *tmp;

    t = time(NULL);
    tmp = localtime(&t);
    strftime(outstr, sizeof(outstr), "%d/%m/%Y %H:%M:%S", tmp);

    fprintf(fres, "STARTED AT: %s\n\n", outstr);
    fputs("------------------------------------------------\n", fres);
    fputs("OPTIMIZATION\n",fres);
    fputs("------------------------------------------------\n\n", fres);
    fputs("Parameters => Min & Max Values of Calculated Results => Quality\n\n", fres);

    machine->setOutput(fres);
    machine->set_CV(cv);
    machine->init();

    machine->setParameters(usep);
    machine->create_random(max_iter);

    machine->save(cond.model.c_str());

    t = time(NULL);
    tmp = localtime(&t);
    strftime(outstr, sizeof(outstr), "%d/%m/%Y %H:%M:%S", tmp);

    fprintf(fres, "\nFINISHED AT: %s\n\n", outstr);
    fprintf(fres, "THANK YOU FOR USING THE PROGRAM!\n");

    if(fres != NULL)
        fclose(fres);

    return 0;
}

int main(int argc, char ** argv)
{

    if(argc == 1)
    {
        help();
        return 0;
    }

    std::string params;
    for(int i=1; i< argc; ++i)
    {
        params += argv[i];
        params += " ";
    }

    std::cerr << "Program launched with command: " << params << std::endl;

    //Launch Yacc to parse command line.
    if(!parse_command_line(params.c_str()))
        return 0;
     
    #ifdef WIN32
   
    const int MAX_MODULE_NAME = 1024;
    char modulename[MAX_MODULE_NAME];

    int len = GetModuleFileName(NULL, modulename, 1024);
    if(len == 0)
    {
        fprintf(stderr, "Cannot determine the module path.\n");
	return EX_USAGE;
    }

    puts("Set env variables for OpenBabel.");

    char * p = strrchr(modulename, '\\');
    if( p == NULL )
    {
        fprintf(stderr, "Cannot analyze module's name.\n");
	return EX_USAGE;
    }

    *++p = '\0';
    
    len = p - modulename;

    char env_babel_lib[MAX_MODULE_NAME + 30];
    sprintf(env_babel_lib, "BABEL_LIBDIR=%s/plugins",modulename);
    putenv(env_babel_lib);

    char env_babel_data[MAX_MODULE_NAME + 30];
    sprintf(env_babel_data, "BABEL_DATADIR=%s/data",modulename);
    putenv(env_babel_data);
   

    #endif

    if(cond.help)
    {
        help();
        return 0;
    }

    return cond.prognosis ? do_prognosis(params)
           : do_model(params);
}

