
/*
 * svr.cpp
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

#include "svr.h"
#include <iostream>


Svr::Svr() : Machine("1-SVM")
{
    std::cout << "SVR is loaded." << std::endl;

    param.svm_type = NU_SVR;
    param.kernel_type = PRECOMPUTED;

    //we dont use these parameter. But the library's inteface expects it.

    param.degree = 3;
    param.gamma = 0;
    param.coef0 = 0;
    param.cache_size = 100;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;

    //?
    dimensionality = 1;

    m_NumParameters = 0;

    Parameters = NULL;
    lp = NULL;
    mp = NULL;

}

void Svr::setCMFA(CMFA *cmfa)
{
    m_cmfa = cmfa;
    m_cmfa->printSelKernels();

    if(lp != NULL)
    {
        free(lp);
        free(mp);
        free(Parameters);
    }

    m_NumParameters = cmfa->count() * 2 + 1;

    Parameters = (double * ) calloc(m_NumParameters , sizeof(double));
    lp = (double * ) calloc(m_NumParameters , sizeof(double));
    mp = (double * ) calloc(m_NumParameters , sizeof(double));

    lp[0] = 0.001, lp[1] = 1e-3;
    for(int i =2; i< m_NumParameters; ++i)
        lp[i] = (i%2) ? 0.001 : 0.0001;

    mp[0] = 0.800, mp[1] = 1e3;
    for(int i =2; i< m_NumParameters; ++i)
        mp[i] = (i%2) ? 1.0 : 0.1 ;

    Parameters[0] = 0.04;
    for(int i =2; i< m_NumParameters; ++i)
        Parameters[i] = (i%2) ? 0.01 : 0.01 ;

}

Svr::~Svr()
{
    if(lp != NULL)
    {
        free(lp);
        free(mp);
        free(Parameters);
    }
}

bool Svr::setData(SEAL *train_mols, SEAL *test_mols)
{
    train = train_mols;
    test = test_mols;

    if(test == NULL || train == NULL)
    {
        fprintf(stderr,"Both test and train sets must be supplied.\n");
        return false;
    }

    return true;
}

bool Svr::build(const double * params, const std::vector<int> &flags, const std::vector<int> &mask)
{

    param.nu = params[0];
    param.C = params[1];

    m_cmfa->setParameters(params + 2);

    int rn = 0;
    const int N = flags.size();

    for(int i = 0; i< N; ++i)
        if(flags[i]) rn++;

    problem.l = rn;

    problem.y = (double *) calloc(rn , sizeof(double));
    problem.x = (struct svm_node **) calloc(rn, sizeof(struct svm_node*));

    for(int i= 0; i< rn; ++i)
        problem.x[i] = (struct svm_node *) calloc(rn +2, sizeof(struct svm_node));

    int i_r = -1;

    for(int i =0; i< N; ++i)
    {
        if(flags[i] == 0) continue;

        i_r++;

        problem.y[i_r] = 1.0;

        problem.x[i_r][i_r+1].value = 1.0;
        problem.x[i_r][i_r+1].index = i_r+1;

        problem.x[i_r][0].value = i_r+1;
        problem.x[i_r][0].index = 0;

        problem.x[i_r][rn+1].value = 0.0;
        problem.x[i_r][rn+1].index = -1;

        int j_r = i_r;
        for(int j=i+1; j< N; ++j)
        {
            if(flags[j] == 0) continue;

            OBMol *mol1 = train->getMolecule(mask[i]);
            OBMol *mol2 = train->getMolecule(mask[j]);

            double K = m_cmfa->calculate(mol1, mol2);

            j_r++;
            problem.x[i_r][j_r+1].value = problem.x[j_r][i_r+1].value = K;

            problem.x[i_r][j_r+1].index = j_r+1;
            problem.x[j_r][i_r+1].index = i_r+1;

        }
    }

    model = svm_train(&problem, &param);

    //predict +1

    struct svm_node * testing = (struct svm_node *) calloc(rn+2 , sizeof(struct svm_node));

    for(int i =0; i< N; ++i)
    {
        if(flags[i] == 1) continue;

        testing[0].value = 1;
        testing[0].index = 0;

        testing[rn+1].value = 0.0;
        testing[rn+1].index = -1;

        OBMol *mol1 = train->getMolecule(mask[i]);

        int i_p = 0;
        for(int j=0; j< N; ++j)
        {
            if(flags[j] == 0) continue;

            OBMol *mol2 = train->getMolecule(mask[j]);
            double K = m_cmfa->calculate(mol1, mol2);

            testing[i_p+1].value = K;
            testing[i_p+1].index = i_p+1;

            i_p++;

        }

        struct result * res = create_result();
        res->y_pred[0] = svm_predict(model, testing);
        res->y_real[0] = 1.0;

        results.push_back(res);

    }

    free(testing);

    for(int i =0; i< rn; ++i)
    {
        free(problem.x[i]);
        problem.x[i] = NULL;
    }

    free(problem.x);
    problem.x = NULL;

    free(problem.y);
    problem.y = NULL;

    svm_free_and_destroy_model(&model);
    svm_destroy_param(&param);

}

bool Svr::load(std::string & filename)
{

}

bool Svr::save(std::string &filename)
{

}

double Svr::statistic()
{
    /*
     One-class classification differs from other CV-methods.
     We calculate here predictions for test decoys structrues;
     */
    predict_decoys();


    return 1.0;
}

void Svr::predict_decoys()
{
    param.nu = Parameters[0];

    int rn = 0;
    const int N = train->getNumberOfMolecules();

    problem.l = N;

    problem.y = (double *) calloc(N , sizeof(double));
    problem.x = (struct svm_node **) calloc(N, sizeof(struct svm_node*));

    for(int i= 0; i< N; ++i)
        problem.x[i] = (struct svm_node *) calloc(N +2, sizeof(struct svm_node));

    int i_r = -1;

    for(int i =0; i< N; ++i)
    {
        i_r++;

        problem.y[i_r] = 1.0;

        problem.x[i_r][i_r+1].value = 1.0;
        problem.x[i_r][i_r+1].index = i_r+1;

        problem.x[i_r][0].value = i_r+1;
        problem.x[i_r][0].index = 0;

        problem.x[i_r][N+1].value = 0.0;
        problem.x[i_r][N+1].index = -1;

        int j_r = i_r;
        for(int j=i+1; j< N; ++j)
        {
            OBMol *mol1 = train->getMolecule(i);
            OBMol *mol2 = train->getMolecule(j);

            double K = m_cmfa->calculate(mol1, mol2);

            j_r++;
            problem.x[i_r][j_r+1].value = problem.x[j_r][i_r+1].value = K;

            problem.x[i_r][j_r+1].index = j_r+1;
            problem.x[j_r][i_r+1].index = i_r+1;

        }
    }

    model = svm_train(&problem, &param);

    //predict +1

    struct svm_node * testing = (struct svm_node *) calloc(N+2 , sizeof(struct svm_node));

    for(int i =0; i< test->getNumberOfMolecules(); ++i)
    {

        testing[0].value = 1;
        testing[0].index = 0;

        testing[N+1].value = 0.0;
        testing[N+1].index = -1;

        OBMol *mol1 = test->getMolecule(i);

        int i_p = 0;
        for(int j=0; j< N; ++j)
        {

            OBMol *mol2 = train->getMolecule(j);
            double K = m_cmfa->calculate(mol1, mol2);

            testing[i_p+1].value = K;
            testing[i_p+1].index = i_p+1;

            i_p++;

        }

        struct result * res = create_result();
        res->y_pred[0] = svm_predict(model, testing);
        res->y_real[0] = 0.0;

        results.push_back(res);

    }

    free(testing);

    for(int i =0; i< rn; ++i)
    {
        free(problem.x[i]);
        problem.x[i] = NULL;
    }

    free(problem.x);
    problem.x = NULL;

    free(problem.y);
    problem.y = NULL;

    svm_free_and_destroy_model(&model);
    svm_destroy_param(&param);

}

bool Svr::predict(double *x, double *y)
{

}
