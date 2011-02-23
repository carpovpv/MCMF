
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
#include <map>

Svr::Svr() : Machine("SVR")
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

    dimensionality = 1; //always 1 !!!

    m_NumParameters = 0;

    Parameters = NULL;
    lp = NULL;
    mp = NULL;

    m_props = NULL;
    NP = 0;

    gp = popen(GNUPLOT,"w"); /* 'gp' is the pipe descriptor */
    if(gp == NULL)
        fprintf(stderr, "Error init gnuplot\n");

}

void Svr::setProps(std::vector<std::string> *props)
{
    m_props = props;
    NP = props->size();

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

    lp[0] = 0.001;
    lp[1] = 1e-3;
    for(int i =2; i< m_NumParameters; ++i)
        lp[i] = (i%2) ? 0 : 0.0001;

    mp[0] = 0.999;
    mp[1] = 10.000;
    for(int i =2; i< m_NumParameters; ++i)
        mp[i] = (i%2) ? 1.0 : 10.0 ;

    Parameters[0] = 0.4;
    Parameters[1] = 0.2;

    for(int i =2; i< m_NumParameters; ++i)
        Parameters[i] = (i%2) ? 0.0001 : 0.0001 ;

}

Svr::~Svr()
{
    if(lp != NULL)
    {
        free(lp);
        free(mp);
        free(Parameters);
    }
    if(gp !=  NULL)
        pclose(gp);
}

bool Svr::setData(SEAL *train_mols, SEAL *test_mols)
{

    if(NP != 1)
    {
        printf("You can model only one property.\n");
        return false;
    }

    train = train_mols;
    test = test_mols;

    if(train == NULL)
    {
        fprintf(stderr,"Train sets must be supplied.\n");
        return false;
    }

    s_data.clear();
    m_data.clear();

    for(int i =0; i< train->getNumberOfMolecules(); ++i)
    {
        std::vector < double > temp;
        for(int j=0; j< NP; ++j)
        {
            std::string data = train->getMolecule(i)->GetData(m_props->at(j))->GetValue();
            temp.push_back(atof(data.c_str()));
        }
        m_data.push_back(temp);
    }

    printf("Dimension: %d\nLoaded: %d\n", NP, m_data.size());
    if(m_data.size() < 10)
    {
        printf("Number of points for SVR must be greater than 10!\n");
        //return false;
    }

    //scale
    std::vector< double > min_p (NP, DBL_MAX);
    std::vector< double > max_p (NP, DBL_MIN);

    for(int i = 0; i < m_data.size(); ++i)
    {
        for(int j =0; j< NP; ++j)
        {
            if(min_p[j] > m_data[i][j])
                min_p[j] = m_data[i][j];
            if(max_p[j] < m_data[i][j])
                max_p[j] = m_data[i][j];
        }
    }

    printf("Min:");
    for(int i =0; i < NP; ++i)
        printf(" %g", min_p[i]);

    printf("\nMax:");
    for(int i =0; i < NP; ++i)
        printf(" %g", max_p[i]);

    printf("\n");

    //enlarge interval
    const double L = 10;

    /*

     ----|-----|-------|-----|--------------------> x
         l1    m1      m2    l2

     m1 and m2 min and max values for the property
     (m1 - l1) / (m2-m1) = (l2-m2)/(m2-m1) = 1.0 / L

     */

    l1.resize(NP);
    l2.resize(NP);

    for(int i =0; i< NP; ++i)
    {
        l1[i] = min_p[i] - (max_p[i] - min_p[i]) / L;
        l2[i] = max_p[i] + (max_p[i] - min_p[i]) / L;
    }

    //scale to [0.1; 0.9]

    for(int i = 0; i< m_data.size() ; ++i)
    {
        std::vector< double > temp(NP);
        for(int j=0; j< NP; ++j)
            temp[j] = 0.9 + 0.8 *(m_data[i][j] - l2[j]) / (l2[j]-l1[j]);

        s_data.push_back(temp);
    }

    //calc SD for q2 metric
    SD = 0.0;
    A = 0.0;

    for(int i =0; i< m_data.size(); ++i)
        A += m_data[i][0];

    A /= m_data.size();
    for(int i =0; i< m_data.size(); ++i)
    {
        double y = m_data[i][0] - A;
        SD += y*y;
    }

    printf("SD = %g\nA= %g\n", SD,A);

    return true;
}

void Svr::clearCache()
{
    m_cmfa->clearCache();
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
        problem.x[i] = (struct svm_node *) calloc(rn+2, sizeof(struct svm_node));

    int i_r = -1;

    for(int i =0; i< N; ++i)
    {
        if(flags[i] == 0) continue;

        i_r++;

        problem.y[i_r] = s_data[mask[i]][0];

        //first element
        problem.x[i_r][0].value = i_r+1.0;
        problem.x[i_r][0].index = 0;

        //last element
        problem.x[i_r][rn+1].value = 0.0;
        problem.x[i_r][rn+1].index = -1;

        int j_r = i_r;
        for(int j=i; j< N; ++j) //j=i+1;
        {

            if(flags[j] == 0) continue;

            OBMol *mol1 = train->getMolecule(mask[i]);
            OBMol *mol2 = train->getMolecule(mask[j]);

            double K = m_cmfa->calculate(mol1, mol2);

            problem.x[i_r][j_r+1].value = K;
            problem.x[j_r][i_r+1].value = K;

            problem.x[j_r][i_r+1].index = i_r+1;
            problem.x[i_r][j_r+1].index = j_r+1;

            j_r++;
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
        res->y_real[0] = m_data[mask[i]][0];
        res->y_pred[0] = l2[0] + (l2[0] - l1[0]) * (res->y_pred[0] - 0.9) / 0.8;

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

bool Svr::save(const char * filename)
{

}

double Svr::statistic()
{

    //calc q2 and RMSE
    double PRESS = 0.0;

    double r = 0.0;
    double xx = 0.0;
    double yy = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    double xy = 0.0;

    const int n = results.size();

    FILE * fp = fopen("temp","w");

    for(int i=0; i< n; ++i)
    {
       double y = (m_data[i][0] - results[i]->y_pred[0]);
       fprintf(fp,"%g %g\n", m_data[i][0], results[i]->y_pred[0] );
       PRESS += y*y;

       xy += m_data[i][0] *  results[i]->y_pred[0];
       xx += m_data[i][0];
       yy += results[i]->y_pred[0];
       x2 += m_data[i][0] * m_data[i][0];
       y2 += results[i]->y_pred[0] *results[i]->y_pred[0];

    }

    fclose(fp);

    fprintf(gp, "plot 'temp', x\n");

    fflush(gp);

    double q2 = 1.0 - PRESS / SD;
    double RMSE = sqrt(PRESS/ results.size());

    r = (n * xy - xx * yy) / ( sqrt(n*x2 - xx*xx) * sqrt(n*y2 - yy*yy));


    printf("Q^2 = %g RMSE = %g R = %g R^2 = %g\n", q2, RMSE, r, r*r);

    //m_cmfa->clearNorms();

    return -RMSE;
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
        for(int j=i; j< N; ++j)
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

bool Svr::predict(OBMol * mol)
{

}
