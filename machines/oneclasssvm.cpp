
/*
 * oneclasssvm.cpp
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

#include "oneclasssvm.h"
#include <iostream>

int OneClassSVM::res_comp(struct res_auc  p1, struct res_auc  p2)
{
    return p1.fpr < p2.fpr;
}

OneClassSVM::OneClassSVM() : Machine("1-SVM")
{
    std::cout << "One Class SVM is loaded." << std::endl;

    param.svm_type = ONE_CLASS;
    param.kernel_type = PRECOMPUTED;

    //we dont use these parameter. But the library's inteface expects it.

    param.degree = 3;
    param.gamma = 0;
    param.coef0 = 0;
    param.cache_size = 100;
    param.C = 1;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;

    dimensionality = 1;

    m_NumParameters = 0;

    Parameters = NULL;
    lp = NULL;
    mp = NULL;

}

void OneClassSVM::setCMFA(CMFA *cmfa)
{
    m_cmfa = cmfa;
    m_cmfa->printSelKernels();

    if(lp != NULL)
    {
        free(lp);
        free(mp);
        free(Parameters);
    }

    m_NumParameters = cmfa->count() *2;

    Parameters = (double * ) calloc(m_NumParameters , sizeof(double));
    lp = (double * ) calloc(m_NumParameters , sizeof(double));
    mp = (double * ) calloc(m_NumParameters , sizeof(double));

    lp[0] = 0.001;
    for(int i =1; i< m_NumParameters; ++i)
        lp[i] = (i%2) ? 0.0001 : 0.001;

    mp[0] = 0.800;
    for(int i =1; i< m_NumParameters; ++i)
        mp[i] = (i%2) ? 0.1 : 1.0 ;

    Parameters[0] = 0.04;
    for(int i =1; i< m_NumParameters; ++i)
        Parameters[i] = (i%2) ? 0.01 : 0.01 ;

}

OneClassSVM::~OneClassSVM()
{
    if(lp != NULL)
    {
        free(lp);
        free(mp);
        free(Parameters);
    }
}

bool OneClassSVM::setData(SEAL *train_mols, SEAL *test_mols)
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

bool OneClassSVM::build(const double * params, const std::vector<int> &flags, const std::vector<int> &mask)
{

    param.nu = params[0];
    m_cmfa->setParameters(params + 1);

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

bool OneClassSVM::load(std::string & filename)
{

}

bool OneClassSVM::save(std::string &filename)
{

}

double OneClassSVM::statistic()
{
    /*
     One-class classification differs from other CV-methods.
     We calculate here predictions for test decoys structrues;
     */
    predict_decoys();

    /*Calculate AUC */

    double max_r = -DBL_MAX;
    double min_r = +DBL_MAX;

    for(int i =0; i < results.size(); ++i)
    {
        double temp = results[i]->y_pred[0];
        if( temp > max_r )
            max_r = temp;
        if(temp < min_r)
            min_r = temp;
    }

    printf("max_r : %g\nmin_r : %g\n", max_r, min_r);

    double step = (max_r - min_r ) / 10000.0;
    std::vector < struct res_auc > auc;

    for(double h = min_r - 10.0 * step; h< max_r + 10.0 * step ; h+= step)
    {

        double tn = 0.0;
        double tp = 0.0;
        double fp = 0.0;
        double fn = 0.0;

        for(int i=0; i< results.size(); ++i)
        {
            if(results[i]->y_real[0] == 1)
            {
                //real positive
                if(results[i]->y_pred[0] >= h) tp++;
                else fn++;
            }
            else
            {
                //real negative
                if(results[i]->y_pred[0] < h) tn++;
                else fp++;
            }
        }

        struct res_auc temp;
        temp.fpr = fp / (tn + fp);
        temp.tpr = tp / (tp + fn);
        temp.threshold = h;

        bool add = true;
        for(int i =0; i<auc.size(); ++i)
            if( fabs(auc[i].tpr - temp.tpr) < 1.0e-4 &&
                    fabs(auc[i].fpr - temp.fpr) < 1.0e-4)
            {
                add = false;
                break;
            }
        if(add)
            auc.push_back(temp);

    }


    std::sort(auc.begin(), auc.end(), res_comp);

    //for(int i =0; i< auc.size(); ++i)
    //	printf("%g %g\n", auc[i].fpr, auc[i].tpr);

    double AUC=0.0;

    double x1, x2, y1, y2;

    x1 = auc[0].fpr;
    y1 = auc[0].tpr;

    for(int i =1; i< auc.size(); ++i)
    {
        x2 = auc[i].fpr;
        y2 = auc[i].tpr;

        AUC+=(x2-x1) * (y1 + y2);
        x1 = x2;
        y1 = y2;
    }

    AUC = AUC * 0.5;
    printf("AUC: %g\n", AUC);

    if(mode)
    {
        for(int i =0; i< auc.size(); i++)
            fprintf(fres, "%g %g\n", auc[i].fpr, auc[i].tpr );

        double ot = 0.0;
        double m_ot1, m1;
        double m_ot2, m2;

        m1 = DBL_MAX;
        m2 = -DBL_MAX;

        m_ot1 = m_ot2 = 0;

        for(int i =0; i<auc.size(); ++i)
        {

            double x = auc[i].fpr;
            double y = auc[i].tpr;
            ot = auc[i].threshold;

            double r = y - (1.0 - x);
            if(r >= 0)
            {
                if( r <= m1)
                {
                    m1 = r;
                    m_ot1 = ot;
                }
            }
            else
            {
                if(r > m2)
                {
                    m2 = r;
                    m_ot2 = ot;
                }
            }

        }
        ot = (m_ot1 + m_ot2) / 2.0;

        double tn = 0.0;
        double tp = 0.0;
        double fp = 0.0;
        double fn = 0.0;

        for(int i=0; i< results.size(); ++i)
        {
            if(results[i]->y_real[0] == 1)
            {
                //real positive
                if(results[i]->y_pred[0] >= ot) tp++;
                else fn++;
            }
            else
            {
                //real negative
                if(results[i]->y_pred[0] < ot) tn++;
                else fp++;
            }
        }

        double fpr = fp / (tn + fp);
        double tpr = tp / (tp + fn);

        fprintf(fres, "\n%.2f %.4f ", AUC, ot);
        fprintf(fres, " %.0f %.0f %.0f %.0f %.2f %.2f ", tn, tp, fn, fp, 100 * tn/(tn+fp),100 * tp/(tp+fn));

        for(int i =0; i< m_NumParameters; ++i)
            fprintf(fres, " %g ", Parameters[i]);

        }

    return AUC;
}

void OneClassSVM::predict_decoys()
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

bool OneClassSVM::predict(double *x, double *y)
{

}
