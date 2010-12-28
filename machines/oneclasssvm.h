
/*
 * oneclasssvm.h
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

#ifndef ONECLASSSVM_H
#define ONECLASSSVM_H

#include "../machine.h"
#include "../cmfa.h"
#include "../svm.h"



class OneClassSVM : public Machine
{
public:
    OneClassSVM();
    ~OneClassSVM();

    bool load(std::string & filename);
    bool save(std::string & filename);

    bool predict(double *x, double *y);
    bool build(const double * params, const std::vector<int> &flags, const std::vector<int> &mask);

    void setCMFA( CMFA * cmfa);
    bool setData(SEAL *train_mols, SEAL *test_mols);

    double statistic();

private:

    struct res_auc
    {
        double fpr;
        double tpr;
        double threshold;
    };
    CMFA * m_cmfa;

    struct svm_parameter param;
    struct svm_problem problem;
    struct svm_model * model;

    void predict_decoys();

    static int res_comp(struct res_auc a1, struct res_auc a2);

};

#endif // ONECLASSSVM_H
