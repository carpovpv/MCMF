
/*
 * machine.cpp
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

#include "machine.h"

#include <vector>
#include <math.h>
#include <algorithm>
#include <nlopt.h>
#include <time.h>
#include "cmfa.h"

Machine::Machine(const std::string & name) : m_name(name)
{
    test = NULL;
    train = NULL;
    fres = NULL;
    m_CV = 10;
    reallast = false;

}

void Machine::setOutput(FILE * f)
{
    fres = f;
}

std::string & Machine::getName()
{
    return m_name;
}

void Machine::setParameters(double * h)
{
    for(int i=0; i< m_NumParameters; ++i)
        Parameters[i] = lp[i] + (mp[i] - lp[i] - 0.1) / (1.00 * RAND_MAX) * rand();

    for(int i =0; i< m_NumParameters; ++i)
        if(h[i] != UNKNOWN_VALUE)
            Parameters[i] = h[i];
}

bool Machine::setData(SEAL *train_mols, SEAL *test_mols)
{
    train = train_mols;
    test = test_mols;
    return true;
}

void Machine::set_CV(int CV)
{
    m_CV = CV;
}

int Machine::get_CV()
{
    return m_CV;
}

void Machine::clearCache()
{
    m_cmfa->clearCache();
}

double Machine::optim(unsigned, const double *m_params, double *, void *ptr)
{

    //printf("Start optim\n");

    Machine * machine = static_cast<Machine *> (ptr);

    if(machine->train == NULL)
    {

        return 0;
    }
    const int N = machine->train->getNumberOfMolecules();

    double CV = machine->get_CV();

    machine->clearCache();

    printf("[");
    fprintf(machine->fres, "[");
    for(int i =0; i< machine->m_NumParameters; ++i)
    {
        printf(" %.6f ", m_params[i]);
        fprintf(machine->fres," %.6f ", m_params[i]);
    }
    printf("] => ");
    fprintf(machine->fres, "] => ");

    for(int i =0; i < machine->results.size(); ++i)
        machine->drop_result(machine->results[i]);

    machine->results.clear();

    std::vector < int > mask(N);
    std::vector < int > flags(N);

    for(int i =0; i< N; ++i)
        mask[i] = i;


    //random_shuffle(mask.begin(), mask.end());

    const int N_CV = ceil(1.0 * N / CV);

    for(int i_cv = 0; i_cv < N; i_cv+= N_CV)
    {
        for(int i=0; i< flags.size(); ++i)
            flags[i] = 1;

        const int l_cv = i_cv + N_CV;
        for(int i= i_cv; i< (l_cv < N ? l_cv : N); ++i)
            flags[i] = 0;

        int rn = 0;
        for(int i = 0; i< N; ++i)
        {
            if(flags[i] == 1) rn++;
            //printf("%d ", flags[i]);
        }

        //printf("Build %d %d\n", i_cv, N);
        machine->build(m_params, flags, mask);
    }

    //printf("Numbers: %d\n", machine->results.size());
    double s = machine->statistic();

    for(int i =0; i < machine->results.size(); ++i)
        machine->drop_result(machine->results[i]);

    machine->results.clear();

    //printf("End optim\n");
    return s;

}

void Machine::setCMFA(CMFA *cmfa)
{

}

struct result * Machine::create_result()
{
    struct result * res = new struct result;
    res->y_real = new double [dimensionality];
    res->y_pred = new double [dimensionality];
    return res;
}

void Machine::drop_result(struct result * res)
{
    delete [] res->y_pred;
    delete [] res->y_real;
    delete res;
}

Machine::~Machine()
{

}

double Machine::create(nlopt_algorithm algo)
{
    //grid search
    /*
        FILE * fp = fopen("grid.search", "w");

        for(double nu = 0.001; nu < 0.8; nu += 0.05)
          for(double c = 1e-5; c< 10; c*=10)
        {
            for(double g = 0.0001; g < 10; g+=10)
            {

                Parameters[0] = nu;
                Parameters[1] = c;
                Parameters[2] = g;

                double q = optim(0, Parameters, NULL, this);

                fprintf(fp, "%g %g %g %g\n", nu,c , g, q);
                fflush(fp);

            }


        }
        fclose(fp);

        return true;
    */

    nlopt_opt opt = nlopt_create(algo, m_NumParameters);
    nlopt_set_max_objective(opt, optim, this);

    nlopt_set_xtol_rel(opt, 1e-3);
    nlopt_set_stopval(opt, 0.96);
    nlopt_set_ftol_rel(opt, 1e-3);

    nlopt_set_lower_bounds( opt, lp);
    nlopt_set_upper_bounds( opt, mp);

    for(int i = 0; i< m_NumParameters; ++i)
        printf("%g <=  %g  <= %g\n", lp[i], Parameters[i], mp[i]);


    double minf = 0.0;
    int err;

    mode = 0;

    if ((err = nlopt_optimize(opt, Parameters, &minf)) < 0)
    {
        fprintf(stderr,"nlopt failed %d!\n", err);
        exit(0);

    }
    else
    {
        printf("found maximum at f(" );
        for(int i =0; i< m_NumParameters; ++i)
            printf(" %g ", Parameters[i]);
        printf (") = %0.10g\n", minf);

        //this cause output statistic to the result file.
        mode = 1;
        optim(0, Parameters, NULL, this);

    }

    nlopt_destroy(opt);

    return minf;

}

double Machine::create_random(int max_iter)
{

    double *best_params = (double *) calloc(m_NumParameters,sizeof(double));
    double best_rmse = -RAND_MAX;

    for(int i=0; i< m_NumParameters; ++i)
        best_params[i] = Parameters[i];

    int iter = 0;

    while(iter++<max_iter)
    {

        //double temp = create();
        double temp = optim(0, Parameters,NULL, this);

        if(temp > best_rmse)
        {
            for(int i=0; i< m_NumParameters; ++i)
                best_params[i] = Parameters[i];
            best_rmse = temp;
        }


        for(int i=0; i< m_NumParameters; ++i)
            Parameters[i] = lp[i] + (mp[i] - lp[i] - 0.1) / (1.00 * RAND_MAX) * rand();

    }

    for(int i =0; i< m_NumParameters; ++i)
        Parameters[i] = best_params[i];

    reallast = true;
    double temp = create();

    if(temp > best_rmse)
    {
        best_rmse = temp;
        for(int i=0; i< m_NumParameters; ++i)
            best_params[i] = Parameters[i];
    }

    printf("Final:\n");
    printf("found maximum at f(" );
    for(int i =0; i< m_NumParameters; ++i)
        printf(" %g ", best_params[i]);
    printf (") = %g\n", best_rmse);

    free(best_params);

}
