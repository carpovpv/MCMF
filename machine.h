
/*
 * machine.h
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

#ifndef MACHINE_H
#define MACHINE_H

#include <string>
#include "seal/seal.h"
#include <nlopt.h>

#define UNKNOWN_VALUE 666
#define GNUPLOT "gnuplot -persist"

struct result
{
    double * y_real;
    double * y_pred;
};

class Machine
{
public:

    Machine(const std::string & name);
    virtual ~Machine();

    virtual bool setData(SEAL * train_mols, SEAL * test_mols);
    virtual double statistic() = 0;

    virtual bool load(const char * filename) {}
    virtual bool save(const char * filename) = 0;

    virtual bool build(const double *,
                       const std::vector< int > & flags,
                       const std::vector<int > & mask) = 0;

    virtual bool predict(OBMol *) = 0;
    virtual void clearCache() {}

    void setParameters(double *);
    void setOutput(FILE *);

    static double optim(unsigned, const double *m_params, double *, void * ptr);

    double create(nlopt_algorithm algo = NLOPT_LN_BOBYQA);
    double create_random(int maxiter = 3);

    virtual void init() {}

    std::string & getName();

protected:

    std::string m_name;
    int m_NumParameters;

    double *Parameters, *lp, *mp;
    int dimensionality;

    SEAL * train, * test;
    bool mode;

    FILE * fres;

    std::vector<struct result * > results;

    struct result * create_result();
    void drop_result(struct result *);

private:

    Machine (const Machine &);
    Machine & operator=(const Machine &);


};

#endif // MACHINE_H
