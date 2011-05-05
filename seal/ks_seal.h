
/*
 * ks_seal.h
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

#ifndef _KS_SEAL
#define _KS_SEAL

#include <openbabel/math/vector3.h>
#include "seal.h"
#include <nlopt.h>
#include <vector>

class KS_Seal : public SEAL
{
public:
    KS_Seal(const char * fp,
            double alpha = 0.5,
            double wE = 1.0,
            double wS = 1.0,
            int probes =1
           );
    virtual ~KS_Seal();
    void set2Mol(OBMol & mol);

private:
    int m_cycles;
    nlopt_opt opt;

    double alpha, wE, wS;
    int probes;

    double lq1, lq2, lq3, lq4;
    vector3 lt;

    double score(unsigned int st1, unsigned int st2);

    double norm( const double *);
    void createQ(const double *, double *);

    void align();
    double optim(double * , double *);

    OBMol mainMol;
    OBMol supMol;
    static double score_function(unsigned , const double *p, double *, void *);

};

#endif

