
/*
 * mcmf.h
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

#ifndef __CMFA_H
#define __CMFA_H

#include <openbabel/mol.h>
#include <boost/utility.hpp>

#include <vector>
#include "kernel.h"

/*!
  Continuous Molecular Field Analysis.
*/

class CMFA :  boost::noncopyable
{
public:
    ~CMFA();
    double calculate(OBMol *mol1, OBMol *mol2, Mode mode = Training);

    void addKernel(CKernel *);
    void delKernel(CKernel *);

    virtual void save(FILE *fp);
    virtual void load(FILE * fp);

    void clear();
    void clearCache();

    void setParameters(const double * h);
    int count() const;

    void printSelKernels();
    const std::string & getKernelName(int i) const;

    std::vector< CKernel * > m_kernels;

private:

    const double * m_h;

    std::map < OBMol *, std::map<OBMol *, double > > gramm;
    //Kernel norms
    std::map < OBMol *, double> norms;

    double calculate_n(OBMol *mol1, OBMol *mol2, Mode mode = Training);
};

#endif

