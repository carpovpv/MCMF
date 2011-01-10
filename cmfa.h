
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

#include <vector>
#include "kernel.h"

/*!
  Continuous Molecular Field Analysis.
*/

class CMFA
{
public:
    CMFA();
    virtual ~CMFA();

    virtual double calculate(OBMol *, OBMol *);

    void addKernel(CKernel *);
    void delKernel(CKernel *);
    void clear();

    void setParameters(const double * h);
    int count() {
        return m_kernels.size();
    }

    void printSelKernels();
    const char * getKernelName(int i);

private:

    std::vector< CKernel * > m_kernels;

    const double * m_h;

    CMFA (const CMFA &);
    CMFA & operator=(const CMFA &);

};

#endif

