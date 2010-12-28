
/*
 * cmfa.cpp
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

#include "cmfa.h"

CMFA::CMFA()
{
}

CMFA::~CMFA()
{

}

double CMFA::calculate(OBMol *mol1, OBMol * mol2)
{

    double s = 0.0;
    double hh = 0.0;
    int i=0;
    for( ; i< m_kernels.size() -1; ++i)
    {
        s+= m_h[i*2 +1] * m_kernels[i]->calculate(mol1, mol2, m_h[i*2]);
        hh+=m_h[i*2+1];
    }

    s+= (1.0- hh) * m_kernels[m_kernels.size() -1]->calculate(mol1, mol2, m_h[i*2]);
    return s;

}

void CMFA::setParameters(const double *h)
{
    m_h = h;
}

void CMFA::addKernel(CKernel * kernel)
{
    bool added = true;
    for(int i=0; i< m_kernels.size(); ++i)
        if(m_kernels[i] == kernel)
        {
            added = false;
            break;
        }

    if(added)
        m_kernels.push_back(kernel);
}

void CMFA::delKernel(CKernel * kernel)
{
    for(int i = m_kernels.size() -1; i>=0; --i)
    {
        if(m_kernels[i] == kernel)
            m_kernels.erase( m_kernels.begin() + i);
    }
}

void CMFA::printSelKernels()
{
    for(int i =0; i<m_kernels.size(); ++i)
        printf("%s is loaded.\n", m_kernels[i]->getName());
}

