
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
    m_norm = true;
}

CMFA::~CMFA()
{

}

void CMFA::setNormalise(bool norm)
{
    m_norm = norm;

}

double CMFA::calculate(OBMol *mol1, bool regime, OBMol * mol2)
{
    struct molkey mol;
    if(regime)
    {

        if( reinterpret_cast<int *> (mol1) < reinterpret_cast<int *> (mol2) )
        {
            mol.first = mol1;
            mol.second = mol2;
        }
        else
        {
            mol.first = mol2;
            mol.second = mol1;
        }

        if(gramm.find(mol) != gramm.end())
            return gramm[mol];
    }

    double s = 0.0;
    double hh = 0.0;

    int i=0;

    for( ; i< m_kernels.size() -1; ++i)
    {
        CKernel * kernel = m_kernels[i];
        double p = kernel->calculate(mol1, regime, mol2, m_h[i*2], m_norm);

        s+= m_h[i*2 +1] * p;
        hh+=m_h[i*2+1];
    }

    CKernel * kernel = m_kernels[m_kernels.size() -1];
    double p = kernel->calculate(mol1, regime, mol2, m_h[i*2], m_norm);
    //printf("-Kernel: %g %d %g\n", s, gramm.size(), p);

    s+= (1.0- hh) * p;

    if(regime)
        gramm[mol] = s;

    //printf("Kernel: %g %d %g\n", s, gramm.size(), hh);
    return s;

}

void CMFA::clearCache()
{
    //printf("Gramm size: %d\n", gramm.size());
    gramm.clear();
    clearNorms();
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

void CMFA::clear()
{
    for(int i = m_kernels.size() -1; i>=0; --i)
        m_kernels.erase( m_kernels.begin() + i);
}

const char * CMFA::getKernelName(int i)
{
    if(i >= m_kernels.size())
        return NULL;
    return m_kernels[i]->getName();
}

void CMFA::clearNorms()
{
    for(int i =0; i< m_kernels.size(); ++i)
        m_kernels[i]->clearNorms();
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

