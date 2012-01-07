
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

double CMFA::calculate_n(OBMol *mol1, OBMol *mol2, Mode mode)
{
    double s = 0.0;
    double hh = 0.0;

    int i;

    for(i=0; i< m_kernels.size() -1; ++i)
    {
        CKernel * kernel = m_kernels[i];
        double p = kernel->calculate(mol1, mol2, m_h[i*2], mode);

        s  += m_h[i*2 +1] * p;
        hh += m_h[i*2+1];
    }

    CKernel * kernel = m_kernels[m_kernels.size() -1];
    double p = kernel->calculate(mol1, mol2, m_h[i*2], mode);

    s+= (1.0 - hh) * p;

    return s;
}

double CMFA::calculate(OBMol *m1, OBMol * m2, Mode regime)
{
    OBMol * mol1 = m1;
    OBMol * mol2 = m2;

    if(regime == Training) //обучение - берем из кэша.
    {
        if( reinterpret_cast<int *> (m1) > reinterpret_cast<int *> (m2) )
        {
            mol1 = m2;
            mol2 = m1;
        }

        if(gramm.find(mol1) != gramm.end())
            if(gramm[mol1].find(mol2) != gramm[mol1].end())
                return gramm[mol1][mol2];
    }

    double s = calculate_n(mol1, mol2, regime);

    //need normalisation

    if(regime == Training)
    {
        std::map < OBMol *, double>::iterator it1 = norms.find(mol1);
        std::map < OBMol *, double>::iterator it2 = norms.find(mol2);

        double s1, s2;
        if(it1 == norms.end())
        {
            s1 = calculate_n(mol1, mol1);
            norms[mol1] = s1;
        }
        else
            s1 = norms[mol1];

        if(it2 == norms.end())
        {
            s2 = calculate_n(mol2, mol2);
            norms[mol2] = s2;
        }
        else
            s2 = norms[mol2];
        s /= sqrt(s1 * s2);
    }
    else
    {
         double s1 = calculate_n(mol1, mol1);
         double s2 = calculate_n(mol2, mol2);

         s /= sqrt(s1 * s2);
    }

    if(regime == Training)
    {
        gramm[mol1][mol2] = s;
    }
    return s;

}

void CMFA::clearCache()
{
    /*printf("Gramm matrix:\n");
    std::map < OBMol *, std::map<OBMol *, double > >::iterator it;
    for(it = gramm.begin(); it!= gramm.end(); ++it)
        printf("%d\n", it->second.size());

    printf("=============\n");*/

    gramm.clear();
    norms.clear();
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

void CMFA::save(FILE *fp)
{
    for(int i =0; i< m_kernels.size(); ++i)
    {
        CKernel * kernel = m_kernels[i];
        kernel->save(fp);
    }
}

void CMFA::load(FILE *fp)
{

}

int CMFA::count() const
{
    return m_kernels.size();
}

const std::string & CMFA::getKernelName(int i) const
{
    return m_kernels[i]->getName();
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
        printf("%s is loaded.\n", m_kernels[i]->getName().c_str());
}

CMFA::~CMFA()
{
    for(int i=0; i< m_kernels.size(); ++i)
        delete m_kernels[i];
}
