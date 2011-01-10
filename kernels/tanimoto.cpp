
/*
 * tanimoto.cpp
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


#include "tanimoto.h"

TanimotoKernel::TanimotoKernel(DescriptorFactory * descr) : CKernel( descr)
{
    name =  "Tanimoto";
    m_remap = descr != NULL ? descr->needMapping() : false;
}

TanimotoKernel::~TanimotoKernel()
{
}

double TanimotoKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma)
{
    if(m_remap)
    {
    }
    else
    {

        const std::vector< struct Descriptor>  &m1 = m_descrfactory->getDescriptors(mol1);
        const std::vector< struct Descriptor>  &m2 = m_descrfactory->getDescriptors(mol2);

        const unsigned n = m1.size();

        double p = 0.0, q=0, r =0;
        unsigned int i;
        for (i=0; i<n; i++)
        {
            if (m1[i].value >=0.5 && m2[i].value >= 0.5) p++;
            if (m1[i].value <0.5 && m2[i].value  >= 0.5) q++;
            if (m1[i].value >= 0.5 && m2[i].value < 0.5) r++;
        }
        return p/ (p+q+r);


    }

    return 0;

}

