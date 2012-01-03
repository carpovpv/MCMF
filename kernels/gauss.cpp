
/*
 * gauss.cpp
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

#include "gauss.h"

GaussKernel::GaussKernel(DescriptorFactory * descr) : CKernel("Gaussian", descr)
{

}

double GaussKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode regime)
{
    const std::vector< double>  m1 = m_descrfactory->getDescriptors(mol1, regime);
    const std::vector< double>  m2 = m_descrfactory->getDescriptors(mol2, Training);

    double s = 0.0;
    const unsigned n = m1.size();

    for(int i=0; i< n; ++i)
        s+= pow((m1[i] - m2[i]), 2);

    return exp(-s * gamma);
}

