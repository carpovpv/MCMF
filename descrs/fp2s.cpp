
/*
 * fp2s.cpp
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

#include "fp2s.h"
#include <openbabel/plugin.h>

FingerPrints2s::FingerPrints2s() :
    DescriptorFactory("FP2"),
    m_ndescr(1024)
{
    m_ob  = OBFingerprint::FindFingerprint("FP2");
    descrs.resize(m_ndescr);
}

const std::vector < double > & FingerPrints2s::getDescriptors(OBMol * mol, Mode regime)
{
    if(regime == Training)
    {
        std::map < OBMol *, std::vector < double > >::iterator it;
        it = m_descrs.find(mol);

        if(it != m_descrs.end())
        {
            return m_descrs[mol];
        }
    }
    //calc descrs...

    std::vector< unsigned int> fp;
    m_ob->GetFingerprint(mol,fp, m_ndescr);

    unsigned int one =  1;
    one = one  << (8 * sizeof(unsigned int) -1) ;

    int N = 0;
    for (int i=0; i< fp.size(); ++i)
    {
        unsigned int temp = fp[i];
        for (int j=0; j< (8 * sizeof(unsigned int )); ++j)
        {
            descrs[N++] = (temp & one) ? 1.0 : 0.0;
            temp = temp << 1  ;
        }
    }
    if(regime == Training)
        m_descrs[mol] = descrs;

    return descrs;
}
