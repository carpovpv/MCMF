
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

FingerPrints2s::FingerPrints2s() : m_ndescr(1024)
{
    m_ob  = OBFingerprint::FindFingerprint("FP2");
    if(m_ob == NULL)
        printf("Problem in FP2\n");
    name = "FP2";
    descrs.resize(m_ndescr);
}

FingerPrints2s::~FingerPrints2s()
{
    m_descrs.clear();
}

bool FingerPrints2s::needMapping() const
{
    return false;
}

const std::vector < struct Descriptor > & FingerPrints2s::getDescriptors(OBMol * mol, bool regime)
{

    if(regime)
    {
        std::map < OBMol *, std::vector < struct Descriptor > >::iterator it;
        it = m_descrs.find(mol);

        if(it != m_descrs.end())
        {
            return m_descrs[mol];
        }
    }
    else
    {

        std::string c = mol->GetData("prognosis")->GetValue();
        long cur = atol(c.c_str());

        if( prev == cur)
            return descrs;
        prev = cur;
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
            descrs[N++].value = (temp & one) ? 1.0 : 0.0;
            temp = temp << 1  ;
        }
    }        
    if(regime)
        m_descrs[mol] = descrs;

    return descrs;
}
