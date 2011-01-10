
/*
 * fp2s.h
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

#ifndef __FP2s_H
#define __FP2s_H

#include "../descfact.h"
#include <openbabel/fingerprint.h>


using namespace OpenBabel;

class FingerPrints2s : public DescriptorFactory
{
public:

    FingerPrints2s();
    ~FingerPrints2s();
    bool needMapping() const;
    const std::vector < struct Descriptor > & getDescriptors(OBMol *);

private:

    std::map < OBMol * , std::vector< struct Descriptor > > m_descrs;
    const unsigned int m_ndescr;
    OBFingerprint * m_ob;
};

#endif
