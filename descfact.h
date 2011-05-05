
/*
 * descfact.h
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

#ifndef __DESCRFACTORY_H
#define __DESCRFACTORY_H

#include <vector>
#include <string>
#include <map>

#include <openbabel/mol.h>

using namespace OpenBabel;

struct Descriptor
{
    double value;
    std::string name;
};

class DescriptorFactory
{
public:
    DescriptorFactory();
    virtual ~DescriptorFactory() ;

    std::string getName() const ;

    virtual const std::vector < struct Descriptor > & getDescriptors(OBMol *, bool) = 0;
    virtual bool needMapping()  const ;

protected:
    std::string name;
    long prev;
    std::vector< struct Descriptor> descrs;

private:
    DescriptorFactory (const DescriptorFactory &);
    DescriptorFactory & operator=(const DescriptorFactory &);


};

#endif

