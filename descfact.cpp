
/*
 * descfact.cpp
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

#include "descfact.h"

DescriptorFactory::DescriptorFactory(const std::string & name)
    :m_name(name)
{
    fprintf(stderr, "Loading %s.\n", name.c_str());
    descrcode = D_UNKNOWNDESCR;
    prev = -1;
    prognosis = false;
}

void DescriptorFactory::setPrognosis(bool ok)
{
    prognosis = ok;
}

void DescriptorFactory::load(FILE *fp)
{

}

void DescriptorFactory::save(FILE *fp) const
{
    fprintf(fp,"Descriptors: %d (%s)\n%d ", descrcode, m_name.c_str(), descrnames.size() );
    for(int i =0; i< descrnames.size(); ++i)
        fprintf(fp, "%s ", descrnames[i].c_str());
    fprintf(fp, "\n");
}

DescriptorFactory::~DescriptorFactory()
{

}

const std::string & DescriptorFactory::getName() const
{
    return m_name;
}

