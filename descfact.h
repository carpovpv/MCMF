
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
#include <boost/utility.hpp>

using namespace OpenBabel;

//Если используется Training, все значения дескрипторов, ядер и норм запоминаются в локальных
//кэшах. Если Prediction, то они всегда вычисляются
//В функциях прогноза, первая молекула -- всегда тестируемая, а вторая принадлежит классу
//обучаемых структур.

enum Mode {Training, Prediction};

class DescriptorFactory : boost::noncopyable
{
public:
    DescriptorFactory(const std::string & name);
    virtual ~DescriptorFactory() ;

    const std::string & getName() const ;

    /// regime true -> training
    virtual const std::vector < double > & getDescriptors(OBMol *, Mode regime) = 0;

    virtual void save(FILE * fp) const;
    virtual void load(FILE * fp);

protected:

    //временный вектор для работы над текущей структурой
    std::vector< double> descrs;

    //кэш дескрипторов для обучения
    std::map < OBMol * , std::vector< double > > m_descrs;

    //названия дескрипторов. Во многих блоках не используется,
    //
    std::vector <std::string> descrnames;

private:

    std::string m_name;

};

#endif

