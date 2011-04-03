
/*
 * seal.cpp
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

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <fstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "seal.h"
#include "../fields.h"

using namespace OpenBabel;
using namespace boost;


SEAL::SEAL(const char * sdf,  const std::vector<std::string> * props)
{

    std::ifstream ifs(sdf);

    if (!ifs)
        throw std::runtime_error("File with structures not found!");

    readStructures(&ifs, props,-1);

    ifs.close();
}

SEAL::SEAL(std::ifstream *ifs, const std::vector<std::string> *props, int mn)
{
    readStructures(ifs, props, mn);
}

void SEAL::readStructures(std::ifstream *ifs, const std::vector<std::string> *props, int mn)
{
    conv = new OBConversion(ifs, &std::cout);
    conv->SetInAndOutFormats("SDF","SDF");

    if (debug)
        std::cerr <<  "Reading the database into memory..." << std::endl;

    OBMol mol;
    int i=0;
    while ( conv->Read(&mol) )
    {
        i++;
        if( mn > 0 && i > mn  )
            break;
        bool add = true;
        if(props)
        {
            for(int j=0; j< props->size(); ++j)
            {
                if(!mol.HasData(props->at(j)))
                        add = false;
                else
                {
                    if(mol.HasData(props->at(j)))
                    {
                        std::string data = mol.GetData(props->at(j))->GetValue();
                        trim(data);
                        if(data.empty())
                            add = false;
                    }
                }
            }
        }

        if(add)
        {
           mol.DeleteHydrogens();
           mols.push_back(mol);
        }

    }

    if (debug)
        std::cerr << "Done reading " << mols.size()  << "." << std::endl;
}

void SEAL::go()
{
    for(int i=0; i< mols.size(); ++i)
        mols[i].Center();

    align();

    for(int i =0; i< m_mols.size(); ++i)
    {
        OBAtom *atom;
        FOR_ATOMS_OF_MOL(atom, m_mols[i])
        {
            Fields * ff = new Fields();
            ff->calcValues(&*atom);
            atom->SetData(ff);
        }
    }
}

OBMol * SEAL::getMolecule(int p)
{
    if ( p < 0 || p >= m_mols.size())
        return NULL;
    return & m_mols[p];
}

void SEAL::align()
{
    for(int i =0; i< mols.size(); ++i)
        m_mols.push_back(mols[i]);
}

int SEAL::getNumberOfMolecules() const
{
    return m_mols.size();
}

SEAL::~SEAL()
{
    delete conv;
}








