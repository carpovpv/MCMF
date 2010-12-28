
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

#include "seal.h"

using namespace OpenBabel;

SEAL::SEAL(const char * sdf)
{

    std::ifstream ifs(sdf);

    if (!ifs)
    {
        throw std::runtime_error("File with structures not found!");
    }

    conv = new OBConversion(&ifs, &std::cout);
    conv->SetInAndOutFormats("SDF","SDF");

    if (debug)
        std::cerr <<  "Reding the database into memory..." << std::endl;

    OBMol mol;
    int i=0;
    while (conv->Read(&mol) )
    {

        mol.DeleteHydrogens();
        mols.push_back(mol);
    }

    if (debug)
        std::cerr << "Done reading." << std::endl;

    ifs.close();
}

void SEAL::go()
{
    for(int i=0; i< mols.size(); ++i)
        mols[i].Center();

    align();
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








