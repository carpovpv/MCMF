
/*
 * seal.h
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

#ifndef _SEAL
#define _SEAL

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <vector>

using namespace OpenBabel;

class SEAL
{
public:
    SEAL(const char * sdf);
    virtual ~SEAL();
    void go();

    int getNumberOfMolecules() const;
    OBMol * getMolecule(int);

private:
    SEAL (const SEAL &);
    SEAL & operator=(const SEAL &);


    static const int  optim = 1000;

    virtual void align();

protected:
    OBConversion * conv;

    static const bool debug = false;

    std::vector <OBMol> mols;
    std::vector <OBMol> m_mols;

};


#endif

