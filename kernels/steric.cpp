
/*
 * steric.cpp
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

#include "steric.h"

StericKernel::StericKernel() : CKernel()
{
    name = "Steric";

    /*
    @ARTICLE{VJGDASMR,
      author = {Vinter, J.G. and Davis, A. and Saunders, M.R.},
      title = {Strategic approaches to drug design, I. An integrated software framework
            for molecular modelling},
      journal = {J. Comp.-Aided Molecular Design},
      year = {1987},
      volume = {1},
      pages = {31-51},
      owner = {carpovpv},
      timestamp = {2010.11.02}
    }
    */

    std::map < int, double > temp;

    temp[Atom_H] = 0.042;
    temp[Atom_C] = 0.067;
    temp[Atom_N] = 0.063;
    temp[Atom_O] = 0.069;
    temp[Atom_F] = 0.068;
    temp[Atom_S] = 0.115;
    temp[Atom_CL]= 0.115;
    temp[Atom_BR]= 0.136;
    temp[Atom_I] = 0.162;

    table[Atom_H] = temp;

    temp.clear();
    temp[Atom_C] = 0.107;
    temp[Atom_N] = 0.100;
    temp[Atom_O] = 0.111;
    temp[Atom_F] = 0.108;
    temp[Atom_S] = 0.183;
    temp[Atom_CL]= 0.183;
    temp[Atom_BR]= 0.215;
    temp[Atom_I] = 0.258;

    table[Atom_C] = temp;

    temp.clear();
    temp[Atom_N] = 0.095;
    temp[Atom_O] = 0.105;
    temp[Atom_F] = 0.102;
    temp[Atom_S] = 0.172;
    temp[Atom_CL]= 0.172;
    temp[Atom_BR]= 0.203;
    temp[Atom_I] = 0.243;

    table[Atom_N] = temp;

    temp.clear();
    temp[Atom_O] = 0.116;
    temp[Atom_F] = 0.112;
    temp[Atom_S] = 0.190;
    temp[Atom_CL]= 0.190;
    temp[Atom_BR]= 0.224;
    temp[Atom_I] = 0.268;

    table[Atom_O] = temp;

    temp.clear();
    temp[Atom_F] = 0.109;
    temp[Atom_S] = 0.185;
    temp[Atom_CL]= 0.185;
    temp[Atom_BR]= 0.217;
    temp[Atom_I] = 0.259;

    table[Atom_F] = temp;

    temp.clear();
    temp[Atom_S] = 0.314;
    temp[Atom_CL]= 0.314;
    temp[Atom_BR]= 0.369;
    temp[Atom_I] = 0.442;

    table[Atom_S] = temp;

    temp.clear();
    temp[Atom_CL]= 0.314;
    temp[Atom_BR]= 0.369;
    temp[Atom_I] = 0.442;

    table[Atom_CL] = temp;

    temp.clear();
    temp[Atom_BR]= 0.434;
    temp[Atom_I] = 0.522;

    table[Atom_BR] = temp;

    temp.clear();
    temp[Atom_I] = 0.623;

    table[Atom_I] = temp;

}

double StericKernel::steric(int atom1, int atom2)
{
    int i1, i2;

    if( atom2 > atom1 )
        i1 = atom1, i2 = atom2;
    else
        i1 = atom2, i2 = atom1;

    if(table.find(i1) != table.end())
    {
        if(table[i1].find(i2) != table[i1].end())
            return table[i1][i2];
        else
            return 0.0;

    }
    return 0.0;
}

StericKernel::~StericKernel()
{
}

double StericKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma)
{
    double s = 0.0;
    double w2 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {
        FOR_ATOMS_OF_MOL(b, mol2)
        {
            int i1 = a->GetAtomicNum();
            int i2 = b->GetAtomicNum();
            double w = steric(i1, i2);
            w2 += w*w;

            double x = a->x() - b->x();
            double y = a->y() - b->y();
            double z = a->z() - b->z();
            s += w * exp ( -gamma/2 * ( x*x + y*y +z*z  ) );
        }
    }
    return s / w2;
}

