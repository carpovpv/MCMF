
/*
 * hydropho.cpp
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

#include "hydropho.h"
#include "../fields.h"

HydrophobicKernel::HydrophobicKernel() : CKernel("Hydrophobic")
{
}

double HydrophobicKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode)
{

    double  s = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {
        Fields * f = dynamic_cast<Fields *>( a->GetData(OBGenericDataType::CustomData0));
        w1 = f->getValue(Fields::Hydrophobic);

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            f = dynamic_cast<Fields *> (b->GetData(OBGenericDataType::CustomData0));
            w2 = f->getValue(Fields::Hydrophobic);

            s += w1 * w2 * exp ( -gamma / 4.0 * ( pow((a->x() - b->x()), 2) + pow((a->y() - b->y()), 2)  + pow((a->z() - b->z()), 2)  ));
        }
    }
    s = s * COEFF(gamma);
    return s;

}

