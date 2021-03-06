
/*
 * abraham.cpp
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

#include "abraham.h"
#include "../fields.h"

AbrahamKernelA::AbrahamKernelA() : CKernel("AbrahamA")
{
    kerncode = D_ABRAHAMA;
}

double AbrahamKernelA::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode)
{
    double  s = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {

        Fields * f = dynamic_cast<Fields *>( a->GetData(OBGenericDataType::CustomData0));
        w1 = f->getValue(Fields::AbrahamA);

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            f = dynamic_cast<Fields *>( b->GetData(OBGenericDataType::CustomData0));
            w2 = f->getValue(Fields::AbrahamA);

            double x = a->x() - b->x();
            double y = a->y() - b->y();
            double z = a->z() - b->z();

            s += w1 * w2 * exp ( -gamma/4.0 * ( x*x + y*y +z*z  ) );
        }
    }

    return COEFF(gamma) * s;
}

AbrahamKernelB::AbrahamKernelB() : CKernel("AbrahamB")
{
    kerncode = D_ABRAHAMB;
}

double AbrahamKernelB::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode)
{
    double  s = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {

        Fields * f = dynamic_cast<Fields *>( a->GetData(OBGenericDataType::CustomData0));
        w1 = f->getValue(Fields::AbrahamB);

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            f = dynamic_cast<Fields *>( b->GetData(OBGenericDataType::CustomData0));
            w2 = f->getValue(Fields::AbrahamB);

            double x = a->x() - b->x();
            double y = a->y() - b->y();
            double z = a->z() - b->z();

            s += w1 * w2 * exp ( -gamma/4.0 * ( x*x + y*y +z*z  ) );
        }
    }

    return COEFF(gamma) * s;
}

AbrahamKernelE::AbrahamKernelE() : CKernel("AbrahamE")
{
    kerncode = D_ABRAHAME;
}

double AbrahamKernelE::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode)
{
    double  s = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {

        Fields * f = dynamic_cast<Fields *>( a->GetData(OBGenericDataType::CustomData0));
        w1 = f->getValue(Fields::AbrahamE);

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            f = dynamic_cast<Fields *>( b->GetData(OBGenericDataType::CustomData0));
            w2 = f->getValue(Fields::AbrahamE);

            double x = a->x() - b->x();
            double y = a->y() - b->y();
            double z = a->z() - b->z();

            s += w1 * w2 * exp ( -gamma/4.0 * ( x*x + y*y +z*z  ) );
        }
    }

    return COEFF(gamma) * s;
}

AbrahamKernelS::AbrahamKernelS() : CKernel("AbrahamS")
{
    kerncode = D_ABRAHAMS;
}

double AbrahamKernelS::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode)
{
    double  s = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {

        Fields * f = dynamic_cast<Fields *>( a->GetData(OBGenericDataType::CustomData0));
        w1 = f->getValue(Fields::AbrahamS);

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            f = dynamic_cast<Fields *>( b->GetData(OBGenericDataType::CustomData0));
            w2 = f->getValue(Fields::AbrahamS);

            double x = a->x() - b->x();
            double y = a->y() - b->y();
            double z = a->z() - b->z();

            s += w1 * w2 * exp ( -gamma/4.0 * ( x*x + y*y +z*z  ) );
        }
    }
    return COEFF(gamma) * s;
}


