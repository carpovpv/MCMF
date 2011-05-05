
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

AbrahamKernelA::AbrahamKernelA() : CKernel()
{
    name = "AbrahamA";

}

double AbrahamKernelA::calculate(OBMol * mol1, bool regime, OBMol * mol2, double gamma, bool norm)
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
    if(norm)
    {

       if(norms.find(mol1) == norms.end())
       {
           norms[mol1] = calculate(mol1,regime, mol1, gamma, false);
       }

       if(norms.find(mol2) == norms.end())
       {
           norms[mol2] = calculate(mol2, regime, mol2, gamma, false);
       }

       s = s/ sqrt(norms[mol1] * norms[mol2]);
    }
    return COEFF * s;
}

AbrahamKernelB::AbrahamKernelB() : CKernel()
{
    name = "AbrahamB";

}

double AbrahamKernelB::calculate(OBMol * mol1, bool regime, OBMol * mol2, double gamma, bool norm)
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
    if(norm)
    {

       if(norms.find(mol1) == norms.end())
       {
           norms[mol1] = calculate(mol1, regime, mol1, gamma, false);
       }

       if(norms.find(mol2) == norms.end())
       {
           norms[mol2] = calculate(mol2, regime, mol2, gamma, false);
       }

       s = s/ sqrt(norms[mol1] * norms[mol2]);
    }
    return COEFF * s;
}

AbrahamKernelE::AbrahamKernelE() : CKernel()
{
    name = "AbrahamE";

}

double AbrahamKernelE::calculate(OBMol * mol1, bool regime, OBMol * mol2, double gamma, bool norm)
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
    if(norm)
    {

       if(norms.find(mol1) == norms.end())
       {
           norms[mol1] = calculate(mol1, regime, mol1, gamma, false);
       }

       if(norms.find(mol2) == norms.end())
       {
           norms[mol2] = calculate(mol2, regime, mol2, gamma, false);
       }

       s = s/ sqrt(norms[mol1] * norms[mol2]);
    }
    return COEFF * s;
}

AbrahamKernelS::AbrahamKernelS() : CKernel()
{
    name = "AbrahamS";

}

double AbrahamKernelS::calculate(OBMol * mol1, bool regime, OBMol * mol2, double gamma, bool norm)
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
    if(norm)
    {

       if(norms.find(mol1) == norms.end())
       {
           norms[mol1] = calculate(mol1, regime, mol1, gamma, false);
       }

       if(norms.find(mol2) == norms.end())
       {
           norms[mol2] = calculate(mol2,regime, mol2, gamma, false);
       }

       s = s/ sqrt(norms[mol1] * norms[mol2]);
    }
    return COEFF * s;
}


