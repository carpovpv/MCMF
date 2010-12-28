
#include "electro.h"

ElectroStaticKernel::ElectroStaticKernel() : CKernel()
{
    name = "Electro-Static kernel";
}

ElectroStaticKernel::~ElectroStaticKernel()
{
}

double ElectroStaticKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma)
{
    double w1 = 0.0;
    double w2 = 0.0;
    double  s = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {
        double q1 = a->GetPartialCharge();

        double x1 = a->x();
        double y1 = a->y();
        double z1 = a->z();

        w1 += q1 * q1;

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            double q2 = b->GetPartialCharge();

            double x2 = b->x();
            double y2 = b->y();
            double z2 = b->z();

            s += q2 * exp( -gamma / 2.0 * ( (x1-x2) * (x1-x2) + (y1-y2)*(y1-y2) + (z1 -z2) * (z1-z2) ));

        }
        s *= q1;
    }

    FOR_ATOMS_OF_MOL(a, mol2)
    {
        double q = a->GetPartialCharge();
        w2 += q*q ;
    }

    w2 = sqrt(w2);
    w1 = sqrt(w1);

    s /= (w1*w2);

    return s;

}

