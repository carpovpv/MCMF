
#include "electro.h"

ElectroStaticKernel::ElectroStaticKernel() : CKernel("ElectroStatic")
{

}

ElectroStaticKernel::~ElectroStaticKernel()
{
}

//ядро зависит от gamma, поэтому его надо все время вычислять
double ElectroStaticKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma, Mode)
{
    double  s = 0.0;
    FOR_ATOMS_OF_MOL(a, mol1)
    {
        double q1 = a->GetPartialCharge();
        FOR_ATOMS_OF_MOL(b, mol2)
        {
            double q2 = b->GetPartialCharge();

            double x = b->x() - a->x();
            double y = b->y() - a->y();
            double z = b->z() - a->z();

            s += q2 * q1 * exp( -gamma / 4.0 * ( x*x + y*y + z*z ));
        }
    }
    s *= COEFF(gamma);
    return s;
}

