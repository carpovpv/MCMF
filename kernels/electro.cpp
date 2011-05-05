
#include "electro.h"

ElectroStaticKernel::ElectroStaticKernel() : CKernel()
{
    name = "Electro-Static";
}

ElectroStaticKernel::~ElectroStaticKernel()
{
}

double ElectroStaticKernel::calculate(OBMol * mol1, bool regime, OBMol * mol2, double gamma, bool norm)
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
    s *= COEFF;

    if(norm)
    {        

        if(norms.find(mol2) == norms.end())
            norms[mol2] = calculate(mol2, regime, mol2, gamma, false);

        if(regime)
        {
            if(norms.find(mol1) == norms.end())
                norms[mol1] = calculate(mol1, regime, mol1, gamma, false);
        }
        else
        {
            std::string c = mol1->GetData("prognosis")->GetValue();
            long cur = atol(c.c_str());

            //printf("MolID=: %ld %ld\n", cur, prev);
            if( prev == cur)
                return s / sqrt(curnorm * norms[mol2]);

            prev = cur;
            curnorm  = calculate(mol1, regime, mol1, gamma, false);
            return s / sqrt(curnorm * norms[mol2]);
        }

        s = s/ sqrt(norms[mol1] * norms[mol2]);
    }
    //printf("Kernel : %d %d %g\n", mol1, mol2, s);

    return s;
}

