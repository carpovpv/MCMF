#include "linear.h"

LinearKernel::LinearKernel(DescriptorFactory * descr) : CKernel( descr)
{
    name = "Linear";
}

double LinearKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma)
{
    const std::vector< struct Descriptor>  m1 = m_descrfactory->getDescriptors(mol1);
    const std::vector< struct Descriptor>  m2 = m_descrfactory->getDescriptors(mol2);

    double s = 0.0;
    const unsigned n = m1.size();

    for(int i=0; i< n; ++i)
        s+= m1[i].value * m2[i].value;

    return s;
}
