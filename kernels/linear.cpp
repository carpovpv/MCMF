#include "linear.h"

LinearKernel::LinearKernel(DescriptorFactory * descr) : CKernel( descr)
{
    name = "Linear";
    if(descr != NULL)
        name += ":" + descr->getName();
}

double LinearKernel::calculate(OBMol * mol1, bool regime, OBMol * mol2, double gamma, bool norm)
{
    const std::vector< struct Descriptor>  m1 = m_descrfactory->getDescriptors(mol1, regime);
    const std::vector< struct Descriptor>  m2 = m_descrfactory->getDescriptors(mol2, false);

    double s = 0.0;
    const unsigned n = m1.size();

    for(int i=0; i< n; ++i)
        s+= m1[i].value * m2[i].value;

    return s;
}
