#include "linear.h"

LinearKernel::LinearKernel(DescriptorFactory * descr) : CKernel("Linear", descr)
{
}

double LinearKernel::calculate(OBMol * mol1, OBMol * mol2, double, Mode regime)
{
    const std::vector< double >  m1 = m_descrfactory->getDescriptors(mol1, regime);
    const std::vector< double >  m2 = m_descrfactory->getDescriptors(mol2, Training);

    double s = 0.0;
    const unsigned n = m1.size();

    for(int i=0; i< n; ++i)
        s+= m1[i] * m2[i];

    return s;
}
