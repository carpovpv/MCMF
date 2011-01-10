#ifndef LINEAR_H
#define LINEAR_H

#include "../kernel.h"

class LinearKernel: public CKernel
{
public:
    LinearKernel(DescriptorFactory *);
    double calculate(OBMol *, OBMol *, double );
};
#endif // LINEAR_H
