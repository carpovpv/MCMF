
#ifndef __ABRAHAM_H
#define __ABRAHAM_H

#include "../kernel.h"
#include <map>

class AbrahamKernel: public CKernel
{
public:
    AbrahamKernel();
    virtual double calculate(OBMol *, OBMol *, double , bool norm = false);
};

#endif

