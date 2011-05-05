
#ifndef __STERICK_H
#define __STERICK_H

#include "../kernel.h"
#include <map>

class StericKernel: public CKernel
{
public:
    StericKernel();  
    virtual double calculate(OBMol *, bool,  OBMol *, double , bool norm = false);
};

#endif

