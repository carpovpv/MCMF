
#ifndef __HYDROPHOBIC_H
#define __HYDROPHOBIC_H

#include "../kernel.h"
#include <map>

class HydrophobicKernel: public CKernel
{
public:
    HydrophobicKernel();
    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training);

};

#endif

