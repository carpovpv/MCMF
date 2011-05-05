
#ifndef __ABRAHAM_H
#define __ABRAHAM_H

#include "../kernel.h"
#include <map>

class AbrahamKernelA: public CKernel
{
public:
    AbrahamKernelA();
    virtual double calculate(OBMol *, bool, OBMol *, double , bool norm = false);
};

class AbrahamKernelB: public CKernel
{
public:
    AbrahamKernelB();
    virtual double calculate(OBMol *, bool, OBMol *, double , bool norm = false);
};

class AbrahamKernelE: public CKernel
{
public:
    AbrahamKernelE();
    virtual double calculate(OBMol *, bool, OBMol *, double , bool norm = false);
};

class AbrahamKernelS: public CKernel
{
public:
    AbrahamKernelS();
    virtual double calculate(OBMol *, bool, OBMol *, double , bool norm = false);
};


#endif

