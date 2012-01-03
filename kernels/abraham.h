
#ifndef __ABRAHAM_H
#define __ABRAHAM_H

#include "../kernel.h"
#include <map>

class AbrahamKernelA: public CKernel
{
public:
    AbrahamKernelA();
    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training);
};

class AbrahamKernelB: public CKernel
{
public:
    AbrahamKernelB();
    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training);
};

class AbrahamKernelE: public CKernel
{
public:
    AbrahamKernelE();
    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training);
};

class AbrahamKernelS: public CKernel
{
public:
    AbrahamKernelS();
    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training);
};


#endif

