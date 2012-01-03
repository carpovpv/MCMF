
#ifndef __HYDROPHOBICV_H
#define __HYDROPHOBICV_H

#include "../kernel.h"
#include <map>

class HydrophobicKernelV: public CKernel
{
public:
    HydrophobicKernelV();
    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training);
private:
    std::map< int, double> hydrophobicity;
    std::map< OBAtom *, int >  cache;

    int isX(OBAtom *);
    int isR(OBAtom *);
    int aromaticHet(OBAtom *);
    int countAl(OBAtom *);
    int countAr(OBAtom *);

    int get_hydrophobicity(OBAtom *);
};

#endif

