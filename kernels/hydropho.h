
#ifndef __HYDROPHOBIC_H
#define __HYDROPHOBIC_H

#include "../kernel.h"
#include <map>

class HydrophobicKernel: public CKernel
{
public:
    HydrophobicKernel();
    virtual ~HydrophobicKernel();
    virtual double calculate(OBMol *, OBMol *, double );
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

