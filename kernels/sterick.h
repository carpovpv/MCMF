
#ifndef __STERIC_H
#define __STERIC_H

#include "../kernel.h"
#include <map>

class StericKernelK: public CKernel
{
public:
    StericKernelK();
    virtual double calculate(OBMol *,  OBMol *, double, Mode mode = Training);
    double steric(int atom1, int atom2);

private:
    /*!
    	map[i][j] -> i,j - atoms' numbers
    */

    enum {Atom_H = 1, Atom_C = 6, Atom_N = 7, Atom_O = 8, Atom_F = 9, Atom_S = 16,
          Atom_CL = 17, Atom_BR = 35, Atom_I = 53
         };

    std::map < int, std::map < int, double > > table;

};

#endif

