
/*
 * kernel.h
 * Copyright (C) Carpov Pavel   2010 <carpovpv@qsar.chem.msu.ru>
                 Baskin Igor I. 2010 <igbaskin@gmail.com>
 *
 * MCMF is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MCMF is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __KERNEL_H
#define __KERNEL_H

#include <openbabel/mol.h>
#include "descfact.h"
#include <string>
#include <math.h>
#include <map>

/*!
  The abstract class for kernel representation. Simple kernels should
  represent the specific molecular fileds. They are  superimposed  in
  CMFA class. The  DescriptorFactory calculetes specific  descriptors
  for the molecule, then the Kernel computes the similarity.
*/

#define COEFF  (sqrt((M_PI * M_PI * M_PI) / gamma / gamma / gamma))

class CKernel
{
public:

    CKernel( DescriptorFactory * descrfactory = NULL) :
        m_descrfactory(descrfactory)
    {
        norms.clear();
         prev = -1;
    }    

    virtual ~CKernel()
    {

    }

    /*!
      Returns the similarity between two molecules within the field.
    */

    virtual double calculate(OBMol *, bool, OBMol *, double, bool norm = false ) = 0;

    /*!
       Return the name of a kernel.
    */
    const char * getName() const {
        return name.c_str();
    }

    void clearNorms() {norms.clear();}

    void setDescriptorFactory(DescriptorFactory *f) { m_descrfactory = f;}

protected:

    std::string name;

    DescriptorFactory * m_descrfactory;

    std::map < OBMol *, double> norms;

    //prognosis
    double curnorm; //for current molecule in the prognosis mode.
    long prev;


};

#endif

