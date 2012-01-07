
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
#include <boost/utility.hpp>

#include <string>
#include <math.h>
#include <map>
#include <stdexcept>

#include "descfact.h"

class KernelFailed : public std::runtime_error
{
public:
    explicit KernelFailed(const std::string& __arg):
            runtime_error(__arg)
    {

    }
};

/*!
  The abstract class for kernel representation. Simple kernels should
  represent the specific molecular fileds. They are  superimposed  in
  CMFA class. The  DescriptorFactory calculetes specific  descriptors
  for the molecule, then the Kernel computes the similarity.
*/

//norms for MCMF approach
inline double COEFF(double gamma)
{
    return M_PI / gamma * sqrt( M_PI / gamma);
}

class CKernel : boost::noncopyable
{
public:


    CKernel(const std::string & name, DescriptorFactory * descrfactory = NULL);
    virtual ~CKernel();

    /*!
      Return the similarity between two molecules within the field.
    */

    virtual double calculate(OBMol *, OBMol *, double, Mode mode = Training ) = 0;

    const std::string & getName() const;

    void save(FILE * fp) const;

protected:

    std::string m_name;
    DescriptorFactory * m_descrfactory;

    int kerncode;

};

#endif

