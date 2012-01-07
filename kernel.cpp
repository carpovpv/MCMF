#include "kernel.h"

CKernel::CKernel(const std::string &name, DescriptorFactory *descrfactory)
    : m_name(name),
      m_descrfactory(descrfactory)
{
    fprintf(stderr, "Loading kernel %s.\n", name.c_str());
    kerncode = D_UNKNOWNKERNEL;
}

void CKernel::save(FILE *fp) const
{    
    if(m_descrfactory != NULL)
        m_descrfactory->save(fp);
    else
        fprintf(fp, "Descriptors: %d 0 (None)\n", D_UNKNOWNDESCR);
    fprintf(fp, "Kernel: %d (%s)\n", kerncode, m_name.c_str());
}

CKernel::~CKernel()
{
    if(m_descrfactory)
        delete m_descrfactory;
}
const std::string & CKernel::getName() const
{
    return m_name;
}

