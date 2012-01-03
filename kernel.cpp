#include "kernel.h"

CKernel::CKernel(const std::string &name, DescriptorFactory *descrfactory)
    : m_name(name),
      m_descrfactory(descrfactory)
{

}

void CKernel::save(FILE *fp) const
{
    fprintf(fp, "Kernel:%s\n", m_name.c_str());
    if(m_descrfactory != NULL)
        m_descrfactory->save(fp);
}

CKernel::~CKernel()
{

}
const std::string & CKernel::getName() const
{
    return m_name;
}

