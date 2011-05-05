#ifndef SPECTROPHORES_H
#define SPECTROPHORES_H


#include "../descfact.h"
#include <openbabel/spectrophore.h>

class Spectrophores : public DescriptorFactory
{
public:

    Spectrophores();
    ~Spectrophores();
    bool needMapping() const;
    const std::vector < struct Descriptor > & getDescriptors(OBMol *, bool);

private:
    OpenBabel::OBSpectrophore s;

    std::map < OBMol * , std::vector< struct Descriptor > > m_descrs;


};

#endif // SPECTROPHORES_H
