#ifndef SPECTROPHORES_H
#define SPECTROPHORES_H


#include "../descfact.h"
#include <openbabel/spectrophore.h>

class Spectrophores : public DescriptorFactory
{
public:

    Spectrophores();
    ~Spectrophores();
    const std::vector < double > & getDescriptors(OBMol *, Mode mode = Training);
    void load(FILE *fp);

private:
    OpenBabel::OBSpectrophore s;

};

#endif // SPECTROPHORES_H
