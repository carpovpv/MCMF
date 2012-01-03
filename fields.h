#ifndef FIELDS_H
#define FIELDS_H

#include <openbabel/generic.h>

using namespace OpenBabel;

class Fields : public OBGenericData
{
public:
    Fields();
    enum FieldType { Hydrophobic,
                     StericR,
                     StericE,
                     AbrahamA,
                     AbrahamB,
                     AbrahamE,
                     AbrahamS
                   };

    double getValue(FieldType type);
    void calcValues(OBAtom *);

protected:

    std::map < Fields::FieldType, double > values;

private:

    void getTripos(const std::string & atom, double * R, double *E);
};

#endif // FIELDS_H
