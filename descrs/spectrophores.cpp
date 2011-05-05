#include "spectrophores.h"

Spectrophores::Spectrophores()
{
    s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
    s.SetResolution(3.0);
    s.SetStereo(OpenBabel::OBSpectrophore::AllStereoSpecificProbes);
    s.SetNormalization(OpenBabel::OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd);

    name = "Spectrophores";

}

Spectrophores::~Spectrophores()
{

}

bool Spectrophores::needMapping() const
{
    return false;
}

const std::vector < struct Descriptor > & Spectrophores::getDescriptors(OBMol * mol, bool regime)
{
    //const bool regime = true;

    if(regime)
    {
         std::map < OBMol *, std::vector < struct Descriptor > >::iterator it;
         it = m_descrs.find(mol);

         if(it != m_descrs.end())
         {
             return m_descrs[mol];
         }

    }
    else
    {

        std::string c = mol->GetData("prognosis")->GetValue();
        long cur = atol(c.c_str());

        //printf("MolID=: %ld %ld\n", cur, prev);
        if( prev == cur && descrs.size())
            return descrs;
        prev = cur;
    }

    std::vector<double> result = s.GetSpectrophore(mol);

    descrs.clear();
    descrs.resize(result.size());

    double ss = 0.0;
    for (unsigned int i(0); i < result.size(); ++i)
    {
        descrs[i].value = result[i];
        ss += result[i] * result[i];

    }

    for(int i = 0; i < result.size(); ++i)
    {
        descrs[i].value /= sqrt(ss);
        //printf("%f ", descrs[i].value);
    }

    //printf("Regime: %d %u\n", regime, mol);

    if(regime)
        m_descrs[mol] = descrs;


    return descrs;

}
