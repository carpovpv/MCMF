#include "spectrophores.h"

Spectrophores::Spectrophores()
    :DescriptorFactory("Spectrophores")
{
    s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
    s.SetResolution(3.0);
    s.SetStereo(OpenBabel::OBSpectrophore::AllStereoSpecificProbes);
    s.SetNormalization(OpenBabel::OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd);

    descrcode = D_SPECTROPHORES;
}

Spectrophores::~Spectrophores()
{

}

const std::vector < double > & Spectrophores::getDescriptors(OBMol * mol, Mode regime)
{
    //const bool regime = true;

    if(regime == Training)
    {
        std::map < OBMol *, std::vector < double > >::iterator it;
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

        if( prev == cur)
            return descrs;
        prev = cur;
    }

    std::vector<double> result = s.GetSpectrophore(mol);

    descrs.clear();
    descrs.resize(result.size());

    double ss = 0.0;
    for (unsigned int i(0); i < result.size(); ++i)
    {
        descrs[i] = result[i];
        ss += result[i] * result[i];

    }

    for(int i = 0; i < result.size(); ++i)
    {
        descrs[i] /= sqrt(ss);
        //printf("%f ", descrs[i].value);
    }

    //printf("Regime: %d %u\n", regime, mol);

    if(regime == Training)
        m_descrs[mol] = descrs;


    return descrs;

}
