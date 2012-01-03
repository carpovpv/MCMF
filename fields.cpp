#include "fields.h"

#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/mol.h>

Fields::Fields() : OBGenericData("mcmf", OBGenericDataType::CustomData0)
{
    ttab.SetFromType("INT");
    ttab.SetToType("SYB");
}

/*
 Parametrisation from triypos.rpm
 */

void Fields::getTripos(const std::string &atom, double *R, double *E)
{
    if(atom == "C.1")   {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "C.2")   {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "C.3")   {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "C.ar")  {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "C.cat") {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "N.1")   {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "N.2")   {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "N.3")   {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "N.4")   {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "N.am")  {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "N.ar")  {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "N.pl3") {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "O.2")   {
        *R = 1.52,   *E = 0.116;
        return;
    }
    if(atom == "O.3")   {
        *R = 1.52,   *E = 0.116;
        return;
    }
    if(atom == "O.co2") {
        *R = 1.52,   *E = 0.116;
        return;
    }
    if(atom == "S.2")   {
        *R = 1.8,    *E = 0.314;
        return;
    }
    if(atom == "S.3")   {
        *R = 1.8,    *E = 0.314;
        return;
    }
    if(atom == "S.o")   {
        *R = 1.7,    *E = 0.314;
        return;
    }
    if(atom == "S.o2")  {
        *R = 1.7,    *E = 0.314;
        return;
    }
    if(atom == "P.3")   {
        *R = 1.8,    *E = 0.314;
        return;
    }
    if(atom == "H")     {
        *R = 1.5,    *E = 0.042;
        return;
    }
    if(atom == "F")     {
        *R = 1.47,   *E = 0.109;
        return;
    }
    if(atom == "Cl")    {
        *R = 1.75,   *E = 0.314;
        return;
    }
    if(atom == "Br")    {
        *R = 1.85,   *E = 0.434;
        return;
    }
    if(atom == "I")     {
        *R = 1.98,   *E = 0.623;
        return;
    }
    if(atom == "Li")    {
        *R = 1.2,    *E = 0.4;
        return;
    }
    if(atom == "Na")    {
        *R = 1.2,    *E = 0.4;
        return;
    }
    if(atom == "K")     {
        *R = 1.2,    *E = 0.4;
        return;
    }
    if(atom == "Ca")    {
        *R = 1.2,    *E = 0.6;
        return;
    }
    if(atom == "Al")    {
        *R = 1.2,    *E = 0.042;
        return;
    }
    if(atom == "Si")    {
        *R = 2.1,    *E = 0.042;
        return;
    }
    if(atom == "O.t3p") {
        *R = 1.76825,*E = 0.15207;
        return;
    }
    if(atom == "H.t3p") {
        *R = 1.008,  *E = 0;
        return;
    }
    if(atom == "O.spc") {
        *R = 1.7766, *E = 0.1554;
        return;
    }
    if(atom == "H.spc") {
        *R = 1.008,  *E = 0;
        return;
    }
    if(atom == "LP")    {
        *R = 1e-05,  *E = 1e-05;
        return;
    }
    if(atom == "Du")    {
        *R = 0,      *E = 0;
        return;
    }
    if(atom == "Du.C")  {
        *R = 0,      *E = 0;
        return;
    }
    if(atom == "ANY")   {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "HEV")   {
        *R = 1.7,    *E = 0.107;
        return;
    }
    if(atom == "HET")   {
        *R = 1.55,   *E = 0.095;
        return;
    }
    if(atom == "HAL")   {
        *R = 1.75,   *E = 0.314;
        return;
    }

    *R = *E = 0.0;
    return ;
}

void Fields::calcValues(OBAtom * atom)
{

    std::string src,dst;

    src = atom->GetType();
    ttab.Translate(dst,src);

    double h = 0.0; //hydrophobicity
    double a = 0.0; //Abraham a
    double b = 0.0; //Abraham b
    double e = 0.0; //Abraham e
    double s = 0.0; //Abraham s

    int ihydr = atom->ImplicitHydrogenCount();

    // the code below was transformed from mol2.c  developed by Baskin I.I.
    if(dst == "C.1")
    {
        h += 0.041270;  /* p1.C__ */
        e += -0.005210; /* p1.C__ */
        s += -0.013906; /* p1.C__ */
        if (ihydr == 1)
        {
            h -= 0.093156; /* p1.CT1 */
            s += 0.045181; /* p1.CT1 */
        }
    }
    else if(dst == "C.2")
    {
        h += 0.041270; /* p1.C__ */
        e += -0.005210; /* p1.C__ */
        s += -0.013906; /* p1.C__ */
        h += 0.006410; /* p.CD_ */
        s += 0.016508; /* p.CD_ */
        if (ihydr == 0)
        {
            h += 0.007629; /* p.CD3 */
            b += 0.026654; /* p.CD3 */
            e += 0.076572; /* p.CD3 */
        }
        else if (ihydr == 1)
        {
            h += 0.202552; /* p.CD2 */
            e += 0.058299; /* p.CD2 */
        }
        else
        {
            h += 0.422564;  /* p.CD1 */
            b += -0.041897; /* p.CD1 */
            e += -0.019698; /* p.CD1 */
            s += -0.098377; /* p.CD1 */
        }
    }
    else if (dst == "C.3")
    {
        h += 0.041270; /* p1.C__ */
        e += -0.005210; /* p1.C__ */
        s += -0.013906; /* p1.C__ */
        a += -0.009684; /* p1.CA_ */
        if (ihydr == 0)
        {
            h += 0.089007; /* p1.CA4 */
            b += 0.096834; /* p1.CA4 */
            e += 0.153114; /* p1.CA4 */
        } else if (ihydr == 1)
        {
            b += 0.041221; /* p1.CA3 */
            e += 0.091218; /* p1.CA3 */
        }
        else if (ihydr == 2)
            h += 0.208288; /* p1.CA2 */
        else
        {
            h += 0.369764;  /* p1.CA1 */
            b += -0.022123; /* p1.CA1 */
            e += -0.102229; /* p1.CA1 */
            s += -0.055649; /* p1.CA1 */
        }
    }
    else if (dst == "C.ar")
    {
        h += 0.041270; /* p1.C__ */
        e += -0.005210; /* p1.C__ */
        s += -0.013906; /* p1.C__ */
        h += 0.194585; /* p1.CB_ */
        e += 0.117666; /* p1.CB_ */
        s += 0.070704; /* p1.CB_ */
        if (ihydr == 0)
        {
            h += 0.085378; /* p1.CB2 */
            b += 0.002903; /* p1.CB2 */
            e += 0.107527; /* p1.CB2 */
            s += 0.046501; /* p1.CB2 */
        } else
        {
            e += -0.044799; /* p1.CB1 */
            s += -0.011864; /* p1.CB1 */
        }
    }
    else if (dst == "C.cat")
    {
        e += -0.005210; /* p1.C__ */
        s += -0.013906; /* p1.C__ */
    }
    else if (dst == "N.1")
    {
        a += 0.077667; /*p1.N__*/
        s += 0.167948; /*p1.N__*/
        b += 0.061824; /*p1.NT_*/
        s += 0.232745; /*p1.NT_*/
    }
    else if (dst == "N.2")
    {
        a += 0.077667; /*p1.N__*/
        s += 0.167948; /*p1.N__*/
        if (ihydr == 1) {
            h += -0.200643; /* p1.ND1 */
        }
        else
        {
            b += -0.763354; /* p1.ND2 */
            e += 0.682859; /* p1.ND2 */
        }
    }
    else if (dst == "N.3")
    {
        a += 0.077667; /*p1.N__*/
        s += 0.167948; /*p1.N__*/
        h += -0.091021; /* p1.NA_ */
        b += 0.367558; /* p1.NA_ */
        e += 0.099602; /* p1.NA_ */
        s += 0.084806; /* p1.NA_ */
        if (ihydr == 0)
        {
            h += -0.495588; /* p.NA3 */
            b += 0.150946; /* p.NA3 */
            e += 0.088502; /* p.NA3 */
        } else if (ihydr == 1)
        {
            h += -0.452520; /* p.NA2 */
            a += 0.163329;
        } else if (ihydr == 2)
        {
            h += -0.763440; /* p.NA1 */
            a += 0.120801; /* p.NA1 */
            s += -0.069707; /* p.NA1 */
        }
    }
    else if (dst == "N.4")
    {
        a += 0.077667; /*p1.N__*/
        s += 0.167948; /*p1.N__*/
        h += -4.382696; /* p1.NC_ */
    }
    else if (dst == "N.am")
    {
        a += 0.077667; /*p1.N__*/
        s += 0.167948; /*p1.N__*/
        h += -0.091021; /* p1.NA_ */
        if (ihydr == 0)
        {
            h += -0.495588; /* p.NA3 */
            a += -0.096006;
        }
        else if (ihydr == 1)
            h += -0.452520; /* p.NA2 */
        else if (ihydr == 2)
            h += -0.763440; /* p.NA1 */
    }
    else if (dst == "N.ar")
    {
        a += 0.077667;  /*p1.N__*/
        s += 0.167948;  /*p1.N__*/
        a += -0.104421; /*p1.NB_*/
        b += 0.282666;  /*p1.NB_*/
    }
    else if (dst == "N.pl3")
    {
        a += 0.077667; /*p1.N__*/
        s += 0.167948; /*p1.N__*/
        h += 0.483597; /* p1.N5_ */
        b += -0.296874; /* p1.N5_ */
        e += 0.060327; /* p1.N5_ */
        s += -0.351281; /* p1.N5_ */
    }
    else if (dst == "O2")
    {
        h += -0.224286; /* p1.O__ */
        h += -0.098496; /* p1.OD_ */
        a += 0.038390; /* p1.OD_ */
        e += -0.011749; /* p1.OD_ */
        s += 0.371065; /* p1.OD_ */
    }
    else if (dst == "O3")
    {
        h += -0.224286; /* p1.O__ */
        h += -0.102250; /* p1.OA_ */
        s += 0.111253; /* p1.OA_ */
        b += 0.147962; /* p1.O__ */
        if (ihydr == 0) {
            h += -0.010867; /* p1.OA2 */
            a += -0.035354; /* p1.OA2 */
            e += -0.020649; /* p1.OA2 */
        } else {
            h += -0.279404; /* p1.OA1 */
            a += 0.445333; /* p1.OA1 */
            b += 0.014953; /* p1.OA1 */
            e += 0.010690; /* p1.OA1 */
            s += 0.027754; /* p1.OA1 */
        }
    }
    else if (dst == "O.co2")
        h = -0.224286; /* p1.O__ */
    else if (dst == "S.2")
    {
        h += -0.077480; /* p1.SD_ */
        h += 0.734821; /* p1.SD1 */
        e += 0.136585; /* p1.SD_ */
    }
    else if (dst == "S.3")
    {
        e += 0.195277; /* p1.S__ */
        s += 0.082278; /* p1.SA_ */
        if (ihydr == 0) {
            h += 0.504755; /* p1.SA2 */
            b += 0.103649; /* p1.SA2 */
            e += 0.071103; /* p1.SA2 */
        }
    }
    else if (dst == "S.o")
    {
        h += -0.077480; /* p1.SD_ */
        h += -1.055756; /* p1.SD3 */
    }
    else if (dst == "S.o2")
    {
        h += -0.077480; /* p1.SD_ */
        h += 0.027576; /* p1.SD4 */
    }
    else if (dst == "P.3")
    {
        b += 0.686142; /* p1.P__ */
        e += 0.103799; /* p1.P__ */
        s += 0.215383; /* p1.P__ */
    }
    else if (dst == "F")
    {
        h += 0.269012; /* p1.F_1 */
        b += -0.077060; /* p1.F_1 */
        e += -0.220887; /* p1.F_1 */
        s += -0.077487; /* p1.F_1 */
    }
    else if (dst == "Cl")
    {
        h += 0.569605; /* p1.Cl1 */
        b += -0.068040; /* p1.Cl1 */
        e += -0.008902; /* p1.Cl1 */
        s += 0.029923; /* p1.Cl1 */
    }
    else if (dst == "Br")
    {
        h += 0.696413; /* p1.Br1 */
        b += -0.074564; /* p1.Br1 */
        e += 0.167323; /* p1.Br1 */
        s += 0.131222; /* p1.Br1 */
    }
    else if (dst == "I")
    {
        h += 0.323219; /* p1.I__ */
        e += 0.509975; /* p1.I__ */
        s += 0.135703; /* p1.I__ */
    }
    else if (dst == "Si")
        h += 0.617895; /* Si4 */

    values[Hydrophobic] = h;
    values[AbrahamA] = a;
    values[AbrahamB] = b;
    values[AbrahamE] = e;
    values[AbrahamS] = s;

    double R, E;
    getTripos(dst, &R,&E);

    values[StericR] = R;
    values[StericE] = E;

}

double Fields::getValue(FieldType type)
{
    return values[type];
}

