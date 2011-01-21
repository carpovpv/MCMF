
/*
 * hydropho.cpp
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

#include "hydropho.h"

HydrophobicKernel::HydrophobicKernel() : CKernel()
{
    name =  "Hydrophobic";
    /*
    @ARTICLE{hydropobicity,
      author = {Viswanadhan, V.N. and Ghose, A.K. and Revankar, G.R., and Robins,
    	R.K.},
      title = {Atomic Physicochemical Parameters for Three Dimensional Structure
    	Directed Quantitative Structure-Activity Relationships. 4. Additional
    	Parameters for Hydrophobic and Dispersive Interactions and Their
    	Application for an Automated Superposition of Certain Naturally Occurring
    	Nucleoside Antibiotics},
      journal = {J. Chem. Inf. Comput. Sci.},
      year = {1989},
      volume = {29},
      pages = {163-272},
      owner = {carpovpv},
      timestamp = {2010.11.06}
    }
    */

    hydrophobicity[1] = -0.6771;
    hydrophobicity[2] = -0.4873;
    hydrophobicity[3] = -0.3633;
    hydrophobicity[4] = -0.1366;
    hydrophobicity[5] = -1.0824;
    hydrophobicity[6] = -0.8370;
    hydrophobicity[7] = -0.6015;
    hydrophobicity[8] = -0.5210;
    hydrophobicity[9] = -0.4042;
    hydrophobicity[10] = 0.3651;
    hydrophobicity[11] = -0.5399;
    hydrophobicity[12] = 0.4011;
    hydrophobicity[13] = 0.2263;
    hydrophobicity[14] =0.8282;
    hydrophobicity[15] = -0.1053;
    hydrophobicity[16] = -0.0681;
    hydrophobicity[17] = -0.2287;
    hydrophobicity[18] = -0.3665;
    hydrophobicity[19] = -0.9188;
    hydrophobicity[20] = -0.0082;
    hydrophobicity[21] = -0.1047;
    hydrophobicity[22] = 0.1513;
    hydrophobicity[24] = 0.0068;
    hydrophobicity[25] = 0.1600;
    hydrophobicity[26] = -0.1033;
    hydrophobicity[27] = 0.0598;
    hydrophobicity[28] = 0.1290;
    hydrophobicity[29] = 0.1652;
    hydrophobicity[30] = 0.2975;
    hydrophobicity[31] = 0.9421;
    hydrophobicity[32] = 0.2074;
    hydrophobicity[33] = -0.1774;
    hydrophobicity[34] = -0.2782;
    hydrophobicity[35] = -0.3630;
    hydrophobicity[36] = -0.0321;
    hydrophobicity[37] = 0.3568;
    hydrophobicity[38] = 0.8255;
    hydrophobicity[39] = -0.1116;
    hydrophobicity[40] = 0.0709;
    hydrophobicity[41] = 0.4571;
    hydrophobicity[42] = -0.1316;
    hydrophobicity[43] = 0.0498;
    hydrophobicity[44] = 0.1847;
    hydrophobicity[46] = 0.4418;
    hydrophobicity[47] = 0.3343;
    hydrophobicity[48] = 0.3161;
    hydrophobicity[49] = -0.1488;
    hydrophobicity[50] = -0.3260;
    hydrophobicity[51] = 0.2099;
    hydrophobicity[52] = 0.3695;
    hydrophobicity[53] = 0.2697;
    hydrophobicity[54] = 0.3647;
    hydrophobicity[56] = 0.1402;
    hydrophobicity[57] = 0.4860;
    hydrophobicity[58] = -0.3514;
    hydrophobicity[59] = 0.1720;
    hydrophobicity[60] = 0.2712;
    hydrophobicity[61] = 1.5810;
    hydrophobicity[64] = 0.1473;
    hydrophobicity[66] = 0.1187;
    hydrophobicity[67] = 0.2805;
    hydrophobicity[68] = 0.3954;
    hydrophobicity[69] = 0.3132;
    hydrophobicity[70] = 0.4238;
    hydrophobicity[71] = 0.8678;
    hydrophobicity[72] = -0.0528;
    hydrophobicity[73] = 0.4198;
    hydrophobicity[74] = 0.1461;
    hydrophobicity[75] = -0.1106;
    hydrophobicity[76] = -2.7640;
    hydrophobicity[77] = -2.7919;
    hydrophobicity[78] = 0.5721;
    hydrophobicity[81] = 0.4174;
    hydrophobicity[82] = 0.2167;
    hydrophobicity[83] = 0.2792;
    hydrophobicity[84] = 0.5839;
    hydrophobicity[85] = 0.3425;
    hydrophobicity[86] = 0.9609;
    hydrophobicity[87] = 0.5594;
    hydrophobicity[88] = 0.4656;
    hydrophobicity[89] = 0.9624;
    hydrophobicity[90] = 0.6345;
    hydrophobicity[91] = 1.0242;
    hydrophobicity[92] = 0.4374;
    hydrophobicity[93] = 0.4332;
    hydrophobicity[94] = 1.2362;
    hydrophobicity[95] = 0.9351;
    hydrophobicity[96] = 1.4350;
    hydrophobicity[99] = 1.7018;
    hydrophobicity[100] = 0.9336;
    hydrophobicity[106] = 0.7268;
    hydrophobicity[107] = 0.6145;
    hydrophobicity[108] = 0.3828;
    hydrophobicity[109] = -0.1708;
    hydrophobicity[110] = 0.3717;
    hydrophobicity[116] = -1.6251;
    hydrophobicity[117] = 0.3308;
    hydrophobicity[120] = 0.0236;

}
/*
bool HydrophobicKernel::isC1sp3(OBAtom * a)
{
    int hybr = a->GetHyb();
    if(hybr!=3) return false;
    int s = 0;
    FOR_NBORS_OF_ATOM(nbr, a)
    {
	int atom = nbr->GetAtomicNum();
	if(atom == 1 || atom == 6 ) continue;
	else s++;
    }

    if(s == 1) return true;
    return false;

}
*/

int HydrophobicKernel::isR(OBAtom * atom)
{
    int r = 0;
    FOR_NBORS_OF_ATOM(nbr, atom)
    {
        if(nbr->IsCarbon())
            r++;
    }
    return r;
}

int HydrophobicKernel::countAl(OBAtom * atom)
{
    int r = 0;
    FOR_NBORS_OF_ATOM(nbr, atom)
    {
        if(!nbr->IsAromatic() && nbr->IsCarbon())
            r++;
    }
    return r;

}

int HydrophobicKernel::countAr(OBAtom * atom)
{
    int r = 0;
    FOR_NBORS_OF_ATOM(nbr, atom)
    {
        if(nbr->IsAromatic())
            r++;
    }
    return r;

}

int HydrophobicKernel::aromaticHet(OBAtom *atom)
{
    int r = 0;
    FOR_NBORS_OF_ATOM(nbr, atom)
    {
        if(nbr->IsAromatic() && !nbr->IsCarbon())
            r++;
    }
    return r;
}

int HydrophobicKernel::isX(OBAtom * atom)
{
    int r = 0;
    FOR_NBORS_OF_ATOM(nbr, atom)
    {
        int atomicNum = nbr->GetAtomicNum();
        if(!(atomicNum == 1 || atomicNum == 6))
            r++;
    }
    return r;
}

int HydrophobicKernel::get_hydrophobicity(OBAtom * a)
{

    const int h = a->ImplicitHydrogenCount();
    const int hyb = a->GetHyb();

    const int r = isR(a);
    const int x = isX(a);
    const int ar = aromaticHet(a);
    const bool rr = a->IsAromatic();

    const int cal = countAl(a);
    const int car = countAr(a);

    const int Z = a->GetAtomicNum();
    int temp = 0;


    switch(Z)
    {
        //carbon
    case 6:
        if(h == 4) return 1;
        if(hyb == 3 && h == 3 && x == 0 && r == 1) return 1;
        if(hyb == 3 && h == 2 && x == 0 && r == 2) return 2;
        if(hyb == 3 && h == 1 && x == 0 && r == 3) return 3;
        if(hyb == 3 && r == 4 ) return 4;
        if(hyb == 3 && h == 3 && x == 1 && r == 0) return 5;
        if(hyb == 3 && h == 2 && r == 1 && x == 1) return 6;
        if(hyb == 3 && h == 2 && x == 2) return 7;
        if(hyb == 3 && h == 1 && r == 2 && x == 1) return 8;
        if(hyb == 3 && h == 1 && r == 1 && x == 2) return 9;
        if(hyb == 3 && h == 1 && x == 3) return 10;
        if(hyb == 3 && h == 0 && r == 3 && x == 1) return 11;
        if(hyb == 3 && r == 2 && x == 2) return 12;
        if(hyb == 3 && r == 1 && x == 3) return 13;
        if(hyb == 3 && x == 4) return 14;
        if(hyb == 2 && h == 2) return 15;
        if(hyb == 2 && h == 1 && r == 1) return 16;
        if(hyb == 2 && r == 2) return 17;
        if(hyb == 2 && h == 1 && x ==1) return 18;
        if(hyb == 2 && r == 1 && x == 1) return 19;
        if(hyb == 2 && x == 2) return 20;
        if(hyb == 1 && h == 1) return 21;
        if(hyb == 1 && r == 1) return 22;
        if(hyb == 1 && x == 1) return 23;
        if(rr == true && h == 1 && ar == 0) return 24;
        if(rr == true && x ==0 && ar == 0 && h ==0) return 25;
        if(rr == true && x == 1 && ar == 0) return 26;
        if(rr == true && h == 1 && ar == 1) return 27;
        if(rr == true && h == 0 && ar == 1 && x == 0) return 28;
        if(rr == true && h == 0 && ar == 1 && x == 1) return 29;
        if(rr == true && h == 1 && ar == 2) return 30;
        if(rr == true && h == 0 && ar == 2 && x == 0) return 31;
        if(rr == true && h == 0 && ar == 2 && x == 1) return 32;
        if(hyb == 2 && h == 1 && x == 1 && cal == 1) return 36;
        if(hyb == 2 && h == 1 && x == 1 && car == 1) return 37;
        if(hyb == 2 && h == 0 && x == 1 && cal == 2) return 38;
        if(hyb == 2 && h == 0 && x == 1 && ar == 1 && r == 2) return 39;
        if(hyb == 2 && h == 0 && x == 2 && r == 1) return 40;
        if(hyb == 1 && h == 0 && r == 1 && x == 1) return 40;
        if(hyb == 2 && h == 0 && x == 3) return 41;
        if(rr == true)
        {
            bool single = false;
            FOR_NBORS_OF_ATOM(nbr, a)
            {
                if(nbr->IsAromatic() && !nbr->IsCarbon() && nbr->CountBondsOfOrder(2) == 0)
                {
                    single = true;
                    break;
                }
            }

            if(single)
            {
                if (h == 1 && ar == 2) return 42;
                if (h == 1 && ar == 1 && r == 1) return 33;
                if (x == 3 && ar == 2) return 44;
                if (r == 1 && ar == 2) return 43;
                if (ar == 1 && h == 0) return 34;
                if(ar == 1 && x == 2) return 35;
            }
        }
        break;
        //hydrogen
    case 1:
        FOR_NBORS_OF_ATOM(nbr, a)
        {
            if(nbr->IsCarbon())
            {
                int nl = 0;
                int nhyb = nbr->GetHyb();
                FOR_NBORS_OF_ATOM(nbr1, &*nbr)
                {
                    if(nbr1->IsCarbon() || nbr1->IsHydrogen()) continue;
                    nl += nbr->GetBond(&*nbr1)->GetBO();
                }
                if(nhyb == 3 && nl == 0 )
                {
                    int nnl = 0;
                    int maxl = 0;

                    bool alpha = false;

                    FOR_NBORS_OF_ATOM(nbr1, &*nbr)
                    {
                        if(!nbr1->IsCarbon()) continue;
                        nnl = 0;
                        FOR_NBORS_OF_ATOM(nbr2, &*nbr1)
                        {
                            if(nbr2->IsCarbon() && nbr2->GetHyb() <3 )
                            {
                                alpha = true;
                                continue;
                            }
                            if(nbr2->IsHydrogen()) continue;
                            nnl++;
                        }
                        if(maxl < nnl) maxl = nnl;
                    }
                    if(alpha == true) return 51;
                    if(maxl == 0) return 46;
                    if(maxl == 1) return 52;
                    if(maxl == 2) return 53;
                    if(maxl == 3) return 54;
                    if(maxl >= 4) return 55;
                }
                else if((nhyb == 3 && nl == 1 )|| (nhyb ==2 && nl == 0)) return 47;
                else if((nhyb == 3 && nl == 2) || (nhyb == 2 && nl == 1 ) || ( nhyb == 1 && nl == 0)) return 48;
                else if((nhyb == 3 && nl == 3) || (nhyb == 2 && nl == 2) || ( nhyb == 1 && nl == 3) || (nhyb == 2 && nl == 3)) return 49;

            }
            else
                return 50;
        }


        break;
        //oxygen
    case 8:

        if(cal == 2) return 59;
        if(cal == 1 && h == 1) return 56;
        if(car == 1 && h == 1) return 57;
        if(car == 1 && cal == 1) return 60;
        if(car == 2) return 60;
        if(rr == true) return 60;

        FOR_NBORS_OF_ATOM(nbr, a)
        if(nbr->IsCarbon() && nbr->GetBond(&*a)->GetBO() > 1) return 57;

        break;
        //selenium
    case 34:
        return 64;
        break;
    case 7: //nitrogen

        if(h==2 && cal == 1) return 66;
        if(h == 1 && cal == 2) return 67;
        if(h== 0 && cal == 3) return 68;
        if(car == 1 && h == 2) return 69;
        if(h ==2 && x ==1 ) return 69;
        if(h==1 && cal ==1 && car ==1 ) return 70;
        if(h==0 && cal == 2 && car ==1 ) return 71;
        if(a->IsAmideNitrogen()) return 72;
        FOR_NBORS_OF_ATOM(nbr, a)
        {
            if( !(nbr->IsCarbon() || nbr->IsHydrogen()))
            {
                FOR_NBORS_OF_ATOM(nbr1, &*nbr)
                if(!(nbr1->IsCarbon() || nbr1->IsHydrogen()))
                    if(nbr->GetBond(&*nbr1)->GetBO()  > 1) return 72;
            }
        }

        if(car >= 2 || (rr == true && a->BOSum() == 2)) return 73;
        if( a->HasDoubleBond() || a->HasBondOfOrder(3))  return 74;
        if(rr == true && !a->HasDoubleBond()) 75;
        if(cal == 1 && a->IsNitroOxygen()) return 77;

        FOR_NBORS_OF_ATOM(nbr, a)
        if(nbr->IsOxygen()) return 76;
        if(a->HasDoubleBond()) return 78;
        return 0;
        break;
        //fluorine
    case 9:
    case 17:
    case 35:
    case 53:
        FOR_NBORS_OF_ATOM(nbr, a)
        {
            if(nbr->IsCarbon())
            {
                int nl = -1;
                int nhyb = nbr->GetHyb();

                FOR_NBORS_OF_ATOM(nbr1, &*nbr)
                {
                    if(nbr1->IsCarbon() || nbr1->IsHydrogen()) continue;
                    nl++;
                }

                switch(Z)
                {
                case 9:
                    if(nhyb == 3 && nl == 1) return 81;
                    if(nhyb == 3 && nl == 2) return 82;
                    if(nhyb == 3 && nl == 3) return 83;
                    if(nhyb == 2 && nl ==1 ) return 84;
                    return 85;
                case 17:
                    if(nhyb == 3 && nl == 1) return 86;
                    if(nhyb == 3 && nl == 2) return 87;
                    if(nhyb == 3 && nl == 3) return 88;
                    if(nhyb == 2 && nl ==1 ) return 89;
                    return 90;
                case 35:
                    if(nhyb == 3 && nl == 1) return 91;
                    if(nhyb == 3 && nl == 2) return 92;
                    if(nhyb == 3 && nl == 3) return 93;
                    if(nhyb == 2 && nl ==1 ) return 94;
                    return 95;
                case 53:
                    if(nhyb == 3 && nl == 1) return 96;
                    if(nhyb == 3 && nl == 2) return 97;
                    if(nhyb == 3 && nl == 3) return 98;
                    return 100;
                }
            }
        }

        if(Z == 9) return 85;
        else if(Z == 17) return 90;
        else if(Z == 35) return 95;
        return 100;

        break;
        //sulphur
    case 16:

        if(h==1) return 106;
        temp = 0;

        FOR_NBORS_OF_ATOM(nbr, a)
        {
            if(nbr->IsOxygen())
                temp++;
        }
        if(temp == 2) return 110;
        if(temp == 1) return 109;

        if(a->HasDoubleBond()) return 108;
        return 107;

        break;
        //phosphorus
    case 15:
        if(x == 1) return 116;
        else if(x == 3) return 120;
        return 117;

        break;
    default :
        return 0;

    }
    return 0;
}
HydrophobicKernel::~HydrophobicKernel()
{
}

double HydrophobicKernel::calculate(OBMol * mol1, OBMol * mol2, double gamma)
{
    double w1 = 0.0;
    double w2 = 0.0;
    double  s = 0.0;

    double w12 = 0.0;
    double w22 = 0.0;

    FOR_ATOMS_OF_MOL(a, mol1)
    {
        int l = 0;
        if(cache.find(&*a) == cache.end())
        {
            l = get_hydrophobicity(&*a);
            cache[&*a] = l;
        }
        else l = cache[&*a];

        if(hydrophobicity.find(l) != hydrophobicity.end()) w1 = hydrophobicity[l];
        else w1 = 0;
        w12 += w1 * w1;

        FOR_ATOMS_OF_MOL(b, mol2)
        {
            int k = 0;
            if(cache.find(&*b) == cache.end())
            {
                k = get_hydrophobicity(&*b);
                cache[&*b] = k;
            }
            else k = cache[&*b];

            if(hydrophobicity.find(k) != hydrophobicity.end()) w2 = hydrophobicity[k];
            else w2 = 0;
            s += w1 * w2 * exp ( -gamma / 4.0 * ( pow((a->x() - b->x()), 2) + pow((a->y() - b->y()), 2)  + pow((a->z() - b->z()), 2)  ));
        }
    }

    FOR_ATOMS_OF_MOL(b, mol2)
    {
        int k = get_hydrophobicity(&*b);
        if(hydrophobicity.find(k) != hydrophobicity.end()) w2 = hydrophobicity[k];
        else w2 = 0;
        w22 += w2 * w2;
    }

    w12 = sqrt(w12);
    w22 = sqrt(w22);

    return s / (w12 * w22);

}

