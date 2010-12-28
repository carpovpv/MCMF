
/*
 * ks_seal.cpp
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

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "ks_seal.h"
#include <nlopt.h>
#include <vector>

using namespace std;

double KS_Seal::score_function(unsigned , const double *p, double *, void * parent)
{
    KS_Seal * cls = static_cast< KS_Seal * > (parent);

    cls->m_cycles++;

    OBMol _supMol( cls->supMol);

    const int N1 = cls->mainMol.NumAtoms();
    const int N2 = cls->supMol.NumAtoms();

    double qv[9];

    vector3 t (p[4], p[5], p[6]);
    cls->createQ(p, qv);

    _supMol.Rotate(qv);
    _supMol.Translate(t);

    double s = 0.0;

    for(int i=1; i<= N1; ++i)
        for(int j=1; j<= N2; ++j)
        {

            OBAtom * atom_i = cls->mainMol.GetAtom(i);
            OBAtom * atom_j = _supMol.GetAtom(j);

            //charges
            double fch_i = atom_i->GetPartialCharge();
            double fch_j = atom_j->GetPartialCharge();

            //vdw
            double v_i = etab.GetVdwRad(atom_i->GetAtomicNum());
            double v_j = etab.GetVdwRad(atom_j->GetAtomicNum());

            double w_ij = cls->wE * fch_i * fch_j + cls->wS * v_i * v_j;

            vector3 coords_i = atom_i->GetVector();
            vector3 coords_j = atom_j->GetVector();

            double distance = coords_i.distSq(coords_j);
            s += w_ij *  exp( - cls->alpha * distance);
        }

    s *= -1.0;
//    fprintf(stderr, "SL %g\n", s);
    return s;
}

KS_Seal::KS_Seal(const char * fp,
                 double _alpha,
                 double _wE,
                 double _wS,
                 int _probes)
    : SEAL(fp),
      alpha(_alpha),
      wE(_wE),
      wS(_wS),
      probes(_probes)
{
    srand(time(NULL));

    opt = nlopt_create(NLOPT_LN_BOBYQA, 7);
    nlopt_set_min_objective(opt, score_function, this);
    nlopt_set_xtol_rel(opt, 1e-4);

    double min_constraint[] = {  -0.5, -0.5 , -0.5, -0.5, -100, -100, -100};
    double max_constraint[] = {  10, 10 , 10, 10, 100, 100, 100};

    nlopt_set_lower_bounds( opt,min_constraint);
    nlopt_set_upper_bounds( opt,max_constraint);

}

double KS_Seal::optim (double *start, double *xmin )
{

    int i;
    double minf;

    int err ;
    if ((err = nlopt_optimize(opt, start, &minf)) < 0)
        fprintf(stderr,"nlopt failed %d!\n", err);

    for(int i=0; i< 7; ++i)
        xmin[i] = start[i];

    fprintf(stderr, "minf: %g\n", minf);
    return minf;

}

double KS_Seal::norm(const double *q)
{
    return (q[0]*q[0] +q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

void KS_Seal::createQ(const double *q, double *Q)
{

//  fprintf(stderr,"qaternion: %g %g %g %g\n", q[0], q[1], q[2], q[3]);

    Q[0] = q[0] * q[0] + q[1]* q[1] - q[2]* q[2] - q[3] *q[3];
    Q[1] = 2.0 * ( q[1]*q[2] +q[0]*q[3]);
    Q[2] = 2.0 * ( q[1]*q[3] - q[0]*q[2]);
    Q[3] = 2.0 * ( q[1]* q[2] - q[0]* q[3]);
    Q[4] = q[0] * q[0] + q[2]* q[2] - q[1]* q[1] - q[3]* q[3];
    Q[5] = 2.0 * (q[2] *q[3] +q[0]*q[1]);
    Q[6] = 2.0 * (q[1] *q[3] +q[0]* q[2]);
    Q[7] = 2.0 * ( q[2] *q[3] -q[0]*q[1] );
    Q[8] = q[0]* q[0] +q[3]*q[3] -q[1]*q[1] - q[2]*q[2];

    double p = norm(q);

    for(int i=0; i< 9; ++i)
        Q[i] = Q[i] / p;

}

double KS_Seal::score(unsigned int st1, unsigned int st2)
{

    mainMol = mols[st1];
    const int N1 = mainMol.NumAtoms();
    const int N2 = mols[st2].NumAtoms();

    double smin = RAND_MAX;
    double sp[7];

    double * coordinates = new double [N2 * 3] ;

    for(int i=0; i< N2*3; ++i)
        coordinates[i] = 0.0;

//   for(int n_mol = 0; n_mol < (st2 - st1); ++n_mol)
    int n_mol = 0;
    {
        vector<int> was1, was2;
        for(int l = 0; l < probes; ++l)
        {

            supMol  = mols[st2];
            double x0 = supMol.GetAtom(1)->x();
            //        std::cerr << "Atom 0 " << x0 << std::endl;

            srand(time(NULL));

            unsigned int for_t1 ;

            const int max_count = 1000;
            int cc = 0;

            bool k = false;
            do {
                for_t1 = 1 + (double) rand() / RAND_MAX * N1;
                k = false;
                if(find(was1.begin(), was1.end(), for_t1) == was1.end())
                {
                    k=true;
                    was1.push_back(for_t1);
                }
                cc++;
            } while(cc < max_count && k==false);

            unsigned int for_t2;
            cc = 0;
            do {
                for_t2 = 1 + (double) rand() / RAND_MAX * N2;
                k = false;
                if(find(was2.begin(), was2.end(), for_t1) == was2.end())
                {
                    k=true;
                    was2.push_back(for_t2);
                }
                cc++;
            } while(cc < max_count && k==false);


            cerr << for_t1 << " " << for_t2 << endl;

            const OBAtom * const atom1 = mainMol.GetAtom(for_t1);
            const OBAtom * const atom2 = supMol.GetAtom(for_t2);

            double q[7];

            q[4] = atom2->GetX() - atom1->GetX();
            q[5] = atom2->GetY() - atom1->GetY();
            q[6] = atom2->GetZ() - atom1->GetZ();

            q[0] = (double) rand() / RAND_MAX  - 0.5;
            q[1] = (double) rand() / RAND_MAX  - 0.5;
            q[2] = (double) rand() / RAND_MAX  - 0.5;
            q[3] = (double) rand() / RAND_MAX  - 0.5;

            double p[7];

            double s = optim(q,p);
            if(s < smin)
            {
                smin = s;
                for(int i=0; i< 7; ++i)
                    sp[i] = p[i];
            }


        }

        double qv[9];

        vector3 t (sp[4], sp[5], sp[6]);
        createQ(sp, qv);

        supMol.Rotate(qv);
        supMol.Translate(t);

        double * n_coordinates =  supMol.GetCoordinates();
        for(int i=0; i< N2*3; ++i)
            coordinates[i] += n_coordinates[i];

    }

    for(int i=0; i< N2*3; ++i)
    {
        coordinates[i] = 0.8;//(-0.5 )*coordinates[i] / (st2 - st1);
//	    fprintf(stderr,"coord %g\n", coordinates[i]);
    }

//   fprintf(stderr, "ST@: %d\n", st2 - st1);

    //supMol.SetCoordinates(coordinates);

    //m_mols.push_back(supMol);
    conv->Write(&supMol);

    delete [] coordinates;

}

KS_Seal::~KS_Seal()
{
    nlopt_destroy(opt);
}

void KS_Seal::align()
{
    if(debug)
        std::cerr << "Starting alignment process...\n";

    m_mols.push_back(mols[0]);
    conv->Write(&mols[0]);

    for(int first = 1; first < mols.size(); ++first)
    {
        std::cerr << "Align " << first << " structure. " ;
        m_cycles = 0;
        score(0, first);
        std::cerr << "Cycles " <<  m_cycles << std::endl;
    }
}

