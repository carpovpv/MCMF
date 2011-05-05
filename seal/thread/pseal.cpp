
#include <string>
#include <stdexcept>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

#include <pthread.h>
#include <nlopt.h>

using namespace OpenBabel;

double wE = 1.0;
double wS = 1.0;
double alpha = 0.5;
int nprobes = 10;

OBMol *templ;

OBConversion conv;

double inline norm(const double *q)
{
    return (q[0]*q[0] +q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

void inline createQ(const double *q, double *Q)
{

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

double score_function(unsigned , const double *p, double *, void * mol)
{
    double s = 0.0;
    OBMol * pmol = static_cast<OBMol *> ( mol);

    const int N1 = templ->NumAtoms();
    const int N2 = pmol->NumAtoms();

    double * coords = new double[ 3 * N2] ;
    double * m = pmol->GetCoordinates();

    for(int i = 0 ; i < 3*N2; ++i)
        coords[i] = m[i];

    double qv[9];

    vector3 t (p[4], p[5], p[6]);
    createQ(p, qv);

    for (int i = 0;i < N2;++i)
    {
            double x = m[i*3  ];
            double y = m[i*3+1];
            double z = m[i*3+2];
            m[i*3  ] = qv[0]*x + qv[1]*y + qv[2]*z + t[0];
            m[i*3+1] = qv[3]*x + qv[4]*y + qv[5]*z + t[1];
            m[i*3+2] = qv[6]*x + qv[7]*y + qv[8]*z + t[2];
    }

    for(int i=1; i<= N1; ++i)
        for(int j=1; j<= N2; ++j)
        {


            OBAtom * atom_i = templ->GetAtom(i);
            OBAtom * atom_j = pmol->GetAtom(j);

            //charges
            double fch_i = atom_i->GetPartialCharge();
            double fch_j = atom_j->GetPartialCharge();

            //vdw
            double v_i = etab.GetVdwRad(atom_i->GetAtomicNum());
            double v_j = etab.GetVdwRad(atom_j->GetAtomicNum());

            double w_ij = wE * fch_i * fch_j + wS * v_i * v_j;

            vector3 coords_i = atom_i->GetVector();
            vector3 coords_j = atom_j->GetVector();

            double distance = coords_i.distSq(coords_j);
            s += w_ij *  exp( - alpha * distance);
        }

    s *= -1.0;

    pmol->SetCoordinates(coords);
//    fprintf(stderr, "SL %g\n", s);

    delete [] coords;
    return s;
}

double optim (double *start, double *xmin , nlopt_opt * opt )
{
    double minf;

    int err ;
    if ((err = nlopt_optimize(*opt, start, &minf)) < 0)
        fprintf(stderr,"nlopt failed %d!\n", err);

    for(int i=0; i< 7; ++i)
        xmin[i] = start[i];

    //fprintf(stderr, "minf: %g %d\n", minf, opt);
    return minf;
}

void *pseal_align(void * pmol)
{
    //fprintf(stderr,"Thread\n");
    OBMol * mol = static_cast<OBMol *> ( pmol );

    //printf("Mol %d\n", mol->NumAtoms());

    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, 7);
    nlopt_set_min_objective(opt, score_function, pmol);
    nlopt_set_xtol_rel(opt, 1e-3);

    double min_constraint[] = {  -0.5, -0.5 , -0.5, -0.5, -100, -100, -100};
    double max_constraint[] = {  10, 10 , 10, 10, 100, 100, 100};

    nlopt_set_lower_bounds( opt,min_constraint);
    nlopt_set_upper_bounds( opt,max_constraint);

////////

    const int N2 = mol->NumAtoms();
    const int N1 = templ->NumAtoms();

    double smin = RAND_MAX;
    double sp[7];

    std::vector<int> was1, was2;
    for(int l = 0; l < nprobes; ++l)
    {

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


        //std::cerr << for_t1 << " " << for_t2 << std::endl;

        const OBAtom * const atom1 = templ->GetAtom(for_t1);
        const OBAtom * const atom2 = mol->GetAtom(for_t2);

        double q[7];

        q[4] = atom2->GetX() - atom1->GetX();
        q[5] = atom2->GetY() - atom1->GetY();
        q[6] = atom2->GetZ() - atom1->GetZ();

        q[0] = (double) rand() / RAND_MAX  - 0.5;
        q[1] = (double) rand() / RAND_MAX  - 0.5;
        q[2] = (double) rand() / RAND_MAX  - 0.5;
        q[3] = (double) rand() / RAND_MAX  - 0.5;

        double p[7];

        double s = optim(q,p,&opt);
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

    double * m = mol->GetCoordinates();
    for (int i = 0;i < N2;++i)
    {
            double x = m[i*3  ];
            double y = m[i*3+1];
            double z = m[i*3+2];
            m[i*3  ] = qv[0]*x + qv[1]*y + qv[2]*z + t[0];
            m[i*3+1] = qv[3]*x + qv[4]*y + qv[5]*z + t[1];
            m[i*3+2] = qv[6]*x + qv[7]*y + qv[8]*z + t[2];
    }

    nlopt_destroy(opt);

    return NULL;
}

int main(int argc , char **argv)
{

    char * fp = NULL;
    int error =0 ;

    pthread_t thread, pt;

    while(argc-- >0)
    {
        if(!strncmp(argv[argc],"--sdf=",6 ))
            fp = argv[argc] + 6;
        else if(!strncmp(argv[argc],"--we=", 5))
            wE = atof(argv[argc] + 5);
        else if(!strncmp(argv[argc],"--ws=", 5))
            wS = atof(argv[argc] + 5);
        else if(!strncmp(argv[argc],"--alpha=", 8))
            alpha = atof(argv[argc] + 8);
        else if(!strncmp(argv[argc], "--np=",5))
            nprobes = atoi(argv[argc] + 5);
    }

    std::cerr << "Filename: " << fp << std::endl;
    std::cerr << "Alpha: " << alpha << std::endl;
    std::cerr << "wE: " << wE << std::endl;
    std::cerr << "wS: " << wS << std::endl;
    std::cerr << "Nprobes: " << nprobes << std::endl;


    conv.SetInFormat("SDF");
    conv.SetOutFormat("SDF");

    std::ifstream ifs(fp);
    if(!ifs)
    {
        error++;
        std::cerr << "Template not found!" << std::endl;
        return error;
    }

    conv.SetInStream(&ifs);
    conv.SetOutStream(&std::cout);

    templ = new OBMol();
    if(!conv.Read(templ))
    {
        error++;
        std::cerr << "Can't read the template!" << std::endl;
        return error;
    }

    templ->DeleteHydrogens();
    templ->Center();

    FOR_ATOMS_OF_MOL(atom, templ)
        atom->GetPartialCharge();

    ifs.close();
    conv.SetInStream(&std::cin);

    OBMol mol, pmol;
    while(true)
    {
        if(!conv.Read(&pmol))
           break;

        pmol.DeleteHydrogens();
        pmol.Center();

        FOR_ATOMS_OF_MOL(atom, pmol)
            atom->GetPartialCharge();

        if(!conv.Read(&mol))
        {
            pthread_create(&pt, NULL, pseal_align, &pmol);
            pthread_join(pt, NULL);
            conv.Write(&pmol);

            break;
        }

        mol.DeleteHydrogens();
        mol.Center();
        FOR_ATOMS_OF_MOL(atom, mol)
            atom->GetPartialCharge();

        pthread_create(&thread, NULL, pseal_align, &pmol);
        pthread_create(&pt, NULL, pseal_align, &mol);

        pthread_join(thread, NULL);
        conv.Write(&pmol);

        pthread_join(pt, NULL);
        conv.Write(&mol);

    }

    delete templ;
    return 0;
}

