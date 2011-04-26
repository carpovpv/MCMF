
#include "../Net.h"
#include <string.h>
#include "../vector_io.h"
#include <iostream>
#include <ctype.h>
#include "../types.h"
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <stdlib.h>

using namespace NeuralNetwork;
using namespace std;
bool debug = true;

double euclid(const double *a, const double  *b, int N)
{
    double p = 0.0, q=0, r =0;
    unsigned int i;
    for (i=0; i<N; i++)
    {
        if (a[i] >=0.5 && b[i]>= 0.5) p++;
        if (a[i] <0.5 && b[i] >= 0.5) q++;
        if (a[i] >= 0.5 && b[i] < 0.5) r++;
    }
    return p/ (p+q+r);
}

double euclid(const vector<double> *a, const vector <double > *b)
{
    double p = 0.0, q=0, r =0;
    unsigned int i;
    for (i=0; i<a->size(); i++)
    {
        if (a->at(i) >=0.5 && b->at(i) >= 0.5) p++;
        if (a->at(i) <0.5 && b->at(i) >= 0.5) q++;
        if (a->at(i) >= 0.5 && b->at(i) < 0.5) r++;
    }
    return p/ (p+q+r);
}

int
main (int argc, char **argv)
{

    char *model = 0;
    char *dsc = 0;

    //if (!argc || argc != 3)
    // {
    //  if (debug)
//	cerr << "Use: ./prognosis -model=file_model -dsc=file_dsc " << endl;
//     return 1;
    //  }
    if (!strncmp (argv[1], "-model=", 7))
        model = (argv[1] + 7);
    if (!strncmp (argv[2], "-dsc=", 5))
        dsc = (argv[2] + 5);
    if (dsc == 0 || model == 0)
    {
        if (debug)
            cerr << "Error in parameters" << endl;
        return 1;
    }
    if (debug)
    {
        cerr << "Model name: " << model << endl;
        cerr << "Dsc file: " << dsc << endl;
    }
    std::ifstream in (model, std::ios::binary);
    if (!in)
    {
        if (debug)
            cerr << "Can't open file " << model << " for reading. " << endl;
        return 1;
    }
    FILE *fp;
    if (*dsc == '\0')
        fp = stdin;
    else
        fp = fopen (dsc, "r");
    if (!fp)
    {
        if (debug)
            cerr << "Can't open dsc file " << dsc << endl;
        return 1;
    }
    if (argc == 4)
    {
        //hierarchical
        FILE *hqsar = fopen(argv[3],"r");
        if (hqsar)
        {
            char hmodel_name[255];
            char hsystem[1024];
            char hmodel[1024];
            sprintf(hmodel, "%s", model);
            char *last = strrchr(hmodel,'/');
            if (last)
            {
                * (last +1 ) = '\0';
                while (fscanf(hqsar,"%s", hmodel_name) == 1)
                {

                    sprintf(last +1 ,"%s.mdl", hmodel_name);
                    sprintf(hsystem, "prognosis -model=%s -dsc=%s > hresults 2>/dev/null", hmodel, dsc );

                    //if(debug) printf("COMMAND: %s\n", hsystem);

                    system(hsystem);

                    FILE *newdescrs = fopen("hresults","r");
                    if (newdescrs)
                    {
                        char descr[1024], value[1024];
                        char buf1[255], buf2[255];
                        descr[0]=value[0]='\0';

                        while (fscanf(newdescrs, "%s : %s :", buf1, buf2)==2)
                        {
                            strcat(descr,buf1);
                            strcat(descr," ");
                            strcat(value, buf2);
                            strcat(value, " ");

                            //if(debug) printf("%s %s\n", buf1, buf2);

                        }
                        fclose(newdescrs);
                        FILE *nndescrs = fopen("ndescrs","w");
                        if (nndescrs)
                        {
                            fprintf(nndescrs,"%s\n%s\n$$$$", descr, value);
                            fclose(nndescrs);
                            system ("paste -d ' ' all.all ndescrs > all1.all 2>/dev/null");
                            system("mv -f all1.all all.all");
                        }

                    }

                }
            }
            fclose(hqsar);
        }
    }

    if (debug)
        cerr << "Reading binary data from model..." << endl;

    vector < string > descriptors, prop, descriptor_blocks;;
    vector < double >dsc_max, dsc_min, prop_max, prop_min;

    vector_read (in, &descriptor_blocks);
    vector_read (in, &descriptors);

    vector_read (in, &dsc_max);
    vector_read (in, &dsc_min);
    vector_read (in, &prop);
    vector_read (in, &prop_max);
    vector_read (in, &prop_min);

    int NN = 1;
    int cls;
    int r =0;

    if (debug)
        cerr << "Descriptor blocks  used in building the model: " <<
             descriptor_blocks.size () << endl;

    if (debug)
        for (int i = 0; i < descriptor_blocks.size (); ++i)
            cerr << descriptor_blocks[i] << endl;


    if (debug)
    {
        cerr << "Descriptors used in building the model: " << descriptors.
             size () << endl;
        cerr << "Reading descriptors from dsc..." << endl;
    }
    int space = 0;
    int q;
    vector < string > dsc_name;


    q = fgetc(fp);
    while (true)
    {
        char buh[100000];
        int h = 0;
        while (isspace(q) && q!='\n')  q = fgetc (fp);

        //if(q=='\n') break;

        space++;
        while (!isspace (q) && q != '\n')
        {
            buh[h++] = q;
            q = fgetc (fp);
        }
        buh[h] = '\0';
        dsc_name.push_back (buh);
        if (q == '\n') break;
    }

    dsc_name.pop_back ();

    if (debug)
        cerr << "Number of descriptors in dsc " << dsc_name.size () << endl;

    if (debug)
    {
        for (int i = 0; i < dsc_name.size (); ++i)
            cerr << "Descr in dsc: " << dsc_name[i] << endl;

        cerr << "Creating a map..." << endl;
    }
    vector < int >dmap (dsc_name.size (), -1);
    for (int i = 0; i < dsc_name.size (); ++i)
    {
        vector < string >::iterator it =
            find (descriptors.begin (), descriptors.end (), dsc_name[i]);
        if (it != descriptors.end ())
            dmap[i] = distance (descriptors.begin (), it);
        if (debug)
            cerr << i << " " << dmap[i] << endl;
    }

    int nn = count (dmap.begin (), dmap.end (), -1);
    if (nn == dmap.size ())
    {
        if (debug)
            cerr <<
                 "Descriptors are incomplete! Model can't be used within this dsc file."
                 << endl;
        return 1;
    }

    NetType nettype = OCCDIABOLO;
//    in.read (reinterpret_cast < char *>(&nettype), sizeof (NetType));

//    nettype = BACKPROP;

    switch (nettype)
    {
    case BACKPROP:
    case OCCDIABOLO:
    {
        if (debug)
            cerr << "Read model..." << endl;

        Net *net = new Net (in);
        if (!net)
        {
            if (debug)
                cerr << "Error in read model." << endl;
            return 1;
        }

        double *x = new double[descriptors.size ()];
        double *y = new double[prop.size ()];

        int NN = 1;

        while (true)
        {


            unsigned char c = fgetc (fp);
            if (c == '$')
                return 0;

            ungetc (c, fp);

//            for (int i = 0; i < descriptors.size (); ++i)
            //              x[i] = 0.9 - 0.8 * dsc_max[i] / (dsc_max[i] - dsc_min[i]);



            for (int i = 0; i < dsc_name.size (); ++i)
            {
                double _x;
                if(fscanf (fp, "%lf", &_x) !=1) exit(0);
                if (dmap[i] != -1)
                {
                    x[dmap[i]] = _x;
                    /*                x[dmap[i]] =
                                        0.8 * x[dmap[i]] / (dsc_max[dmap[i]] -
                                                            dsc_min[dmap[i]]) + 0.9 -
                                        0.8 * dsc_max[dmap[i]] / (dsc_max[dmap[i]] -
                                                                  dsc_min[dmap[i]]);
                      */
                }
            }


//    for(int i=0;i<descriptors.size();++i)
//    cout << x[i] << " ";
//    cout << endl;

            net->run (x, y);


            if (debug)
                cerr << NN++ << " ";

            /*            if (nettype == BACKPROP)
                        {
                            for (int i = 0; i < prop.size (); ++i)
                            {
                                cout << prop[i] << " : " ;
                                y[i] =
                                    (y[i] - 0.9 +
                                     0.8 * prop_max[i] / (prop_max[i] -
                                                          prop_min[i])) * (prop_max[i] -
                                                                           prop_min[i]) / 0.8;

                                cout << y[i] << " : ";
                            }

                            cout << endl;
            		continue;
                        }
                        else
            */
            if (nettype == OCCDIABOLO)
            {

                double dist = euclid(x,y,descriptors.size());
                printf("%s %g\n", argv[3],dist);

//	    	    std::vector < double> constraints;
//		    vector_read(in, &constraints);
//		    printf("feature : %s :", (dist >= constraints[0] && dist <= constraints[1]) ? "yes" : "no");
            }

//            break;		// only for web version
        }
//
        delete net;
        break;
    }
    case OCCK:
    {
        int K;
        in.read (reinterpret_cast < char *>(&K), sizeof (int));
        //      printf("One-class: %d\n",K);
        std::vector< std::vector< double>  > kvectors, kconstraints;
        std::vector< double > temp;

        for (int i=0; i< K; ++i)
        {
            vector_read(in, &temp);
            kvectors.push_back(temp);
            temp.clear();
        }

        //      printf("One-class: %d vectors read\n",kvectors.size());

        for (int i=0; i< K; ++i)
        {
            vector_read(in, &temp);
            kconstraints.push_back(temp);
//		  printf("Size: %d\n", temp.size());
            temp.clear();
        }
        //    printf("One-class: %d constraints vectors read\n",kconstraints[0].size());

        temp.clear();
        temp.resize(kvectors[0].size(),0.0);

        while (true)
        {

            unsigned char c = fgetc (fp);
            if (c == '$')
                break;
            ungetc (c, fp);


            for (int i = 0; i < dsc_name.size (); ++i)
            {
                double _x;
                if (fscanf (fp, "%lf", &_x)!=1) goto oc2;
                if (dmap[i]!=-1) temp[dmap[i]] = _x;
            }


            int mc = 0;
            double md = DBL_MAX;

            for (int i=0; i< K; ++i)
            {
                double dist = euclid(&kvectors[i], &temp);
//		      printf("Dist to %d = %g\n", i+1, dist);
                if ( md > dist) {
                    md = dist;
                    mc = i;
                }
            }
//	      printf("Winner class: %d %g %g\n", mc,kconstraints[mc][0], kconstraints[mc][1]);

            if (md >= kconstraints[mc][0] && md <= kconstraints[mc][1])
                printf("feature : yes :\n");
            else
                printf("feature : no :\n");


            break;
        }
oc2:
        ;
        break;
    }

    }
f2:

    fclose (fp);
    return 0;

}
