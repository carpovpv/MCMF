
#include "Net.h"

#include <time.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <signal.h>
#include "vector_io.h"

using namespace NeuralNetwork;

std::vector < std::vector < double > >
  prop;

std::vector < int >
  training;			// only "pointers" to training data
std::vector < int >
  validation;			//    ----------//----------
std::vector < int >
  prediction;			//    ----------//----------

std::vector < std::vector < double > >
  data;

// edges for descriptors
std::vector < double >
  normalization_max;
std::vector < double >
  normalization_min;
//edges for properties
std::vector < double >
  norm_max;
std::vector < double >
  norm_min;
double
quad (double val, double i)
{
  return val + i * i;
}

Net * net = NULL;

void sig_handler(int sig)
{
 if(sig == SIGINT)
        {
                if(net==NULL) exit(0);
                net->stop = true;
                printf("Terminated by user!\n");
        }
}


int
count (int i)
{
  int
    s = 0;
  for (int j = 0; j < prop.size (); j++)
    {
      if (prop.at (j).at (i) != UNKNOWN_VALUE)
        ++s;
    }

  return s;
}

void
find_min_max (double M)
{
  for (int i = 0; i < data.at (0).size (); ++i)
    {
      normalization_max.push_back (data.at (0).at (i));
      normalization_min.push_back (data.at (0).at (i));
    }
  for (int i = 0; i < data.size (); ++i)
    for (int j = 0; j < data.at (0).size (); ++j)
      if (data.at (i).at (j) > normalization_max.at (j))
        normalization_max.at (j) = data.at (i).at (j);
      else if (data.at (i).at (j) < normalization_min.at (j))
        normalization_min.at (j) = data.at (i).at (j);

  //enlatge interval
  for (int j = 0; j < data.at (0).size (); ++j){

        double x0 = normalization_max.at (j) ;
        double y0 = normalization_min.at (j);

        normalization_max.at (j) = (x0 + y0 + M*(x0-y0))/2;
        normalization_min.at (j) = -normalization_max[j] +x0 + y0;

  }

  for (int i = 0; i < prop.at(0).size (); ++i)
    {
      norm_max.push_back ( prop.at(0).at (i));
      norm_min.push_back ( prop.at(0).at (i));
    }
  for (int i = 0; i < prop.size (); ++i)
    for (int j = 0; j < prop.at(i).size(); ++j)
      if (prop.at (i).at (j) != UNKNOWN_VALUE && (prop.at (i).at (j) > norm_max.at (j)))
        norm_max.at (j) = prop.at (i).at (j);
      else if (prop.at (i).at (j) != UNKNOWN_VALUE
               && (prop.at (i).at (j) < norm_min.at (j)))
        norm_min.at (j) = prop.at (i).at (j);

 //enlarge intrerval
 for (int j = 0; j < prop.at(0).size(); ++j)
 {

        double x0 = norm_max[j];
        double y0 = norm_min[j];

        norm_max.at (j) = (x0 + y0 + M*(x0-y0))/2.0;
        norm_min.at (j) = -norm_max[j] + x0 +y0;

//	std::cout << norm_max[j] << " " << norm_min[j] << std::endl;

 }

  std::cout << "Normalization ... ";

  for (int i = 0; i < data.size (); ++i)
    for (int j = 0; j < data.at (0).size (); ++j)
      data.at (i).at (j) =
        0.8 * data.at (i).at (j) / (normalization_max.at (j) -
                                    normalization_min.at (j)) + 0.9 -
        0.8 * normalization_max.at (j) / (normalization_max.at (j) -
                                          normalization_min.at (j));

  for (int i = 0; i < prop.size (); ++i)
    {
      for (int j = 0; j < prop.at (0).size (); ++j)
        if (prop.at (i).at (j) != UNKNOWN_VALUE)
          {
            prop.at (i).at (j) =
              0.8 * prop.at (i).at (j) / (norm_max.at (j) - norm_min.at (j)) +
              0.9 - 0.8 * norm_max.at (j) / (norm_max.at (j) -
                                             norm_min.at (j));
          }
    }
  std::cout << "... fininshed." << std::endl;
}

class ArrayRangeExampleFactory:
public ExampleFactory
{
public:
  ArrayRangeExampleFactory (int source, int initLower, int initUpper):
  currentSource (source),
  currentExample (initLower),
  lower (initLower),
  upper (initUpper)
  {
  }
  void
  getExample (int inputSize, real * input, int outputSize, real * output)
  {
    if (currentSource == 1)
      {

        for (int i = 0; i < inputSize; ++i)
        input[i] =  data.at (training.at (currentExample)).at (i);
        for (int i = 0; i < outputSize; ++i)
        output[i] = prop.at (training.at (currentExample)).at (i);
      }
    else if(currentSource == 2)
      {

        for (int i = 0; i < inputSize; ++i)
          input[i] = data.at (prediction.at (currentExample)).at (i);
        for(int i=0;i<outputSize;++i)
        output[i] = prop.at (prediction.at (currentExample)).at (i);

      }
    else
      {
        for (int i = 0; i < inputSize; ++i)
          input[i] = data.at (validation.at (currentExample)).at (i);
        for(int i=0;i<outputSize;++i)
        output[i] = prop.at (validation.at (currentExample)).at (i);

      }

    currentExample++;
    if (currentExample >= upper)
      currentExample = lower;
  }
  int
  numExamples ()
  {
    return upper - lower + 1;
  }
private:
  int
    currentSource;
  int
    currentExample;
  int
    lower,
    upper;
};


int
main (int argc, char *argv[])
{

  signal(SIGINT,sig_handler);
  FILE *
    set = NULL, *dsc = NULL;
  char
    str[255];
  char
    buf[20];

  int N = 5,NN=1000;
  double M = 4.0 / 3.0;
  int epochs = 100;

  bool logarithm = false;
  bool novalid = false;

  std::vector < std::string > names;
  std::vector < std::string> descriptor_blocks;
  std::vector < int >
    properties;
  std::vector < int > only_valid;
  std::vector < int >
    sliding;
  std::vector < int > _layers;
  int NLayers;

  for (int i = 1; i < argc; ++i)
    {
      if (!strncmp(argv[i],"-N=",3))
      {
          char *p = argv[i] + 3;
          while (true)
            {
              int
                k;
              sscanf (p, "%d", &k);
              _layers.push_back (k);
              while (isdigit (*p))
                p++;
              if (*p == ',')
                p++;
              else if (*p == '\0')
                break;
            }
            NLayers= _layers.size();
      }

      if (!strncmp(argv[i],"-epochs=",8)) epochs = atoi(argv[i]+8);

      if (!strncmp(argv[i],"-cycles=",8)) NN = atoi(argv[i]+8);
      if (!strncmp(argv[i],"-M=",3)) M = 1.0 / atof(argv[i]+3);
      if (!strncmp(argv[i],"-log",4)) logarithm = true;

      if (!strncmp (argv[i], "-set=", 5))
        set = fopen (argv[i] + 5, "r");
      if (!strncmp(argv[i], "-novalid", 8))
        novalid = true;

      if( !strncmp(argv[i],"--valid=",8))
      {
          char *p = argv[i] + 8;
          if(!isdigit(*p))
          {
                  //validation set is stored in a file
                  FILE *fpValid = fopen(p,"r");
                  if(fpValid == NULL)
                  {
                          fprintf(stderr,"File with the validation set cannot be opened!\n");
                          return 0;
                  }
                  int val;
                  while(fscanf(fpValid,"%d",&val)==1)
                          only_valid.push_back(val);
                  fclose(fpValid);
          }
          else
          while (true)
            {
              int
                k;
              sscanf (p, "%d", &k);
              only_valid.push_back (k);
              while (isdigit (*p))
                p++;
              if (*p == ',')
                p++;
              else if (*p == '\0')
                break;
            }
      }

      if (!strncmp (argv[i], "-prop=", 6))
        {
          char *
            p = argv[i] + 6;
          while (true)
            {
              int
                k;
              sscanf (p, "%d", &k);
              properties.push_back (k);
              while (isdigit (*p))
                p++;
              if (*p == ',')
                p++;
              else if (*p == '\0')
                break;
            }
        }
      if (!strncmp (argv[i], "-sliding=", 9))
        {
          char *
            p = argv[i] + 9;
          while (true)
            {
              int
                k;
              sscanf (p, "%d", &k);
              sliding.push_back (k);
              while (isdigit (*p))
                p++;
              if (*p == ',')
                p++;
              else if (*p == '\0')
                break;
            }
          if (sliding.size () != 3 || (sliding.at (0) == sliding.at (1)))
            {
              std::
                cout <<
                "Sliding syntax: -sliding=first_validation,first_prediction,step"
                << std::endl;
              std::
                cout <<
                "Note: first validation and prediction points are not to be equal."
                << std::endl;
              if (set)
                fclose (set);
              if (dsc)
                fclose (dsc);
              return 1;
            }
        }
      if (!strncmp (argv[i], "-dsc=", 5))
        dsc = fopen (argv[i] + 5, "r");
    }
  if (!set || !dsc)
    {
      printf ("Unable to open set/dsc file\n");
      return 1;
    }

  printf ("Set opened. Parse set...\n");
  do
      fgets (str, 254, set);
  while (strncmp (str, "AA",2));

  fgets (str, 254, set);

  while (strncmp (str, "IN", 2))
    {
      char *
        p = strchr (str, '=');
      strcpy (buf, p + 1);
      buf[strlen (buf) - 1] = '\0';
      names.push_back (buf);
//      printf("%s\n",buf );
      fgets (str, 254, set);

    }

  int
    num_properties = properties.size ()? properties.size () : names.size ();
  printf ("Number of selected properties %d\n", num_properties);

  if (properties.size ())
    if (*std::max_element (properties.begin (), properties.end ()) >
        names.size ())
      {
        std::cout << "Check selected properties...\n" << std::endl;
        fclose (set);
        fclose (dsc);
        return 1;
      }

  std::vector < double >
  buffer (num_properties);

  int
    record = 0;
  fgets (str, 254, set);
  while (true)
    {
      for (int i = 0; i < buffer.size (); ++i)
        buffer.at (i) = 666;
      while (strncmp (str, "A", 1))
        fgets (str, 254, set);
      while (strncmp (str, "IN", 2))
        {
//	  printf("S: %s\n", str);
          char *
            p = strchr (str, 'A');
          p++;
          int
            k;
          double
            pr;
          sscanf (p, "%d", &k);
          p = strchr (str, ' ');
          p++;
          sscanf (p, "%lf", &pr);
          if (properties.size ())
            {
              std::vector < int >::iterator
                it = std::find (properties.begin (), properties.end (), k);
              if (it != properties.end ())
                {
                  int
                    i;
                  for (i = 0; *it != properties.at (i); i++);
                  buffer.at (i) = pr;
                }
            }
          else
            {
              buffer.at (k - 1) = pr;
            }
          if (!fgets (str, 254, set))
            {
              prop.push_back (buffer);
              goto fine;
            }
        }
      prop.push_back (buffer);

  //    printf("Hello %d %s\n", buffer.size(), str);

      record++;
      fgets (str, 254, set);
      if (feof (set))
        break;
    }
fine:
  fclose (set);


printf("Properties: \n");

// if(logarithm)
//  for(int i=0;i<prop.size();++i)
  //{
//	  for(int j=0;j<prop.at(0).size();++j)
//		  prop[i][j] = log(prop[i][j]+1.0);
//		  printf(" %g",prop.at(i).at(j));
//	  printf("\n");
  //}


  std::cout << "Number of points (all):  " << prop.size () << std::endl;
  std::cout << "Data read successfully. Set file closed." << std::endl;
  std::cout << "List of properties: " << std::endl;

  std::vector <string > prop_model;
  if (properties.size ())
    for (int i = 0; i < properties.size (); ++i)
    {
      prop_model.push_back(names.at(properties.at(i)-1));
      std::cout << "\t" << names.at (properties.at (i) -
                                     1) << "(" << count (i) << ")" << std::
        endl;
    }
  else
    for (int i = 0; i < names.size (); ++i){
      std::cout << "\t" << names.
        at (i) << "(" << count (i) << ")" << std::endl;
      prop_model.push_back(names.at(i));
    }

  std::cout << "Dsc opend. Parse dsc..." << std::endl;

  int
    space = 0;
  int
    q;

  std::vector < std::string > dsc_name;

  while ((q = fgetc (dsc)) != '\n')
    if (q == ' ' || q =='\t')
      space++;

std::cout << "Number of descriptors: " << space << std::endl;

rewind(dsc);

for(int a=0;a<space;++a)
{

char buh[1000];
fscanf(dsc,"%s",buh);
dsc_name.push_back(buh);

}

  while (true)
    {
      std::vector < double >
      buffer (space);
      for (int i = 0; i < space; ++i)
        {
          double
            k;
          if (fscanf (dsc, "%lf", &k) != 1)
            goto fine_1;
          buffer.at (i) = k;
//	  std::cout << "Descr: " << k << std::endl;
        }
      data.push_back (buffer);
    }
fine_1:
  std::cout << "Loaded " << data.size () << " descriptors." << std::endl;
  std::cout << "Number of descriptors equal to number of data. Ok." << std::
    endl;

// define descriptor blocks

while(fscanf(dsc,"%s",str)==1)
{
        char *p = str;

        bool t = true;
        bool dig = true;

        while(*p!='\0')
        {
                if(*p!='$') t = false;
                if(!isdigit(*p)) dig = false;
                p++;
        }

        if ( !t && !dig)
                descriptor_blocks.push_back(str);

}


  fclose (dsc);
  std::cout << "Dsc closed." << novalid << std::endl;

  if(!novalid)
{
  if (num_properties == 1 && sliding.size () == 3 || only_valid.size() > 0)
          std::cout << "Ordinary sliding performed." << std::endl;
  else std::cout << "Sliding only on data not properties. " << std::endl;

if(only_valid.size() == 2)
{
        //I mean sliding
        int first = only_valid[0];
        int step = only_valid[1];

        only_valid.clear();
        for(int i=0;i<data.size();++i)
                if(((i-first)%step)==0) only_valid.push_back(i);
}


  if(sliding.size())
{

      std::cout << "Perform simple sliding on data... " << std::endl;
      std::cout << "First validation point: " << sliding.at (0) << std::endl;
      std::cout << "First prediction point: " << sliding.at (1) << std::endl;
      std::cout << "Sliding step: " << sliding.at (2) << std::endl;

      for (int i = 0; i < data.size (); ++i)
        {
          if ((i - sliding.at (0)) % sliding.at (2) == 0)
            {
              validation.push_back (i);
              continue;
            }
          if ((i - sliding.at (1)) % sliding.at (2) == 0)
            prediction.push_back (i);
          else
            training.push_back (i);
        }
}
else
{
  for (int i = 0; i < data.size (); ++i)
        {
                if(find(only_valid.begin(), only_valid.end(),i+1)==only_valid.end())

                        training.push_back(i);
                else
                        validation.push_back(i);
        }

}
}
else
{
        // novalidation one-class classification problems
         for (int i = 0; i < data.size (); ++i)
                training.push_back(i);
}

      std::cout << "Training data size:   " << training.size () << std::endl;
      std::cout << "Validation data size: " << validation.
        size () << std::endl;
      std::cout << "Prediction data size: " << prediction.
        size () << std::endl;

      printf("Number of input neurons: %d\n",N);

 // find_min_max (M);

  time_t _time;
  time(&_time);
  printf("Start at: %s\n",ctime(&_time));

  srand ((int)time(NULL));

  int *  layers = (int *) calloc( NLayers + 2, sizeof(int));
  layers[0] = data.at(0).size();
  layers[NLayers+1] = prop.at(0).size();
  for(int y=1; y< NLayers+1; y++)
        layers[y] = _layers[y-1];

  for(int y=0; y<NLayers +2; ++y)
        printf("Layer: %d %d\n", y, layers[y]);

   net = new Net (NLayers + 2 , layers, 0.05, 0.5, 1.0,NN);
  net->randomizeWeights ();

  ArrayRangeExampleFactory train (1, 0, training.size ());

  ArrayRangeExampleFactory valid (0, 0, validation.size ());
  ArrayRangeExampleFactory pred  (2, 0, prediction.size ());

  real error ;
  if(novalid)
         error = net->autotrain (train, train, epochs, 1.05f);
  else
        error = net->autotrain (train, valid, epochs, 1.05f);

  std::cout << "Final validation error: " << error << std::endl;
 /*
  double * yy = new double[1024];
  double * rr = new double[1024];

  fprintf(stderr,"Training:\n");

  for(int i=0;i<training.size();++i)
  {
        for (int j = 0; j < data.at (0).size (); ++j)
        yy[j] = data.at (training.at (i)).at (j);
        net->run (yy, rr);

        fprintf(stderr,"\n");

  }

  fprintf(stderr,"Validation:\n");

  for(int i=0;i<validation.size();++i)
  {
        for (int j = 0; j < data.at (0).size (); ++j)
        yy[j] = data.at (validation.at (i)).at (j);
        net->run (yy, rr);

        fprintf(stderr,"\n");

  }
   fprintf(stderr,"Prediction:\n");


  for(int i=0;i<prediction.size();++i)
  {
        for (int j = 0; j < data.at (0).size (); ++j)
        yy[j] = data.at (prediction.at (i)).at (j);

        net->run (yy, rr);
        double tpp = 0, tqq =0 , trr = 0;

        for(int q=0; q< 1024; ++q)
        {
                if(yy[q] == 1 && rr[q] >= 0.5) tpp++;
                if(yy[q] == 0 && rr[q] >= 0.5) trr++;
                if(yy[q] == 1 && rr[q] < 0.5) tqq++;
        }

        //fprintf(stderr,"1 %g\n",  tpp /( tpp + tqq + trr) );
        fprintf(stderr,"%g %g\n", yy[0], rr[0]);

  }
  */
//  delete yy;
//  delete rr;

time(&_time);
printf("Finished at: %s\n",ctime(&_time));

char *mm;
bool mod_yes=false;
char model[3000];

for(int i=0;i<argc;++i)
{
if(!strncmp(argv[i],"-model=",7)) mm = argv[i] + 7;
mod_yes = true;
}

    model[0]=0;
    for(int i=0;i<argc;++i){
         strcat(model,argv[i]);
         strcat(model," ");
}
if(!mod_yes) mm = model;
printf("Command: %s\n %s\n",model, mm);

/*
if(!novalid)
{
printf("Perform sensitivity analysis:\n");

net->_train_sens(1,train);
net->_train_sens(1,valid);
if(prediction.size())
net->_train_sens(1,pred);

int yp = training.size() + validation.size();
if(prediction.size()) yp +=  prediction.size();

for(int i=0;i<net->sensitivity.size();++i)
{
std::cout << dsc_name.at(i);
printf("\t%12.5f\n",(net->sensitivity.at(i) / yp-0.9 + 0.8*normalization_max.at(i)/ (normalization_max.at(i) - normalization_min.at(i))) * (normalization_max.at(i)-normalization_min.at(i))/0.8);
}
}
*/

    net->doneTraining ();
    std::ofstream out(mm, std::ios::binary);

    vector_write(out,&descriptor_blocks);
    vector_write(out,&dsc_name);
    vector_write(out,&normalization_max);
    vector_write(out,&normalization_min);
    vector_write(out,&prop_model);
    vector_write(out,&norm_max);
    vector_write(out,&norm_min);


  net->save(out);
//free(layers);

return 0;

}
