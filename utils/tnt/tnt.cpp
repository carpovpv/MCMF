
// Calculate true and false predictions

#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

double threshold;

void calc(vector < vector < double > > *training)
{
	double tn = 0.0;
	double tp = 0.0;
	double _fp = 0.0;
	double fn = 0.0;

	for(int i=0; i< training->size(); ++i)
	{
		if(training->at(i)[0]>1.0e-3)
		{
			//real positive
			if(training->at(i)[1] >= threshold) tp++;
			else fn++;
		}
		else
		{
			//real negative
			if(training->at(i)[1] < threshold) tn++;
			else 
			{
				_fp++;
			}
		}
	}
	double s = training->size();

//	printf("Size: %5g\n", s);
//	printf("TN: = %g (%g) TP: = %g (%g) FP: = %g (%g) FN: = %g (%g)\n", tn, tn/s, tp,tp /s,  _fp, _fp/s,  fn, fn/s);

//	printf(" _FP: %g\n TN: %g\n TP: %g\n FN: %g\n", _fp, tn, tp, fn);
	printf(" %.0f %.0f %.0f %.0f %.1f %.1f\n", tn, tp, fn, _fp, 100 * tn/(tn+_fp),100 * tp/(tp+fn));


}


int main(int argc, char **argv)
{
	if(argc != 2)
	{
		fprintf(stderr,"Using ./tn  threshold\n");
		return 1;
	}
	threshold = atof(argv[1]);

	vector< vector < double> > training;
	vector < double > temp(2,0);
	int cur = 0;

	double x,y;
	while(scanf("%lf %lf", &x, &y)==2)
	{

		temp[0]=x;
		temp[1]=y;
		training.push_back(temp);
	}
	calc(&training);
	return 0;
}
