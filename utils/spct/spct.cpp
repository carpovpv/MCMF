
// Calculate true and false predictions

#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <limits.h>

using namespace std;

void calc(vector < vector < double > > *training, double threshold)
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
	printf(" %.6f %.6f %.6f\n",_fp/(tn+_fp),tp/(tp+fn), threshold);


}


int main(int argc, char **argv)
{

	char str[255];

	vector< vector < double> > training;
	vector < double > temp(2,0);
	int cur = 0;
	double x,y;

	double maxh = -100000;
	double minh = 100000;

	while(scanf("%lf %lf", &x, &y)== 2)
	{

 		temp[0]=x;
		temp[1]=y;

		if(y>maxh) maxh = y;
		if(y<minh) minh = y;

		training.push_back(temp);
	}
        const double step = (maxh - minh) / 10000;

	for(double h=minh-0.1*step; h<=maxh+0.1*step; h+=step)
		calc(&training, h);


	return 0;
}
