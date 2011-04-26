
#include <stdio.h>

int main(int argc, char **argv)
{
	FILE *fp = stdin;
	double S=0.0;

	double x1, x2, y1, y2;

	fscanf(fp,"%lf %lf",&x1, &y1);

	while(fscanf(fp,"%lf %lf",&x2, &y2)==2)
	{
		S+=(x2-x1) * (y1 + y2);
		x1 = x2;
		y1 = y2;
	}

	S = S * 0.5;

	fclose(fp);

	printf("AUC: = %g\n", S);
	return 0;
}
