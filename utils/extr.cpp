
#include <stdio.h>
#include <string.h>

int main(int argc, char ** argv)
{
	
	FILE * fp = fopen(argv[1], "r");
	long pos = 0, prevpos=0;
	long maxpos = 0;
	double auc = 0.0;
	double prevauc = 0.0;

	char str[BUFSIZ+1];
	while(fgets(str, BUFSIZ, fp))
	{		
		int k = strlen(str);
		if(k < 2)
		{
//			pos = ftell(fp);
		}
		else if(k > 20)
		{
			sscanf(str, "%lf", &auc);
			
			//printf("%f %ld\n", auc, pos);
			if(auc > prevauc)
			{
				prevauc = auc;
				prevpos = pos;
				
			}
			pos = ftell(fp);
		}

	}

//	printf("%f %ld\n", prevauc, prevpos);
	
	fseek(fp, prevpos, SEEK_SET);
	while(true && !feof(fp))
	{
		fgets(str, BUFSIZ, fp);
		int k = strlen(str);
		printf("%s", str);

		if(k < 2) 
		{
			fgets(str, BUFSIZ, fp);
			printf("%s", str);
			break;
		}
	}
	
	
}

