
#include <stdio.h>

int main( int argc, char **argv)
{
	FILE * fp = fopen(argv[1], "r");
        unsigned int c;
	
	while((c=fgetc(fp)) != '\n')
		;

	printf("AA\n");
	for(int y=1; y<= 1024; y++)
		printf("A%d=%d\n", y,y);
	
	printf("IN 1\n");

	int i=0, u, d = 2;
	fscanf(fp,"%d", &u);
	do
	{
		printf("A%d  %d\n", i+1, u); 
		if(i==1023)
		{
			printf("IN %d\n", d++);
			i=0;			
		}
		else i++;		
	}
	while(fscanf(fp, "%d", &u) == 1);
 
	return 0;
}

