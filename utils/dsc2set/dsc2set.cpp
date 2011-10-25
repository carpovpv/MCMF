
#include <stdio.h>
#include <ctype.h>

int main( int argc, char **argv)
{
	FILE * fp = fopen(argv[1], "r");
        unsigned int c;
	
        int N = 1;
	while((c=fgetc(fp)) != '\n')
            if(c == ' ') N++;

        //printf("%d\n", N);
        //return 0;

	printf("AA\n");
        for(int y=1; y< N; y++)
		printf("A%d=%d\n", y,y);
	
	printf("IN 1\n");

	int i=0, u, d = 2;
        while(fscanf(fp, "%d", &u) == 1)
	{
		printf("A%d  %d\n", i+1, u); 
                if(i==N-2)
		{
                    char c;
                    while( isspace(c = fgetc(fp)) );

                    ungetc(c,fp);

                    //fprintf(stderr,"%c\n", c);
                    if(c != '$')
                        printf("IN %d\n", d++);

                    i=0;
		}
		else i++;		
	}
        char s[255];
        fgets(s, 254, fp);
        fprintf(stderr, "%s", s);
	return 0;
}

