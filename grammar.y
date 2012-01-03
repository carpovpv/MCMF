%{
#include <stdio.h>
#include <vector>
#include <string>

#include "parser.h"

int yyparse(void);
int yylex(void);  

void yyerror( char * s )
{
    fprintf( stderr ,"%s\n" , s );
}

ConditionParams cond;
ConditionKernel kern;


%}
%union {double dval; 
       char * cval;
       int ival;}
%token <dval> LEX_NUMBER <ival> LEX_POSINT
%token LEX_MACHINE LEX_KERNEL LEX_MODEL LEX_RESULTS LEX_CV LEX_SDF_TRAIN LEX_SDF_TEST LEX_PROGNOSIS LEX_HELP LEX_PARAMS
%token LEX_DESCR  <ival> D_DESCR  D_MNA D_FP2 D_SPECTR K_GAUSS <ival> K_KERNEL K_LINEAR K_TANIMOTO K_ELECTROST K_STERIC      
%token K_HYDROPHOBIC <cval> M_1SVM <cval>M_SVR <cval> LEX_ID LEX_OPER LEX_MAX_ITER   
%%
	
	cmd    : cmd stmt 
	       | stmt
	       ;
	stmt   : cmfa
	       | machine
	       | model
	       | maxiter
	       | cv
	       | sdftrain
	       | sdftest   
	       | parameters			
               | results
               | LEX_HELP               { cond.help = true;      }
	       | LEX_PROGNOSIS    	{ cond.prognosis = true; }
	       ;
	cmfa   : LEX_KERNEL '=' kernels 	
	       ;
	kernels: kernels ',' kernel
	       | kernel
	       ;
	kernel : K_KERNEL '(' D_DESCR ')'	
		{
			kern.kernel = $1; 
			kern.descr = $3; 
			cond.kernels.push_back(kern);
		}
               | K_KERNEL			
		{ 
			kern.kernel = $1; 
			kern.descr = 0; 
			cond.kernels.push_back(kern); 
		}
	       ;
	machine: LEX_MACHINE '=' M_1SVM	        
		{ 
			if(!cond.machine.empty())
			{
				fprintf(stderr,"Machine has already been set! Error in syntax, check it!\n");
				YYABORT;
			}
			else
				cond.machine = $3; 
		}
	       | LEX_MACHINE '=' M_SVR
		{ 
			if(!cond.machine.empty())
			{
				fprintf(stderr,"Machine has already been set! Error in syntax, check it!\n");
				YYABORT;
			}
			else
				cond.machine = $3; 
		}
		;
        model  : LEX_MODEL '=' LEX_ID           	
               	{
			if(cond.model.empty()) 
				cond.model = $3;
			else
			{
				fprintf(stderr, "Model parameter has already been set! Error in syntax, check it!\n");
				YYABORT;
			}
		}
		;
        results: LEX_RESULTS '=' LEX_ID
		{
			if(cond.results.empty()) 
				cond.results = $3;
			else
			{
				fprintf(stderr, "Results parameter has already be set! Error in syntax, check it!\n");
				YYABORT;
			}
		}
		;
	maxiter: LEX_MAX_ITER '=' LEX_POSINT
		{
			cond.max_iter = $3;
		}
		;
	cv     : LEX_CV '=' LEX_POSINT
		{
			cond.cv = $3;
		}
		;
	sdftrain : LEX_SDF_TRAIN '=' LEX_ID
		{
			if(cond.sdf_train.empty()) 
				cond.sdf_train = $3;
			else
			{
				fprintf(stderr, "Sdf-train parameter has already be set! Error in syntax, check it!\n");
				YYABORT;
			}
		}
		;
	sdftest: LEX_SDF_TEST '=' LEX_ID      
		{
			if(cond.sdf_test.empty()) 
				cond.sdf_test = $3;
			else
			{
				fprintf(stderr, "Sdf-test parameter has already be set! Error in syntax, check it!\n");
				YYABORT;
			}
		}
		;
	parameters : LEX_PARAMS '=' params;
	params : params ',' LEX_NUMBER        { cond.params.push_back($3); }
	       | params ',' LEX_POSINT	      { cond.params.push_back($3); }
	       | LEX_NUMBER 		      { cond.params.push_back($1); }
               | LEX_POSINT                   { cond.params.push_back($1); }
               ;
	

%%

#include "lex.yy.cpp"

bool parse_command_line( char *str)
{
    yy_delete_buffer(YY_CURRENT_BUFFER);
    YY_BUFFER_STATE buf = yy_scan_string("help sdf-test = asdf11.dat model=12.056.dat kernel=gauss(fp2),electrostatic cv=10 sdf-train=q h=2.3,4,0.156 results=qw machine=1-svm max-iter=2");

    yyparse();
    printf("%s\n", "OK");
}
