%{
  //Команды для лекса. Благодаря нему мы поддерживаем 
  //весьма произвольный порядок следования аргументов.
extern "C"  int yywrap()
{
return 1;
}
%}
delim [ \t\n]
ws {delim}+
digit [0-9]
letter [a-zA-Z\.\-_]
number {digit}+(\.{digit}+)?(E[+-]?{digit}+)?
id {digit}*{letter}({letter}|{digit})*
posint {digit}+
%%
ws 		{}
1-svm  		{yylval.cval = yytext; return M_1SVM;}
svr             {yylval.cval = yytext; return M_SVR;}
max-iter     	{return LEX_MAX_ITER;}
machine 	{return LEX_MACHINE;}
kernel 		{return LEX_KERNEL;}
model		{return LEX_MODEL;}
results 	{return LEX_RESULTS;}
cv 		{return LEX_CV;}
sdf-test 	{return LEX_SDF_TEST;}
sdf-train 	{return LEX_SDF_TRAIN;}
prognosis 	{return LEX_PROGNOSIS;}
"=" 		{return '=';}
"," 		{return ',';}
"(" 		{return '(';}
")" 		{return ')';}
electrostatic 	{yylval.ival = K_ELECTROST; return K_KERNEL;}
steric 		{yylval.ival = K_STERIC; return K_KERNEL;}
hydrophobic 	{yylval.ival = K_HYDROPHOBIC; return K_KERNEL;}
gauss 		{yylval.ival = K_GAUSS; return K_KERNEL;}
linear 		{yylval.ival = K_LINEAR; return K_KERNEL;}
tanimoto 	{yylval.ival = K_TANIMOTO; return K_KERNEL;}
fp2 		{yylval.ival = D_FP2; return D_DESCR;}
spectrophores 	{yylval.ival = D_SPECTR; return D_DESCR;}
mna		{yylval.ival = D_MNA; return D_DESCR;}
h 		{return LEX_PARAMS;}
help 		{return LEX_HELP; }
{posint}        {yylval.ival = atoi(yytext); return LEX_POSINT;}
{number} 	{yylval.dval = atof(yytext); return LEX_NUMBER;}
{id} 		{yylval.cval = yytext; return LEX_ID;}
.		{}
%%


