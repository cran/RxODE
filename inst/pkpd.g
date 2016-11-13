{
#include "dparse_tables.h"
}

start: (assign | deriv)+;


inputs: 'Inputs' LC NAME (',' NAME)* ',' 'file' '=' '"' FILE '"' RC;
outputs: 'Outputs' LC NAME (',' NAME)* ',' 'file' '=' '"' FILE '"' RC;
parameters_para: 'Parameters' LC assign+ RC;
initialize_para: 'Initialize' LC assign+ RC;
estimation_para: 'Estimation' LC assign+ RC;
dynamics_para: 'Dynamics' LC (assign | deriv)+ RC;
model_para: 'Model' LC (assign | distr)+ RC;
iiv_para: 'Inter_subj_vars' LC (assign | distr | paras)+ RC;
para_inits_para: 'Para_inits' LC paras+ RC;


assign: NAME '=' expr ';';	
deriv: 'd/dt' LP NAME RP '=' expr ';';
paras: vars '=' expr ';';
distr: (NAME | vars) '~' norml ';';
norml: 'N' '(' expr ',' expr ')';
expr: term (('+'|'-') term)*;
term: factor (('*'|'/') factor)*;
factor: ('+'|'-') factor | power | atom;
power: atom '^' atom;
atom: LP expr RP | function | sub_atom;
function: NAME LP exprargs RP;
exprargs: expr (',' expr)*;
sub_atom: NAME | vars | NUMBER;
vars: ('theta' | 'eta' | 'eps' | 'om' | 'sg') LP digit+ RP;


FILE ::= (letter | digit | '_' | '.')+;
NAME ::= (letter|'_') (letter | digit | '_')*;
letter ::= "[a-zA-Z]";
digit ::= "[0-9]";
NUMBER ::= integer | floatnumber;
integer ::= nonzerodigit digit* | '0';
floatnumber ::= pointfloat | exponentfloat;
pointfloat ::= intpart? fraction | intpart '.';
exponentfloat ::= (intpart | pointfloat) exponent;
intpart ::= digit+;
fraction ::= "." digit+;
exponent ::= ("e" | "E") ("+" | "-")? digit+;
nonzerodigit ::= "[1-9]";
digit ::= "[0-9]";
whitespace: ( "[ \t\r\n]+" | singleLineComment )*;
singleLineComment: '#' "[^\n]*" '\n';


LP ::= '(' ;
RP ::= ')' ;
LC ::= '{' ;
RC ::= '}' ;
QM ::= '?' ;

