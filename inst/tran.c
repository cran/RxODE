#include <stdio.h>
#include <string.h>
#include "dparse_tree.h"
#define max(a,b) (a)>(b) ? (a):(b)
#define MXSYM 5000
#define MXDER 500
#define MXLEN 1200
#define MXBUF 2400
#define SBPTR sb.s+sb.o


char *sbuf_read(char *pathname);  /* defined in util.h */
extern D_ParserTables parser_tables_gram;


typedef struct symtab {
  char *ss;			/* symbol string: all vars*/
  int vo[MXSYM];	/* offset of symbols */
  int lh[MXSYM];	/* lhs symbols? */
  int di[MXDER];	/* ith of state vars */
  int nv;			/* nbr of symbols */
  int ix;			/* ith of curr symbol */
  int fn;			/* curr symbol a fn?*/
  int nd;			/* nbr of dydt */
} symtab;
symtab tb;

typedef struct sbuf {
  char s[MXBUF];	/* curr print buffer */
  int o;			/* offset of print buffer */
} sbuf;
sbuf sb;			/* buffer w/ current parsed & translated line */
        			/* to be stored in a temp file */

static FILE *fpIO, *fp_inits;


/* new symbol? if no, find it's ith */
int new_or_ith(const char *s) {
  int i, len, len_s=strlen(s);

  if (!tb.nv) return 1;
  if (tb.fn) return 0;
  if (!strcmp("t", s)) return 0;
  if (!strcmp("podo", s)) return 0;
  if (!strcmp("tlast", s)) return 0;

  for (i=0; i<tb.nv; i++) {
    len = tb.vo[i+1] - tb.vo[i] - 1;  /* -1 for added ',' */
    if (!strncmp(tb.ss+tb.vo[i], s, max(len, len_s))) {	/* note we need take the max in order not to match a sub-string */
      tb.ix = i;
      return 0;
    }
  }
  return 1;
}

void wprint_node(int depth, char *name, char *value, void *client_data) {
  sprintf(SBPTR, " %s", value);
  sb.o += strlen(value)+1;
}

void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
#ifndef __RODE__
  if (
      !strcmp("estimation_para", name) ||
      !strcmp("initialize_para", name) ||
      !strcmp("parameters_para", name) ||
      !strcmp("model_para", name) ||
      !strcmp("iiv_para", name) ||
      !strcmp("outputs", name) ||
      !strcmp("inputs", name)
     ) return;
#endif
  int nch = d_get_number_of_children(pn), i;
  char *value = (char*)dup_str(pn->start_loc.s, pn->end);
  char pexpr[80];


  if (!strcmp("NAME", name) && new_or_ith(value)) {
    static int pos=0;
    sprintf(tb.ss+pos, "%s,", value);
    pos += strlen(value)+1;
    tb.vo[++tb.nv] = pos;
  }

  if (!strcmp("LP", name)) {sprintf(SBPTR, "("); sb.o++;}
  if (!strcmp("RP", name)) {sprintf(SBPTR, ")"); sb.o++;}
  if (!strcmp(",",  name)) {sprintf(SBPTR, ","); sb.o++;}

  if (
      !strcmp("sub_atom", name) ||
      !strcmp("comp_op", name) ||
      !strcmp("or", name) ||
      !strcmp("and", name) ||
      !strcmp("not", name) ||
      !strcmp("+", name) ||
      !strcmp("-", name) ||
      !strcmp("*", name) ||
      !strcmp("/", name) ||
      !strcmp("=", name)
     )
    fn(depth, name, value, client_data);
  free(value);

  depth++;
  if (nch != 0) {

    if (!strcmp("power", name)) {
      sprintf(SBPTR, " pow(");
      sb.o+=5;
    }

    for (i = 0; i < nch; i++) {
      if (!strcmp("deriv", name) && i< 2) continue;
      if (!strcmp("deriv", name) && i==3) continue;
      if (!strcmp("deriv", name) && i==4) continue;

      tb.fn = (!strcmp("function", name) && i==0) ? 1 : 0;
      D_ParseNode *xpn = d_get_child(pn,i);
      wprint_parsetree(pt, xpn, depth, fn, client_data);

      if (!strcmp("power", name) && i==0) {
        sprintf(SBPTR, ",");
        sb.o++;
      }

      if (!strcmp("function", name) && i==0) {
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
        sprintf(SBPTR, " %s", v);
        sb.o += strlen(v)+1;
        free(v);
      }

      if (!strcmp("deriv", name) && i==2) {
        sprintf(sb.s, "DADT[%d] = InfusionRate[%d] +", tb.nd, tb.nd);
        sb.o = strlen(sb.s);

        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
        new_or_ith(v);
        tb.lh[tb.ix] = 1;
        tb.di[tb.nd] = tb.ix;
        tb.nd++;
        free(v);
        continue;
      }

      if (!strcmp("assign", name) && i==0) {
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
        sprintf(sb.s, "%s", v);
        sb.o = strlen(v);

        new_or_ith(v);
        tb.lh[tb.ix] = 1;
        free(v);
      }
    }

    if (!strcmp("assign", name) || !strcmp("deriv", name))
      fprintf(fpIO, "%s;\n", sb.s);
    if (!strcmp("paras", name))
      fprintf(fp_inits, "%s;\n", sb.s);

    if (!strcmp("power", name)) {
      sprintf(SBPTR, ")");
      sb.o++;
    }
  }

}

void retieve_var(int i, char *buf) {
  int len;

  len = tb.vo[i+1] - tb.vo[i] - 1;
  strncpy(buf, tb.ss+tb.vo[i], len);
  buf[len] = 0;
}

void err_msg(int chk, const char *msg, int code)
{
  if(!chk) {
    if (tb.ss) free(tb.ss);
    fprintf(stderr, msg);
    exit(code);
  }
}

/* when prnt_vars() is called, user defines the behavior in "case" */
void prnt_vars(int lhs, const char *pre_str, const char *post_str) {
  int i, j;
  static int scenario=0;
  char buf[64];

  printf(pre_str);
  for (i=0, j=0; i<tb.nv; i++) {
    if (lhs && tb.lh[i]>0) continue;
    j++;
    retieve_var(i, buf);
    switch(scenario) {
      case 0: printf(i<tb.nv-1 ? "\t%s,\n" : "\t%s;\n", buf); break;
      case 1: printf("\t%s = par_ptr[%d];\n", buf, j-1); break;
      default: break;
    }
  }
  printf(post_str);
  scenario++;
}

void codegen() {
  int i, j;
  char sLine[MXLEN+1];
  char buf[64];

  char *hdft[]=
    {
      "#include <math.h>\nextern long dadt_counter;\nextern double InfusionRate[99];\nextern double *par_ptr;\nextern double podo;\nextern double tlast;\n\n// prj-specific differential eqns\nvoid dydt(unsigned int neq, double t, double *A, double *DADT)\n{\n",
      "    dadt_counter++;\n}\n\n"
    };

  fpIO = fopen( "out2.txt", "r" );
  err_msg((int) fpIO, "Coudln't access out2.txt.\n", -1);

  printf("%s", hdft[0]);
  prnt_vars(0, "double\n", "\n");	/* declare all used vars */
  prnt_vars(1, "", "\n");			/* pass system pars */

  for (i=0; i<tb.nd; i++) {			/* name state vars */
    retieve_var(tb.di[i], buf);
    printf("\t%s = A[%d];\n", buf, i);
  }
  printf("\n");

  while(fgets(sLine, MXLEN, fpIO))	/* parsed eqns */
    printf("\t%s", sLine);
  printf("%s", hdft[1]);

#ifdef __RODE__
  //output nder & system pars -- to be read by dvode() in R.
  fprintf(stderr, "%d ", tb.nd);
  for (i=0; i<tb.nv; i++) {
    if (tb.lh[i]>0) continue;
	retieve_var(i, buf);
	fprintf(stderr, "%s ", buf);
  }
#endif

  fclose(fpIO);
}

void inits() {
  tb.ss = (char *) malloc(64*MXSYM);
  err_msg((int) tb.ss, "error allocating vars", 1);

  tb.vo[0]=0;
  memset(tb.lh, 0, MXSYM);
  tb.nv=0;
  tb.nd=0;
  tb.fn=0;
}

int main(int argc, char *argv[]) {
  char *buf;
  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_gram, 1024);
  p->save_parse_tree = 1;

  if (argc!=2) {
    fprintf(stderr,"Usage: %s FILE_to_parse\n",argv[0]);
    return -1;
  } else {
    buf = sbuf_read(argv[1]);
    err_msg((int) buf, "error: empty buf\n", -2);
  }

  if ((pn=dparse(p, buf, strlen(buf))) && !p->syntax_errors) {
    inits();
    fpIO = fopen( "out2.txt", "w" );
#ifndef __RODE__
    fp_inits = fopen( "inits.txt", "w" );
#endif
    err_msg((int) fpIO, "error opening out2.txt\n", -2);
    wprint_parsetree(parser_tables_gram, pn, 0, wprint_node, NULL);
    fclose(fpIO);
    if (fp_inits) fclose(fp_inits);

    codegen();
    remove("out2.txt");
    free(tb.ss);
  } else {
    printf("\nfailure\n");
  }
  return 0;
}

