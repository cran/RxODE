#define USE_FC_LEN_T
#define STRICT_R_HEADER
#include "sbuf.h"
#include "tran.h"

// Taken from dparser and changed to use Calloc
int rc_buf_read(const char *pathname, char **buf, int *len) {
  struct stat sb;
  int fd;
  *buf = 0;
  *len = 0;
  fd = open(pathname, O_RDONLY);
  if (fd <= 0)
    return -1;
  memset(&sb, 0, sizeof(sb));
  fstat(fd, &sb);
  *len = sb.st_size;
  *buf = Calloc(*len + 3,char);
  // MINGW likes to convert cr lf => lf which messes with the size
  size_t real_size = read(fd, *buf, *len);
  (*buf)[real_size] = 0;
  (*buf)[real_size + 1] = 0;
  *len = real_size;
  close(fd);
  return *len;
}

// Taken from dparser and changed to use Calloc
char * rc_sbuf_read(const char *pathname) {
  char *buf;
  int len;
  if (rc_buf_read(pathname, &buf, &len) < 0)
    return NULL;
  return buf;
}



void sIniTo(sbuf *sbb, int to) {
  if (sbb->s != NULL) Free(sbb->s);
  sbb->s    = Calloc(to, char);
  sbb->sN   = to;
  sbb->s[0] = '\0';
  sbb->o    = 0;
}

void sIni(sbuf *sbb) {
  sIniTo(sbb, SBUF_MXBUF);
}

void sFree(sbuf *sbb) {
  if (sbb->s != NULL) Free(sbb->s);
  sNull(sbb);
}

void sFreeIni(sbuf *sbb) {
  sFree(sbb);
  sIni(sbb);
}

void sAppendN(sbuf *sbb, const char *what, int n) {
  if (sbb->sN == 0) sIni(sbb);
  if (sbb->sN <= 2 + n + sbb->o){
    int mx = sbb->o + 2 + n + SBUF_MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    sbb->sN = mx;
  }
  sprintf(sbb->s+sbb->o, "%s", what);
  sbb->o +=n;
}

void sAppend(sbuf *sbb, const char *format, ...) {
  if (sbb->sN == 0) sIni(sbb);
  if (format == NULL) return;
  int n = 0;
  va_list argptr, copy;
  va_start(argptr, format);
  va_copy(copy, argptr);
#if defined(_WIN32) || defined(WIN32)
  n = vsnprintf(NULL, 0, format, copy) + 1;
#else
  char zero[2];
  n = vsnprintf(zero, 0, format, copy) + 1;
#endif
  va_end(copy);
  if (sbb->sN <= sbb->o + n + 1) {
    int mx = sbb->o + n + 1 + SBUF_MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    sbb->sN = mx;
  }
  vsnprintf(sbb->s+ sbb->o, sbb->sN - sbb->o, format, argptr);
  va_end(argptr);
  sbb->o += n-1;
}

void sPrint(sbuf *sbb, const char *format, ...) {
  if (sbb->sN == 0) sIni(sbb);
  sClear(sbb);
  if (format == NULL) return;
  int n = 0;
  va_list argptr, copy;
  va_start(argptr, format);
  va_copy(copy, argptr);
#if defined(_WIN32) || defined(WIN32)
  n = vsnprintf(NULL, 0, format, copy) + 1;
#else
  char zero[2];
  n = vsnprintf(zero, 0, format, copy) + 1;
#endif
  va_end(copy);
  if (sbb->sN <= sbb->o + n + 1){
    int mx = sbb->o + n + 1 + SBUF_MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    sbb->sN = mx;
  }
  vsnprintf(sbb->s+ sbb->o, sbb->sN - sbb->o, format, argptr);
  va_end(argptr);
  sbb->o += n-1;
}

void lineIni(vLines *sbb) {
  if (sbb->s != NULL) Free(sbb->s);
  sbb->s = Calloc(SBUF_MXBUF, char);
  sbb->sN = SBUF_MXBUF;
  sbb->s[0]='\0';
  sbb->o = 0;
  if (sbb->lProp != NULL) Free(sbb->lProp);
  if (sbb->line != NULL) Free(sbb->line);
  if (sbb->lType != NULL) Free(sbb->lType);
  if (sbb->os != NULL) Free(sbb->os);
  sbb->lProp = Calloc(SBUF_MXLINE, int);
  sbb->lType = Calloc(SBUF_MXLINE, int);
  sbb->line = Calloc(SBUF_MXLINE, char*);
  sbb->os = Calloc(SBUF_MXLINE, int);
  sbb->nL=SBUF_MXLINE;
  sbb->lProp[0] = -1;
  sbb->lType[0] = 0;
  sbb->n = 0;
}

void lineFree(vLines *sbb) {
  if (sbb->s != NULL) Free(sbb->s);
  if (sbb->lProp != NULL) Free(sbb->lProp);
  if (sbb->line != NULL) Free(sbb->line);
  if (sbb->lType != NULL) Free(sbb->lType);
  if (sbb->os != NULL) Free(sbb->os);
  lineNull(sbb);
}

void addLine(vLines *sbb, const char *format, ...) {
  if (sbb->sN == 0) lineIni(sbb);
  if (format == NULL) return;
  int n = 0;
  va_list argptr, copy;
  va_start(argptr, format);
  va_copy(copy, argptr);
  errno = 0;
  // Try first.
#if defined(_WIN32) || defined(WIN32)
  n = vsnprintf(NULL, 0, format, copy);
#else
  char zero[2];
  n = vsnprintf(zero, 0, format, copy);
#endif
  if (n < 0){
    parseFree(0);
    Rf_errorcall(R_NilValue, _("encoding error in 'addLine' format: '%s' n: %d; errno: %d"), format, n, errno);
  }
  va_end(copy);
  if (sbb->sN <= sbb->o + n){
    int mx = sbb->sN + n + 2 + SBUF_MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    // The sbb->line are not correct any longer because the pointer for sbb->s has been updated;
    // Fix them
    for (int i = sbb->n; i--;){
      sbb->line[i] = &(sbb->s[sbb->os[i]]);
    }
    sbb->sN = mx;
  }
  vsnprintf(sbb->s + sbb->o, sbb->sN - sbb->o, format, argptr);
  va_end(argptr);
  if (sbb->n + 2 >= sbb->nL){
    int mx = sbb->nL + n + 2 + SBUF_MXLINE;
    sbb->lProp = Realloc(sbb->lProp, mx, int);
    sbb->lType = Realloc(sbb->lType, mx, int);
    sbb->line = Realloc(sbb->line, mx, char*);
    sbb->os = Realloc(sbb->os, mx, int);
    sbb->nL = mx;
  }
  sbb->line[sbb->n]=&(sbb->s[sbb->o]);
  sbb->os[sbb->n]= sbb->o;
  sbb->o += n + 1; // n should include the \0 character
  sbb->n = sbb->n+1;
  sbb->lProp[sbb->n] = -1;
  sbb->lType[sbb->n] = 0;
  sbb->os[sbb->n]= sbb->o;
}

void curLineProp(vLines *sbb, int propId){
  sbb->lProp[sbb->n] = propId;
}

void curLineType(vLines *sbb, int propId) {
  sbb->lType[sbb->n] = propId;
}

void doDot(sbuf *out, char *buf) {
  for (int k = 0; k < (int)strlen(buf); k++){
    if (buf[k] == '.'){
      sAppend(out,"_DoT_");
      if (rx_syntax_allow_dots == 0){
	updateSyntaxCol();
	trans_syntax_error_report_fn(NODOT);
      }
    } else {
      sPut(out,buf[k]);
    }
  }
}

void doDot2(sbuf *sb, sbuf *sbDt, char *buf) {
  for (int k = 0; k < (int)strlen(buf); k++) {
    if (buf[k] == '.') {
      sAppend(sb,"_DoT_");
      sAppend(sbDt,"_DoT_");
      if (rx_syntax_allow_dots == 0){
	updateSyntaxCol();
	trans_syntax_error_report_fn(NODOT);
      }
    } else {
      sPut(sb,buf[k]);
      sPut(sbDt,buf[k]);
    }
  }
}


