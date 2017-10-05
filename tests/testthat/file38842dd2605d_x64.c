#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if (_jac_counter_val() == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if (_dadt_counter_val() == 0) Rprintf
#define LHS_Rprintf Rprintf
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define R_pow Rx_pow
#define R_pow_di Rx_pow_di

// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t);
typedef double (*RxODE_transit3)(double t, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef double (*RxODE_transit4)(double t, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val);
typedef long (*RxODE_cnt) ();
typedef void (*RxODE_inc) ();
typedef double (*RxODE_val) ();
typedef SEXP (*RxODE_ode_solver) (SEXP sexp_theta, SEXP sexp_inits, SEXP sexp_lhs, SEXP sexp_time, SEXP sexp_evid,SEXP sexp_dose, SEXP sexp_pcov, SEXP sexp_cov, SEXP sexp_locf, SEXP sexp_atol, SEXP sexp_rtol, SEXP sexp_hmin, SEXP sexp_hmax, SEXP sexp_h0, SEXP sexp_mxordn, SEXP sexp_mxords, SEXP sexp_mx,SEXP sexp_stiff, SEXP sexp_transit_abs, SEXP sexp_object, SEXP sexp_extra_args, SEXP sexp_matrix, SEXP sexp_add_cov);
typedef void (*RxODE_assign_fn_pointers)(void (*fun_dydt)(unsigned int, double, double *, double *),void (*fun_calc_lhs)(double, double *, double *),void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),void (*fun_update_inis)(SEXP _ini_sexp),int fun_jt,int fun_mf, int fun_debug);

typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef void (*RxODE_ode_solver_0_6_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc,double hmin, double hmax,double h0,int mxordn,int mxords,int mxstep);
typedef double (*RxODE_solveLinB)(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
typedef double (*RxODE_sum_prod)(double *input, int n);
// Give par pointers
RxODE_vec _par_ptr, _InfusionRate;
RxODE_update_par_ptr _update_par_ptr;
RxODE_cnt _dadt_counter_val, _jac_counter_val;
RxODE_inc _dadt_counter_inc, _jac_counter_inc;
RxODE_val podo, tlast;
RxODE_transit4 _transit4;
RxODE_transit3 _transit3;
RxODE_fn _safe_log, safe_zero, factorial, _as_zero, abs_log, abs_log1p;
RxODE_fn2 sign_exp, Rx_pow;
RxODE_fn2i Rx_pow_di;
RxODE_assign_fn_pointers _assign_fn_pointers;
RxODE_ode_solver_old_c _old_c;
RxODE_ode_solver_0_6_c _c_0_6;
RxODE_solveLinB solveLinB;
RxODE_sum_prod _sum1, _prod1;

double _sum(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _sum1(p, n);
  Free(p);
  return s;
}

double _prod(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _prod1(p, n);
  Free(p);
  return s;
}

extern double _sign(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double s = 1;
  for (unsigned int i = 0; i < n; i++){
    s = sign(va_arg(valist, double))*s;
    if (s == 0){
      break;
    }
  }
  va_end(valist);
  return s;
}

extern void file38842dd2605d_x64_ode_solver_ptr();



// prj-specific differential eqns
void file38842dd2605d_x64_dydt(unsigned int _neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double 
		C2,
	centr,
	V2,
	C3,
	peri,
	V3,
	depot,
	KA,
	CL,
	Q,
	eff,
	Kin,
	Kout,
	EC50;

	(void)t;	(void)C2;
	(void)centr;
	(void)V2;
	(void)C3;
	(void)peri;
	(void)V3;
	(void)depot;
	(void)KA;
	(void)CL;
	(void)Q;
	(void)eff;
	(void)Kin;
	(void)Kout;
	(void)EC50;

	_update_par_ptr(t);
	V2 = _par_ptr(0);
	V3 = _par_ptr(1);
	KA = _par_ptr(2);
	CL = _par_ptr(3);
	Q = _par_ptr(4);
	Kin = _par_ptr(5);
	Kout = _par_ptr(6);
	EC50 = _par_ptr(7);

	depot = __zzStateVar__[0];
	centr = __zzStateVar__[1];
	peri = __zzStateVar__[2];
	eff = __zzStateVar__[3];

	C2=centr/safe_zero(V2);
	C3=peri/safe_zero(V3);
	__DDtStateVar__[0]=_InfusionRate(0)+-KA*depot;
	__DDtStateVar__[1]=_InfusionRate(1)+KA*depot-CL*C2-Q*C2+Q*C3;
	__DDtStateVar__[2]=_InfusionRate(2)+Q*C2-Q*C3;
	__DDtStateVar__[3]=_InfusionRate(3)+Kin-Kout*(1-C2/safe_zero((EC50+C2)))*eff;
    _dadt_counter_inc();
}

// Jacobian derived vars
void file38842dd2605d_x64_calc_jac(unsigned int _neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  _jac_counter_inc();
}
// Functional based initial conditions.
void file38842dd2605d_x64_inis(SEXP _ini_sexp){
	double *__zzStateVar__ = REAL(_ini_sexp);
	double t=0;
double 
		C2,
	centr,
	V2,
	C3,
	peri,
	V3,
	depot,
	KA,
	CL,
	Q,
	eff,
	Kin,
	Kout,
	EC50;

	(void)t;	(void)C2;
	(void)centr;
	(void)V2;
	(void)C3;
	(void)peri;
	(void)V3;
	(void)depot;
	(void)KA;
	(void)CL;
	(void)Q;
	(void)eff;
	(void)Kin;
	(void)Kout;
	(void)EC50;

	_update_par_ptr(0.0);
	V2 = _par_ptr(0);
	V3 = _par_ptr(1);
	KA = _par_ptr(2);
	CL = _par_ptr(3);
	Q = _par_ptr(4);
	Kin = _par_ptr(5);
	Kout = _par_ptr(6);
	EC50 = _par_ptr(7);

	depot = __zzStateVar__[0];
	centr = __zzStateVar__[1];
	peri = __zzStateVar__[2];
	eff = __zzStateVar__[3];

	C2=centr/safe_zero(V2);
	C3=peri/safe_zero(V3);
  __zzStateVar__[0]=depot;
  __zzStateVar__[1]=centr;
  __zzStateVar__[2]=peri;
  __zzStateVar__[3]=eff;
}
// prj-specific derived vars
void file38842dd2605d_x64_calc_lhs(double t, double *__zzStateVar__, double *_lhs) {
double 
		__DDtStateVar_0__,
	__DDtStateVar_1__,
	__DDtStateVar_2__,
	__DDtStateVar_3__,
	C2,
	centr,
	V2,
	C3,
	peri,
	V3,
	depot,
	KA,
	CL,
	Q,
	eff,
	Kin,
	Kout,
	EC50;

	(void)t;	(void)__DDtStateVar_0__;
	(void)__DDtStateVar_1__;
	(void)__DDtStateVar_2__;
	(void)__DDtStateVar_3__;
	(void)C2;
	(void)centr;
	(void)V2;
	(void)C3;
	(void)peri;
	(void)V3;
	(void)depot;
	(void)KA;
	(void)CL;
	(void)Q;
	(void)eff;
	(void)Kin;
	(void)Kout;
	(void)EC50;

	_update_par_ptr(t);
	V2 = _par_ptr(0);
	V3 = _par_ptr(1);
	KA = _par_ptr(2);
	CL = _par_ptr(3);
	Q = _par_ptr(4);
	Kin = _par_ptr(5);
	Kout = _par_ptr(6);
	EC50 = _par_ptr(7);

	depot = __zzStateVar__[0];
	centr = __zzStateVar__[1];
	peri = __zzStateVar__[2];
	eff = __zzStateVar__[3];

	C2=centr/safe_zero(V2);
	C3=peri/safe_zero(V3);
	__DDtStateVar_0__=_InfusionRate(0)+-KA*depot;
	__DDtStateVar_1__=_InfusionRate(1)+KA*depot-CL*C2-Q*C2+Q*C3;
	__DDtStateVar_2__=_InfusionRate(2)+Q*C2-Q*C3;
	__DDtStateVar_3__=_InfusionRate(3)+Kin-Kout*(1-C2/safe_zero((EC50+C2)))*eff;

	_lhs[0]=C2;
	_lhs[1]=C3;
}
extern SEXP file38842dd2605d_x64_model_vars(){
	SEXP lst      = PROTECT(allocVector(VECSXP, 12));
	SEXP names    = PROTECT(allocVector(STRSXP, 12));
	SEXP params   = PROTECT(allocVector(STRSXP, 8));
	SEXP lhs      = PROTECT(allocVector(STRSXP, 2));
	SEXP state    = PROTECT(allocVector(STRSXP, 4));
	SEXP stateRmS = PROTECT(allocVector(INTSXP, 4));
	int *stateRm  = INTEGER(stateRmS);
	SEXP sens     = PROTECT(allocVector(STRSXP, 0));
	SEXP fn_ini   = PROTECT(allocVector(STRSXP, 0));
	SEXP dfdy     = PROTECT(allocVector(STRSXP, 0));
	SEXP tran     = PROTECT(allocVector(STRSXP, 12));
	SEXP trann    = PROTECT(allocVector(STRSXP, 12));
	SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));
	SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));
	SEXP model    = PROTECT(allocVector(STRSXP, 4));
	SEXP modeln   = PROTECT(allocVector(STRSXP, 4));
	SET_STRING_ELT(lhs,0,mkChar("C2"));
	SET_STRING_ELT(params,0,mkChar("V2"));
	SET_STRING_ELT(lhs,1,mkChar("C3"));
	SET_STRING_ELT(params,1,mkChar("V3"));
	SET_STRING_ELT(params,2,mkChar("KA"));
	SET_STRING_ELT(params,3,mkChar("CL"));
	SET_STRING_ELT(params,4,mkChar("Q"));
	SET_STRING_ELT(params,5,mkChar("Kin"));
	SET_STRING_ELT(params,6,mkChar("Kout"));
	SET_STRING_ELT(params,7,mkChar("EC50"));
	SET_STRING_ELT(state,0,mkChar("depot"));
	stateRm[0] = 0;
	SET_STRING_ELT(state,1,mkChar("centr"));
	stateRm[1] = 0;
	SET_STRING_ELT(state,2,mkChar("peri"));
	stateRm[2] = 0;
	SET_STRING_ELT(state,3,mkChar("eff"));
	stateRm[3] = 0;
	SET_STRING_ELT(modeln,0,mkChar("model"));
	SET_STRING_ELT(model,0,mkChar("\n         C2 = centr/V2;\n         C3 = peri/V3;\n         d/dt(depot) =-KA*depot;\n         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;\n         d/dt(peri)  =                    Q*C2 - Q*C3;\n         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;"));
	SET_STRING_ELT(modeln,1,mkChar("normModel"));
	SET_STRING_ELT(model,1,mkChar("C2=centr/V2;\nC3=peri/V3;\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n"));
	SET_STRING_ELT(modeln,2,mkChar("parseModel"));
	SET_STRING_ELT(model,2,mkChar("C2=centr/safe_zero(V2);\nC3=peri/safe_zero(V3);\n__DDtStateVar__[0] = _InfusionRate(0) + -KA*depot;\n__DDtStateVar__[1] = _InfusionRate(1) + KA*depot-CL*C2-Q*C2+Q*C3;\n__DDtStateVar__[2] = _InfusionRate(2) + Q*C2-Q*C3;\n__DDtStateVar__[3] = _InfusionRate(3) + Kin-Kout*(1-C2/safe_zero((EC50+C2)))*eff;\n"));
	SET_STRING_ELT(modeln,3,mkChar("expandModel"));
	SET_STRING_ELT(model,3,mkChar("\n         C2 = centr/V2;\n         C3 = peri/V3;\n         d/dt(depot) =-KA*depot;\n         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;\n         d/dt(peri)  =                    Q*C2 - Q*C3;\n         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;"));
	SEXP ini    = PROTECT(allocVector(REALSXP,0));
	SEXP inin   = PROTECT(allocVector(STRSXP, 0));
	SET_STRING_ELT(names,0,mkChar("params"));
	SET_VECTOR_ELT(lst,  0,params);
	SET_STRING_ELT(names,1,mkChar("lhs"));
	SET_VECTOR_ELT(lst,  1,lhs);
	SET_STRING_ELT(names,2,mkChar("state"));
	SET_VECTOR_ELT(lst,  2,state);
	SET_STRING_ELT(names,3,mkChar("trans"));
	SET_VECTOR_ELT(lst,  3,tran);
	SET_STRING_ELT(names,5,mkChar("model"));
	SET_VECTOR_ELT(lst,  5,model);
	SET_STRING_ELT(names,4,mkChar("ini"));
	SET_VECTOR_ELT(lst,  4,ini);
	SET_STRING_ELT(names,6,mkChar("md5"));
	SET_VECTOR_ELT(lst,  6,mmd5);
	SET_STRING_ELT(names,7,mkChar("podo"));
	SET_VECTOR_ELT(lst,  7,ScalarLogical(0));
	SET_STRING_ELT(names,8,mkChar("dfdy"));
	SET_VECTOR_ELT(lst,  8,dfdy);
	SET_STRING_ELT(names,9,mkChar("sens"));
	SET_VECTOR_ELT(lst,  9,sens);
	SET_STRING_ELT(names,10,mkChar("fn.ini"));
	SET_VECTOR_ELT(lst,  10,fn_ini);
	SET_STRING_ELT(names,11,mkChar("state.ignore"));
	SET_VECTOR_ELT(lst,  11,stateRmS);
	SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
	SET_STRING_ELT(mmd5,0,mkChar("4e0bffaa4479a301beedbb697218d9dc"));
	SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
	SET_STRING_ELT(mmd5,1,mkChar("af410a14969a98b2b69123f188a415ba"));
	SET_STRING_ELT(trann,0,mkChar("jac"));
	SET_STRING_ELT(tran,0,mkChar("fullint"));
	SET_STRING_ELT(trann,1,mkChar("prefix"));
	SET_STRING_ELT(tran, 1,mkChar("file38842dd2605d_x64_"));
	SET_STRING_ELT(trann,2,mkChar("dydt"));
	SET_STRING_ELT(tran, 2,mkChar("file38842dd2605d_x64_dydt"));
	SET_STRING_ELT(trann,3,mkChar("calc_jac"));
	SET_STRING_ELT(tran, 3,mkChar("file38842dd2605d_x64_calc_jac"));
	SET_STRING_ELT(trann,4,mkChar("calc_lhs"));
	SET_STRING_ELT(tran, 4,mkChar("file38842dd2605d_x64_calc_lhs"));
	SET_STRING_ELT(trann,5,mkChar("model_vars"));
	SET_STRING_ELT(tran, 5,mkChar("file38842dd2605d_x64_model_vars"));
	SET_STRING_ELT(trann,6,mkChar("ode_solver"));
	SET_STRING_ELT(tran, 6,mkChar("file38842dd2605d_x64_ode_solver"));
	SET_STRING_ELT(trann,7,mkChar("ode_solver_sexp"));
	SET_STRING_ELT(tran, 7,mkChar("file38842dd2605d_x64_ode_solver_sexp"));
	SET_STRING_ELT(trann,8,mkChar("ode_solver_0_6"));
	SET_STRING_ELT(tran, 8,mkChar("file38842dd2605d_x64_ode_solver_0_6"));
	SET_STRING_ELT(trann,9,mkChar("ode_solver_focei_eta"));
	SET_STRING_ELT(tran, 9,mkChar("file38842dd2605d_x64_ode_solver_focei_eta"));
	SET_STRING_ELT(trann,10,mkChar("ode_solver_ptr"));
	SET_STRING_ELT(tran, 10,mkChar("file38842dd2605d_x64_ode_solver_ptr"));
	SET_STRING_ELT(trann,11,mkChar("inis"));
	SET_STRING_ELT(tran, 11,mkChar("file38842dd2605d_x64_inis"));
	setAttrib(tran, R_NamesSymbol, trann);
	setAttrib(mmd5, R_NamesSymbol, mmd5n);
	setAttrib(model, R_NamesSymbol, modeln);
	setAttrib(ini, R_NamesSymbol, inin);
	setAttrib(lst, R_NamesSymbol, names);
	UNPROTECT(17);
	return lst;
}
extern void file38842dd2605d_x64_ode_solver(
                    int *neq,
                    double *theta,      //order:
                    double *time,
                    int *evid,
                    int *ntime,
                    double *inits,
                    double *dose,
                    double *ret,
                    double *atol,
                    double *rtol,
                    int *stiff,
                    int *transit_abs,
                    int *nlhs,
                    double *lhs,
                    int *rc
                    ){
  // Backward compatible ode solver for 0.5* C interface
  file38842dd2605d_x64_ode_solver_ptr();
  _old_c(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

void file38842dd2605d_x64_ode_solver_0_6(int *neq,
                        double *theta,  //order:
                        double *time,
                        int *evid,
                        int *ntime,
                        double *inits,
                        double *dose,
                        double *ret,
                        double *atol,
                        double *rtol,
                        int *stiff,
                        int *transit_abs,
                        int *nlhs,
                        double *lhs,
                        int *rc,
                        double hmin,
                        double hmax,
                        double h0,
                        int mxordn,
                        int mxords,
                        int mxstep) {
  // Backward compatible ode solver for 0.5* C interface
  file38842dd2605d_x64_ode_solver_ptr();
  _c_0_6(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc,
	hmin, hmax, h0, mxordn, mxords, mxstep);
}

extern void file38842dd2605d_x64_ode_solver_ptr  (){
  _assign_fn_pointers(file38842dd2605d_x64_dydt , file38842dd2605d_x64_calc_lhs , file38842dd2605d_x64_calc_jac, file38842dd2605d_x64_inis, 2 , 22,
#ifdef __DEBUG__
                      1
#else
                      0
#endif
                      );
}

extern SEXP file38842dd2605d_x64_ode_solver_sexp (// Parameters
                                 SEXP sexp_theta,
                                 SEXP sexp_inits,
                                 SEXP sexp_lhs,
				 // Events
				 SEXP sexp_time,
				 SEXP sexp_evid,
				 SEXP sexp_dose,
				 // Covariates
				 SEXP sexp_pcov,
				 SEXP sexp_cov,
				 SEXP sexp_locf,
				 // Solver Options
				 SEXP sexp_atol,
				 SEXP sexp_rtol,
				 SEXP sexp_hmin,
				 SEXP sexp_hmax,
				 SEXP sexp_h0,
				 SEXP sexp_mxordn,
				 SEXP sexp_mxords,
				 SEXP sexp_mx,
				 SEXP sexp_stiff,
				 SEXP sexp_transit_abs,
				 // Object Creation
				 SEXP sexp_object,
				 SEXP sexp_extra_args,
				 SEXP sexp_matrix,
				 SEXP sexp_add_cov){
  file38842dd2605d_x64_ode_solver_ptr();
  RxODE_ode_solver ode_solver = (RxODE_ode_solver) R_GetCCallable("RxODE","RxODE_ode_solver");
  return ode_solver(sexp_theta,sexp_inits,sexp_lhs,sexp_time,sexp_evid,sexp_dose,sexp_pcov,sexp_cov,sexp_locf,sexp_atol,
		    sexp_rtol,sexp_hmin,sexp_hmax,sexp_h0,sexp_mxordn,sexp_mxords,sexp_mx,sexp_stiff,sexp_transit_abs,
		    sexp_object,sexp_extra_args,sexp_matrix,sexp_add_cov);
}

static R_NativePrimitiveArgType file38842dd2605d_x64_ode_solverrx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType  file38842dd2605d_x64_ode_solver_0_6rx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

//Initilize the dll to match RxODE's calls
void R_init_file38842dd2605d_x64 (DllInfo *info){
  // Get the RxODE calling interfaces
  _InfusionRate   = (RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRate");
  _update_par_ptr = (RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptr");
  _par_ptr = (RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptr");
  _dadt_counter_val = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_val");
  _jac_counter_val  = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_val");
  _dadt_counter_inc = (RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_inc");
  _jac_counter_inc  = (RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_inc");
  podo  = (RxODE_val) R_GetCCallable("RxODE","RxODE_podo");
  tlast = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlast");
  factorial=(RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _transit3 = (RxODE_transit3) R_GetCCallable("RxODE","RxODE_transit3");
  _transit4 = (RxODE_transit4) R_GetCCallable("RxODE","RxODE_transit4");
  _safe_log =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  safe_zero =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  _as_zero = (RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _assign_fn_pointers=(RxODE_assign_fn_pointers) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","RxODE_ode_solver_old_c");
  _c_0_6 = (RxODE_ode_solver_0_6_c)R_GetCCallable("RxODE","RxODE_ode_solver_0_6_c");
  _sum1   = (RxODE_sum_prod)R_GetCCallable("RxODE","RxODE_sum");
  _prod1 = (RxODE_sum_prod) R_GetCCallable("RxODE","RxODE_prod");
  sign_exp = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  Rx_pow = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_pow");
  Rx_pow_di = (RxODE_fn2i) R_GetCCallable("RxODE","RxODE_pow_di");
  abs_log = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  abs_log1p = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  solveLinB = (RxODE_solveLinB) R_GetCCallable("RxODE","RxODE_solveLinB");
  // Register the outside functions
  R_RegisterCCallable("file38842dd2605d_x64","file38842dd2605d_x64_ode_solver",       (DL_FUNC) file38842dd2605d_x64_ode_solver);
  R_RegisterCCallable("file38842dd2605d_x64","file38842dd2605d_x64_ode_solver_sexp",  (DL_FUNC) file38842dd2605d_x64_ode_solver_sexp);
  R_RegisterCCallable("file38842dd2605d_x64","file38842dd2605d_x64_ode_solver_0_6",   (DL_FUNC) file38842dd2605d_x64_ode_solver_0_6);
  R_RegisterCCallable("file38842dd2605d_x64","file38842dd2605d_x64_ode_solver_ptr",   (DL_FUNC) file38842dd2605d_x64_ode_solver_ptr);

  /* R_registerRoutines(info, NULL, NULL, NULL, NULL); */
  /* R_useDynamicSymbols(info,TRUE); */

  static const R_CMethodDef cMethods[] = {
    {"file38842dd2605d_x64_ode_solver", (DL_FUNC) &file38842dd2605d_x64_ode_solver, 15, file38842dd2605d_x64_ode_solverrx_t},
    {"file38842dd2605d_x64_ode_solver_0_6", (DL_FUNC) &file38842dd2605d_x64_ode_solver_0_6, 21, file38842dd2605d_x64_ode_solver_0_6rx_t},
    {NULL, NULL, 0, NULL}
  };
  
  R_CallMethodDef callMethods[]  = {
    {"file38842dd2605d_x64_ode_solver_ptr", (DL_FUNC) &file38842dd2605d_x64_ode_solver_ptr, 0},
    {"file38842dd2605d_x64_ode_solver_sexp", (DL_FUNC) &file38842dd2605d_x64_ode_solver_sexp, 23},
    {"file38842dd2605d_x64_model_vars", (DL_FUNC) &file38842dd2605d_x64_model_vars, 0},
    {"file38842dd2605d_x64_inis", (DL_FUNC) &file38842dd2605d_x64_inis, 1},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
  // Register the function pointers so if someone directly calls the
  // ode solvers directly, they use the last loaded RxODE model.
  file38842dd2605d_x64_ode_solver_ptr();
}
