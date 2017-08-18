/* File: dfitpackmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * See http://cens.ioc.ee/projects/f2py2e/
 * Generation date: Fri Aug 18 19:08:32 2017
 * $Revision:$
 * $Date:$
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *dfitpack_error;
static PyObject *dfitpack_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
  PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
  fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#define CHECKARRAY(check,tcheck,name) \
  if (!(check)) {\
    PyErr_SetString(dfitpack_error,"("tcheck") failed for "name);\
    /*goto capi_fail;*/\
  } else 
#define CHECKSCALAR(check,tcheck,name,show,var)\
  if (!(check)) {\
    char errstring[256];\
    sprintf(errstring, "%s: "show, "("tcheck") failed for "name, var);\
    PyErr_SetString(dfitpack_error,errstring);\
    /*goto capi_fail;*/\
  } else 

/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int int_from_pyobj(int* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyInt_Check(obj)) {
    *v = (int)PyInt_AS_LONG(obj);
    return 1;
  }
  tmp = PyNumber_Int(obj);
  if (tmp) {
    *v = PyInt_AS_LONG(tmp);
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (int_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = dfitpack_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}

static int double_from_pyobj(double* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyFloat_Check(obj)) {
#ifdef __sgi
    *v = PyFloat_AsDouble(obj);
#else
    *v = PyFloat_AS_DOUBLE(obj);
#endif
    return 1;
  }
  tmp = PyNumber_Float(obj);
  if (tmp) {
#ifdef __sgi
    *v = PyFloat_AsDouble(tmp);
#else
    *v = PyFloat_AS_DOUBLE(tmp);
#endif
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = dfitpack_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/
  /* start usercode multiline (0) */


static double dmax(double* seq,int len) {
  double val;
  int i;
  if (len<1)
    return -1e308;
  val = seq[0];
  for(i=1;i<len;++i)
    if (seq[i]>val) val = seq[i];
  return val;
}
static double dmin(double* seq,int len) {
  double val;
  int i;
  if (len<1)
    return 1e308;
  val = seq[0];
  for(i=1;i<len;++i)
    if (seq[i]<val) val = seq[i];
  return val;
}

static int calc_regrid_lwrk(int mx, int my, int kx, int ky,
                            int nxest, int nyest) {
 int u = MAX(my, nxest);
 return 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+u;
}

  /* end multiline (0)*/

/* See f2py2e/rules.py */
extern void F_FUNC(fpchec,FPCHEC)(double*,int*,double*,int*,int*,int*);
extern void F_FUNC(bispeu,BISPEU)(double*,int*,double*,int*,double*,int*,int*,double*,double*,double*,int*,double*,int*,int*);
extern void F_FUNC(regrid,REGRID)(int*,int*,double*,int*,double*,double*,double*,double*,double*,double*,int*,int*,double*,int*,int*,int*,double*,int*,double*,double*,double*,double*,int*,int*,int*,int*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/*********************************** fpchec ***********************************/
static char doc_f2py_rout_dfitpack_fpchec[] = "\
ier = fpchec(x,t,k)\n\nWrapper for ``fpchec``.\
\n\nParameters\n----------\n"
"x : input rank-1 array('d') with bounds (m)\n"
"t : input rank-1 array('d') with bounds (n)\n"
"k : input int\n"
"\nReturns\n-------\n"
"ier : int";
/* extern void F_FUNC(fpchec,FPCHEC)(double*,int*,double*,int*,int*,int*); */
static PyObject *f2py_rout_dfitpack_fpchec(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,int*,double*,int*,int*,int*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double *x = NULL;
  npy_intp x_Dims[1] = {-1};
  const int x_Rank = 1;
  PyArrayObject *capi_x_tmp = NULL;
  int capi_x_intent = 0;
  PyObject *x_capi = Py_None;
  int m = 0;
  double *t = NULL;
  npy_intp t_Dims[1] = {-1};
  const int t_Rank = 1;
  PyArrayObject *capi_t_tmp = NULL;
  int capi_t_intent = 0;
  PyObject *t_capi = Py_None;
  int n = 0;
  int k = 0;
  PyObject *k_capi = Py_None;
  int ier = 0;
  static char *capi_kwlist[] = {"x","t","k",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOO:dfitpack.fpchec",\
    capi_kwlist,&x_capi,&t_capi,&k_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable x */
  ;
  capi_x_intent |= F2PY_INTENT_IN;
  capi_x_tmp = array_from_pyobj(NPY_DOUBLE,x_Dims,x_Rank,capi_x_intent,x_capi);
  if (capi_x_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 1st argument `x' of dfitpack.fpchec to C/Fortran array" );
  } else {
    x = (double *)(PyArray_DATA(capi_x_tmp));

  /* Processing variable t */
  ;
  capi_t_intent |= F2PY_INTENT_IN;
  capi_t_tmp = array_from_pyobj(NPY_DOUBLE,t_Dims,t_Rank,capi_t_intent,t_capi);
  if (capi_t_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 2nd argument `t' of dfitpack.fpchec to C/Fortran array" );
  } else {
    t = (double *)(PyArray_DATA(capi_t_tmp));

  /* Processing variable k */
    f2py_success = int_from_pyobj(&k,k_capi,"dfitpack.fpchec() 3rd argument (k) can't be converted to int");
  if (f2py_success) {
  /* Processing variable ier */
  /* Processing variable m */
  m = len(x);
  /* Processing variable n */
  n = len(t);
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
      Py_BEGIN_ALLOW_THREADS
        (*f2py_func)(x,&m,t,&n,&k,&ier);
      Py_END_ALLOW_THREADS
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("i",ier);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  /* End of cleaning variable n */
  /* End of cleaning variable m */
  /* End of cleaning variable ier */
  } /*if (f2py_success) of k*/
  /* End of cleaning variable k */
  if((PyObject *)capi_t_tmp!=t_capi) {
    Py_XDECREF(capi_t_tmp); }
  }  /*if (capi_t_tmp == NULL) ... else of t*/
  /* End of cleaning variable t */
  if((PyObject *)capi_x_tmp!=x_capi) {
    Py_XDECREF(capi_x_tmp); }
  }  /*if (capi_x_tmp == NULL) ... else of x*/
  /* End of cleaning variable x */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************* end of fpchec *******************************/

/*********************************** bispeu ***********************************/
static char doc_f2py_rout_dfitpack_bispeu[] = "\
z,ier = bispeu(tx,ty,c,kx,ky,x,y)\n\nWrapper for ``bispeu``.\
\n\nParameters\n----------\n"
"tx : input rank-1 array('d') with bounds (nx)\n"
"ty : input rank-1 array('d') with bounds (ny)\n"
"c : input rank-1 array('d') with bounds ((nx-kx-1)*(ny-ky-1))\n"
"kx : input int\n"
"ky : input int\n"
"x : input rank-1 array('d') with bounds (m)\n"
"y : input rank-1 array('d') with bounds (m)\n"
"\nReturns\n-------\n"
"z : rank-1 array('d') with bounds (m)\n"
"ier : int";
/* extern void F_FUNC(bispeu,BISPEU)(double*,int*,double*,int*,double*,int*,int*,double*,double*,double*,int*,double*,int*,int*); */
static PyObject *f2py_rout_dfitpack_bispeu(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,int*,double*,int*,double*,int*,int*,double*,double*,double*,int*,double*,int*,int*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double *tx = NULL;
  npy_intp tx_Dims[1] = {-1};
  const int tx_Rank = 1;
  PyArrayObject *capi_tx_tmp = NULL;
  int capi_tx_intent = 0;
  PyObject *tx_capi = Py_None;
  int nx = 0;
  double *ty = NULL;
  npy_intp ty_Dims[1] = {-1};
  const int ty_Rank = 1;
  PyArrayObject *capi_ty_tmp = NULL;
  int capi_ty_intent = 0;
  PyObject *ty_capi = Py_None;
  int ny = 0;
  double *c = NULL;
  npy_intp c_Dims[1] = {-1};
  const int c_Rank = 1;
  PyArrayObject *capi_c_tmp = NULL;
  int capi_c_intent = 0;
  PyObject *c_capi = Py_None;
  int kx = 0;
  PyObject *kx_capi = Py_None;
  int ky = 0;
  PyObject *ky_capi = Py_None;
  double *x = NULL;
  npy_intp x_Dims[1] = {-1};
  const int x_Rank = 1;
  PyArrayObject *capi_x_tmp = NULL;
  int capi_x_intent = 0;
  PyObject *x_capi = Py_None;
  double *y = NULL;
  npy_intp y_Dims[1] = {-1};
  const int y_Rank = 1;
  PyArrayObject *capi_y_tmp = NULL;
  int capi_y_intent = 0;
  PyObject *y_capi = Py_None;
  double *z = NULL;
  npy_intp z_Dims[1] = {-1};
  const int z_Rank = 1;
  PyArrayObject *capi_z_tmp = NULL;
  int capi_z_intent = 0;
  int m = 0;
  double *wrk = NULL;
  npy_intp wrk_Dims[1] = {-1};
  const int wrk_Rank = 1;
  PyArrayObject *capi_wrk_tmp = NULL;
  int capi_wrk_intent = 0;
  int lwrk = 0;
  int ier = 0;
  static char *capi_kwlist[] = {"tx","ty","c","kx","ky","x","y",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOO:dfitpack.bispeu",\
    capi_kwlist,&tx_capi,&ty_capi,&c_capi,&kx_capi,&ky_capi,&x_capi,&y_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable tx */
  ;
  capi_tx_intent |= F2PY_INTENT_IN;
  capi_tx_tmp = array_from_pyobj(NPY_DOUBLE,tx_Dims,tx_Rank,capi_tx_intent,tx_capi);
  if (capi_tx_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 1st argument `tx' of dfitpack.bispeu to C/Fortran array" );
  } else {
    tx = (double *)(PyArray_DATA(capi_tx_tmp));

  /* Processing variable ty */
  ;
  capi_ty_intent |= F2PY_INTENT_IN;
  capi_ty_tmp = array_from_pyobj(NPY_DOUBLE,ty_Dims,ty_Rank,capi_ty_intent,ty_capi);
  if (capi_ty_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 2nd argument `ty' of dfitpack.bispeu to C/Fortran array" );
  } else {
    ty = (double *)(PyArray_DATA(capi_ty_tmp));

  /* Processing variable kx */
    f2py_success = int_from_pyobj(&kx,kx_capi,"dfitpack.bispeu() 4th argument (kx) can't be converted to int");
  if (f2py_success) {
  /* Processing variable ky */
    f2py_success = int_from_pyobj(&ky,ky_capi,"dfitpack.bispeu() 5th argument (ky) can't be converted to int");
  if (f2py_success) {
  /* Processing variable x */
  ;
  capi_x_intent |= F2PY_INTENT_IN;
  capi_x_tmp = array_from_pyobj(NPY_DOUBLE,x_Dims,x_Rank,capi_x_intent,x_capi);
  if (capi_x_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 6th argument `x' of dfitpack.bispeu to C/Fortran array" );
  } else {
    x = (double *)(PyArray_DATA(capi_x_tmp));

  /* Processing variable ier */
  /* Processing variable nx */
  nx = len(tx);
  /* Processing variable ny */
  ny = len(ty);
  /* Processing variable c */
  c_Dims[0]=(nx-kx-1)*(ny-ky-1);
  capi_c_intent |= F2PY_INTENT_IN;
  capi_c_tmp = array_from_pyobj(NPY_DOUBLE,c_Dims,c_Rank,capi_c_intent,c_capi);
  if (capi_c_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 3rd argument `c' of dfitpack.bispeu to C/Fortran array" );
  } else {
    c = (double *)(PyArray_DATA(capi_c_tmp));

  CHECKARRAY(len(c)==(nx-kx-1)*(ny-ky-1),"len(c)==(nx-kx-1)*(ny-ky-1)","3rd argument c") {
  /* Processing variable m */
  m = len(x);
  /* Processing variable z */
  z_Dims[0]=m;
  capi_z_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE|F2PY_INTENT_C;
  capi_z_tmp = array_from_pyobj(NPY_DOUBLE,z_Dims,z_Rank,capi_z_intent,Py_None);
  if (capi_z_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `z' of dfitpack.bispeu to C/Fortran array" );
  } else {
    z = (double *)(PyArray_DATA(capi_z_tmp));

  /* Processing variable lwrk */
  lwrk = kx+ky+2;
  /* Processing variable y */
  y_Dims[0]=m;
  capi_y_intent |= F2PY_INTENT_IN;
  capi_y_tmp = array_from_pyobj(NPY_DOUBLE,y_Dims,y_Rank,capi_y_intent,y_capi);
  if (capi_y_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 7th argument `y' of dfitpack.bispeu to C/Fortran array" );
  } else {
    y = (double *)(PyArray_DATA(capi_y_tmp));

  /* Processing variable wrk */
  wrk_Dims[0]=lwrk;
  capi_wrk_intent |= F2PY_INTENT_HIDE|F2PY_INTENT_CACHE;
  capi_wrk_tmp = array_from_pyobj(NPY_DOUBLE,wrk_Dims,wrk_Rank,capi_wrk_intent,Py_None);
  if (capi_wrk_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `wrk' of dfitpack.bispeu to C/Fortran array" );
  } else {
    wrk = (double *)(PyArray_DATA(capi_wrk_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
      Py_BEGIN_ALLOW_THREADS
        (*f2py_func)(tx,&nx,ty,&ny,c,&kx,&ky,x,y,z,&m,wrk,&lwrk,&ier);
      Py_END_ALLOW_THREADS
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("Ni",capi_z_tmp,ier);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    Py_XDECREF(capi_wrk_tmp);
  }  /*if (capi_wrk_tmp == NULL) ... else of wrk*/
  /* End of cleaning variable wrk */
  if((PyObject *)capi_y_tmp!=y_capi) {
    Py_XDECREF(capi_y_tmp); }
  }  /*if (capi_y_tmp == NULL) ... else of y*/
  /* End of cleaning variable y */
  /* End of cleaning variable lwrk */
  }  /*if (capi_z_tmp == NULL) ... else of z*/
  /* End of cleaning variable z */
  /* End of cleaning variable m */
  } /*CHECKARRAY(len(c)==(nx-kx-1)*(ny-ky-1))*/
  if((PyObject *)capi_c_tmp!=c_capi) {
    Py_XDECREF(capi_c_tmp); }
  }  /*if (capi_c_tmp == NULL) ... else of c*/
  /* End of cleaning variable c */
  /* End of cleaning variable ny */
  /* End of cleaning variable nx */
  /* End of cleaning variable ier */
  if((PyObject *)capi_x_tmp!=x_capi) {
    Py_XDECREF(capi_x_tmp); }
  }  /*if (capi_x_tmp == NULL) ... else of x*/
  /* End of cleaning variable x */
  } /*if (f2py_success) of ky*/
  /* End of cleaning variable ky */
  } /*if (f2py_success) of kx*/
  /* End of cleaning variable kx */
  if((PyObject *)capi_ty_tmp!=ty_capi) {
    Py_XDECREF(capi_ty_tmp); }
  }  /*if (capi_ty_tmp == NULL) ... else of ty*/
  /* End of cleaning variable ty */
  if((PyObject *)capi_tx_tmp!=tx_capi) {
    Py_XDECREF(capi_tx_tmp); }
  }  /*if (capi_tx_tmp == NULL) ... else of tx*/
  /* End of cleaning variable tx */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************* end of bispeu *******************************/

/******************************** regrid_smth ********************************/
static char doc_f2py_rout_dfitpack_regrid_smth[] = "\
nx,tx,ny,ty,c,fp,ier = regrid_smth(x,y,z,[xb,xe,yb,ye,kx,ky,s])\n\nWrapper for ``regrid_smth``.\
\n\nParameters\n----------\n"
"x : input rank-1 array('d') with bounds (mx)\n"
"y : input rank-1 array('d') with bounds (my)\n"
"z : input rank-1 array('d') with bounds (mx*my)\n"
"\nOther Parameters\n----------------\n"
"xb : input float, optional\n    Default: dmin(x,mx)\n"
"xe : input float, optional\n    Default: dmax(x,mx)\n"
"yb : input float, optional\n    Default: dmin(y,my)\n"
"ye : input float, optional\n    Default: dmax(y,my)\n"
"kx : input int, optional\n    Default: 3\n"
"ky : input int, optional\n    Default: 3\n"
"s : input float, optional\n    Default: 0.0\n"
"\nReturns\n-------\n"
"nx : int\n"
"tx : rank-1 array('d') with bounds (nxest)\n"
"ny : int\n"
"ty : rank-1 array('d') with bounds (nyest)\n"
"c : rank-1 array('d') with bounds ((nxest-kx-1)*(nyest-ky-1))\n"
"fp : float\n"
"ier : int";
/* extern void F_FUNC(regrid,REGRID)(int*,int*,double*,int*,double*,double*,double*,double*,double*,double*,int*,int*,double*,int*,int*,int*,double*,int*,double*,double*,double*,double*,int*,int*,int*,int*); */
static PyObject *f2py_rout_dfitpack_regrid_smth(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,int*,double*,int*,double*,double*,double*,double*,double*,double*,int*,int*,double*,int*,int*,int*,double*,int*,double*,double*,double*,double*,int*,int*,int*,int*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  int iopt = 0;
  int mx = 0;
  double *x = NULL;
  npy_intp x_Dims[1] = {-1};
  const int x_Rank = 1;
  PyArrayObject *capi_x_tmp = NULL;
  int capi_x_intent = 0;
  PyObject *x_capi = Py_None;
  int my = 0;
  double *y = NULL;
  npy_intp y_Dims[1] = {-1};
  const int y_Rank = 1;
  PyArrayObject *capi_y_tmp = NULL;
  int capi_y_intent = 0;
  PyObject *y_capi = Py_None;
  double *z = NULL;
  npy_intp z_Dims[1] = {-1};
  const int z_Rank = 1;
  PyArrayObject *capi_z_tmp = NULL;
  int capi_z_intent = 0;
  PyObject *z_capi = Py_None;
  double xb = 0;
  PyObject *xb_capi = Py_None;
  double xe = 0;
  PyObject *xe_capi = Py_None;
  double yb = 0;
  PyObject *yb_capi = Py_None;
  double ye = 0;
  PyObject *ye_capi = Py_None;
  int kx = 0;
  PyObject *kx_capi = Py_None;
  int ky = 0;
  PyObject *ky_capi = Py_None;
  double s = 0;
  PyObject *s_capi = Py_None;
  int nxest = 0;
  int nyest = 0;
  int nx = 0;
  double *tx = NULL;
  npy_intp tx_Dims[1] = {-1};
  const int tx_Rank = 1;
  PyArrayObject *capi_tx_tmp = NULL;
  int capi_tx_intent = 0;
  int ny = 0;
  double *ty = NULL;
  npy_intp ty_Dims[1] = {-1};
  const int ty_Rank = 1;
  PyArrayObject *capi_ty_tmp = NULL;
  int capi_ty_intent = 0;
  double *c = NULL;
  npy_intp c_Dims[1] = {-1};
  const int c_Rank = 1;
  PyArrayObject *capi_c_tmp = NULL;
  int capi_c_intent = 0;
  double fp = 0;
  double *wrk = NULL;
  npy_intp wrk_Dims[1] = {-1};
  const int wrk_Rank = 1;
  PyArrayObject *capi_wrk_tmp = NULL;
  int capi_wrk_intent = 0;
  int lwrk = 0;
  int *iwrk = NULL;
  npy_intp iwrk_Dims[1] = {-1};
  const int iwrk_Rank = 1;
  PyArrayObject *capi_iwrk_tmp = NULL;
  int capi_iwrk_intent = 0;
  int kwrk = 0;
  int ier = 0;
  static char *capi_kwlist[] = {"x","y","z","xb","xe","yb","ye","kx","ky","s",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOO|OOOOOOO:dfitpack.regrid_smth",\
    capi_kwlist,&x_capi,&y_capi,&z_capi,&xb_capi,&xe_capi,&yb_capi,&ye_capi,&kx_capi,&ky_capi,&s_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable iopt */
  iopt = 0;
  /* Processing variable x */
  ;
  capi_x_intent |= F2PY_INTENT_IN;
  capi_x_tmp = array_from_pyobj(NPY_DOUBLE,x_Dims,x_Rank,capi_x_intent,x_capi);
  if (capi_x_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 1st argument `x' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    x = (double *)(PyArray_DATA(capi_x_tmp));

  /* Processing variable y */
  ;
  capi_y_intent |= F2PY_INTENT_IN;
  capi_y_tmp = array_from_pyobj(NPY_DOUBLE,y_Dims,y_Rank,capi_y_intent,y_capi);
  if (capi_y_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 2nd argument `y' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    y = (double *)(PyArray_DATA(capi_y_tmp));

  /* Processing variable kx */
  if (kx_capi == Py_None) kx = 3; else
    f2py_success = int_from_pyobj(&kx,kx_capi,"dfitpack.regrid_smth() 5th keyword (kx) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(1<=kx && kx<=5,"1<=kx && kx<=5","5th keyword kx","regrid_smth:kx=%d",kx) {
  /* Processing variable ky */
  if (ky_capi == Py_None) ky = 3; else
    f2py_success = int_from_pyobj(&ky,ky_capi,"dfitpack.regrid_smth() 6th keyword (ky) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(1<=ky && ky<=5,"1<=ky && ky<=5","6th keyword ky","regrid_smth:ky=%d",ky) {
  /* Processing variable s */
  if (s_capi == Py_None) s = 0.0; else
    f2py_success = double_from_pyobj(&s,s_capi,"dfitpack.regrid_smth() 7th keyword (s) can't be converted to double");
  if (f2py_success) {
  CHECKSCALAR(0.0<=s,"0.0<=s","7th keyword s","regrid_smth:s=%g",s) {
  /* Processing variable nx */
  /* Processing variable ny */
  /* Processing variable fp */
  /* Processing variable ier */
  /* Processing variable mx */
  mx = len(x);
  CHECKSCALAR(mx>kx,"mx>kx","hidden mx","regrid_smth:mx=%d",mx) {
  /* Processing variable my */
  my = len(y);
  CHECKSCALAR(my>ky,"my>ky","hidden my","regrid_smth:my=%d",my) {
  /* Processing variable z */
  z_Dims[0]=mx*my;
  capi_z_intent |= F2PY_INTENT_IN;
  capi_z_tmp = array_from_pyobj(NPY_DOUBLE,z_Dims,z_Rank,capi_z_intent,z_capi);
  if (capi_z_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting 3rd argument `z' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    z = (double *)(PyArray_DATA(capi_z_tmp));

  CHECKARRAY(len(z)==mx*my,"len(z)==mx*my","3rd argument z") {
  /* Processing variable xb */
  if (xb_capi == Py_None) xb = dmin(x,mx); else
    f2py_success = double_from_pyobj(&xb,xb_capi,"dfitpack.regrid_smth() 1st keyword (xb) can't be converted to double");
  if (f2py_success) {
  /* Processing variable xe */
  if (xe_capi == Py_None) xe = dmax(x,mx); else
    f2py_success = double_from_pyobj(&xe,xe_capi,"dfitpack.regrid_smth() 2nd keyword (xe) can't be converted to double");
  if (f2py_success) {
  /* Processing variable yb */
  if (yb_capi == Py_None) yb = dmin(y,my); else
    f2py_success = double_from_pyobj(&yb,yb_capi,"dfitpack.regrid_smth() 3rd keyword (yb) can't be converted to double");
  if (f2py_success) {
  /* Processing variable ye */
  if (ye_capi == Py_None) ye = dmax(y,my); else
    f2py_success = double_from_pyobj(&ye,ye_capi,"dfitpack.regrid_smth() 4th keyword (ye) can't be converted to double");
  if (f2py_success) {
  /* Processing variable nxest */
  nxest = mx+kx+1;
  CHECKSCALAR(nxest>=2*(kx+1),"nxest>=2*(kx+1)","hidden nxest","regrid_smth:nxest=%d",nxest) {
  /* Processing variable nyest */
  nyest = my+ky+1;
  CHECKSCALAR(nyest>=2*(ky+1),"nyest>=2*(ky+1)","hidden nyest","regrid_smth:nyest=%d",nyest) {
  /* Processing variable tx */
  tx_Dims[0]=nxest;
  capi_tx_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_tx_tmp = array_from_pyobj(NPY_DOUBLE,tx_Dims,tx_Rank,capi_tx_intent,Py_None);
  if (capi_tx_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `tx' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    tx = (double *)(PyArray_DATA(capi_tx_tmp));

  /* Processing variable ty */
  ty_Dims[0]=nyest;
  capi_ty_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_ty_tmp = array_from_pyobj(NPY_DOUBLE,ty_Dims,ty_Rank,capi_ty_intent,Py_None);
  if (capi_ty_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `ty' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    ty = (double *)(PyArray_DATA(capi_ty_tmp));

  /* Processing variable c */
  c_Dims[0]=(nxest-kx-1)*(nyest-ky-1);
  capi_c_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_c_tmp = array_from_pyobj(NPY_DOUBLE,c_Dims,c_Rank,capi_c_intent,Py_None);
  if (capi_c_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `c' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    c = (double *)(PyArray_DATA(capi_c_tmp));

  /* Processing variable lwrk */
  lwrk = calc_regrid_lwrk(mx,my,kx,ky,nxest,nyest);
  /* Processing variable kwrk */
  kwrk = 3+mx+my+nxest+nyest;
  /* Processing variable wrk */
  wrk_Dims[0]=lwrk;
  capi_wrk_intent |= F2PY_INTENT_HIDE|F2PY_INTENT_CACHE;
  capi_wrk_tmp = array_from_pyobj(NPY_DOUBLE,wrk_Dims,wrk_Rank,capi_wrk_intent,Py_None);
  if (capi_wrk_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `wrk' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    wrk = (double *)(PyArray_DATA(capi_wrk_tmp));

  /* Processing variable iwrk */
  iwrk_Dims[0]=kwrk;
  capi_iwrk_intent |= F2PY_INTENT_HIDE|F2PY_INTENT_CACHE;
  capi_iwrk_tmp = array_from_pyobj(NPY_INT,iwrk_Dims,iwrk_Rank,capi_iwrk_intent,Py_None);
  if (capi_iwrk_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(dfitpack_error,"failed in converting hidden `iwrk' of dfitpack.regrid_smth to C/Fortran array" );
  } else {
    iwrk = (int *)(PyArray_DATA(capi_iwrk_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
      Py_BEGIN_ALLOW_THREADS
        (*f2py_func)(&iopt,&mx,x,&my,y,z,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nx,tx,&ny,ty,c,&fp,wrk,&lwrk,iwrk,&kwrk,&ier);
      Py_END_ALLOW_THREADS
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("iNiNNdi",nx,capi_tx_tmp,ny,capi_ty_tmp,capi_c_tmp,fp,ier);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    Py_XDECREF(capi_iwrk_tmp);
  }  /*if (capi_iwrk_tmp == NULL) ... else of iwrk*/
  /* End of cleaning variable iwrk */
    Py_XDECREF(capi_wrk_tmp);
  }  /*if (capi_wrk_tmp == NULL) ... else of wrk*/
  /* End of cleaning variable wrk */
  /* End of cleaning variable kwrk */
  /* End of cleaning variable lwrk */
  }  /*if (capi_c_tmp == NULL) ... else of c*/
  /* End of cleaning variable c */
  }  /*if (capi_ty_tmp == NULL) ... else of ty*/
  /* End of cleaning variable ty */
  }  /*if (capi_tx_tmp == NULL) ... else of tx*/
  /* End of cleaning variable tx */
  } /*CHECKSCALAR(nyest>=2*(ky+1))*/
  /* End of cleaning variable nyest */
  } /*CHECKSCALAR(nxest>=2*(kx+1))*/
  /* End of cleaning variable nxest */
  } /*if (f2py_success) of ye*/
  /* End of cleaning variable ye */
  } /*if (f2py_success) of yb*/
  /* End of cleaning variable yb */
  } /*if (f2py_success) of xe*/
  /* End of cleaning variable xe */
  } /*if (f2py_success) of xb*/
  /* End of cleaning variable xb */
  } /*CHECKARRAY(len(z)==mx*my)*/
  if((PyObject *)capi_z_tmp!=z_capi) {
    Py_XDECREF(capi_z_tmp); }
  }  /*if (capi_z_tmp == NULL) ... else of z*/
  /* End of cleaning variable z */
  } /*CHECKSCALAR(my>ky)*/
  /* End of cleaning variable my */
  } /*CHECKSCALAR(mx>kx)*/
  /* End of cleaning variable mx */
  /* End of cleaning variable ier */
  /* End of cleaning variable fp */
  /* End of cleaning variable ny */
  /* End of cleaning variable nx */
  } /*CHECKSCALAR(0.0<=s)*/
  } /*if (f2py_success) of s*/
  /* End of cleaning variable s */
  } /*CHECKSCALAR(1<=ky && ky<=5)*/
  } /*if (f2py_success) of ky*/
  /* End of cleaning variable ky */
  } /*CHECKSCALAR(1<=kx && kx<=5)*/
  } /*if (f2py_success) of kx*/
  /* End of cleaning variable kx */
  if((PyObject *)capi_y_tmp!=y_capi) {
    Py_XDECREF(capi_y_tmp); }
  }  /*if (capi_y_tmp == NULL) ... else of y*/
  /* End of cleaning variable y */
  if((PyObject *)capi_x_tmp!=x_capi) {
    Py_XDECREF(capi_x_tmp); }
  }  /*if (capi_x_tmp == NULL) ... else of x*/
  /* End of cleaning variable x */
  /* End of cleaning variable iopt */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/***************************** end of regrid_smth *****************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"fpchec",-1,{{-1}},0,(char *)F_FUNC(fpchec,FPCHEC),(f2py_init_func)f2py_rout_dfitpack_fpchec,doc_f2py_rout_dfitpack_fpchec},
  {"bispeu",-1,{{-1}},0,(char *)F_FUNC(bispeu,BISPEU),(f2py_init_func)f2py_rout_dfitpack_bispeu,doc_f2py_rout_dfitpack_bispeu},
  {"regrid_smth",-1,{{-1}},0,(char *)F_FUNC(regrid,REGRID),(f2py_init_func)f2py_rout_dfitpack_regrid_smth,doc_f2py_rout_dfitpack_regrid_smth},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "dfitpack",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif

#if PY_VERSION_HEX >= 0x03000000
#define RETVAL m
PyMODINIT_FUNC PyInit_dfitpack(void) {
#else
#define RETVAL
PyMODINIT_FUNC initdfitpack(void) {
#endif
  int i;
  PyObject *m,*d, *s;
#if PY_VERSION_HEX >= 0x03000000
  m = dfitpack_module = PyModule_Create(&moduledef);
#else
  m = dfitpack_module = Py_InitModule("dfitpack", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module dfitpack (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module 'dfitpack' is auto-generated with f2py (version:2).\nFunctions:\n"
"  ier = fpchec(x,t,k)\n"
"  z,ier = bispeu(tx,ty,c,kx,ky,x,y)\n"
"  nx,tx,ny,ty,c,fp,ier = regrid_smth(x,y,z,xb=dmin(x,mx),xe=dmax(x,mx),yb=dmin(y,my),ye=dmax(y,my),kx=3,ky=3,s=0.0)\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  dfitpack_error = PyErr_NewException ("dfitpack.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));



/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"dfitpack");
#endif

  return RETVAL;
}
#ifdef __cplusplus
}
#endif
