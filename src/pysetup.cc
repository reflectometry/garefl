/* This program is public domain. */

/* Implementation and wrapping of example C fit functions */

#include <Python.h>
/* Support for pre-2.5 python versions */
#if (PY_VERSION_HEX < 0x02050000)
typedef int Py_ssize_t;
#endif

extern "C" {
#include "refl.h"
};

#include <stdio.h>
#ifdef HAVE_DLOPEN
#include <dlfcn.h>
#include<iostream>
#endif

#define VECTOR(obj,buf,len) \
    do { \
        if (PyObject_AsWriteBuffer(obj, (void **)&buf, &len) == 0) \
          len /= sizeof(*buf); \
        else \
          len = -1; \
    } while (0)

double Inf = 0.;

/* Application control parameters */
// FIXME: shouldn't be using globals here
int MODELS;
int approximate_roughness = 0;
int weighted = 1;
fitinfo *fit;
fit_constraints *constraints = NULL;

/* Copy parameters between models. */
void tied_parameters(fitinfo fit[])
{
  int i,k;

  /* Rescue the free parameters from the model. */
  for (i=0; i < fit[1].pars.n; i++)
    fit[1].pars.value[i] = *(fit[1].pars.address[i]);

  /* Go through all layers copying parameters from model 0
   * to the other models. This is more work than we strictly
   * need to do (we only really need to update the parameters
   * which are varying in the fit), but the extra cost is
   * trivial compared to the cost of calculating the profile
   * and the reflectivity, so let the code be simpler.
   */
  for (k=0; k < fit[0].m.n; k++) {
    for (i=1; i < MODELS; i++) {
      fit[i].m.d[k] = fit[0].m.d[k];
      fit[i].m.rho[k] = fit[0].m.rho[k];
      fit[i].m.mu[k] = fit[0].m.mu[k];
      fit[i].m.rough[k] = fit[0].m.rough[k];
#ifdef HAVE_MAGNETIC
      fit[i].m.theta[k] = fit[0].m.theta[k];
      fit[i].m.thetarough[k] = fit[0].m.thetarough[k];
      fit[i].m.P[k] = fit[0].m.P[k];
      fit[i].m.Prough[k] = fit[0].m.Prough[k];
#endif
    }
  }

  /* Restore the free parameters to the model. */
  for (i=0; i < fit[1].pars.n; i++)
    *(fit[1].pars.address[i]) = fit[1].pars.value[i];
}

void save_models(fitinfo *fit)
{
  char name[100];
  int i;

  for (i=0; i < MODELS; i++) {
    sprintf(name,"model%d.dat",i);
    model_print(&fit[i].m, name);
    sprintf(name,"profile%d.dat",i);
    profile_print(&fit[i].p, name);
    sprintf(name,"fit%d.dat",i);
    /*{int i; for (i=0;i<fit.nQ;i++) printf("Q[%d]: %g\n",i,fit.fitQ[i]);}*/
    if (fit[i].datatype == FIT_POLARIZED) {
      int n=strlen(name);
      name[n+1]='\0';
      name[n]='A'; 
      data_printfit(name, &fit[i].dataA, fit[i].fitA);
      name[n]='B'; 
      data_printfit(name, &fit[i].dataB, fit[i].fitB);
      name[n]='C'; 
      data_printfit(name, &fit[i].dataC, fit[i].fitC);
      name[n]='D'; 
      data_printfit(name, &fit[i].dataD, fit[i].fitD);
    } else {
      data_printfit(name, &fit[i].dataA, fit[i].fitA);
    }
  }
}


double update_models(fitinfo *fit)
{
  int n = 0;
  double sumsq = 0.;
  int i;

  if (*constraints) (*constraints)(fit);
  for (i=0; i < MODELS; i++) {
    int n_i = 0;
    double sumsq_i = 0.;
    fit_update(&fit[i], approximate_roughness);
    if (weighted) fit_wsumsq(&fit[i],&n_i,&sumsq_i);
    else fit_sumsq(&fit[i],&n_i,&sumsq_i);
    fit[i].chisq_est = sumsq_i/n_i;
    n += n_i;
    sumsq += sumsq_i;
  }		   
  /* printf("sumsq=%10g, n=%4d, pars=%d\n",sumsq,n,fit[0].pars.n); */
  return n < fit[0].pars.n ? sumsq : sumsq / (n - fit[0].pars.n) ;
}

double step_fn(int n, const double p[], void *user_data)
{
  fitinfo* fit = (fitinfo*)user_data;
  double chisq;

  pars_set(&fit[0].pars, p);
  chisq = update_models(fit);
#if 0
  printf("function step");
  pars_print(&fit[0].pars);
  printf("chisq=%g\n",chisq);
#endif
  return chisq;
}


typedef fitinfo* (*setupfn)(int *pmodels, void**);
static char reflmodel_doc[]=
"Reflectometry objective function";
static PyObject* py_reflmodel(PyObject*obj,PyObject*args)
{
  const char *modelfile;
  int ok = PyArg_ParseTuple(args, "s:reflmodel.reflmodel",&modelfile);
  if (!ok) return NULL;

#ifdef HAVE_DLOPEN
  void *libhandle = dlopen (modelfile, RTLD_LAZY);
  if (libhandle == NULL) {
    const char *err1 = dlerror();
    std::cout << err1 << std::endl;
    perror(modelfile);
    PyErr_Format(PyExc_ValueError, "Unable to open library %s", modelfile);
    return NULL;
  }
  setupfn fn = (setupfn)dlsym(libhandle, "setup_models");
  const char *err = dlerror();
  if (err != NULL) {
    PyErr_Format(PyExc_ValueError, "Unable to load setup_models: %s",err);
    return NULL;
  }
  fit = (*fn)(&MODELS, (void**)&constraints);
  if (fit == NULL) {
    PyErr_Format(PyExc_ValueError, "Unable to run setup_models: %s",err);
    return NULL;
  }

  return Py_BuildValue("OO",
		       PyCObject_FromVoidPtr((void *)step_fn, NULL),
		       PyCObject_FromVoidPtr((void *)fit, NULL));
//#else
//  PyErr_SetString(PyExc_NotImplementedError,
  //                 "This platform doesn't support dlopen");
  //  return NULL;
#else
  fit = (setup_models)(&MODELS)
#endif
}

static char savemodel_doc[]=
"Save reflectometry model parameters";
static PyObject* py_savemodel(PyObject*obj,PyObject*args)
{
  PyObject* pobj;
  Py_ssize_t np;
  double* p;

  int ok = PyArg_ParseTuple(args, "O:reflmodel.savemodel", &pobj);
  if (!ok) return NULL;
  
  VECTOR(pobj, p, np);

  if (PyErr_Occurred()) return NULL;
  step_fn(np, p, fit);
  save_models(fit);

  return Py_BuildValue("");
}

static char numpars_doc[]=
"Number of parameters to the model";
static PyObject* py_numpars(PyObject*obj,PyObject*args)
{
  int np = fit[0].pars.n;

  PyArg_ParseTuple(args, ":reflmodel.numpars");
  return Py_BuildValue("i", np);

}
static char get_parname_doc[]=
"returns paramater name at given index";
static PyObject* py_get_parname(PyObject* obj, PyObject* args)
{
  const char* name;
  int i;
  if (!PyArg_ParseTuple(args, "i:reflmodel.get_parname", &i))
    return NULL;
  
  
  name = pars_name(&fit[0].pars,i);
  return Py_BuildValue("s", name);

}

static char save_staj_doc[]=
"save fit as .sta file";
static PyObject* py_save_staj(PyObject* obj, PyObject* args)
{
  char* filename_base;
  if(!PyArg_ParseTuple(args, "s:reflmodel.save_staj", &filename_base))
     return NULL;

  for(int i = 0; i < MODELS; i++)
    {
      char filename[FILENAME_MAX];
      sprintf(filename, "%s%d%s", filename_base, i, ".sta");
      fit_save_staj(&fit[i], filename);
    }
  
  return Py_BuildValue("");
}

static char getpars_doc[]=
"Returns model parameters";
static PyObject* py_getpars(PyObject*obj,PyObject*args)
{
  PyObject *pobj,*loobj,*hiobj;
  Py_ssize_t np,nlo,nhi;
  double *p, *lo, *hi;
  if (!PyArg_ParseTuple(args, "OOO:reflmodel.savemodel", &pobj, &loobj, &hiobj))
    return NULL;
  
  VECTOR(pobj, p, np);
  VECTOR(loobj, lo, nlo);
  VECTOR(hiobj, hi, nhi);

  if (PyErr_Occurred()) return NULL;


  if (!(fit[0].pars.n == np && np == nlo && nlo == nhi))
    return NULL;

  for (int i = 0; i < fit[0].pars.n; i++)
    {
      p[i] = pars_peek(&fit[0].pars, i);
      lo[i] = pars_min(&fit[0].pars, i);
      hi[i] = pars_max(&fit[0].pars, i);
    }
  
  return Py_BuildValue("");
}



static PyMethodDef methods[] = {
    {"reflmodel", py_reflmodel, METH_VARARGS, reflmodel_doc},
    {"savemodel", py_savemodel, METH_VARARGS, savemodel_doc},
    {"numpars", py_numpars, METH_VARARGS, numpars_doc},
    {"getpars", py_getpars, METH_VARARGS, getpars_doc},
    {"get_parname", py_get_parname, METH_VARARGS, get_parname_doc},
    {"save_staj", py_save_staj, METH_VARARGS, save_staj_doc},
    {0}
};


#if defined(WIN32) && !defined(__MINGW32__)
__declspec(dllexport)
#endif
extern "C" void initreflmodel(void) 
{
  Inf = 1./Inf;
  Py_InitModule3("reflmodel",methods,"Reflectometry fitting objective function");
  
}

/* $Id: _cfitfn_example.cc 11 2006-03-30 08:44:57Z pkienzle $ */
