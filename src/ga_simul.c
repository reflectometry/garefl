/* This program is public domain */

#include "reflconfig.h"
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <signal.h>
#include <assert.h>
#include <math.h>
#ifdef HAVE_SETRLIMIT
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include "refl.h"
#include "ga.h"
#include "amoeba.h"
#undef USE_QUAD_FIT
#ifdef USE_QUAD_FIT
#include "newuoa/newuoa.h"
#endif

#if defined(__WIN32__) && !defined(_POSIX_VERSION)
#define USE_WIN32_SIGINT 1
#endif

#define USE_NLLS_FIT /* True if user should have levenberg-marquardt option */

/* CVS information */
static char cvsid[] = "$Id$";

/* Fake fortran main required to keep some linkers happy. */
#ifdef F77_DUMMY_MAIN
# ifdef __cplusplus
extern "C"
# endif
int F77_DUMMY_MAIN() { return 1; }
#endif


/* Signal handler */
int halt = 0;

#ifdef USE_WIN32_SIGINT
#include <windows.h>
static BOOL CALLBACK
win32_sigint_handler (DWORD sig)
{
  switch(sig) {
    case CTRL_C_EVENT: halt = 1; break;
    case CTRL_BREAK_EVENT:   
    case CTRL_CLOSE_EVENT:
    case CTRL_LOGOFF_EVENT:
    case CTRL_SHUTDOWN_EVENT: halt = -1; break;
  }

  // Return TRUE if the event was handled, or FALSE if another handler 
  // should be called.
  return TRUE;
}

void set_signal_handlers(void)
{
  // Intercept windows console control events.
  // Note that the windows console signal handlers chain, so if 
  // install_signal_handlers is called more than once in the same program,
  // then first call the following to avoid duplicates:
  //     SetConsoleCtrlHandler(win32_sigint_handler, FALSE);
  if (!SetConsoleCtrlHandler(win32_sigint_handler, TRUE)) {
    fprintf(stderr,"SetConsoleCtrlHandler failed with %ld\n", GetLastError());
  }
 
  // Let the user close the console window or shutdown without the
  // pesky dialog. XXX FIXME XXX should this be user configurable?
  SetProcessShutdownParameters(0x280, SHUTDOWN_NORETRY);
#ifdef HAVE_NICE
  nice(10);
#endif
}

#else /* !USE_WIN32_SIGINT */

void  INThandler(int sig) { 
 /* Use halt=1 rather than halt_ga() directly so that
  * Ctrl-C is responsive during a randomize operation.
  */
  if (halt!=0) {
      fprintf(stderr,"You're insisting...\n");
      exit(1);
  }
  halt=1; 
  signal(sig, INThandler);
}
void  TERMhandler(int sig) { 
  /* Stop within ten minutes of being told to. */
  if (halt<0) --halt;
  else halt=-1;
  if (halt < -600) {
    fprintf(stderr,"Aborting...\n");
    exit(1);
  }
  signal(sig, TERMhandler);
}

void set_signal_handlers(void)
{
#ifdef HAVE_NICE
  nice(10);
#endif
#ifdef SIGXCPU
  signal(SIGXCPU, TERMhandler);
#endif
  signal(SIGTERM, TERMhandler);
  signal(SIGINT, INThandler);
}
#endif /* !USE_WIN32_SIGINT */

void cpulimit(double hours)
{
#ifdef HAVE_SETRLIMIT
  struct rlimit r;
  r.rlim_cur = (int)(hours*3600.);
  r.rlim_max = RLIM_INFINITY;
  setrlimit(RLIMIT_CPU,&r);
#endif
}




/* Application control parameters */
int MODELS;
double CHI_CUT = 500.0;
double bestchi=1e308;
double *bestpars=NULL;
int weighted = 1;
int approximate_roughness = 0;
double portion = 1.;
char  parFile[64];
int log_improvement = 1;
FILE *parFD;
FILE *searchFD;
fitinfo *fit;
Settings set; /* GA fit parameters */
fit_constraints *constraints = NULL;
int step_num = 0;
int start_step = 0;
double trust_region = 0.1;

void tic(void)
{
  start_step = step_num;
}

void toc(void)
{
  printf("#steps = %d\n",step_num-start_step);
}

/* GA is a minimizer, so convert chisq to GA fitness value. */
double chisq2fitness(double chisq)
{
  double fitness;
  
  /* ga maximizes so invert chisq */
  fitness = CHI_CUT*MODELS - chisq;

  /* Squeeze poor models into the range 0-1 */
  if (fitness < 1.) fitness = 1.+2.*atan(fitness-1.)/M_PI;

  return fitness;	   
}

double fitness2chisq(double fitness)
{
  double chisq;
  
  /* Unsqueeze poor models from the range 0-1 */
  if (fitness < 1.) fitness = tan((fitness + 1.)*M_PI/2.) + 1.;
  
  chisq = CHI_CUT*MODELS - fitness;

  return chisq;
}

void write_ga_pop(const char filename[])
{
  FILE *file = stdout;
  fitpars *pars = &fit[0].pars;
  int p,i;

  if (filename != NULL) {
  	file = fopen(filename, "w");
  	if (file == NULL) return;
  }
  fprintf(file,"    %25s: ","fitness");
  for (i=0; i < set.np; i++) {
  	double chisq = fitness2chisq(set.pop->indiv[i].fitness);
  	fprintf(file,"%15g ",chisq);
  }
  fprintf(file,"\n");
  for (p=0; p < set.nParams; p++) {
  	fprintf(file,"%3d %25s: ",p, (pars->name[p]?pars->name[p]:""));
	for (i=0; i < set.np; i++) {
	  pars_set01(pars, set.pop->indiv[i].value);
	  if (*constraints) (*constraints)(fit);
	  fprintf(file,"%15g ",pars_peek(pars,p));
	}
	fprintf(file,"\n");
  }
  if (filename != NULL) fclose(file);
}


void check_halt(void);
void record_improvement(void);
double do_step(fitinfo *fit, double portion);
double update_models(fitinfo *fit);
double partial_update_models(fitinfo *fit,double portion,double best);

static void chisq_plot(fitinfo *fit, int nth, double portion, int steps)
{
  FILE *file;
  file = fopen("chisq.dat","w");
  if (file == NULL) {
    printf("Could not open file chisq.dat\n");
    return;
  }

  /* chisq plot for the nth parameter */
  if (nth >= 0 && nth < fit[0].pars.n) {
    int i,k;
    double d = pars_min(&fit[0].pars,nth);
    double step = (pars_max(&fit[0].pars,nth)-d)/steps;

    /* Print headers */
    printf("Chisq plot for %s [%i]\n%12s %12s",
	   pars_name(&fit[0].pars,nth), nth,"value","chisq");
    fprintf(file, "# Chisq plot for %s [%i]\n# %10s %12s",
	    pars_name(&fit[0].pars,nth), nth,"value","chisq");
    if (MODELS>1) for (k=0; k < MODELS; k++) {
      printf(" %8s %3d","model",k);
      fprintf(file," %8s %3d","model",k);
    }
    printf("\n");
    fprintf(file,"\n");	

    /* Print data */
    for (i=0; i <= steps; i++) {
      double chisq;

      /* Set parameter in model */
      pars_poke(&fit[0].pars,nth,d);

      /* Compute full or partial chisq, depending on portion */
      if (portion < 1.)	chisq = partial_update_models(fit,portion,1e308);
      else chisq = update_models(fit);

      /* Report chisq and sumsq/df for data set */
      printf("%12.4g %12.2f",d,chisq);
      fprintf(file,"%.15g %.15g",d,chisq);
      if (MODELS>1) for (k=0; k < MODELS; k++) {
	printf(" %12.2f",fit[k].chisq_est);
	fprintf(file," %.15g",fit[k].chisq_est);
      }
      printf("\n");
      fprintf(file,"\n");	

      /* Move parameter. */
      d += step;
    }
  } else {
    pars_print(&fit[0].pars);
  }
  fclose(file);
}

static void chisq_surf(fitinfo *fit, int n1, int n2, double portion, int steps)
{
  FILE *file;
  file = fopen("chisqsurf.dat","w");
  if (file == NULL) {
    printf("Could not open file chisqsurf.dat\n");
    return;
  }

  /* chisq plot for the nth parameter */
  if (n1 < 0 || n1 >= fit[0].pars.n || n2 < 0 || n2 >= fit[0].pars.n) {
    pars_print(&fit[0].pars);
  } else {
    int i,j,k;
    double d1 = pars_min(&fit[0].pars,n1);
    double d2 = pars_min(&fit[0].pars,n2);
    double step1 = (pars_max(&fit[0].pars,n1)-d1)/steps;
    double step2 = (pars_max(&fit[0].pars,n2)-d2)/steps;


    /* Print headers */
    printf("Chisq plot for %s x %s [%dx%d]\n%12s %12s",
	   pars_name(&fit[0].pars,n1), pars_name(&fit[0].pars,n2),
	   n1,n2,"value","chisq");
    fprintf(file, "# Chisq plot for %s x %s [%dx%d]\n# %10s %12s",
	   pars_name(&fit[0].pars,n1), pars_name(&fit[0].pars,n2),
	   n1,n2,"value","chisq");
    if (MODELS>1) for (k=0; k < MODELS; k++) {
      printf(" %8s %3d","model",k);
      fprintf(file," %8s %3d","model",k);
    }
    printf("\n");
    fprintf(file,"\n");	

    /* Print data */
    for (i=0; i <= steps; i++) {
      double d = d2;
      for (j=0; j <= steps; j++) {
        double chisq;

        /* Set parameter in model */
        pars_poke(&fit[0].pars,n1,d1);
        pars_poke(&fit[0].pars,n2,d);

        /* Compute full or partial chisq, depending on portion */
        if (portion < 1.) chisq = partial_update_models(fit,portion,1e308);
        else chisq = update_models(fit);

        /* Report chisq and sumsq/df for data set */
        printf("%12.4g %12.4g %12.2f",d1,d,chisq);
        fprintf(file,"%.15g %.15g %.15g",d1,d,chisq);
        if (MODELS>1) for (k=0; k < MODELS; k++) {
  	  printf(" %12.2f",fit[k].chisq_est);
	  fprintf(file," %.15g",fit[k].chisq_est);
        }
        printf("\n");
        fprintf(file,"\n");	

        /* Move parameter. */
	d += step2;
      }
      printf("\n");
      fprintf(file,"\n");	
      d1 += step1;
    }
  }
  fclose(file);
}

double step_fn(int n, const double p[], void *user_data)
{
  fitinfo* fit = (fitinfo*)user_data;
  double chisq;

  check_halt();
  pars_set(&fit[0].pars, p);
  chisq = do_step(fit,portion);
#if 0
  printf("function step");
  pars_print(&fit[0].pars);
  printf("chisq=%g\n",chisq);
#endif
  return chisq;
}


double step_fn01(const int n, const double p[], void *user_data)
{
  fitinfo* fit = (fitinfo*)user_data;
  double chisq;

  check_halt();
  pars_set01(&fit[0].pars, p);
  chisq = do_step(fit,portion);
#if 0
  printf("function step 01");
  pars_print(&fit[0].pars);
  printf("chisq=%g\n",chisq);
#endif
  return chisq;
}
 
double fortran_step_fn01(const int *n, const double p[], void *user_data)
{
  return step_fn01(*n, p, user_data);
}

#ifdef USE_NLLS_FIT
#define INCLUDING_CMD_NLLS
#include "cmd_nlls.c"
#endif

#ifdef USE_QUAD_FIT
void cmd_quadfit(void)
{
  const int n = fit[0].pars.n;
  const int npnt = 2*n;
  double *p = (double *)malloc(sizeof(double)*(NEWUOA_WORKSIZE(n,npnt)+n));
  double *work = p+n;
  double rhobeg,rhoend;
  int print, maxiter;

  assert(work != NULL);

  /* Find starting point in 0-1 box */
  pars_set(&fit[0].pars, bestpars);
  pars_get01(&fit[0].pars, p);

  rhobeg = trust_region;
  rhoend = rhobeg*1e-4;
  print = 1;
  maxiter = 500;
  tic();
  newuoa(fortran_step_fn01, &n, &npnt, p, &rhobeg, &rhoend, 
	 &print, &maxiter, work, fit);
  toc();

  /* Add the best back to population; */
  setChromosome(&set, 0, p);
  free(p);
}
#endif

void cmd_amoeba(void)
{
  simplex s;
  double *pbest;
  int ndim = fit[0].pars.n;
  double *work = (double *)malloc((2*ndim+amoeba_worksize(ndim))*sizeof(double));
  double *bounds = work+amoeba_worksize(ndim);
  double *po;
  int i,j;

  write_pop_backup(&set); /* In case fit crashes */
  /* Record what we are doing */
  if (parFD != NULL) {
    fprintf(parFD,"# %15d   Starting amoeba\n", GetGen(&set));
    fflush(parFD); 
  }

  /* Suppress logging during amoeba */
  log_improvement = 0;

  /* Make sure memory allocation was successful */
  assert(work != NULL);

  /* Bounded in a [0,1] box */
  for (i=0; i < ndim; i++) {
    bounds[i] = 0.;
    bounds[i+ndim] = 1.;
  }

  /* Initialize amoeba structure */
  amoeba_init(&s, ndim, bounds, work, step_fn01, fit);

  /* Start fit from the current population best, in 0-1 notation */
  po = amoeba_VERTEX(&s,0);
  getChromosome(&set,fittest(&set),po);
  amoeba_VALUE(&s,0) = update_models(fit);

  /* Random restart */
  /* Note: consider using a normal ball instead of a complete space search */
  for (i=1; i <= ndim; i++) {
    double *pi = amoeba_VERTEX(&s,i);
    for (j=0; j < ndim; j++) pi[j] = fmod(po[j]+frandom()*trust_region, 1.);
    pars_set01(&fit[0].pars,pi);
    amoeba_VALUE(&s,i) = update_models(fit);
  }

#if 0
  printf("Initial simplex\n");
  for (i=0; i <=ndim; i++) {
    pars_set01(&fit[0].pars,amoeba_VERTEX(&s,i));
    pars_print(&fit[0].pars);
    printf("chisq=%g\n",y[i]);
  }
#endif

  /* Pull toward the nearest minimum */
  tic(); 
  pbest = amoeba(&s,1e-8,500);
  toc();

  /* Restore logging, recording amoeba results */
  log_improvement = 1;
  pars_set(&fit[0].pars,bestpars);
  update_models(fit);
  record_improvement();

  /* Add the best back to population; */
  setChromosome(&set, 0, pbest);
  free(work);
}

void save_all_staj(void)
{
  char filename[20];
  int i;

  for (i=0; i < MODELS; i++) {
    if (fit[i].datatype == FIT_POLARIZED) {
      sprintf(filename,"reflfit%d.sta",i);
    } else {
      sprintf(filename,"reflfit%d.staj",i);
    }
    fit_save_staj(&fit[i],filename);
  }
}

void cmd_save_staj(void)
{
  /* Use best parameter set */
  pars_set(&fit[0].pars, bestpars);
  /* Apply all constraints */
  if (*constraints) (*constraints)(fit);
  /* Print to file */
  save_all_staj();
}


void cmd_accelerate(void)
{
  char buffer[100];
  printf("Enter new portion (current is %g): ",portion*100.);
  fflush(stdout);
  if (fgets(buffer,sizeof(buffer),stdin)!=NULL) {
    int p;
    if (sscanf(buffer,"%d",&p)>0) {
      if (p<1) p = 1;
      else if (p>100) p = 100;
      if (floor(portion*100.+0.5) != p) {
	if (parFD != NULL)
	  fprintf(parFD,"# %15d   Changing portion to %d%%\n", GetGen(&set),p);
	portion = p/100.;
      }
    }
  }
}

void cmd_change_range(void)
{
  char buffer[100];
  const char *name;
  double pmin, pmax;
  int ipar = -1;

  ipar = pars_select(&fit[0].pars);
  if (ipar < 0) return;
  name = pars_name(&fit[0].pars,ipar);
  pmin = pars_min(&fit[0].pars,ipar);
  pmax = pars_max(&fit[0].pars,ipar);

#if 0
  /* Print the value for this parameter from the best set. */
  printf(" Fittest %17s [%3i]= %12g in [%g,%g]\n", name, ipar, bestpars[ipar]);
#endif

  printf("Enter new range (current is [%g %g]) for par %i (%s): ",
	 pmin, pmax, ipar, name);
  fflush(stdout);
  if (fgets(buffer,sizeof(buffer),stdin)!=NULL) {
    float minval, maxval;
    if (sscanf(buffer,"%f %f",&minval,&maxval)==2) { 
      printf("Setting range [%g %g] for par %i (%s)\n", 
	     minval, maxval, ipar, name);
      if (parFD != NULL)
	fprintf(parFD,"# %15d   Setting range [%g %g] for par %i (%s)\n", 
		GetGen(&set), minval, maxval, ipar, name);
      // First we take care of the population
      setRange( &set, ipar, pmin, pmax-pmin, minval, maxval-minval );
      // Then we fix the parameters in memory
      pars_set_range(&fit[0].pars, ipar, minval, maxval);
    } else {
      printf("Range must be entered as a low-high pair\n");
    }
  }
}

void cmd_change_parameter(void)
{
  int ipar;
  double val;

  bestchi = 1e308;
  pars_constrained(fit[0].pars.n, bestpars, &fit[0].pars);
  ipar = pars_select(&fit[0].pars);
  if (ipar >= 0) ipar = pars_enter_value(&fit[0].pars,ipar,&val);
  if (ipar >= 0) {
    const char *name = pars_name(&fit[0].pars,ipar);
    printf("Setting par %i (%s) to %g\n", ipar, name, val);
    if (parFD != NULL)
      fprintf(parFD,"# %15d   Setting par %i (%s) to %g\n", 
	      GetGen(&set), ipar, name, val);
    setParValue(&set,ipar, pars_to_01(&fit[0].pars,ipar,val));
  }
}

void cmd_chisq_surf(void)
{
  int ipar=-1, jpar=-1;
  printf("1 ");
  ipar = pars_select(&fit[0].pars);
  if (ipar < 0) return;
  printf("2 ");
  jpar = pars_select(&fit[0].pars);
  if (jpar < 0) return;
  pars_set(&fit[0].pars,bestpars);
  chisq_surf(fit, ipar, jpar, 1., 20);
}

void cmd_chisq_plot(void)
{
  int ipar = -1;
  ipar = pars_select(&fit[0].pars);
  if (ipar >= 0) {
    pars_set(&fit[0].pars,bestpars);
    chisq_plot(fit, ipar, 1., 20);
    printf("Current best = %g, chi^2 = %g]\n",
	   pars_peek(&fit[0].pars,ipar), bestchi);
 } 
}

void cmd_print_best(void)
{
  int i;
  double chisq;

  /* Compute partial chisq of best */
  pars_set(&fit[0].pars, bestpars);
  chisq = do_step(fit,1.); /* Force full chisq calc. */

  /* Print parameters. */
  printf("Generation %d: Best parameters\n",GetGen(&set));
  pars_print(&fit[0].pars);
  printf("Chi^2 = %g\n",chisq);
  for (i=0; i<MODELS; i++)
    printf("Model %d: df=%d, sumsq/df=%g\n",i,fit[i].nQ,fit[i].chisq_est);
}

void cmd_randomize(void)
{
  int ipar = -1;

  bestchi = 1e308;
  ipar = pars_select(&fit[0].pars);
  printf("Randomizing (%i)\n", ipar);
  if (parFD != NULL)
    fprintf(parFD,"# %15d   Randomizing %i\n", GetGen(&set), ipar);
  randomizePop(&set,ipar);
}

void cmd_approximate(void)
{
  approximate_roughness = !approximate_roughness;
  printf("Setting approximate roughness to %d\n", approximate_roughness);
  if (parFD != NULL)
    fprintf(parFD,"# %15d   Setting approximate roughness to %d\n",
	    GetGen(&set), approximate_roughness);
}

void cmd_quit(void)
{
  /* exit */
  printf("Bye bye\n");
  if (parFD != NULL) fclose(parFD);
  if (searchFD != NULL) fclose(searchFD);
  exit(0);
}


void cmd_write_population(void)
{
  write_pop(&set);
}

void cmd_change_trust_region(void)
{
  char buffer[100];
  printf("Enter new trust region in range (0,1] (current is %g): ",
	 trust_region);
  fflush(stdout);
  if (fgets(buffer,sizeof(buffer),stdin)!=NULL) {
    if (sscanf(buffer,"%lg",&trust_region)>0) {
      if (trust_region<1e-5) trust_region = 1e-5;
      else if (trust_region>1.) trust_region = 1.;
      if (parFD != NULL)
	fprintf(parFD,"# %15d   Changing trust region power to %g\n", 
		GetGen(&set),trust_region);
    }
  }
}

void halt_ga (void) {
  char buffer[100];
  time_t now;
  now = time(NULL);
  printf("Generation number %i  -- %s\n",GetGen(&set),ctime(&now));
  printf("Select an option:\n");
  printf("   q (quit)\n");
  printf("   w (write population)\n");
  printf("   b (write and quit)\n");
  printf("   a (run amoeba from current best)\n");
  printf("   t (change trust region for amoeba)\n");
#ifdef USE_NLLS_FIT
  printf("   l (run Levenberg-Marquardt from current best)\n");
#endif
#ifdef USE_QUAD_FIT
  printf("   Q (run Powell's NEWUOA from current best)\n");
#endif
  printf("   A (accelerate)\n");
  printf("   S (approximate roughness)\n");
  printf("   r (randomize)\n");
  printf("   x (plot chi^2)\n");
  printf("   X (plot chi^2 surface)\n");
  printf("   c (change a parameter)\n");
  printf("   p (print current best values)\n");
  printf("   R (Change range for a parameter)\n");
  printf("   j (write staj files for current best values)\n");
  printf("   any other character to continue\n");
  if (fgets(buffer,sizeof(buffer),stdin) == NULL) return;
  switch (buffer[0]) {
  case 'b': cmd_write_population(); /* Fall through to 'q' */
  case 'q': cmd_quit();             break;
  case 'w': cmd_write_population(); break;
  case 'A': cmd_accelerate();       break;
  case 'S': cmd_approximate();      break;
  case 'R': cmd_change_range();     break;
    /* Notice that the pop.dat will be saved with the NEW range,
     * which might cause problems if the program is stopped and restarted */
  case 'p': cmd_print_best();       break;
  case 'x': cmd_chisq_plot();       break;
  case 'X': cmd_chisq_surf();       break;
  case 'c': cmd_change_parameter(); break;
  case 't': cmd_change_trust_region(); break;
  case 'r': cmd_randomize();        break;
#ifdef USE_NLLS_FIT
  case 'l': cmd_nlls();             break;
#endif
#ifdef USE_QUAD_FIT
  case 'Q': cmd_quadfit();          break;
#endif
  case 'a': cmd_amoeba();           break;
  case 'j': cmd_save_staj();        break;
  }
  printf("Continuing...\n");
  fflush(stdout);
}

/* Check for Ctrl-C or other interupts */
void check_halt(void) 
{
  if (!halt) return;

  if (halt < 0) {
    printf("Halting with halt<0\n"); fflush(stdout);
    write_pop_backup(&set);
    if (parFD) fclose(parFD);
    exit(0);
  } else {
    printf("Halting with halt>0\n"); fflush(stdout);
    halt = 0;
    halt_ga();
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

double partial_update_models(fitinfo *fit, double portion, double best)
{
  int n = 0;
  double sumsq = 0.;
  int i;
  
  if (*constraints) (*constraints)(fit);
  fit_partial(&fit[0],approximate_roughness,portion,best,weighted,&n,&sumsq);
  for (i=1; i < MODELS; i++) {
    fit_partial(&fit[i],approximate_roughness,portion,best,weighted,&n,&sumsq);
    /* if (sumsq/n > 10.*best) break; */
  }
  return n < fit[0].pars.n ? sumsq : sumsq / (n - fit[0].pars.n) ;  
}

void save_models(fitinfo *fit)
{
  char name[100];
  int i, k;
  int incoh = MODELS; /* Number the incoherent models starting at N */

  for (i=0; i < MODELS; i++) {
  	/* Save the model and profile */
    sprintf(name,"model%d.dat",i);
    model_print(&fit[i].m, name);
    sprintf(name,"profile%d.dat",i);
    profile_print(&fit[i].p, name);

	  /* Save the supplementary incoherent models and profiles. */
    for (k=0; k < fit[i].number_incoherent; k++) {
    	profile p;
    	
    	/* Save the model */
    	sprintf(name,"model%d.dat",incoh);
    	model_print(fit[i].incoherent_models[k], name);
    	
    	/* Need to regenerate the profile since it was tossed. */
    	profile_init(&p);
    	model_profile(fit[i].incoherent_models[k], &p);
    	sprintf(name,"profile%d.dat",incoh);
    	profile_print(&p, name);
    	profile_destroy(&p);
    	
    	incoh++; /* Next incoherent profile goes to next available number */
    }
    
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

/* Results improved; record it */
void record_improvement(void)
{
  int i;
  time_t now;

  /* Print to screen */
  now = time(NULL);
  printf("Generation %d: best chisq=%g  --  %s\n", 
	GetGen(&set), 
	bestchi, ctime(&now));

  if (!log_improvement) return;

  pars_print(&fit[0].pars);
  printf("\n");
  fflush(stdout);

  /* Print to file */
  if (parFD != NULL) {
    fprintf(parFD,"%17d %17g", GetGen(&set), bestchi);
    for (i=0; i < fit->pars.n; i++) 
      fprintf(parFD," %17g",fit[0].pars.value[i]);
    fprintf(parFD,"\n");
    fflush(parFD);
  }

  /* Print plots */
  save_models(fit);
}

/* Evaluate model with current parameters */
double do_step(fitinfo *fit, double portion)
{
  int i;
  double chisq;

  step_num++;

  /* Compute reflectivity and calculate chisq */
  if (portion < 0.6)  {
    chisq = partial_update_models(fit, portion, bestchi);
    if (chisq < bestchi) chisq = update_models(fit);
  } else {
    chisq = update_models(fit);
  }
  

  /* Log the parameter set being searched */
  // printf("done eval with fit[0].m.d[1] = %g\n",fit[0].m.d[1]);
  // printf("printing %g\n",pars_peek(&fit[0].pars,0));
  if (searchFD != NULL) {
    fprintf(searchFD, "%g", chisq);
    for (i=0; i < fit[0].pars.n; i++)
      fprintf(searchFD, " %g", pars_peek(&fit[0].pars,i));
    fprintf(searchFD, "\n");
    fflush(searchFD);
  }

  /* Notify user of an improved fitness, including
   * logging the models in profile#.dat and fit#.dat
   */
  if (chisq < bestchi) {
    bestchi = chisq;
    pars_get(&fit[0].pars,bestpars);
    record_improvement();
  }

  return chisq;
}

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

double step_ga(int n, const double *p, void *user_data)
{
  double chisq;

  chisq = step_fn01(n,p,user_data);
  return chisq2fitness(chisq);
}

void init_ga_fit(fitinfo *fit, int argc, char *argv[])
{
  char cwd[128];
  int k;

  // Number of fit parameters
  set.nParams   = fit->pars.n;
  // Fitness function
  set.function  = step_ga;
  set.funcParms = fit;
  // Backup population file
  strncpy(set.popFile,"pop.dat",sizeof(set.popFile));
  set.popFile[sizeof(set.popFile)-1]='\0';

  /* squeeze pars into [0,1] */
  pars_set_zero_one(&fit[0].pars); 
  printf("Starting parameter values in model:\n");
  pars_print(&fit[0].pars);

  /* Prepare output file */
  if (set.popOption) { 
    /* Continuing a fit --- append to output */
    parFD  = fopen("par.dat", "a");
    fprintf(parFD, "# Restarting from stored population\n");
  } else {
    /* Starting a new fit --- write headers */
    parFD  = fopen("par.dat", "w");
    for (k=0; k<MODELS; k++)
      fprintf(parFD, "#Datafile %i: %s ; ", k, fit[k].dataA.file);
    fprintf(parFD,  "\n#%16s %17s","Generation","Chisq_best");
    for (k=0; k < fit->pars.n; k++) 
      fprintf(parFD, " %17.17s", fit->pars.name[k]);
    fprintf(parFD,"\n#%16s %17s","","parameter range:");
    for (k=0; k < fit->pars.n; k++)
      fprintf(parFD, " %8g,%-8g", fit->pars.min[k],fit->pars.min[k]+fit->pars.range[k]);
    fprintf(parFD,"\n# --------------- -----------------");
    for (k=0; k < fit->pars.n; k++) fprintf(parFD, " -----------------");
    fprintf(parFD,"\n"); fflush(parFD);
  }

  /* Save the startup command to the parameter file */
  fprintf(parFD,"# Startup command (%s):", getcwd(cwd,sizeof(cwd))?cwd:"?");
  for (k=0; k < argc; k++) fprintf(parFD, " %s", argv[k]);
  fprintf(parFD,"\n");

  /* Initialize the data structures */
  ga_init(&set,fit->pars.value);
}

void final_ga_fit(void)
{
  /* Write best population */
  write_pop_backup(&set);
  fclose(parFD);
  parFD = NULL;
}

void print_usage(void)
{
  printf("usage: fit [-flags]\n");
  printf("  -A <%%>    accelerated by computing at most %% partial chisq\n");
  printf("  -a         amoeba optimizer\n");
  printf("  -c <value> define the chisq value cut-off\n");
  printf("  -e         elitism --- keep the best model in the population\n");
  printf("  -F         fit the parameters, writing to fit.log each step (default)\n");
  printf("  -g         write profile.dat and fit.dat[ABCD]\n");
#ifdef HAVE_SETRLIMIT
  printf("  -H         max CPU hours (default 18)\n");
#endif
  printf("  -i         keep initial parameters in initial population\n");
  printf("  -j         save models as staj file\n");
#ifdef USE_NLLS_FIT
  printf("  -l         Levenberg-Marquardt optimizer\n");
#endif
  printf("  -L         log all parameter sets and associated chisq to fit.log\n");
  printf("  -m         print the initial model and profile to the screen\n");
  printf("  -n         limit the number of generations for the GA\n");
  printf("  -N         create new pop_##.dat file each trace period\n");
  printf("  -p         use saved population from pop.dat\n");
#ifdef USE_QUAD_FIT
  printf("  -Q         use Powell's NEWUOA unconstrained fit\n");
#endif
  printf("  -r <n:lo:hi> set the range for parameter n\n");
  printf("  -r <n:val> set the value for parameter n\n");
  printf("  -S         use approximate roughness\n");
  printf("  -s         seed value for random number generator\n");
  printf("  -T         trace period for writing pop_bak.dat (0 for none)\n");
  printf("  -w         force unweighted fit\n");
  printf("  -W         force weighted fit\n");
  printf("  -x <n:lo:hi:steps> print chisq landscape of parameter Pn\n");
  printf("             or of the Pm-Pn surface if -x is repeated\n");
  printf("Output\n");
  printf("  model#.dat    : best model   (d rho mu P theta)\n");
  printf("  fit#.dat[ABCD]: best fit     (Q R dR fit)\n");
  printf("  profile#.dat  : best profile (z rho mu P theta)\n");
  printf("  chisq.dat     : last chisq plot (par chisq chisq0 chisq1...\n");
  printf("  pop.dat       : current population\n");
  printf("  par.dat       : parameter improvement and command history\n");
  printf("For plots use gaplot rho|mu|P|theta|fit|chisq [file.ps]\n"); 
}

/* Make the space arbitrarily bigger, and randomize the start condition */
static void pars_explode(fitpars *p, double g)
{
  int n = pars_count(p);
  int i;
  for (i=0; i < n; i++) {
    double min = pars_min(p,i);
    double max = pars_max(p,i);
    pars_set_range(p,i,min,min+(max-min)*g);
    pars_poke01(p,i,frandom());
  }
}

int main(int argc, char *argv[])
{
  enum { 
    GA, AMOEBA, NLLS, QUADFIT, 
    CHISQ, PRINT_MODEL, PRINT_PROFILE, SAVE_STAJ 
  } action;
  int nth=0, n1=0, xdims=0, k, steps=20;
  int generations = 1000000000;
  int ch;
  double hours = 18.;

  //~ fitinfo *fit;
  srand(time(NULL));

  /* Preload model so that args can adjust parameters */
  fit = setup_models(&MODELS);
  if (fit == NULL) exit(1);
  /* Make sure the user specified a resolution */
  for (k=0; k < MODELS; k++) assert(fit[k].dataA.have_resolution==1);

  // Initialize settings
  portion = 1.0;
  initSettings(&set);
  set.initOption = 0;
  set.iElite = 0;
  action = GA;
  while((ch = getopt(argc, argv, "A:f:x:n:T:r:s:t:c:v:X:H:peSFLmgijwWNalQ?")) != -1) {
    switch(ch) {
    case 'e':
      set.iElite = 1;
      break;
      
    case 'F':
      action = GA;
      break;

    case 'm':
      // Print model and exit
      action = PRINT_MODEL;
      break;

    case 'g':
      // Print profile and exit
      action = PRINT_PROFILE;
      break;

    case 'x':
      // Chi^2 plot of nth paramater
      do {
	double lo,hi;
	if (xdims++) n1=nth;
        int n = sscanf(optarg,"%d:%lg:%lg:%d",&nth,&lo,&hi,&steps);
        if (n>=3 && nth >= 0 && nth < pars_count(&fit[0].pars))
	    pars_set_range(&fit[0].pars,nth,lo,hi);
      } while(0);
      action = CHISQ;
      break;

    case 'j':
      action = SAVE_STAJ;
      break;

    case 'a':
      action = AMOEBA;
      break;

    case 'l':
      action = NLLS;
      break;

    case 'Q':
      action = QUADFIT;
      break;

    case 'A':
      k = atoi(optarg);
      if (k < 1) portion = 0.01;
      else if (k > 100) portion = 1.0;
      else portion = k/100.;
      break;

    case 'S':
      approximate_roughness = 1;
      break;

    case 's':
      k = atoi(optarg);
      srand(k);
      break;

    case 'n':
      // Stop after then nth generation
      generations = atoi(optarg);
      break;

    case 'v':
      set.np = atoi(optarg);
      printf("Setting the number of individuals in a population to %i\n", set.np);
      break;

    case 'p':
      // Use preexisting population file
      set.popOption = 1;
      break;

    case 'i':
      // Keep initial parameters in initial population
      set.initOption = 1;
      break;

    case 't':
      // Algo type
      set.algo_type = atoi(optarg);
      break;

    case 'c':
      // Chi^2 cut
      CHI_CUT = atof(optarg);
      break;

    case 'L':
      // Keep a log of every parameter set tested
      searchFD = fopen("fit.log","w");
      break;

    case 'N':
      // Generate new pop files each trace
      set.trace_overwrite = 0;
      break;

    case 'T':
      // Population trace period
      set.trace_period = atoi(optarg);
      break;

    case 'w':
      weighted = 1;
      break;

    case 'H':
      hours = atof(optarg);
      break;

    case 'W':
      weighted = 0;
      break;

    case '?':
      print_usage();
      exit(0);
      break;

    case 'r':
      do {
	double lo=0., hi=0.;
	int k,n;
	n = sscanf(optarg,"%d:%lg:%lg",&k,&lo,&hi);
        if (n==3 && k >= 0 && k < pars_count(&fit[0].pars))
	    pars_set_range(&fit[0].pars,k,lo,hi);
        else if (n==2 && k >= 0 && k < pars_count(&fit[0].pars))
	    pars_poke(&fit[0].pars,k,lo);
      } while (0);
      break;

    case 'X':
      do {
        double explode = atof(optarg);
        if (explode > 0. && explode != 1.) pars_explode(&fit[0].pars,explode);
      } while (0);
      break;

    default:
      exit(1);
      break;
    }
  }

  set_signal_handlers();
  cpulimit(hours);

  /* Apply constraints given new limits */
  if (constraints != NULL) (*constraints)(fit);

  /* Some place to save the best parameters */
  bestpars = (double *)malloc(sizeof(*bestpars)*pars_count(&fit[0].pars));
  assert(bestpars != NULL);

#if 0
  /* Ignore best value given by the model */
  pars_get(&fit[0].pars,bestpars);
  bestchi = update_models(fit);
#endif

  // Testing ground
  // Large population Ula-200
  //~ set.np=200; 
  // Large population Ula-100
  // set.np=100;
  // Fitness slope fslope-100
  //~ set.fitnessSlope=1.0;
  // Fitness slope fslope-75
  //~ set.fitnessSlope=0.75;
  // crossing prob cross-95
  //~ set.pCross=0.95;
  // crossing prob cross-60
  set.pCross=0.60;
  
  
  switch (action) {
  case GA:
    init_ga_fit(fit,argc,argv);
    ga_fit(&set,generations);
    final_ga_fit();

    break;

  case AMOEBA:
    // Need to separate amoeba from GA; currently faking ga init to set up data
    //set.np = 1;
    init_ga_fit(fit,argc,argv);
    cmd_amoeba();
    cmd_amoeba();
    cmd_amoeba();
    final_ga_fit();
    break;

  case NLLS:
#ifdef USE_NLLS_FIT
    // Need to separate nlls from GA; currently faking ga init to set up data
    //set.np = 1;
    init_ga_fit(fit,argc,argv);
    cmd_nlls();
    printf("LM completed...\n");fflush(stdout);
    final_ga_fit();
#endif
    break;

  case QUADFIT:
#ifdef USE_QUAD_FIT
    //set.np = 1;
    init_ga_fit(fit,argc,argv);
    cmd_quadfit();
    final_ga_fit();
#endif
    break;

  case PRINT_MODEL:
    // Print model and exit
    for (k=0; k < MODELS; k++) {
      printf("== Model %d ==\n",k);
      model_profile(&fit[k].m, &fit[k].p);
      profile_print(&fit[k].p,NULL);
      model_print(&fit[k].m,NULL);
    }
    break; 

  case PRINT_PROFILE:
    /* Print the initial profile and fit to file */
    update_models(fit);
    save_models(fit);
    break;

  case CHISQ:
    bestchi = 1e308;
    if (xdims > 1) 
      chisq_surf(fit,n1,nth,portion,steps);
    else
      chisq_plot(fit,nth,portion,steps);
    break;

  case SAVE_STAJ:
    save_all_staj();
    break;
  }

  if (searchFD != NULL) fclose(searchFD);
  for (k=0; k < MODELS; k++) fit_destroy(&fit[k]);
  exit(0);
  return 0;
}

/* $Id$ */
