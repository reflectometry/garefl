/**
 
  Genetic fitting algorithm 
  Based on PIKAIA, Charbonneau and Knapp (1995).
  Described in NCAR/TN-418+IA [National Center for Atmospheric Research Technical Note]
  
  @author   M.Doucet, UMd/NIST <doucet@nist.gov>
  @version  March 31 , 2004

*/
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ga.h"

/// CVS information
static char cvsid[] = "$Id$";


/** Comparison function used for sorting */
static int compare(const void *e1, const void *e2);

/** Ranks the population */
static int rankPop( Settings *set );

/** Select parents from the population */
static int selectPar( Settings *set );

/** Encode parameter values */
static void encode( Individual *indiv );

/** Cross two parents to create offsprings */
static void cross( Real pCross, Individual *indiv1, Individual *indiv2 );

/** Uniform mutation */
static void mutate( Real pmut, Individual *indiv );

/** Decode chromos... */
static void decode( Individual *indiv );

/** Replace the population with the newly generated one */
static void newpop( Settings *set);

/** Adjust mutation probability */
static Real adjmut(Settings *set);



// Random number generator
Real frandom(void) { return (Real)rand()/RAND_MAX; }

//-------------------------------------------------------------
// Element table used for sorting 
struct element {int i; Real data;} Table[PMAX];

//-------------------------------------------------------------
/** Find the fittest individual */
int fittest( Settings *set)
{
  int k, best;

  best = 0;
  for (k=1; k < set->np; k++) {
    if (set->pop->indiv[k].fitness > set->pop->indiv[best].fitness) 
      best = k;
  }
  return best;
}


//-------------------------------------------------------------
// Comparison function used for sorting
static int compare(const void *e1, const void *e2) {
  Real v1 = ((struct element *)e1)->data;
  Real v2 = ((struct element *)e2)->data;
  return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

//-------------------------------------------------------------
/** Ranks the population */
static int rankPop( Settings *set ) {
  int i;
  for(i=0;i<set->np;i++) {
    Table[i].data = set->pop->indiv[i].fitness;
    Table[i].i    = i;
  }
  qsort( (void *) Table, (size_t) set->np, sizeof(struct element), compare);
  // index gives each entry in increasing order of fitness [0 to np-1]
  for(i=0;i<set->np;i++) set->pop->index[i] = Table[i].i;
  // rank gives the ranking for a given entry (ip) [1 to 10] 
  for(i=0;i<set->np;i++) set->pop->indiv[set->pop->index[i]].rank = set->np-i;
  return set->pop->index[set->np-1];
}

//-------------------------------------------------------------
/** Select parents from the population */
static int selectPar( Settings *set ) { 
  int i;
  Real dice, rtfit;
  
  dice = frandom()*set->np*(set->np+1);
  rtfit=0;
  for(i=0;i<set->np;i++) {
    rtfit = rtfit+(set->np+1)+set->fitnessSlope*((set->np+1)-2*set->pop->indiv[i].rank);
    if (rtfit>=dice)  return i;
  }
  return -1;
}

//-------------------------------------------------------------
/** Encode parameter values */
static void encode( Individual *indiv ) { 
  Real z;
  int ip, i, j, k;

  //~ if(indiv->encoded==1) return;

  z = pow(10.0,indiv->nPrec);
  k=0;
  for(i=0;i<indiv->nParams;i++) {
    ip=(int)(indiv->value[i]*z);
    for(j=indiv->nPrec-1;j>=0;j--) {
      indiv->chromo[k+j]=(int)(fmod(ip,10.0));
      ip=ip/10;
    }
    k+=indiv->nPrec;
  }
  indiv->encoded=1;
}

//-------------------------------------------------------------
/** Cross two parents to create offsprings */
static void cross( Real pCross, Individual *indiv1, Individual *indiv2 ) {

  if(frandom()<pCross) {
    // cross-over point
    // Allow complete exchange of parameters
    int isp1 = (int)(frandom()*indiv1->nParams*indiv1->nPrec);
    int i;
  
    for(i=isp1;i<indiv1->nParams*indiv1->nPrec;i++){
      int t = indiv2->chromo[i];
      indiv2->chromo[i]=indiv1->chromo[i];
      indiv1->chromo[i]=t;
    }
  }
}

//-------------------------------------------------------------
/** Uniform mutation */
static void mutate( Real pmut, Individual *indiv ) {
  int i;
  for(i=0;i<indiv->nParams*indiv->nPrec;i++) 
    if(frandom()<pmut) indiv->chromo[i] = (int)(frandom()*10.0);
}

//-------------------------------------------------------------
/** Decode chromos... */
static void decode( Individual *indiv ){
  Real z;
  int ip, i, j, k;
  z = pow(10.0,-(indiv->nPrec));
  k=0;

  for(i=0;i<indiv->nParams;i++) {
    ip=0;
    for(j=0;j<indiv->nPrec;j++) ip=10*ip+indiv->chromo[k+j];
    indiv->value[i]=ip*z;
    k+=indiv->nPrec;
  }
}

//-------------------------------------------------------------
/** Replace the population with the newly generated one */
static void newpop( Settings *set ) {
  int i, k, best = -1;
  
  // If elitist option, keep the fittest
  if(set->iElite==1) best = fittest(set);
  
  for(i=0;i<set->np;i++) {
    if (i != best) {
      for(k=0;k<set->nParams;k++) 
	set->pop->indiv[i].value[k] = set->pop->newph[i].value[k];
      set->pop->indiv[i].fitness
	= set->function(set->nParams,set->pop->newph[i].value,set->funcParms);
    }
  }
  
  // Rank population
  rankPop(set);
}


//-------------------------------------------------------------
void setChromosome(Settings *set, int n, Real *d)
{
  int i;

  assert(n>=0 && n<set->np);
  for (i=0; i < set->nParams; i++) set->pop->indiv[n].value[i]=d[i];
  set->pop->indiv[n].fitness 
    = set->function(set->nParams,set->pop->indiv[n].value,set->funcParms);
}

//-------------------------------------------------------------
Real getChromosome(Settings *set, int n, Real *d)
{
  int i;

  assert(n>=0 && n<set->np);
  for (i=0; i < set->nParams; i++) d[i] = set->pop->indiv[n].value[i];
  return set->pop->indiv[n].fitness;
}  

//-------------------------------------------------------------
//~ Real adjmut(Settings *set, int *ifit, Real *fitness){
static Real adjmut( Settings *set ){
  Real rdif=0, rdiflo=0.05, rdifhi=0.25, delta=1.5, old_value;
  int i;
  
  // Adjustment based on fitness differential
  if(set->iMutate==2 || set->iMutate==5) {
    rdif=fabs( set->pop->indiv[set->pop->index[set->nParams-1]].fitness-
                 set->pop->indiv[set->pop->index[set->nParams/2-1]].fitness )/
             ( set->pop->indiv[set->pop->index[set->nParams-1]].fitness+
                 set->pop->indiv[set->pop->index[set->nParams/2-1]].fitness );
  }
  
  // Adjustment based on normalized metric distance
  if(set->iMutate==3 || set->iMutate==6) {
    rdif=0;
    for(i=0;i<set->nParams;i++) 
        rdif += (set->pop->indiv[set->pop->index[set->nParams-1]].value[i]
                  -set->pop->indiv[set->pop->index[set->nParams/2-1]].value[i])*
                (set->pop->indiv[set->pop->index[set->nParams-1]].value[i]
                  -set->pop->indiv[set->pop->index[set->nParams/2-1]].value[i]);
    rdif=sqrt(rdif)/set->nParams;
  }
  
  old_value = set->pMutate;
  if(rdif<=rdiflo) {
    if(set->pMutateMax<set->pMutate*delta) set->pMutate = set->pMutateMax;
      else set->pMutate = set->pMutate*delta;
  } else if(rdif>=rdifhi) {
    if(set->pMutateMin>set->pMutate/delta) set->pMutate = set->pMutateMin;
      else set->pMutate = set->pMutate/delta;
  }

  if (set->pMutate != old_value) 
    printf("Mutation prob. changed from %g to %g\n", old_value, set->pMutate); 
  return set->pMutate;
}

//-------------------------------------------------------------
// Internal function to write the population to a file
static int _write_pop( const char *filename, Settings *set) {
  int ip,k;
  FILE *file;

  printf("(%i) Writing population [%i x %i] in %s\n", 
     GetGen(set), set->np, set->nParams, filename);

  file = fopen(filename,"w");
  if (file == NULL) {
    perror(filename);
    return -1;
  }

  fprintf(file,"%i\n",GetGen(set));
  for(ip=0; ip<set->np; ip++) {
    for(k=0; k<set->nParams; k++) 
      fprintf(file,"%15g ",set->pop->indiv[ip].value[k]);
    fprintf(file,"\n");
  }
  fclose(file);
  return 0;
}

// Write a population snapshot to pop_*.dat
int write_pop_backup( Settings *set) {
	char popFile[80];

    if (set->trace_overwrite)
      strncpy(popFile,"pop_bak.dat",sizeof(popFile));
    else
      snprintf(popFile,sizeof(popFile),"pop_%i.dat",set->genNumber);
    popFile[sizeof(set->popFile)-1]='\0';
    
    return _write_pop(popFile, set);
}

// Write the final population to popFile (probably pop.dat)
int write_pop( Settings *set) {
  return _write_pop(set->popFile, set);
}

int read_pop( Settings *set, const char *filename ) {
  float buffer;
  int k, ip;

  FILE *file = fopen(filename,"r");
  if (file == NULL) {
    perror(set->popFile);
    return -1;
  }

  if (fscanf(file,"%i",&set->genNumber) != 1) set->genNumber=0;
  k = -1; ip = 0;
  while(fscanf(file,"%f",&buffer)==1) {
    if (k>=set->nParams-1) {
      k=0; 
      if (++ip>=set->np) {
      	printf("Population is too large.\n");
        fclose(file);
      	return -1;
      }
    } else {
      k++;
    }
    set->pop->indiv[ip].value[k] = (Real)buffer;
  }
  fclose(file);
  if (ip != set->np-1 || k != set->nParams-1) {
  	printf("Population is too small (stopped at %d,%d)\n",ip,k);
  	return -1;
  }
  // _write_pop("read.dat",set);  // ...so we can check that reading works...
  return 0;
}
//-------------------------------------------------------------
void statReportHeader( Settings *set ) {
  FILE *file;
  time_t now;
  if (!set->print_stats) return;
  now = time(NULL);
  file = fopen("report.dat","w");
  fprintf(file,"Starting fit -- %s\n",ctime(&now));
  fclose(file);
}
//-------------------------------------------------------------
void statReport( Settings *set ) {
  int i;
  Real fitmin=-1, fitmax=-1, fitsum=0;
  FILE *file;
  if (!set->print_stats) return;
  for(i=0; i<set->np; i++) {
    if (fitmin<0 || set->pop->indiv[i].fitness<fitmin) 
      fitmin=set->pop->indiv[i].fitness;
    if (set->pop->indiv[i].fitness>fitmax) 
      fitmax=set->pop->indiv[i].fitness;
    fitsum += set->pop->indiv[i].fitness;
  }
  
  file = fopen("report.dat","a");
  fprintf(file,"%i %g %g %g %g\n",GetGen(set), fitmin, fitmax, fitsum, fitsum/set->np);
  fclose(file);
}


//-------------------------------------------------------------
void setParValue( Settings *set, int ipar, Real value ) {
    int ip;
    for(ip=0; ip<set->np; ip++) {
        set->pop->indiv[ip].value[ipar] = value;
        set->pop->indiv[ip].fitness = 
          set->function(set->nParams,set->pop->indiv[ip].value,set->funcParms);
    }
}
//-------------------------------------------------------------
void setRange( Settings *set, int ipar, Real oldmin, Real oldrange,
                Real newmin, Real newrange ) {
    int ip;
    Real oldvalue;
    for(ip=0; ip<set->np; ip++) {
        oldvalue = set->pop->indiv[ip].value[ipar];
        set->pop->indiv[ip].value[ipar] = (oldmin+oldrange*oldvalue-newmin)/newrange;
    }
}

//-------------------------------------------------------------
void randomizePop( Settings *set, int par ) {
  int k, ip, best;

  best = fittest(set);
  if(par>=0 && par<set->nParams) {
    printf("Randomizing parameter %i\n", par);
    for(ip=0; ip<set->np; ip++) {
      if(ip != best) {
        set->pop->indiv[ip].value[par] = frandom();
        set->pop->indiv[ip].fitness 
	  = set->function(set->nParams,set->pop->indiv[ip].value,set->funcParms);
      }
    }
  } else if(par<0) {
    printf("Randomizing the whole set of parameters\n");
    for(ip=0; ip<set->np; ip++) {
      if (ip != best) {
	for(k=0; k<set->nParams; k++) set->pop->indiv[ip].value[k] = frandom();
	set->pop->indiv[ip].fitness 
	  = set->function(set->nParams,set->pop->indiv[ip].value,set->funcParms);
      }
    }
  } else {
    printf("Bad randomize parameter: nothing done\n");
  }
}

//-------------------------------------------------------------
void initSettings( Settings *set ) {
  // Number of fit parameters
  set->nParams   = 0;
  // Number of individuals in population
  set->np        = 50;
  // Encoding constant
  set->nd        = 4;
  // Cross probability
  set->pCross    = 0.85;
  // Mutation probability
  set->pMutate   = 0.2;
  // Mutation probability adjustment parameter
  set->iMutate   = 3;
  // Mutation probability min and max
  set->pMutateMin=.01;
  set->pMutateMax=.35;
  set->fitnessSlope=0.5;
  set->iElite = 1;
  // Load from pop datafile
  set->popOption = 0;
  // Algorithm type
  set->algo_type = 0;
  // Keep initial parameters in the population
  set->initOption = 0;
  set->genNumber=0;
  set->trace_period = 20;
  set->trace_overwrite = 1;
  set->print_stats = 0;
}

//-------------------------------------------------------------
int GetGen( Settings *set ) { return set->genNumber; }

//-------------------------------------------------------------
void initPop( Settings *set ) {
  int i;
  set->pop = (Population *) malloc(sizeof(Population));
  for(i=0;i<set->np;i++) {
    set->pop->indiv[i].fitness = 0.;
      set->pop->indiv[i].value = (Real *)malloc(set->nParams*sizeof(Real));
      set->pop->indiv[i].chromo = (int *)malloc(set->nParams*set->nd*sizeof(int));
      set->pop->newph[i].value = (Real *)malloc(set->nParams*sizeof(Real));
      set->pop->newph[i].chromo = (int *)malloc(set->nParams*set->nd*sizeof(int));
      set->pop->indiv[i].nParams = set->nParams;
      set->pop->indiv[i].nPrec = set->nd;
      set->pop->newph[i].nParams = set->nParams;
      set->pop->newph[i].nPrec = set->nd;
  }
}

//-------------------------------------------------------------
void initIndividual( Settings *set, Individual *indiv ) {
    indiv->value = (Real *)malloc(set->nParams*sizeof(Real));
    indiv->chromo = (int *)malloc(set->nParams*set->nd*sizeof(int));
    indiv->nParams = set->nParams;
    indiv->nPrec = set->nd;
}

void freeIndividual( Individual *indiv ) {
	free(indiv->value); indiv->value = NULL;
    free(indiv->chromo); indiv->chromo = NULL;
}

//-------------------------------------------------------------
void copyValues( Individual *in, Individual *out ) {
  int i;
  for(i=0;i<in->nParams;i++) out->value[i] = in->value[i];
};

//-------------------------------------------------------------
//Population initialization
void ga_init( Settings *set, Real *params ) {
  int n             = set->nParams;
  int np            = set->np;
  int ip, k;

  assert(np <= PMAX);
  assert(n <= NMAX);
  printf("Starting ga_fit...\n");
  
  // Allocate data memory
  initPop(set);

  // Initialize (seed first, then random)
  if (set->popOption==1) {
    printf("Reading population [%i x %i] from %s\n", 
	       set->np, set->nParams, set->popFile);
    if (read_pop(set, set->popFile) != 0) exit(1);
  } else {
    printf("Generating random population [%i x %i]\n", 
	       set->np, set->nParams);
    ip = 0;
    if(set->initOption==1) {
      printf("Keeping current parameters in the population\n");
      for(k=0; k<n; k++) set->pop->indiv[0].value[k] = params[k];
      ip++;
    }
    while (ip < np) {
      for(k=0; k<n; k++) set->pop->indiv[ip].value[k] = frandom();
      ip++;
    }
  }

  // Compute initial fitness of each member of the pop.
  for(ip=0; ip<np; ip++) {
    set->pop->indiv[ip].fitness 
      = set->function(n,set->pop->indiv[ip].value,set->funcParms);
  }

  // Rank population
  rankPop(set);

  // Prepare status reports
  statReportHeader(set);
}

// Main fit loop
Real ga_fit( Settings *set, int ngenerations ) {
  int np            = set->np;
  Real pcross     = set->pCross;
  Real pmut       = set->pMutate;
  int ip, ig;
  int ip1, ip2;
  Individual indiv1, indiv2;

  initIndividual(set, &indiv1);
  initIndividual(set, &indiv2);
  
  // Main generation loop
  for(ig=0;ig<ngenerations;ig++) {
    set->genNumber++;
    
    /* write pop_bak.dat every X generations */
    if (set->genNumber && set->trace_period &&
        (set->genNumber % set->trace_period)==0){
       write_pop_backup(set);
    }
  
    // Main population loop
    for(ip=0;ip<np/2;ip++) {
      
      // Find two parents
      ip1 = selectPar(set);
      ip2 = ip1;
      while(ip2==ip1) ip2 = selectPar(set);
      copyValues(&set->pop->indiv[ip1], &indiv1);
      copyValues(&set->pop->indiv[ip2], &indiv2);
      
      // Encode
      encode(&indiv1);
      encode(&indiv2);
      
      // breed
      cross(pcross, &indiv1, &indiv2);
      mutate( pmut, &indiv1 ); 
      mutate( pmut, &indiv2 ); 
    
      // Decode offsprings
      decode(&indiv1);
      decode(&indiv2);
      
      // Insert into population
      copyValues(&indiv1, &set->pop->newph[2*ip]);
      copyValues(&indiv2, &set->pop->newph[2*ip+1]);
    }

    newpop(set);
    adjmut(set);
    statReport(set);
  }
  
  freeIndividual(&indiv1);
  freeIndividual(&indiv2);

#ifdef DEBUG
  for(k=0;k<5;k++) {
    int i;
    printf("------------------------------------\n");
    printf("Solution %d:     Fitness: %10.5e\n", 
	   k+1,set->pop->indiv[index[set->np-1-k]].fitness);
    for(i=0;i<(set->nParams);i++) 
      printf("%d: %f\n",i,set->pop[set->pop->indiv[index[set->np-1-k]]][i]);
  }
  printf("------------------------------------\n");
#endif
  
  return set->pop->indiv[fittest(set)].fitness;
}


#ifdef TEST

//-------------------------------------------------------------
// Fake fitness function for testing
static Real fakefunk( int n, Real *vars, void *dummy ) {
  Real r;
  r = (vars[0]-0.5)*(vars[0]-0.5) + (vars[1]-0.5)*(vars[1]-0.5);
  return cos(9.0*acos(-1)*r)*cos(9.0*acos(-1)*r)*exp(-r*r/0.15);
}

int main(int argc, char *argv[]) {

  // Fit parameters
  Real parms[2];
  Real f;
  int i;
  
  Settings set;
  set.iElite    = 1;
  // Number of fit parameters
  set.nParams   = 2;
  // Number of individuals in population
  set.np        = 100;
  set.nd        = 5;
  // Cross probability
  set.pCross    = 0.85;
  // Mutation probability
  set.pMutate   = 0.05;
  set.iMutate   = 2;
  set.pMutateMin=.0005;
  set.pMutateMax=.25;
  // Fitness function
  set.function  = fakefunk;
  set.funcParms = 0;

  ga_init( &set, parms );  
  f = ga_fit( &set, 100 );
  printf("Fitness: %f\n", f);
  for(i=0;i<set.nParams;i++){
    printf("%d: %f\n",i,parms[i]);
  }
  return 0;
}
#endif

/*
 * $Log$
 * Revision 1.4  2005/11/17 20:48:35  pkienzle
 * Make sure amoeba is initialized from best individual of the population;
 * distinguish automatic and user requested state backup files;
 * refactor parameter display; add 'trust region' command to menu; hide
 * some ga internals; partial code cleanup in ga.
 *
 * Revision 1.3  2005/10/07 16:59:02  pkienzle
 * Suppress debugging info.
 *
 * Revision 1.2  2005/08/24 18:36:57  pkienzle
 * Don't create report.dat when running---file gets too big.
 *
 * Revision 1.1  2005/08/02 00:18:19  pkienzle
 * initial release
 *
 * Revision 1.29  2005/04/15 17:07:54  pkienzle
 * error checking on read/write pop; fix generation number on startup
 *
 * Revision 1.28  2005/04/14 21:46:28  pkienzle
 * Update the whole population before calculating fitness to make it robust
 * against interrupt and change of parameter values midgeneration.
 *
 * Revision 1.27  2005/03/13 17:55:54  pkienzle
 * Remove spurious warnings
 *
 * Revision 1.26  2005/03/10 23:34:53  pkienzle
 * -NT# trace every # in pop_##.dat; -T# trace every # in pop_bak.dat
 *
 * Revision 1.25  2005/03/10 20:36:54  pkienzle
 * Use autoconf to build the system
 * Use HAVE_MAGNETIC as flag
 * Put configuration flags in reflconfig.h
 * Change ga_refl examples to ga_simul
 * Add trace option to ga for writing periodic pop files
 * Add dummy main for fortran in case some system needs it
 * Remove preprocessor directives from fortran (use r4xa exclusively)
 *
 * Revision 1.24  2005/02/01 17:17:54  pkienzle
 * getChromosone returns fitness; objective function takes const double *
 *
 * Revision 1.23  2004/10/21 13:57:54  doucet
 * Rewrote parts of ga.c
 *
 * Revision 1.22  2004/09/22 21:09:32  doucet
 * Fixed bug in pars_zero_one. The [0,1] values must be stored if the solution is to be used later.
 *
 * Revision 1.21  2004/09/22 20:12:56  doucet
 * Put 'Data' in 'Settings'
 *
 * Revision 1.20  2004/09/17 21:55:46  doucet
 * We can now change the range of a parameter using ctrl-c
 * BEWARE: the new ranges are used when writing pop.dat!
 *
 */
