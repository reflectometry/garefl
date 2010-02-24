#ifndef _GA_H
#define _GA_H 1

#ifndef Real
# define Real double
#endif

#define NMAX 250
#define PMAX 1024

// Individual
typedef struct {
  // Chromosome representation
  Real *value;
  int    *chromo;
  int    nParams;
  int    nPrec;
  Real fitness;
  int    rank;
  int    encoded;
} Individual;

// Population
typedef struct {
  Individual indiv[PMAX];
  Individual newph[PMAX];
  // index gives each entry in increasing order of fitness [0 to np-1]
  int index[PMAX];
} Population;

typedef struct {
  // Elite option
  int iElite;
  // Number of fit parameters
  int nParams;
  //Number of individuals in population
  int np;
  // Encoding constant
  int nd;
  // Cross probability
  Real pCross;
  // Mutation probability
  Real pMutate;
  // Number of generations
  int nGenerations;
  // Mutation probability adjustment parameter
  int iMutate;
  Real pMutateMin;
  Real pMutateMax;
  Real fitnessSlope;
  // Fitness function
  Real (*function)(int n, const Real*v, void*p);
  void *funcParms;
  int popOption;
  int algo_type;
  int initOption;
  char popFile[64];
  // Data
  Population *pop;
  int genNumber;
  // Write population ever X generations.
  int trace_period;
  int trace_overwrite;
  // Add to the report every generation. 
  int print_stats;
} Settings;

/** Random number generator */
Real frandom();

/** get/set chromosomes */
/** Find the fittest individual */
Real getChromosome(Settings *set, int n, Real *d);
void setChromosome(Settings *set, int n, Real *d);
int fittest(Settings *set);

/** Perform fit */
void ga_init( Settings *set, Real *params );
Real ga_fit( Settings *set, int nGenerations );

/** Read and write population */
int write_pop_backup( Settings *set);
int write_pop( Settings *set);
int read_pop( Settings *set, const char *file);

void initSettings( Settings *set );
void randomizePop( Settings *set, int par );
int GetGen( Settings *set );

/** Sets a parameter to a given value for all members of the current population */
void setParValue( Settings *set, int ipar, Real value);
void setRange( Settings *set, int ipar, Real oldmin, Real oldrange,
                Real newmin, Real newrange );
#endif
