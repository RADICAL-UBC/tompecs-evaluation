#include <stdlib.h> /* srand, rand */
#include <stdio.h> /* printf, NULL */
#include <math.h>
#include <string.h> // memcpy
#include <assert.h> // assert
#include <locale.h>
#include <float.h> /* for machine epsilons; to do floating point comparisons */
#include <stdbool.h>
#include <time.h>
#include <string.h> /* memset */
#include "gurobi_c.h"

#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "ransampl.h"

#define NO_ACTION -1 

char essentially_equal(double a, double b) {
  return fabs(a - b) <=
    ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * DBL_EPSILON );
}

char approximately_equal(double a, double b) {
  return fabs(a - b) <=
    ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * DBL_EPSILON );
}

char definitely_greater_than(double a, double b) {
  return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * DBL_EPSILON );
}

char definitely_less_than(double a, double b) {
  return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * DBL_EPSILON );
}

char definitely_less_than_or_equal(double a, double b) {
  return (definitely_less_than(a, b) || essentially_equal(a, b));
}

double *uunifast(int num_components, double vector_sum) {
  
  srand(time(NULL));
  
  double *uniform_random_vector = calloc(num_components, sizeof(double));
  
  /* for(int index = 0; index < num_components; ++index) { */
    
  /*   uniform_random_vector[index] = (double)1/((double)num_components); */
  /*   printf( "pmf element: %lf\n", uniform_random_vector[index]); */
    
  /* } */
  
  /* return uniform_random_vector; */
  
  double sum = vector_sum;
  double next_sum;
  
  double uniform_real_number = 0;
  int index = 0;
  for(index = 0; index < num_components - 1; ++index) {
    
    uniform_real_number = (double)rand() / (double)RAND_MAX;
    next_sum = sum * pow(uniform_real_number,
			 1/(double)(num_components - index));
    uniform_random_vector[index] = sum - next_sum;
    sum = next_sum;
  }
  uniform_random_vector[index] = sum;
  return uniform_random_vector;
}

/* Returns: a probability mass function (pmf) on {1, ..., max_execution_time} */
/* as an array of probabilities of max_execution_time elements. */
double *generate_execution_time_pmf(int max_execution_time) {
  return uunifast(max_execution_time, 1);
}

/* To find the number of an array defined using []. Does not work for */
/* dynamically allocated arrays */
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

#define NUM_CRITICALITY_LEVELS 2
#define NUM_JOBS 3

typedef unsigned long ul;

char arr_num_jobs_of_criticality[NUM_CRITICALITY_LEVELS];
ul num_states = 0;

/* Values that state->finish_signals[i] assume */
/* #define FINISHED     0 */
/* #define NOT_FINISHED 1 */

enum CRITICALITY {LO=0, HI, UNKNOWN};

enum FINISH_SIGNALS {
  FINISHED = 0,
  NOT_FINISHED,
  MAX_FINISH_SIGNAL};

/* Values that state->error_flag assumes */
enum ERROR_FLAGS {
  NOT_ERROR=0,
  ERROR,
  POTENTIAL_ERROR,
  ERROR_ABSORB,
  MAX_ERROR_FLAG
};

ul *initial_allocation_vector;
char *initial_finish_signal_vector; 
int initial_error_flag;

typedef struct state_action_pair state_action_pair;
typedef struct state_action_pair_pointer state_action_pair_pointer;
/* typedef struct out_state_pointer out_state_pointer; */

typedef struct state_tag {

  /* those are pointer to elements in state_actions_pairs_list */
  state_action_pair *admissible_state_action_pairs_list;
  int num_admissible_saps;
  
  state_action_pair_pointer *in_state_action_pairs_list;

  /* But field of size = 1 bit */
  char is_initial : 1;
  
  double constraint_cost;
  /* double objective_cost; */
  
  ul time;
  struct state_tag *next;
  
} state;

state *state_tree = NULL;
state *state_tree_tail = NULL;

state **arr_states; /* Multidimensional state pointers */

/* The state of the simulation */
typedef struct sim_state_tag {
  
  ul   *allocation_vector;
  char *finish_signal_vector;
  int error_flag;
  
} sim_state_t;

/* The state of the simulation */
typedef struct sim_stats_tag {
  
  int *arr_num_deadline_misses_per_job;
  /* int *arr_num_errors_per_criticality; */
  int num_errors;
  
} sim_stats_t;

struct state_action_pair_pointer {
  

  state_action_pair_pointer *next;
  ul sap_index;
  /* specific to the parent state that has this sap as an in_state_action_pair. 
     Irrelevant if this object is an admissible state_action pair */
  double prob_parent_state_given_sap;
  
};
/* /\* This is contained in a state_action_pair that that lead to the state *\/ */
/* struct out_state_pointer { */
  
/*   state *s; /\* A pointer to an element in state_tree *\/ */
/*   out_state_pointer *next; */

/*   double prob_state_given_containing_sap; */
/*   double cost_state_given_sap; */
  
/* }; */


state *create_state() {
  
  state *s = malloc( sizeof(state) ); 
  s->next				= NULL;
  s->time				= 0;
  s->admissible_state_action_pairs_list	= NULL;
  s->num_admissible_saps		= 0;
  
  s->in_state_action_pairs_list		= NULL;

  /* s->objective_cost			= 0; */
  s->constraint_cost			= 0;
  s->is_initial				= false;
  return s;
  
}

struct state_action_pair {
  
  ul index;
  /* state *s; /\* A pointer to an element in state_tree *\/ */
  int action;
  state_action_pair* next;
  double objective_cost;
  /* out_state_pointer* out_states_list; */
};

ul num_state_action_pairs = 0;

/* Those will be the optimization variables */
state_action_pair* state_action_pairs_list;

/* Used to add a state_action pair to the end of state_action_pairs_list */
state_action_pair* state_action_pairs_list_tail;


int horizon = 0;

/* typedef struct state_pointer_per_allocation_tag { */
  
/*   state* s; /\* A pointer to an element in state_tree *\/ */
/*   struct state_pointer_per_allocation_tag *next; */
  
/* } state_pointer_per_allocation; */

/* state_pointer_per_allocation **arr_allocation_state_pointers; */

typedef struct job_tag {
  int deadline;
  int criticality;
  int arr_wcets_at_criticalities[NUM_CRITICALITY_LEVELS];
  
  /* The execution pmf is a probability mass function on  
     { 1, ..., arr_wcets_at_criticalities[criticality] } */
  double *execution_time_pmf;
  
  /* Might be assigned by a fixed-priority algorithm. */
  int priority;
  
} job ;

job arr_jobs[NUM_JOBS];
int *arr_job_indexes_per_criticality[NUM_CRITICALITY_LEVELS];


ul dims[2 * NUM_JOBS + 1];

#define ARR2(x1,x2,y1,y2,r) (arr_states					\
			     [dims[1]*dims[2]*dims[3]*dims[4]*x1 +	\
			      dims[2]*dims[3]*dims[4]*x2+		\
			      dims[3]*dims[4]*y1+			\
			      dims[4]*y2+				\
			      r])

#define ARR3(x1,x2,x3,y1,y2,y3,r) (arr_states	\
  [dims[1]*dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*x1 + \
   dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*x2+\
   dims[3]*dims[4]*dims[5]*dims[6]*x3+\
   dims[4]*dims[5]*dims[6]*y1+\
   dims[5]*dims[6]*y2+\
   dims[6]*y3+\
   r])

#define ARR4(x1,x2,x3,x4,y1,y2,y3,y4,r) (arr_states	\
  [dims[1]*dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x1 + \
   dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x2+\
   dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x3+\
   dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x4+\
   dims[5]*dims[6]*dims[7]*dims[8]*y1+\
   dims[6]*dims[7]*dims[8]*y2+\
   dims[7]*dims[8]*y3+\
   dims[8]*y4+\
   r])

#define ARR5(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,r) (arr_states	\
  [dims[1]*dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*dims[9]*dims[10]*x1 + \
   dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*dims[9]*dims[10]*x2+\
   dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*dims[9]*dims[10]*x3+\
   dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*dims[9]*dims[10]*x4+\
   dims[5]*dims[6]*dims[7]*dims[8]*dims[9]*dims[10]*x5+\
   dims[6]*dims[7]*dims[8]*dims[9]*dims[10]*y1+\
   dims[7]*dims[8]*dims[9]*dims[10]*y2+\
   dims[8]*dims[9]*dims[10]*y3+\
   dims[9]*dims[10]*y4+\
   dims[10]*y5+\
   r])

state *lookup_state_ptr_by_params(ul *allocation_vector,
				  char *finish_signal_vector,
				  int error_flag) {
  state *s = NULL;
  
#if NUM_JOBS == 2
  
  s = ARR2(allocation_vector[0],
	   allocation_vector[1],
	   finish_signal_vector[0],
	   finish_signal_vector[1],
	   error_flag);
  
#elif NUM_JOBS == 3
  
  s = ARR3(allocation_vector[0],
	   allocation_vector[1],
	   allocation_vector[2],
	   finish_signal_vector[0],
	   finish_signal_vector[1],
	   finish_signal_vector[2],
	   error_flag);
  
#elif NUM_JOBS == 4
  
  s = ARR4(allocation_vector[0],
	   allocation_vector[1],
	   allocation_vector[2],
	   allocation_vector[3],
	   finish_signal_vector[0],
	   finish_signal_vector[1],
	   finish_signal_vector[2],
	   finish_signal_vector[3],
	   error_flag);
  
#elif NUM_JOBS == 5
  
  s = ARR5(allocation_vector[0],
	   allocation_vector[1],
	   allocation_vector[2],
	   allocation_vector[3],
	   allocation_vector[4],
	   finish_signal_vector[0],
	   finish_signal_vector[1],
	   finish_signal_vector[2],
	   finish_signal_vector[3],
	   finish_signal_vector[4],
	   error_flag);
#endif  
  return s;
}


/* Once a state_action pair is created, it cannot be destroyed and MUST be 
   linked because this function updates the global variable 
   state_action_pairs_count  */
state_action_pair* create_state_action_pair() {
  
  state_action_pair *sap = malloc(sizeof(state_action_pair));
  
  /* sap->s = NULL; */
  sap->action	      = 0;
  sap->next	      = NULL;
  sap->objective_cost = 0;
  /* sap->out_states_list = NULL; */
  return sap;
  
}

ul *clone_ul_array(ul *array, int num_elements) {
 
  ul *arr_copy = malloc(num_elements * sizeof(ul));
  memcpy(arr_copy, array, num_elements * sizeof(ul));
  
  return arr_copy;
  
}

char *clone_char_array(char *array, int num_elements) {
 
  char *arr_copy = malloc(num_elements * sizeof(char));
  memcpy(arr_copy, array, num_elements * sizeof(char));
  
  return arr_copy;
  
}

state* add_state(ul *allocation_vector,
		 char *finish_signal_vector,
		 int error_flag) {
  
  state* s = create_state();

#if NUM_JOBS == 2
  
  ARR2(allocation_vector[0],
       allocation_vector[1],
       finish_signal_vector[0],
       finish_signal_vector[1],
       error_flag) = s;

#elif NUM_JOBS == 3
  
  ARR3(allocation_vector[0],
       allocation_vector[1],
       allocation_vector[2],
       finish_signal_vector[0],
       finish_signal_vector[1],
       finish_signal_vector[2],
       error_flag) = s;

#elif NUM_JOBS == 4
  
  ARR4(allocation_vector[0],
       allocation_vector[1],
       allocation_vector[2],
       allocation_vector[3],
       finish_signal_vector[0],
       finish_signal_vector[1],
       finish_signal_vector[2],
       finish_signal_vector[3],
       error_flag) = s;

#elif NUM_JOBS == 5
  ARR5(allocation_vector[0],
       allocation_vector[1],
       allocation_vector[2],
       allocation_vector[3],
       allocation_vector[4],
       finish_signal_vector[0],
       finish_signal_vector[1],
       finish_signal_vector[2],
       finish_signal_vector[3],
       finish_signal_vector[4],
       error_flag) = s;
  
#endif
  
  /* add to state_tree */ 
  
  /* Link next_state to the list of states */
  if(state_tree == NULL) {
    state_tree = s;
    state_tree_tail = state_tree;
  } else {
    state_tree_tail->next = s;
    state_tree_tail = state_tree_tail->next;
  }
  
  num_states++;
  
  /* if((num_states % 10000000) == 0) { */
  /*   setlocale(LC_NUMERIC, ""); */
    
  /*   printf( "%'lu\n", num_states); */
  /* } */
  
  return s;
  
}

void add_state_action_pair(state *s, state_action_pair* sap) {
  
  /* assert(state_action_pairs_list_tail != NULL); */
  sap->index = num_state_action_pairs++;
  
  /* add sap to state s's admissible_state_action_pairs_list */
  /* state_action_pair_pointer *sap_pointer = */
  /*   malloc(sizeof(state_action_pair_pointer)); */
  
  /* sap_pointer->sap = sap; */
  sap->next = s->admissible_state_action_pairs_list;
  s->admissible_state_action_pairs_list = sap;
  s->num_admissible_saps++;
  
}

/* void add_admissible_state_action_pair(state* s, state_action_pair* sap) { */
/*   state_action_pair_pointer *sap_pointer = */
/*     malloc(sizeof(state_action_pair_pointer)); */
  
/*   sap_pointer->sap = sap; */
/*   sap_pointer->next = s->admissible_state_action_pairs_list; */
/*   s->admissible_state_action_pairs_list = sap_pointer; */
  
/* } */

void add_in_state_action_pair_to_state(state *s,
				       ul sap_index,
				       double prob_s_given_sap) {
  
  state_action_pair_pointer *in_state_action_pair_ptr =
    malloc(sizeof(state_action_pair_pointer));

  in_state_action_pair_ptr->sap_index = sap_index;
  
  in_state_action_pair_ptr->prob_parent_state_given_sap = prob_s_given_sap;
  
  in_state_action_pair_ptr->next = s->in_state_action_pairs_list;
  
  s->in_state_action_pairs_list = in_state_action_pair_ptr;
  /* s->num_total_saps_with_nonzero_coefficients++; */
  
}

/* void add_out_state_to_state_action_pair(state_action_pair *sap, */
/* 					state* out_state, */
/* 					double transition_cost, */
/* 					double prob_out_state_given_sap) { */
  
/*   out_state_pointer *out_state_container = malloc(sizeof(out_state_pointer)); */
/*   out_state_container->s = out_state; */
/*   out_state_container->prob_state_given_containing_sap = */
/*     prob_out_state_given_sap; */
  
/*   out_state_container->cost_state_given_sap = transition_cost; */
  
/*   out_state_container->next = sap->out_states_list; */
/*   sap->out_states_list = out_state_container;  */
  
/* } */

bool any_job_just_missed_deadline(int criticality,
				  int time,
				  char *arr_finish_signal_vector) {
  
  int num_jobs = arr_num_jobs_of_criticality[criticality];
  
  for(int index = 0; index < num_jobs; ++index) {
    int job_index = arr_job_indexes_per_criticality[criticality][index];
    job *j = &arr_jobs[job_index];
    if((time == j->deadline) &&
       (arr_finish_signal_vector[job_index] == NOT_FINISHED)){ 
      
      return true;
    }
  }
  return false;
}

ul get_sum_allocations(int criticality, ul *allocation_vector) {
  
  ul sum_allocations = 0;
  int num_jobs = arr_num_jobs_of_criticality[criticality];
  
  for(int index = 0; index < num_jobs; ++index) {
    int job_index = arr_job_indexes_per_criticality[criticality][index];
    sum_allocations += allocation_vector[job_index];
  }
  
  return sum_allocations;
}

char is_admissible_action(ul *allocation_vector,
			  char *finish_signal_vector,
			  int action) {
  
  if(finish_signal_vector[action] == FINISHED) {return false;}
  
  int criticality = arr_jobs[action].criticality;
  if(criticality == HI) { return true; }
  
  char num_hi_jobs = arr_num_jobs_of_criticality[HI];
  
  for(int index = 0; index < num_hi_jobs; ++index) {
    int hi_job_index = arr_job_indexes_per_criticality[HI][index];
    job *hi_job = &arr_jobs[hi_job_index];
    if ((allocation_vector[hi_job_index] >=
	 hi_job->arr_wcets_at_criticalities[LO]) &&
	(finish_signal_vector[hi_job_index] == NOT_FINISHED)) {
      return false;
    }
  }
  return true;
}

/* int states_equal(state *s1, state *s2) { */
/*   /\* printf( "begin state equals\n"); *\/ */
/*   for(int index = 0; index < NUM_JOBS; ++index) { */
/*     if(s1->total_allocations[index] != s2->total_allocations[index] || */
/*        s1->finish_signals[index] != s2->finish_signals[index]) { */
/*       /\* printf( "end state equals\n"); *\/ */
/*       return false; */
/*     } */
/*   } */
/*   /\* printf( "end state equals\n"); *\/ */
/*   return (s1->error_flag == s2->error_flag);  */
/* } */


/* returns the criticality level realization of the system as */
/* given by the state. If the criticality level realization cannot */
/* be determined yet from the state, it will return UNKNOWN  */
int get_system_criticality_level_realization(ul   *arr_allocation,
					     char *arr_finish_signals) {
  
  int all_hi_jobs_finished_at_lo_demands = true;
  
  int num_hi_jobs = arr_num_jobs_of_criticality[HI];
  
  for(int index = 0; index < num_hi_jobs; ++index) {
    
    int hi_job_index = arr_job_indexes_per_criticality[HI][index];
    
    
    job *hi_job = &arr_jobs[hi_job_index];
    
    if(arr_allocation[hi_job_index] >
       hi_job->arr_wcets_at_criticalities[LO]) {
      return HI;
    }
    
    if ((arr_allocation[hi_job_index] >=
	 hi_job->arr_wcets_at_criticalities[LO]) &&
	(arr_finish_signals[hi_job_index] == NOT_FINISHED)) {
      return HI;
    }
    
    
    if(!((arr_allocation[hi_job_index] <=
	  hi_job->arr_wcets_at_criticalities[LO]) &&
	 (arr_finish_signals[hi_job_index] == FINISHED))) {
      
      all_hi_jobs_finished_at_lo_demands = false;
    }
    
  }
  
  if(all_hi_jobs_finished_at_lo_demands) {
    return LO;
  }
  
  return UNKNOWN;
  
}

bool have_all_jobs_finished(char *finish_signal_vector) {
  for(int index = 0; index < NUM_JOBS; ++index) {
    if(finish_signal_vector[index] == NOT_FINISHED) { return false; }
  }
  return true;
}

/* Computes Q(next_state | state, action), where next_state is specified by 
   giving next_allocation_vector, next_finish_signal_vector,
   and next_error_flag, and (previous) state is given by allocation_vector,
   finish_signal_vector, and error_flag */
/* Assumes: tha state transition is valid and action != NO_ACTION */
double compute_transition_probability(ul *prev_allocation_vector,
				      char* prev_finish_signal_vector,
				      int prev_error_flag,
				      ul *next_allocation_vector,
				      char* next_finish_signal_vector,
				      int next_error_flag,
				      int action) {
  
  assert(action != NO_ACTION);
  
  double prob = 0;
  
  char prev_finish_signal = prev_finish_signal_vector[action];
  assert(prev_finish_signal == NOT_FINISHED);

  char next_finish_signal = next_finish_signal_vector[action];
  ul prev_allocation = prev_allocation_vector[action];
  ul next_allocation  = next_allocation_vector[action];
  
  assert((next_allocation == prev_allocation + 1));
  
  job *j = &arr_jobs[action];
  int max_wcet = j->arr_wcets_at_criticalities[j->criticality];
  
  assert(next_allocation <= max_wcet);
  
  double denom = 0;
  int allocation;
  for(allocation = next_allocation; allocation <= max_wcet; ++allocation) {
    denom += j->execution_time_pmf[allocation - 1];
  }
  double prob_finished = j->execution_time_pmf[next_allocation - 1]/denom;
  
  if(next_finish_signal == FINISHED) {
    
    return prob_finished;
    
  } else {
    
    return ((double)1 - prob_finished);
    
  }
  
}


/* double get_transition_cost(ul *prev_allocation_vector, */
/* 			   char* prev_finish_signal_vector, */
/* 			   char prev_error_flag, */
/* 			   ul *next_allocation_vector, */
/* 			   char* next_finish_signal_vector, */
/* 			   char next_error_flag, */
/* 			   int action) { */



/* } */

int get_next_error_flag(int previous_error_flag,
			ul next_time,
			ul *next_allocation_vector,
			char *next_finish_signal_vector) {
  
  int sys_criticality =
    get_system_criticality_level_realization
    (next_allocation_vector, next_finish_signal_vector);
  
  bool lo_job_just_missed_deadline =
    any_job_just_missed_deadline
    (LO, next_time, next_finish_signal_vector);
  
  bool hi_job_just_missed_deadline =
    any_job_just_missed_deadline
    (HI, next_time, next_finish_signal_vector);

  /* if(hi_job_just_missed_deadline) { */
  /*   printf( "Hi job just missed deadline!\n"); */
  /* } */
  
  /* There can only be one next_error_flag per 
     state-action-finish_signal*/
  switch (previous_error_flag) {
    
  case ERROR:
    /* Next error flag can only be ERROR_ABSORB */
    return ERROR_ABSORB;
    break;
  case ERROR_ABSORB:
    /* Next error flag can only be ERROR_ABSORB */
    return ERROR_ABSORB;
    break;
  case NOT_ERROR:
    
    /* Next error flag can be any {NOT_ERROR, POTENTIAL_ERROR, ERROR} */
    
    /* Check for POTENTIAL_ERROR */
    if(sys_criticality == UNKNOWN &&
       lo_job_just_missed_deadline &&
       !hi_job_just_missed_deadline) {
      
      return POTENTIAL_ERROR;
      
    }
    
    /* Check for ERROR */
    /* NOT_ERROR -> ERROR happens if either a HI criticality job just
       misses its deadline (i.e., at next_time), or a LO criticality job 
       misses its deadline and the system criticality level realization 
       at next_time is LO */
    
    if(hi_job_just_missed_deadline ||
       (sys_criticality == LO && lo_job_just_missed_deadline)) {
      
      /* ERROR */
      return ERROR;
      
    }
    
    /* Neither POTENTIAL_ERROR nor ERROR. Then only NOT_ERROR */
    return NOT_ERROR;
    
    break;
    
  case POTENTIAL_ERROR:
    {
      /* Next error flag can be any {NOT_ERROR, POTENTIAL_ERROR, ERROR} */
      
      /* Check for ERROR */
      /* POTENTIAL_ERROR -> ERROR happens if either the system 
	 criticality level is revealed as LO at next_time, or a HI 
	 criticality job just misses its deadline at next_time */
      
      if(hi_job_just_missed_deadline || (sys_criticality == LO)) {
	/* ERROR */
	return ERROR;
      }
      
      /* Check for POTENTIAL_ERROR */
      if((sys_criticality == UNKNOWN) && !hi_job_just_missed_deadline) {
	return POTENTIAL_ERROR;
      }
      
      /* Check for NOT_ERROR */
      /* POTENTIAL_ERROR -> NOT_ERROR happens if the system criticality 
	 level is realized at next_time (= time + 1) as HI and no 
	 HI criticality job misses its deadline at (time + 1) */
      if(sys_criticality == HI && !hi_job_just_missed_deadline) {
	
	return NOT_ERROR;
	
      }
      
      break;
    }
  }
  
}



/* root_state is added already whenever this function is called  */
/* root state is such that at least one job still requires execution
   (root_finish_signal_vector[i] = NOT_FINISHED for some job i) */
void expand_state_tree_at(state *root_state_ptr,
			  ul *root_allocation_vector,
			  char *root_finish_signal_vector,
			  char root_error_flag) {
  
  
  ul next_time = root_state_ptr->time + 1;
  
  /* for(int index = 0; index < NUM_JOBS; ++index) { */
  /*   next_time += root_allocation_vector[index]; */
  /* } */
  
  /* printf( "%lu\n", next_time); */
  
  assert(root_state_ptr != 0);
  
  ul *next_allocation_vector      = NULL;
  char *next_finish_signal_vector = NULL;
  
  for(int action_index = 0; action_index < NUM_JOBS; ++action_index) {
    
    if(!is_admissible_action(root_allocation_vector,
			     root_finish_signal_vector,
			     action_index)) {
      continue;
    }

    if(next_allocation_vector == NULL) {
      next_allocation_vector = clone_ul_array(root_allocation_vector, NUM_JOBS);
    }
    
    /* Then action "action_index" is a possible next action */
    
    state_action_pair *root_state_action_pair = create_state_action_pair();
    
    root_state_action_pair->action = action_index;
    /* root_state_action_pair->s = root_state_ptr; */
    
    /* Add (root_state, index) as a state-action pair 
       After this operation, and index is assigned to root_state_action_pair */
    add_state_action_pair(root_state_ptr, root_state_action_pair);
    
    next_allocation_vector[action_index]++;    
    
    job *current_job = &arr_jobs[action_index]; 
    
    int current_job_max_wcet =
      current_job->arr_wcets_at_criticalities[current_job->criticality];
    
    int num_possible_finish_signals_for_current_job;
    char arr_possible_finish_signals_for_current_job[2];
    
    if(current_job_max_wcet == next_allocation_vector[action_index]){
      num_possible_finish_signals_for_current_job = 1;
      arr_possible_finish_signals_for_current_job[0] = FINISHED;
    } else {
      assert(next_allocation_vector[action_index] < current_job_max_wcet);
      num_possible_finish_signals_for_current_job = 2;
      arr_possible_finish_signals_for_current_job[0] = NOT_FINISHED;
      arr_possible_finish_signals_for_current_job[1] = FINISHED;
    }
    
    
    for(int finish_signal_index = 0;
	finish_signal_index < num_possible_finish_signals_for_current_job;
	++finish_signal_index) {
      
      /* The system criticality level realization as determined by the 
	 vectors next_allocation_vector, next_finish_signal_vector */

      int next_finish_signal =
	arr_possible_finish_signals_for_current_job[finish_signal_index];
      
      if(next_finish_signal_vector == NULL) {
	next_finish_signal_vector =
	  clone_char_array(root_finish_signal_vector, NUM_JOBS);
      }
      
      next_finish_signal_vector[action_index] = next_finish_signal;
      
      int next_error_flag = get_next_error_flag(root_error_flag,
						 next_time,
						 next_allocation_vector,
						 next_finish_signal_vector);	
      state* next_state =
	lookup_state_ptr_by_params(next_allocation_vector,
				   next_finish_signal_vector,
				   next_error_flag);
      
      if(next_state == 0) { /* next_state does not exist yet */
	next_state = add_state(next_allocation_vector,
			       next_finish_signal_vector,
			       next_error_flag);
	
	next_state->time = next_time;
	
	next_state->constraint_cost = (next_error_flag == ERROR) ? 1 : 0;
	
	bool all_jobs_finished =
	  have_all_jobs_finished(next_finish_signal_vector);
	
	if(all_jobs_finished) {
	  /* next_state is terminal, do not expand it */
	  
	  state_action_pair *terminal_state_action_pair =
	    create_state_action_pair();
	  
	  terminal_state_action_pair->action = NO_ACTION;
	  /* terminal_state_action_pair->s = next_state; */
	  
	  /* Add (next_state, NO_ACTION) as a state-action pair */
	  add_state_action_pair(next_state, terminal_state_action_pair);
	} else {
	  
	  if(next_time < horizon) {
	    expand_state_tree_at(next_state,
				 next_allocation_vector,
				 next_finish_signal_vector,
				 next_error_flag);
	  }
	}
	
      }
      /* root_state_action_pair appears for the first time, so add it 
	 as an in_state_action_pair to next_state after computing Q */
      
      double prob_next_state_given_root_state_action_pair =
	compute_transition_probability(root_allocation_vector,
				       root_finish_signal_vector,
				       root_error_flag,
				       next_allocation_vector,
				       next_finish_signal_vector,
				       next_error_flag,
				       action_index);
      
      
      /* Add it only if Q(next_state | root_state_action_pair) > 0
	 (which is the coefficient of rho(root_state_action_pair) 
	 in the constraint corresponding to next_state)
	 and root_state is such that not all jobs finished execution (i.e., 
	 root_state is in the transient set S'). Since we only expand 
	 at a state if not all jobs finished execution, root_state is 
	 readily in S' so no need to check for that condition. */
      if(definitely_greater_than
	 (prob_next_state_given_root_state_action_pair, 0)){
	
	add_in_state_action_pair_to_state
	  (next_state,
	   root_state_action_pair->index,
	   prob_next_state_given_root_state_action_pair);
	
	/* even if next_state exists, root_state_action_pair appears for 
	   the first time here, so add next_state as an out_state to 
	   root_state_action_pair */
	
	int root_state_sys_criticality =
	  get_system_criticality_level_realization
	  (root_allocation_vector, root_finish_signal_vector);
	
	int next_state_sys_criticality =
	  get_system_criticality_level_realization
	  (next_allocation_vector, next_finish_signal_vector);
	
	if(root_state_sys_criticality == UNKNOWN &&
	   next_state_sys_criticality == HI) {
	  
	  root_state_action_pair->objective_cost +=
	    (double)get_sum_allocations(LO, next_allocation_vector) *
	    prob_next_state_given_root_state_action_pair;
	}
	
      } else {
	/* printf( "Q: %lf\n", prob_next_state_given_root_state_action_pair); */
      }
      
      
      /* double transition_cost = get_transition_cost(root_allocation_vector, */
      /* 					     root_finish_signal_vector, */
      /* 					     root_error_flag, */
	/* 					     next_allocation_vector, */
	/* 					     next_finish_signal_vector, */
	/* 					     next_error_flag, */
	/* 					     action_index); */
      
      
      
	/* add_out_state_to_state_action_pair */
	/*   (root_state_action_pair, */
	/*    next_state, */
	/*    transition_cost, */
	/*    prob_next_state_given_root_state_action_pair); */
      
      next_finish_signal_vector[action_index] =
	root_finish_signal_vector[action_index];
    }
    next_allocation_vector[action_index] =
      root_allocation_vector[action_index];
  }
  
  if(next_finish_signal_vector != NULL){ free(next_finish_signal_vector); }
  if(next_allocation_vector    != NULL){ free(next_allocation_vector); }
}
 

void build_state_tree() {
  
  state *initial_state_ptr = add_state(initial_allocation_vector,
				       initial_finish_signal_vector,
				       initial_error_flag);

  assert(initial_state_ptr != NULL);

  initial_state_ptr->is_initial = true;
  
  expand_state_tree_at(initial_state_ptr,
		       initial_allocation_vector,
		       initial_finish_signal_vector,
		       initial_error_flag);
  
}
double min(double a, double b) {
  if(definitely_greater_than(a,b)) { return b; }
  return a;
}
double max(double a, double b) {
  if(definitely_greater_than(a,b)) { return a; }
  return b;
}
double get_prob_job_demands_at_most(job *j, int demand_upper_bound) {
  
  double prob = 0;
  int max_demand = j->arr_wcets_at_criticalities[j->criticality];
  demand_upper_bound = min(demand_upper_bound, max_demand);
  for(int demand = 1; demand <= demand_upper_bound; ++demand) {
    prob += j->execution_time_pmf[demand - 1];
  }
  return prob;
  
}

typedef struct int_list_tag {
  int index;
  struct int_list_tag* next;
  struct int_list_tag* prev; 
} int_list;

/* returns whether or not the instance is OCBP-Schedulable, and if so,
   sets arr_jobs[i].job->priority according to OCBP */
bool compute_OCBP_priorities() {
  
  int priority = 0;
  int_list* job_indexes_list = NULL;
  int_list* job_index = NULL;
  
  for (int index = 0; index < NUM_JOBS; ++index) {
    job_index = malloc(sizeof(int_list));
    job_index->index = index;
    job_index->next = job_indexes_list;
    job_index->prev = NULL;
    if(job_indexes_list != NULL) {
      job_indexes_list->prev = job_index;
    } 
    job_indexes_list = job_index;
    
  }
  
  char is_lowest_priority_job_found = false;
  int num_remaining_jobs = NUM_JOBS;
  
  while(num_remaining_jobs > 0) {
    
    is_lowest_priority_job_found = false;

    int_list *job_indexes_list_tail = job_indexes_list;
    
    while(job_indexes_list_tail != NULL) {
    
      int sum_exec_times = 0;
      job *lowest_priority_job = &arr_jobs[job_indexes_list_tail->index];
      
      int_list *job_indexes_list_tail2 = job_indexes_list;
      
      while(job_indexes_list_tail2 != NULL) {
	
	job *j = &arr_jobs[job_indexes_list_tail2->index];
	if(lowest_priority_job->criticality < j->criticality) {
	  sum_exec_times +=
	    j->arr_wcets_at_criticalities[lowest_priority_job->criticality];
	} else {
	  sum_exec_times += j->arr_wcets_at_criticalities[j->criticality];
	}
	job_indexes_list_tail2 = job_indexes_list_tail2->next;
      }
      
      if(sum_exec_times <= lowest_priority_job->deadline) {
	lowest_priority_job->priority = priority;
	++priority;
	// delete the job_indexes_list_tail
	if(job_indexes_list_tail->prev != NULL) { /* middle or last */
	  (job_indexes_list_tail->prev)->next = job_indexes_list_tail->next;
	} else {
	  /* job_indexes_list_tail is the head  */
	  assert(job_indexes_list == job_indexes_list_tail);
	  job_indexes_list = job_indexes_list_tail->next;
	}
	if(job_indexes_list_tail->next != NULL) { /* not last element */
	  (job_indexes_list_tail->next)->prev = job_indexes_list_tail->prev;
	}
	job_indexes_list_tail->next = NULL;
	job_indexes_list_tail->prev = NULL;
	free(job_indexes_list_tail);
	job_indexes_list_tail = NULL;
	is_lowest_priority_job_found = true;
	--num_remaining_jobs;
	break;
      } else {
	is_lowest_priority_job_found = false;
	job_indexes_list_tail = job_indexes_list_tail->next;
      }
    }
   
    if(!is_lowest_priority_job_found) {
      break;
    }
    
  }
  
  int_list *job_indexes_list_tail = job_indexes_list;
  while(job_indexes_list_tail != NULL) {
    job_indexes_list_tail = job_indexes_list_tail->next;
    free(job_indexes_list);
    job_indexes_list = job_indexes_list_tail;
  }
  
  return is_lowest_priority_job_found;
}

void print_exec_time_sample(int *arr_exec_time_sample) {

  printf("Job demands sample:     ");
  for(int index = 0; index < NUM_JOBS; ++index) {
    printf("J%d: %d", index, arr_exec_time_sample[index]);
  }
  printf("\n");
}


/* void print_state(ul time, */
/* 		 ul *allocation_vector, */
/* 		 char *finish_signal_vector, */
/* 		 int error_flag) { */
  
/*   printf("time: %lu,    x:[", time); */
/*   int index; */
/*   for(index = 0; index < NUM_JOBS; ++index) { */
/*     if(index < NUM_JOBS - 1) { */
/*       printf("%d   ", allocation_vector[index]); */
/*     } else { */
/*       printf("%d],    y:[", allocation_vector[index]); */
/*     } */
    
/*   } */
  
  
/*   for(index = 0; index < NUM_JOBS; ++index) { */
/*     if(index < NUM_JOBS - 1) { */
/*       printf( "%d  ", finish_signal_vector[index]); */
/*     } else { */
/*       printf( "%d],    r: ", finish_signal_vector[index]); */
/*     } */
    
/*   } */

/*   char *error_str = "hello"; */
/*   if(error_flag == ERROR) { */
/*     error_str = "error"; */
/*   } else if(error_flag == NOT_ERROR) { */
/*     error_str = "not error"; */
/*   } else if(error_flag == POTENTIAL_ERROR) { */
/*     error_str = "potential error"; */
/*   } else if(error_flag == ERROR_ABSORB) { */
/*     error_str = "error\'"; */
/*   } */
  
/*   printf( "%s, ", error_str); */

/*   int crit_level = */
/*     get_system_criticality_level_realization(allocation_vector, */
/* 					     finish_signal_vector); */
/*   printf( "crit: %d\n", crit_level); */
/* } */

char *get_criticality_string(int criticality) {
  if(criticality == LO) {
    return "LO";
  } else if(criticality == HI) {
    return "HI";
  } else {
    assert(false);
  }
}

void print_job_set() {
  
  printf( "\nNumber of jobs: %d\n", NUM_JOBS);
  fflush(stdout);
  
  int index = 0;
  for(index = 0; index < NUM_JOBS; ++index) {
    job *j = &arr_jobs[index];
    printf( "J%d   deadline: %d, criticality: %s, WCET(LO): %d, WCET(HI): %d\n",
	   index,
	   j->deadline,
	   get_criticality_string(j->criticality),
	   j->arr_wcets_at_criticalities[LO],
	   j->arr_wcets_at_criticalities[j->criticality]);
    fflush(stdout);
  }
  
  
}

int main(int argc, char *argv[]){
  
  num_state_action_pairs = 0;
  
  double epsilon_lo = 0.4;
  double epsilon_hi = 0.3;
  
  /* If job is HI criticality, must always maintain:
     arr_wcets_at_criticalities[LO] < arr_wcets_at_criticalities[HI] */

  /* Input jobset */
  
  job *j = &arr_jobs[0];
  j->deadline    = 60;
  j->criticality = LO;
  j->arr_wcets_at_criticalities[LO] = 50;
  /* j->arr_wcets_at_criticalities[HI] = 50; */
  j->execution_time_pmf =
    generate_execution_time_pmf(j->arr_wcets_at_criticalities[j->criticality]);
  
  j = &arr_jobs[1];
  j->deadline    = 140;
  j->criticality = LO;
  j->arr_wcets_at_criticalities[LO] = 75;
  /* No need to set arr_wcets_at_criticalities[HI] in case job is LO 
     criticality. Any other entries of arr_wcets_at_criticalities will not be
     looked at in the code in this case. */
  j->execution_time_pmf =
    generate_execution_time_pmf(j->arr_wcets_at_criticalities[j->criticality]);
  
  j = &arr_jobs[2];
  j->deadline    = 190;
  j->criticality = HI;
  j->arr_wcets_at_criticalities[LO] = 80;
  j->arr_wcets_at_criticalities[HI] = 81;
  j->execution_time_pmf =
    generate_execution_time_pmf(j->arr_wcets_at_criticalities[j->criticality]);
  
  j = &arr_jobs[3];
  j->deadline    = 70;
  j->criticality = HI;
  j->arr_wcets_at_criticalities[LO]        = 8;
  j->arr_wcets_at_criticalities[HI]        = 25;
  j->execution_time_pmf =
    generate_execution_time_pmf(j->arr_wcets_at_criticalities[j->criticality]);
  

  
  j = &arr_jobs[4];
  j->deadline    = 300;
  j->criticality = LO;
  j->arr_wcets_at_criticalities[LO]        = 20;
  j->execution_time_pmf =
    generate_execution_time_pmf(j->arr_wcets_at_criticalities[j->criticality]);

  print_job_set(); 
  
  bool is_ocbp_schedulable = compute_OCBP_priorities();
  
  /* Allocating memory for an array of state structs, and NOT pointers to 
     state */
  
  initial_allocation_vector     = calloc(NUM_JOBS, sizeof(ul));
  initial_finish_signal_vector  = calloc(NUM_JOBS, sizeof(char));
  initial_error_flag            = NOT_ERROR;
  
  for(int crit_index = 0; crit_index < NUM_CRITICALITY_LEVELS; ++crit_index){
    arr_num_jobs_of_criticality[crit_index] = 0;
  }
  
  ul state_array_size = 1;
  
  horizon = 0;
  int index = 0;

  // Initialize random number generator:
  gsl_rng_env_setup();
  gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
  gsl_rng_set(rng, time(NULL));
  // Allocate workspace, and precompute tables:
  /* printf(  "Precomputing Alias sampling tables ...\n" ); */

  ransampl_ws *arr_exec_times_ws[NUM_JOBS];
  
  for(index = 0; index < NUM_JOBS; ++index) {
    
    initial_finish_signal_vector[index] = NOT_FINISHED;
    
    job * j = &arr_jobs[index];
    ul max_wcet = (ul)j->arr_wcets_at_criticalities[j->criticality];
    horizon += max_wcet;
    dims[index] = max_wcet + 1;
    state_array_size *= dims[index];
    arr_num_jobs_of_criticality[j->criticality] += 1;

    arr_exec_times_ws[index] = ransampl_alloc( max_wcet );
    ransampl_set( arr_exec_times_ws[index], j->execution_time_pmf );
    
  }
  
  for(index = NUM_JOBS; index < 2*NUM_JOBS; ++index) {
    dims[index] = MAX_FINISH_SIGNAL;
    state_array_size *= dims[index];
  }
  dims[2 * NUM_JOBS] = MAX_ERROR_FLAG;
  
  state_array_size *= dims[2 * NUM_JOBS];
  
  arr_states = calloc(state_array_size, sizeof(state*));  
    
  /* arr_jobs_per_criticality = */
  /*   malloc(NUM_CRITICALITY_LEVELS * sizeof(job *)); */
  
  for(int crit_index = 0; crit_index < NUM_CRITICALITY_LEVELS; ++crit_index){
    
    arr_job_indexes_per_criticality[crit_index] =
      malloc(arr_num_jobs_of_criticality[crit_index] * sizeof(int));
  }
  
  char *arr_num_jobs_of_criticality_copy =
    clone_char_array(arr_num_jobs_of_criticality, 2);
  
  /* printf( "%d\n", arr_num_jobs_of_criticality_copy[HI]);  */
  
  for(int job_index = 0; job_index < NUM_JOBS; ++job_index){
    
    job *j = &arr_jobs[job_index];
    int criticality = j->criticality;
    arr_job_indexes_per_criticality[criticality]
      [arr_num_jobs_of_criticality_copy[criticality] - 1] = job_index;
    --arr_num_jobs_of_criticality_copy[criticality];
    
  }
  
  free(arr_num_jobs_of_criticality_copy);
  
  printf( "Creating state tables ...\n");
  fflush(stdout);
  build_state_tree();

  setlocale(LC_NUMERIC, "");
  printf( "Building LP\nNumber of variables: %'lu\nNumber of constraints: %'lu\n", num_state_action_pairs, (num_states + 1));
  fflush(stdout);
  
  /**--------------------Begin Building and Solving LP----------------------**/
  GRBenv   *master_env   = NULL;
  GRBmodel *model = NULL;
  int       error = 0;
  
  double    *arr_solution = calloc(num_state_action_pairs, sizeof(double));

  // <NOTE>
  // By BNA
  // ======
  // state_action_pair->index corresponds to its index in arr_solution
  // and arr_objective_coefficients. arr_variable_indexes contains the
  // state_action_pair indexes that have nonzero coefficients in the
  // constraint being constructed.
  // </NOTE>
  int *arr_variable_indexes = calloc(num_state_action_pairs, sizeof(int));
  
  double *arr_constraint_coefficients = calloc(num_state_action_pairs,
					       sizeof(int));
  double *arr_objective_coefficients = calloc(num_state_action_pairs,
					      sizeof(double));
  int       optimization_status;
  double    objective_value;
  
  /* Create environment */
  
  error = GRBloadenv(&master_env, "lp.log");
  if (error) goto QUIT;

  /* Create an empty model */

  error = GRBnewmodel(master_env,&model, "lp", 0, NULL, NULL, NULL, NULL, NULL);
  if (error) goto QUIT;
  
  printf( "Creating objective coefficients array ... \n");
  fflush(stdout);
  
  /* Add variables */
  state_tree_tail = state_tree;
  while(state_tree_tail != NULL) {
    
    state_action_pair *sap_list =
      state_tree_tail->admissible_state_action_pairs_list;
    
    while(sap_list != NULL) {
      arr_objective_coefficients[sap_list->index] =
	sap_list->objective_cost;

      /* error = GRBaddvar(model, 0, NULL, NULL, sap_list->objective_cost, */
      /* 			0.0, GRB_INFINITY, GRB_CONTINUOUS, NULL); */
      /* if (error) goto QUIT; */
      
      /* if(definitely_greater_than(sap_list->objective_cost, 0)) { */
      /* 	printf( "objective cost: %lf\n", sap_list->objective_cost); */
      /* } */
      sap_list = sap_list->next;
    } 
    
    state_tree_tail = state_tree_tail->next;
  } 

  printf( "Adding variables ... \n");
  fflush(stdout);
  
  error = GRBXaddvars(model, num_state_action_pairs, 0, NULL, NULL, NULL,
  		     arr_objective_coefficients, NULL, NULL, NULL, NULL);
  if (error) goto QUIT;
  
  /* Change objective sense to minimization */
  
  /* error = GRBsetdblattr(model, GRB_DBL_ATTR_MODELSENSE, GRB_MINIMIZE); */
  /* if (error) goto QUIT; */
  
  GRBenv *model_env = GRBgetenv(model);
  
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_THREADS, 10); */
  /* if (error) goto QUIT; */
  
  error = GRBsetintparam(model_env, GRB_INT_PAR_CROSSOVER, 0);
  if (error) goto QUIT;

  /* GRB_INT_PAR_METHOD = 0: Primal Simplex */
  /*                      1: Dual Simplex */
  /*                      2: Barrier */
  /*                      3: Nondeterministic Concurrent */
  /*                      4: Dterministic Concurrent */
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_METHOD, 0); */
			 
  if (error) goto QUIT;
  
  /* Integrate new variables */
  
  error = GRBupdatemodel(model);
  if (error) goto QUIT;

  printf( "Creating constraints ... \n");
  fflush(stdout);
  /* Create the constraints */

  ul num_nonzero_constraint_coefficients = 0;
  double constraint_rhs = 0;
  
  state_tree_tail = state_tree;
  /* int num_constraints = 0; */

  /* Add the constraints per state, for all states */
  while(state_tree_tail != NULL) {
    
    state_action_pair *admissible_saps_list_tail =
      state_tree_tail->admissible_state_action_pairs_list;
    
    num_nonzero_constraint_coefficients = 0;
    
    while(admissible_saps_list_tail != NULL) {
      
      arr_variable_indexes[num_nonzero_constraint_coefficients] =
	admissible_saps_list_tail->index;
      
      arr_constraint_coefficients[num_nonzero_constraint_coefficients] = 1;
      ++num_nonzero_constraint_coefficients;
      admissible_saps_list_tail = admissible_saps_list_tail->next;
      
      
    } 
    
    
    state_action_pair_pointer *in_sap_ptrs_list_tail =
      state_tree_tail->in_state_action_pairs_list;
    
    while(in_sap_ptrs_list_tail != NULL) {
      
      arr_variable_indexes[num_nonzero_constraint_coefficients] =
	in_sap_ptrs_list_tail->sap_index;
      
      arr_constraint_coefficients[num_nonzero_constraint_coefficients] =
	-1 * in_sap_ptrs_list_tail->prob_parent_state_given_sap;
      
      ++num_nonzero_constraint_coefficients;
      in_sap_ptrs_list_tail = in_sap_ptrs_list_tail->next;
      
    } 

    constraint_rhs = state_tree_tail->is_initial ? 1 : 0; 
    
    error = GRBaddconstr(model,
			 num_nonzero_constraint_coefficients,
			 arr_variable_indexes,
			 arr_constraint_coefficients,
			 GRB_EQUAL,
			 constraint_rhs,
			 NULL);
    if (error) goto QUIT;

    /* ++num_constraints; */
    /* printf( "Added constraint %d\n", num_constraints); */
    
    state_tree_tail = state_tree_tail->next;
  } 
  
  
  /* Compute delta */
  double prob_sys_criticality_lo = 1;
  double prob_sys_criticality_hi = 1;
  
  for(int job_index = 0; job_index < NUM_JOBS; ++job_index) {
    
    job *j = &arr_jobs[job_index];
    int lo_demand  = j->arr_wcets_at_criticalities[LO];
    
    prob_sys_criticality_lo *= get_prob_job_demands_at_most(j, lo_demand);
    
  }
  
  prob_sys_criticality_hi = 1 - prob_sys_criticality_lo;
  
  double delta = min((prob_sys_criticality_lo * epsilon_lo),
		     (prob_sys_criticality_hi * epsilon_hi));

  printf( "epsilon_lo %lf\n", epsilon_lo); 
  printf( "epsilon_hi %lf\n", epsilon_hi); 
  
  printf( "prob_sys_criticality_lo %lf\n", prob_sys_criticality_lo); 
  printf( "prob_sys_criticality_hi %lf\n", prob_sys_criticality_hi); 
  printf( "Delta %lf\n", delta); 
  fflush(stdout);
  
  /* Add the Delta constraint */
  num_nonzero_constraint_coefficients = 0;
  state_tree_tail = state_tree;
  while(state_tree_tail != NULL) {
    
    double sap_constraint_cost = state_tree_tail->constraint_cost;
    
    if(!essentially_equal(sap_constraint_cost, 0)) {
      
      state_action_pair *admissible_saps_list_tail =
	state_tree_tail->admissible_state_action_pairs_list;
      
      while(admissible_saps_list_tail != NULL) {
	
	arr_variable_indexes[num_nonzero_constraint_coefficients] =
	  admissible_saps_list_tail->index;
      
	arr_constraint_coefficients[num_nonzero_constraint_coefficients] =
	  sap_constraint_cost;
	++num_nonzero_constraint_coefficients;
	admissible_saps_list_tail = admissible_saps_list_tail->next;
	
      } 
      
      
    }
    state_tree_tail = state_tree_tail->next;
  }
  /* printf( "%lu\n", num_nonzero_constraint_coefficients); */
  error = GRBaddconstr(model,
  		       num_nonzero_constraint_coefficients,
  		       arr_variable_indexes,
  		       arr_constraint_coefficients,
  		       GRB_LESS_EQUAL,
  		       delta,
  		       NULL);
  
  /* Optimize model */

  printf( "Solving LP ... \n");
  fflush(stdout);
  
  error = GRBoptimize(model);
  if (error) goto QUIT;
  
  /* Write model to 'lp.lp' */
  
  error = GRBwrite(model, "lp.lp");
  if (error) goto QUIT;
  
  /* Capture solution information */
  
  error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimization_status);
  if (error) goto QUIT;

  if (optimization_status == GRB_OPTIMAL) {
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objective_value);
    if (error) goto QUIT;
    printf( "Optimal objective: %.4e\n", objective_value);
    fflush(stdout);
    
    error = GRBgetdblattrarray(model,
			       GRB_DBL_ATTR_X,
			       0,
			       num_state_action_pairs,
			       arr_solution);
    if (error) goto QUIT;
    
    /* printf( "  x=%.0f, y=%.0f, z=%.0f\n", sol[0], sol[1], sol[2]); */
  } else if (optimization_status == GRB_INF_OR_UNBD) {
    printf( "Model is infeasible or unbounded\n");
    fflush(stdout);
    goto QUIT;
    
  } else {
    printf( "Optimization was stopped early\n");
    fflush(stdout);
    goto QUIT;
  }
  
  /**--------------------End Building and Solving LP------------------------**/

    
  /**---------------------------Begin Simulation------------------------**/


  const int num_sim_samples = 100000;

  sim_state_t *sim_state = malloc(sizeof(sim_state_t));
  sim_stats_t *sim_stats  = malloc(sizeof(sim_stats_t));

  sim_stats->arr_num_deadline_misses_per_job =
    calloc(NUM_JOBS, sizeof(int));
  
  /* sim_stats->arr_num_errors_per_criticality = */
  /*   calloc(NUM_CRITICALITY_LEVELS, sizeof(int)); */
  
  sim_stats->num_errors = 0;
  
  printf( "\n\nSimulating execution over %d demand samples (vectors) \n",
	 num_sim_samples);
  fflush(stdout);
  
  int sample_index;
  
  for(sample_index = 0; sample_index < num_sim_samples; ++sample_index) {
    
    int arr_exec_time_samples[NUM_JOBS];
    for(int index = 0; index < NUM_JOBS; ++index) {
      /* ransampl_draw(...) generates an index in {0, ..., max_wcet - 1} */
      arr_exec_time_samples[index] =
	ransampl_draw( arr_exec_times_ws[index],
		       gsl_rng_uniform(rng),
		       gsl_rng_uniform(rng) ) + 1;
      
    }

    /* print_exec_time_sample(arr_exec_time_samples); */
    
    /* printf( "sample ok ... \n"); */
    
    sim_state->error_flag	    = initial_error_flag;
    sim_state->allocation_vector    =
      clone_ul_array(initial_allocation_vector, NUM_JOBS);
    sim_state->finish_signal_vector =
      clone_char_array(initial_finish_signal_vector, NUM_JOBS);
    
    state *s;
    ul sim_time = 0;
    
    bool *arr_job_missed_deadlines = calloc(NUM_JOBS, sizeof(bool));
    
    /* simulation loop */
    while(!have_all_jobs_finished(sim_state->finish_signal_vector)) {
      
      
      
      /* print_state(sim_time, */
      /* 		  sim_state->allocation_vector, */
      /* 		  sim_state->finish_signal_vector, */
      /* 		  sim_state->error_flag); */
      
      /* Lookup all the admissible state-action pairs for state s */
      s = lookup_state_ptr_by_params(sim_state->allocation_vector,
				     sim_state->finish_signal_vector,
				     sim_state->error_flag);
      
      assert(s != NULL);
      
      state_action_pair *admissible_saps_list_tail =
	s->admissible_state_action_pairs_list;
      
      double *sap_pmf = malloc((s->num_admissible_saps) * sizeof(double));
      state_action_pair **admissible_sap_pointers =
	malloc((s->num_admissible_saps) * sizeof(state_action_pair*));
      
      int index = 0;
      double sap_lp_solution = 0;
      double sum_admissible_saps_lp_solutions = 0;
      
      while(admissible_saps_list_tail != NULL) {
	admissible_sap_pointers[index] = admissible_saps_list_tail;
	
	sap_lp_solution = arr_solution[admissible_saps_list_tail->index];
	
	/* printf( "sap solution: %.*lf\n", 200, sap_lp_solution); */
	
	sap_pmf[index] = sap_lp_solution;
	sum_admissible_saps_lp_solutions += sap_lp_solution;
	++index;
	admissible_saps_list_tail = admissible_saps_list_tail->next;
      }

      /* sum_admissible_saps_lp_solutions can be zero. If so, choose the action 
      /* arbitrarily. Choose the first action */

      int action = 0;
      if(definitely_greater_than(sum_admissible_saps_lp_solutions, 0)) {
	
	/* Normalize sap_pmf */
	for(index = 0; index < s->num_admissible_saps; ++index) {
	  /* printf( "sap: %lf \n", sap_pmf[index]); */
	  sap_pmf[index] = (sap_pmf[index]/sum_admissible_saps_lp_solutions);
	  /* printf( "sap normalized: %lf \n", sap_pmf[index]); */
	  assert(definitely_less_than_or_equal(sap_pmf[index], 1));
	}
	/* printf( "ok\n"); */
	/* Now sample a state-action pair (sap) pointer from sap_pmf */
	/* action is the index of a job  */
	ransampl_ws *sap_index_ws = ransampl_alloc( s->num_admissible_saps );
	ransampl_set( sap_index_ws, sap_pmf );
	int sap_index_decision = ransampl_draw( sap_index_ws,
						gsl_rng_uniform(rng),
						gsl_rng_uniform(rng) );
	
	/* printf( "number of admissible actions: %d\n", s->num_admissible_saps); */
	
	action = admissible_sap_pointers[sap_index_decision]->action;

	ransampl_free( sap_index_ws );
	free(sap_pmf);
	free(admissible_sap_pointers);
	
      } else {
	
	action = admissible_sap_pointers[0]->action;
	
      }
      
      /* printf( "Job (%d) chosen at time (%lu) with probability %lf\n\n", */
      /* 	     action, sim_time, sap_pmf[sap_index_decision]); */
      
      /* printf( "after acting\n"); */
      
      assert(action != NO_ACTION);
      assert(sim_state->finish_signal_vector[action] != FINISHED);
      
      /* Update the state */
      sim_state->allocation_vector[action]++;
      
      if(sim_state->allocation_vector[action] ==
	 arr_exec_time_samples[action]) {
	sim_state->finish_signal_vector[action] = FINISHED;
      }
      
      ++sim_time;
      
      
	int e =
	  get_next_error_flag(sim_state->error_flag,
			      sim_time,
			      sim_state->allocation_vector,
			      sim_state->finish_signal_vector);
	
	sim_state->error_flag = e;

      /* Check if a job just missed its deadline */
      for(int job_index = 0; job_index < NUM_JOBS; ++job_index) {
	job *j = &arr_jobs[job_index];
	bool job_not_finished =
	  (sim_state->finish_signal_vector[job_index] == NOT_FINISHED);
	
	if((sim_time == j->deadline) && job_not_finished) {
	  assert(arr_job_missed_deadlines[job_index] == false);
	  sim_stats->arr_num_deadline_misses_per_job[job_index]++;
	  arr_job_missed_deadlines[job_index] = true;
	  if(j->criticality == HI) {
	    assert(sim_state->error_flag == ERROR);
	  }
	}
      }
      
      if(sim_state->error_flag == ERROR) {
	sim_stats->num_errors++;
      }
      
      
      
    } // end simulation per sample while loop
    
    free(sim_state->allocation_vector);
    free(sim_state->finish_signal_vector);
    free(arr_job_missed_deadlines);
    
  } //end looping over samples for loop
  
  
  /* Display simulation statistics */
  printf( "+=======================Simulation Stats=======================\n");
  printf( "+ Number of samples: %d\n", num_sim_samples);
  printf( "+ Number of deadline misses per job:\n     ", num_sim_samples);
  for(int job_index = 0; job_index < NUM_JOBS; ++job_index) {
    
    printf( "J%d: (%d/%d)   ",
	   job_index,
	   sim_stats->arr_num_deadline_misses_per_job[job_index],
	   num_sim_samples);
    
  }
  
  printf( "\n+ Total number of error executions: (%d/%d)\n",
	 sim_stats->num_errors,
	 num_sim_samples);
  printf( "+==============================================================\n");
  fflush(stdout);




  

    
  /*   const int M=1000000; */
  /*   int i, m; */
  
/*   // Discrete probability distribution example: */
/*   const int n = 9; */
/*   // states of Austria */
/*   const char* names[] = { */
/*     "Wien", "Niederoesterreich", "Oberoesterreich", "Tirol", */
/*     "Kaernten", "Salzburg", "Vorarlberg", "Burgenland", "Steiermark" }; */
/*   // inhabitants in millions as of 2011 according to www.statistik.at */
/*   double p[] = { 1.721573, 1.614661, 1.415020, .711161, */
/* 		 .558056, .532713, .370833, .285377, .1211506 }; */
  
  
    
/* #ifdef RANSAMPL_SHALL_REVEAL_ITS_INTERNALS */
/*     // Inspect tables: */
/*     printf(  "  %-3s  %-3s  %-9s\n", "i", "alias", "prob" ); */
/*     for ( int i=0; i<n; ++i ) */
/*       printf(  "  %3i  %3i  %9.7f\n", i, ws->alias[i], ws->prob[i] ); */
/* #endif */
    
/*     // Draw M random samples; accumulate statistic in histogram 'cumul': */
/*     printf(  "Drawing %i samples ...\n", M ); */
/*     double cumul[n]; */
/*     for ( i=0; i<n; ++i ) */
/*       cumul[i] = 0; */
/*     for ( m=0; m<M; ++m ) { */
/*       i = ransampl_draw( ws, gsl_rng_uniform(rng), gsl_rng_uniform(rng) ); */
/*       cumul[i] += 1; */
/*     } */
    
/*     // Print given probability and obtained frequency: */
/*     printf(  "Result (input->output):\n"); */
/*     double sum = 0; */
/*     for ( int i=0; i<n; ++i ) */
/*         sum += p[i]; */
/*     printf(  "  %-18s  %-9s  %-9s  %-9s\n", "state", "N (Mio.)", "rel", "sim" ); */
/*     for ( int i=0; i<n; ++i ) */
/*         printf(  "  %-18s  %9.7f  %9.7f  %9.7f\n", */
/*                 names[i], p[i], p[i]/sum, ((double)cumul[i])/M ); */

/*     // Free workspace and terminate: */
/*     ransampl_free( ws ); */
/*     return 0; */
  
  
  /* setlocale(LC_NUMERIC, ""); */
  /* printf( "Number of states: %'lu \n", num_states); */
  
  QUIT:
  
  /* Error reporting */
  
  if (error) {
    printf( "ERROR: %s\n", GRBgeterrormsg(master_env));
    exit(1);
  }
  
  /* Free model */
  
  GRBfreemodel(model);
  
  /* Free environment */
  
  GRBfreeenv(master_env);
  
  char *sched_str = is_ocbp_schedulable ? "YES" : "NO";
  printf( "OCBP Schdulability: %s\n", sched_str);
  fflush(stdout);
  
  /* char str[20]; */
  /* printf( "Done. Press any key to continue ...\n"); */
  /* scanf("%s", str); */
  return 0;
  
  /* ul num_invlid_states = 0; */
  /* for(int x1=0; x1 < dims[0]; ++x1) { */
  /*   for(int x2=0; x2 < dims[1]; ++x2) { */
  /*     for(int x3=0; x3 < dims[2]; ++x3) { */
  /* 	for(int x4=0; x4 < dims[3]; ++x4) { */
  /* 	  for(int y1=0; y1 < dims[4]; ++y1) { */
  /* 	    for(int y2=0; y2 < dims[5]; ++y2) { */
  /* 	      for(int y3=0; y3 < dims[6]; ++y3) { */
  /* 		for(int y4=0; y4 < dims[7]; ++y4) { */
  /* 		  for(int r=0;   r < dims[8];   ++r) { */
		    
  /* 		    ul index = dims[1]*dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x1 +  */
  /* 		      dims[2]*dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x2+ */
  /* 		      dims[3]*dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x3+ */
  /* 		      dims[4]*dims[5]*dims[6]*dims[7]*dims[8]*x4+ */
  /* 		      dims[5]*dims[6]*dims[7]*dims[8]*y1+ */
  /* 		      dims[6]*dims[7]*dims[8]*y2+ */
  /* 		      dims[7]*dims[8]*y3+ */
  /* 		      dims[8]*y4+ */
  /* 		      r; */
		    
  /* 		    if(arr_states[index] != 0 ){ */
  /* 		      ++num_invlid_states; */
  /* 		    } */
  /* 		    /\* if(ARR4(x1, x2,x3,x4,y1, y2,y3,y4,r) != 0) { *\/ */
		    
  /* 		    /\*   ++num_invlid_states; *\/ */
  /* 		    /\* } *\/ */
  /* 		  } */
  /* 		} */
  /* 	      } */
  /* 	    } */
  /* 	  } */
  /* 	} */
  /*     } */
  /*   } */
  /* } */
  
  /* /\* setlocale(LC_NUMERIC, ""); *\/ */
  /* printf( "Number of invalid states: %'lu\n", num_invlid_states); */
  



  
/*   int a; */
/*   int b; */
/*   int c; */
/*   int d; */
/*   int e; */

/*   int num_states = 0; */
/*   for(int sum = limit; sum >= 1; --sum) { */
/*     for(a = 0; a <= fmin(sum, arr_jobs[0].wcet); ++a) { */
/*       for(b = 0; b <= fmin(sum - a, arr_jobs[1].wcet); ++b) { */
	
/* /\* #if NUM_JOBS > 2 *\/ */
	
/* 	for(c = 0; c <= fmin(sum - a - b, arr_jobs[2].wcet); ++c) { */
/* 	  for(d = 0; d <= fmin(sum - a - b - c, arr_jobs[3].wcet); ++d) { */
/* 	    e = fmin(sum - a - b - c - d, arr_jobs[4].wcet); */
/* 	    if (a+b+c+d+e != sum) { */
/* 	      /\* printf( "Reject.\n"); *\/ */
/* 	      continue; */
/* 	    } */
/* 	    conductor->next = (state *) malloc( sizeof(state) ); */
	    
/* 	    conductor->next->time = sum; */
	  
/* 	    conductor->next->total_allocations[0] = a; */
/* 	    conductor->next->total_allocations[1] = b; */
/* 	    conductor->next->total_allocations[2] = c; */
/* 	    conductor->next->total_allocations[3] = d; */
/* 	    conductor->next->total_allocations[4] = e; */
/* 	    ++num_states; */
/* 	    /\* printf( "[%d %d %d %d]\n", w,x,y,z); *\/ */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*     printf( "%d\n", sum); */
    
/*   } */
/*   printf( "Number of states: %d\n", num_states); */

  /* fclose(fp); */
  
}
