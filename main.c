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

FILE *fdebug;

  
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

double min(double a, double b) {
  if(definitely_greater_than(a,b)) { return b; }
  return a;
}
double max(double a, double b) {
  if(definitely_greater_than(a,b)) { return a; }
  return b;
}

int imin(int a, int b) {
  if(a < b) { return a; }
  return b;
}


int imax(int a, int b) {
  if(a > b) { return a; }
  return b;
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

#define NUM_JOBS 3
#define NUM_JOBS_SUBMDP 2

typedef unsigned long ul;

/* A jobs criticality can be any of { 0, ..., NUM_CRITICALITY_LEVELS - 1 } */
#define NUM_CRITICALITY_LEVELS 3
char arr_num_jobs_of_criticality[NUM_CRITICALITY_LEVELS];


enum FINISH_SIGNALS {
  FINISHED = 0,
  NOT_FINISHED,
  MAX_FINISH_SIGNAL
};

/* Values that state->error_flag assumes */
enum ERROR_FLAGS {
  NOT_ERROR=0,
  ERROR,
  ERROR_ABSORB,
  MAX_ERROR_FLAG
};

ul *initial_allocation_vector;
char *initial_finish_signal_vector; 
int initial_error_flag;

typedef struct state_action_pair state_action_pair;
typedef struct state_action_pair_pointer state_action_pair_pointer;
/* typedef struct out_state_pointer out_state_pointer; */

enum STATE_TYPE {
  ABSORBING = 0,
  NON_ABSORBING
};

typedef struct state_tag {

  /* those are pointer to elements in state_actions_pairs_list */
  state_action_pair *admissible_state_action_pairs_list;
  int num_admissible_saps;
  
  state_action_pair_pointer *in_state_action_pairs_list;

  /* But field of size = 1 bit */
  char is_initial : 1;
  
  ul time;
  struct state_tag *next;

  int type;
  ul index; /* used to index into arr_primal_solution */
  
} state;

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
  state_action_pair *sap;
  /* specific to the parent state that has this sap as an in_state_action_pair. 
     Irrelevant if this object is an admissible state_action pair */
  double prob_parent_state_given_sap;
  
};

/* Each state_action_pair object contains a list of out_state_pointers 
   of all valid out states */
/* struct out_state_pointer { */
  
/*   out_state_pointer *next; */
/*   ul out_state_index; */
  
/*   double prob_out_state_given_sap; */
  
/* }; */

state *create_state(int type, ul time) {

  if(type == ABSORBING) {
    assert(time != NULL);
  }
  
  state *s = malloc( sizeof(state) );
  s->type                               = type;
  s->next				= NULL;
  s->time				= time;
  s->admissible_state_action_pairs_list	= NULL;
  s->num_admissible_saps		= 0;
  
  s->in_state_action_pairs_list		= NULL;
  
  s->is_initial				= false;
  
  return s;
  
}

struct state_action_pair {
  
  ul index;
  state *s; /* A pointer to an element in state_tree */
  int action;
  state_action_pair* next;
  double objective_cost;
  /* out_state_pointer* out_state_ptrs_list; */
};

/* void add_out_state_to_state_action_pair(ul out_state_index, */
/* 					state_action_pair *sap, */
/* 					double prob_out_state_given_sap); */

/* typedef struct state_pointer_per_allocation_tag { */
  
/*   state* s; /\* A pointer to an element in state_tree *\/ */
/*   struct state_pointer_per_allocation_tag *next; */
  
/* } state_pointer_per_allocation; */

/* state_pointer_per_allocation **arr_allocation_state_pointers; */

typedef struct job_tag {
  int deadline;
  int criticality;
  int wcet;
  
  /* The execution pmf is a probability mass function on  
     { 1, ..., arr_wcets_at_criticalities[criticality] } */
  double *execution_time_pmf;
  
  /* Might be assigned by a fixed-priority algorithm. */
  int priority;
  
} job_t;

job_t arr_jobs[NUM_JOBS];
int arr_deadlines[NUM_JOBS];
int *arr_job_indexes_per_criticality[NUM_CRITICALITY_LEVELS];

#define NUM_JOB_PAIRS (NUM_JOBS * (NUM_JOBS - 1) / 2)

/* Number of dimensions of the multidimensional array hodling the x,y and 
   vectors: sum of number of compoenents of all x,y*/
#define NUM_DIMS (2 * NUM_JOBS_SUBMDP)


typedef struct submdp_tag {

  /* indexes of the two job structs in arr_jobs */
  int arr_job_indexes[NUM_JOBS_SUBMDP]; 

  int horizon;
  
  state *state_tree;
  state *state_tree_tail;

  ul num_absorbing_states;
  ul num_nonabsorbing_states;
  
  state **arr_nonabsorbing_states; /* Multidimensional state pointers */
  state **arr_absorbing_states; 
  
  ul num_state_action_pairs;
  
  /* Those will be used in the optimization */
  /* state_action_pair_pointer *state_action_pair_ptrs_list; */
  /* state_action_pair_pointer *state_action_pair_ptrs_list_tail; */
  
  ul dims[NUM_DIMS];
  
  /* /\* Add also all the solution to the subMDP comprised of these two jobs *\/ */
  /* double    *arr_dual_solution; */
  /* double    *arr_primal_solution; */

  GRBmodel *grb_model;
  
} submdp_t;

submdp_t *arr_submdps[NUM_JOB_PAIRS];

void expand_state_tree_at(submdp_t *, state *, ul *, char *);



#define SUBMDP_STATE_INDEX(x1,x2,y1,y2, dims) (dims[1]*dims[2]*dims[3]*x1 + \
					       dims[2]*dims[3]*x2+	\
					       dims[3]*y1+		\
					       y2)

/* #define SUBMDP_STATE_INDEX(x, y, dims) (dims[1]*dims[2]*dims[3]*x[0] +	\ */
/* 					dims[2]*dims[3]*x[1]+		\ */
/* 					dims[3]*y[0]+			\ */
/* 					y[1])  */

#define SUBMDP_NONABSORBING_STATE_PTR(submdp, x1,x2,y1,y2)		\
  (submdp->arr_nonabsorbing_states[SUBMDP_STATE_INDEX(x1,x2,y1,y2,	\
						      submdp->dims)])


/* #define SUBMDP_NONABSORBING_STATE_PTR(submdp, x, y)			\ */
/*   (submdp->arr_nonabsorbing_states[SUBMDP_STATE_INDEX(x, y, submdp->dims)]) */


ul get_total_submdp_state_count(submdp_t *submdp) {
  
  return submdp->num_absorbing_states + submdp->num_nonabsorbing_states;
  
}

job_t *lookup_job_by_submdp_job_index(submdp_t * submdp, int submdp_job_index) {

  assert(0 <= submdp_job_index && submdp_job_index < NUM_JOBS_SUBMDP);

  return &arr_jobs[submdp->arr_job_indexes[submdp_job_index]];
  
  
}
state *lookup_absorbing_state_ptr(submdp_t *submdp, ul time) {
  
  
  int index;
  /* state *s = NULL; */
  
  for(index = 0; index < submdp->num_absorbing_states; ++index) {
    
    state *potential_state = submdp->arr_absorbing_states[index];
    
    if(potential_state->time == time) {
      return potential_state;
    }
    
  }
  return NULL;
}
/* If type == ABSORBING, then allocation_vector and finish_signal_vector 
   can be NULLs. If type == NON_ABSORBING, then allocation_vector and 
   finish_signal_vector cannot be NULL, but time can be NULL. */ 
state *lookup_state_ptr_by_params(submdp_t *submdp,
				  int type,
				  ul time,
				  ul *allocation_vector,
				  char *finish_signal_vector) {
  state *s = NULL;
  
  if(type == ABSORBING) {
    
    s = lookup_absorbing_state_ptr(submdp, time);
    
  } else if (type == NON_ABSORBING) {
    
    s = SUBMDP_NONABSORBING_STATE_PTR(submdp,
    				      allocation_vector[0],
    				      allocation_vector[1],
    				      finish_signal_vector[0],
    				      finish_signal_vector[1]);
    
    /* s = SUBMDP_NONABSORBING_STATE_PTR(submdp, */
    /* 				      allocation_vector, */
    /* 				      finish_signal_vector); */

    
    
  }
  
  return s;
  
  /* #elif NUM_JOBS == 3 */
  
  /*   s = ARR3(allocation_vector[0], */
/* 	   allocation_vector[1], */
/* 	   allocation_vector[2], */
/* 	   finish_signal_vector[0], */
/* 	   finish_signal_vector[1], */
/* 	   finish_signal_vector[2], */
/* 	   error_flag); */
  
/* #elif NUM_JOBS == 4 */
  
/*   s = ARR4(allocation_vector[0], */
/* 	   allocation_vector[1], */
/* 	   allocation_vector[2], */
/* 	   allocation_vector[3], */
/* 	   finish_signal_vector[0], */
/* 	   finish_signal_vector[1], */
/* 	   finish_signal_vector[2], */
/* 	   finish_signal_vector[3], */
/* 	   error_flag); */
  
/* #elif NUM_JOBS == 5 */
  
/*   s = ARR5(allocation_vector[0], */
/* 	   allocation_vector[1], */
/* 	   allocation_vector[2], */
/* 	   allocation_vector[3], */
/* 	   allocation_vector[4], */
/* 	   finish_signal_vector[0], */
/* 	   finish_signal_vector[1], */
/* 	   finish_signal_vector[2], */
/* 	   finish_signal_vector[3], */
/* 	   finish_signal_vector[4], */
/* 	   error_flag); */
/* #endif   */
  /* return s; */
}


/* Once a state_action pair is created, it cannot be destroyed and MUST be 
   linked because this function updates the global variable 
   state_action_pairs_count  */
state_action_pair* create_state_action_pair(submdp_t *submdp) {
  
  state_action_pair *sap = malloc(sizeof(state_action_pair));
  
  sap->action	      = 0;
  sap->next	      = NULL;
  sap->objective_cost = 0.0;
  /* sap->out_states_list = NULL; */
  
  sap->index = submdp->num_state_action_pairs;
  submdp->num_state_action_pairs++;
  
  /* state_action_pair_pointer *sap_pair_ptr = */
  /*   malloc(sizeof(state_action_pair_pointer)); */
  
  /* sap_pair_ptr->sap = sap; */

  /* /\* Link sap_pair_ptr to the list of sap_pointers *\/ */
  /* if(submdp->state_action_pair_ptrs_list == NULL) { */
  /*   submdp->state_action_pair_ptrs_list = sap_pair_ptr; */
  /*   submdp->state_action_pair_ptrs_list_tail = */
  /*     submdp->state_action_pair_ptrs_list; */
  /* } else { */
  /*   submdp->state_action_pair_ptrs_list_tail->next = sap_pair_ptr; */
  /*   submdp->state_action_pair_ptrs_list_tail = */
  /*     submdp->state_action_pair_ptrs_list_tail->next; */
  /* }   */
  
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

/* lookup_state_ptr_by_params(...) should be called prior to calling 
add_state to make sure that a state pointer has not been created already. */
state* add_state(submdp_t *submdp,
		 int type,
		 ul time,
		 ul *allocation_vector,
		 char *finish_signal_vector) {
  
  state *s = create_state(type, time);
  s->index = get_total_submdp_state_count(submdp);
  
  if(type == ABSORBING) {
    
    assert( (NUM_JOBS_SUBMDP + 1 <= time) && (time <= submdp->horizon + 1) );
    
    /* printf("num absorbing states: %d, horizon: %d\n", submdp->num_absorbing_states, submdp->horizon); */
    
    assert(submdp->num_absorbing_states <
    	   submdp->horizon - NUM_JOBS_SUBMDP + 1);
    
    submdp->arr_absorbing_states[submdp->num_absorbing_states] = s;
    
    submdp->num_absorbing_states++;
    
  } else if (type == NON_ABSORBING) {
    
    SUBMDP_NONABSORBING_STATE_PTR(submdp,
    				  allocation_vector[0],
    				  allocation_vector[1],
    				  finish_signal_vector[0],
    				  finish_signal_vector[1]) = s;

    /* SUBMDP_NONABSORBING_STATE_PTR(submdp, */
    /* 				  allocation_vector, */
    /* 				  finish_signal_vector) = s; */
    
    submdp->num_nonabsorbing_states++;
  }
  
  
  /* add to state_tree 
  
  /* Link next_state to the list of states */
  if(submdp->state_tree == NULL) {
    submdp->state_tree = s;
    submdp->state_tree_tail = submdp->state_tree;
  } else {
    submdp->state_tree_tail->next = s;
    submdp->state_tree_tail = submdp->state_tree_tail->next;
  }  
  /* if((num_states % 10000000) == 0) { */
  /*   setlocale(LC_NUMERIC, ""); */
  
  /*   printf( "%'lu\n", num_states); */
  /* } */
  
  return s;
  
}
  
void add_state_action_pair(submdp_t *submdp,
			   state *s,
			   state_action_pair* sap) {
  
  /* assert(state_action_pairs_list_tail != NULL); */
  
  /* add sap to state s's admissible_state_action_pairs_list */
  /* state_action_pair_pointer *sap_pointer = */
  /*   malloc(sizeof(state_action_pair_pointer)); */
  
  /* sap_pointer->sap = sap; */
  sap->s = s;
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
				       state_action_pair *sap,
				       double prob_s_given_sap) {
  
  state_action_pair_pointer *in_state_action_pair_ptr =
    malloc(sizeof(state_action_pair_pointer));

  in_state_action_pair_ptr->sap = sap;
  
  in_state_action_pair_ptr->prob_parent_state_given_sap = prob_s_given_sap;
  
  in_state_action_pair_ptr->next = s->in_state_action_pairs_list;
  
  s->in_state_action_pairs_list = in_state_action_pair_ptr;
  /* s->num_total_saps_with_nonzero_coefficients++; */
  
}

/* void add_out_state_to_state_action_pair(ul out_state_index, */
/* 					state_action_pair *sap, */
/* 					double prob_out_state_given_sap) { */

/*   out_state_pointer *out_state_ptr = malloc(sizeof(out_state_pointer)); */
/*   out_state_ptr->out_state_index = out_state_index; */
  
/*   out_state_ptr->prob_out_state_given_sap = prob_out_state_given_sap; */
  
/*   /\* out_state_container->cost_state_given_sap = transition_cost; *\/ */
  
/*   out_state_ptr->next = sap->out_state_ptrs_list; */
/*   sap->out_state_ptrs_list = out_state_ptr; */
  
/* } */

/* bool any_job_just_missed_deadline(int criticality, */
/* 				  int time, */
/* 				  char *arr_finish_signal_vector) { */
  
/*   int num_jobs = arr_num_jobs_of_criticality[criticality]; */
  
/*   for(int index = 0; index < num_jobs; ++index) { */
/*     int job_index = arr_job_indexes_per_criticality[criticality][index]; */
/*     job_t *j = &arr_jobs[job_index]; */
/*     if((time == j->deadline) && */
/*        (arr_finish_signal_vector[job_index] == NOT_FINISHED)){  */
      
/*       return true; */
/*     } */
/*   } */
/*   return false; */
/* } */

/* ul get_sum_allocations(int criticality, ul *allocation_vector) { */
  
/*   ul sum_allocations = 0; */
/*   int num_jobs = arr_num_jobs_of_criticality[criticality]; */
  
/*   for(int index = 0; index < num_jobs; ++index) { */
/*     int job_index = arr_job_indexes_per_criticality[criticality][index]; */
/*     sum_allocations += allocation_vector[job_index]; */
/*   } */
  
/*   return sum_allocations; */
/* } */

/* allocation_vector and finish_signal_vector are never NULL */
char is_admissible_action(ul *allocation_vector,
			  char *finish_signal_vector,
			  int action) {
  
  return (finish_signal_vector[action] == NOT_FINISHED);
  
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
/* Assumes: tha state transition is valid */
double compute_transition_probability(submdp_t *submdp,
				      ul *prev_allocation_vector,
				      char* prev_finish_signal_vector,
				      ul *next_allocation_vector,
				      char* next_finish_signal_vector,
				      int action) {
  
  
  
  double prob = 0.0;

  if(action == NO_ACTION) {
    assert(next_allocation_vector   == NULL);
    assert(next_finish_signal_vector == NULL);
    
    assert((prev_finish_signal_vector == NULL) ||
	   have_all_jobs_finished(prev_finish_signal_vector));
    
    prob = 1.0;
    return prob;
    
  } 

  assert(prev_finish_signal_vector != NULL);
  assert(prev_allocation_vector    != NULL);

  assert(next_finish_signal_vector != NULL);
  assert(next_allocation_vector    != NULL);

  
  
  char prev_finish_signal = prev_finish_signal_vector[action];
  assert(prev_finish_signal == NOT_FINISHED);
  
  char next_finish_signal = next_finish_signal_vector[action];
  ul prev_allocation = prev_allocation_vector[action];
  ul next_allocation  = next_allocation_vector[action];
  
  assert((next_allocation == prev_allocation + 1));
  
  job_t *j = lookup_job_by_submdp_job_index(submdp, action);
  int max_wcet = j->wcet;
  
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
    
    return ((double)1.0 - prob_finished);
    
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
  
  /* int sys_criticality = */
  /*   get_system_criticality_level_realization */
  /*   (next_allocation_vector, next_finish_signal_vector); */
  
  /* bool lo_job_just_missed_deadline = */
  /*   any_job_just_missed_deadline */
  /*   (LO, next_time, next_finish_signal_vector); */
  
  /* bool hi_job_just_missed_deadline = */
  /*   any_job_just_missed_deadline */
  /*   (HI, next_time, next_finish_signal_vector); */

  /* /\* if(hi_job_just_missed_deadline) { *\/ */
  /* /\*   printf( "Hi job just missed deadline!\n"); *\/ */
  /* /\* } *\/ */
  
  /* /\* There can only be one next_error_flag per  */
  /*    state-action-finish_signal*\/ */
  /* switch (previous_error_flag) { */
    
  /* case ERROR: */
  /*   /\* Next error flag can only be ERROR_ABSORB *\/ */
  /*   return ERROR_ABSORB; */
  /*   break; */
  /* case ERROR_ABSORB: */
  /*   /\* Next error flag can only be ERROR_ABSORB *\/ */
  /*   return ERROR_ABSORB; */
  /*   break; */
  /* case NOT_ERROR: */
    
  /*   /\* Next error flag can be any {NOT_ERROR, POTENTIAL_ERROR, ERROR} *\/ */
    
  /*   /\* Check for POTENTIAL_ERROR *\/ */
  /*   if(sys_criticality == UNKNOWN && */
  /*      lo_job_just_missed_deadline && */
  /*      !hi_job_just_missed_deadline) { */
      
  /*     return POTENTIAL_ERROR; */
      
  /*   } */
    
  /*   /\* Check for ERROR *\/ */
  /*   /\* NOT_ERROR -> ERROR happens if either a HI criticality job just */
  /*      misses its deadline (i.e., at next_time), or a LO criticality job  */
  /*      misses its deadline and the system criticality level realization  */
  /*      at next_time is LO *\/ */
    
  /*   if(hi_job_just_missed_deadline || */
  /*      (sys_criticality == LO && lo_job_just_missed_deadline)) { */
      
  /*     /\* ERROR *\/ */
  /*     return ERROR; */
      
  /*   } */
    
  /*   /\* Neither POTENTIAL_ERROR nor ERROR. Then only NOT_ERROR *\/ */
  /*   return NOT_ERROR; */
    
  /*   break; */
    
  /* case POTENTIAL_ERROR: */
  /*   { */
  /*     /\* Next error flag can be any {NOT_ERROR, POTENTIAL_ERROR, ERROR} *\/ */
      
  /*     /\* Check for ERROR *\/ */
  /*     /\* POTENTIAL_ERROR -> ERROR happens if either the system  */
  /* 	 criticality level is revealed as LO at next_time, or a HI  */
  /* 	 criticality job just misses its deadline at next_time *\/ */
      
  /*     if(hi_job_just_missed_deadline || (sys_criticality == LO)) { */
  /* 	/\* ERROR *\/ */
  /* 	return ERROR; */
  /*     } */
      
  /*     /\* Check for POTENTIAL_ERROR *\/ */
  /*     if((sys_criticality == UNKNOWN) && !hi_job_just_missed_deadline) { */
  /* 	return POTENTIAL_ERROR; */
  /*     } */
      
  /*     /\* Check for NOT_ERROR *\/ */
  /*     /\* POTENTIAL_ERROR -> NOT_ERROR happens if the system criticality  */
  /* 	 level is realized at next_time (= time + 1) as HI and no  */
  /* 	 HI criticality job misses its deadline at (time + 1) *\/ */
  /*     if(sys_criticality == HI && !hi_job_just_missed_deadline) { */
	
  /* 	return NOT_ERROR; */
	
  /*     } */
      
  /*     break; */
  /*   } */
  /* } */
  return ERROR;
}

void set_state_action_pair_dual_cost(submdp_t *submdp,
				     state_action_pair *sap,
				     char *finish_signal_vector) {
  
  if(sap->s->type == ABSORBING) {
    sap->objective_cost = 0.0;
  } else {
    int job_index;
    /* printf("type: %d, index:%d\n", root_state_ptr->type, */
    /* 	   root_state_action_pair->index);  */
    assert(essentially_equal(sap->objective_cost, 0.0));
    
    for (job_index = 0; job_index < NUM_JOBS_SUBMDP; ++job_index) {
      job_t *job = lookup_job_by_submdp_job_index(submdp, job_index);
      
      if((finish_signal_vector[job_index] == NOT_FINISHED) &&
	 (sap->s->time >= job->deadline)) {
	++sap->objective_cost;
      }
      
    }
    /* printf("%f\n\n", root_state_action_pair->dual_objective_cost);  */
  }
  
  
}

/* If the state to be expanded (root) is such that all jobs finished,
   then next_allocation_vector = next_finish_signal_vector  = NULL.
   An absorbing state is never expanded, so root_allocation_vector and
   root_finish_signal_vector are never NULL */
void _expand(submdp_t *submdp,
	     state *root_state_ptr,
	     state_action_pair *root_state_action_pair,
	     ul *root_allocation_vector,
	     char *root_finish_signal_vector,
	     ul *next_allocation_vector,
	     char *next_finish_signal_vector) {
  
  if(root_state_ptr->time == submdp->horizon + 1) {
    assert(root_state_ptr->type == ABSORBING);
    return;
  }
  
  ul next_time = root_state_ptr->time + 1;
  
  int action = root_state_action_pair->action;
  
  state *next_state;
  
  char next_state_exists = true; 
    
  if(action == NO_ACTION) {
    
    /* Here next state is absorbing. Add it if it does not exist */
    
    next_state = lookup_state_ptr_by_params(submdp,
					    ABSORBING,
					    next_time,
					    NULL,
					    NULL);
    
    if(next_state == NULL) { /* next_state not added yet */
      
      next_state = add_state(submdp, ABSORBING, next_time, NULL, NULL);
      
      next_state_exists = false;
      
    }
    
  } else {
    next_state = lookup_state_ptr_by_params(submdp,
					    NON_ABSORBING,
					    next_time,
					    next_allocation_vector,
					    next_finish_signal_vector);
    
    if(next_state == 0 || next_state == NULL) {

      
      /* next_state does not exist yet */
      next_state_exists = false;
      
      next_state = add_state(submdp,
			     NON_ABSORBING,
			     next_time,
			     next_allocation_vector,
			     next_finish_signal_vector);
      
      /* next_state->constraint_cost = (next_error_flag == ERROR) ? 1 : 0; */
      
    }
    
  }
  
  if((!next_state_exists) && (next_time <= submdp->horizon + 1)) {
    
    expand_state_tree_at(submdp,
			 next_state,
			 next_allocation_vector,
			 next_finish_signal_vector);
  }
  
  /* root_state_action_pair appears for the first time, so add it 
     as an in_state_action_pair to next_state after computing Q */
  
  double prob_next_state_given_root_state_action_pair =
    compute_transition_probability(submdp,
				   root_allocation_vector,
				   root_finish_signal_vector,
				   next_allocation_vector,
				   next_finish_signal_vector,
				   action);
  
  /*, Add it only if Q(next_state | root_state_action_pair) > 0
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
       root_state_action_pair,
       prob_next_state_given_root_state_action_pair);

    /* add_out_state_to_state_action_pair */
    /*   (next_state->index, */
    /*    root_state_action_pair, */
    /*    prob_next_state_given_root_state_action_pair); */

  }
    
  /* even if next_state exists, root_state_action_pair appears for 
     the first time here, so add next_state as an out_state to 
     root_state_action_pair */
  
  /* Compute the immediate dual objective cost */
  
  
    
}


/* root_state is added already whenever this function is called  */
/* Invairant: root_state_ptr =! NULL and 
   root_state_ptr->time <= submdp->horizon */
// <NOTE>
// By BNA
// state_action_pair objects are stored as as elements of
// admissible_state_action_pairs_list, ad they can be retrieved by taking
// the union of the elements of admissible_state_action_pairs_list across all
// states, and the sap->next is the next sap in this list.
// To keep a reference to a sap, create a state_action_pair_pointer object
// and include a reference to sap in it. This is what we do to keep a global
// array of all saps in submdp->state_action_pair_ptrs_list.
// </NOTE>
void expand_state_tree_at(submdp_t *submdp,
			  state *root_state_ptr,
			  ul *root_allocation_vector,
			  char *root_finish_signal_vector) {
  
  assert(root_state_ptr != 0);

  state_action_pair *root_state_action_pair;
  
  ul *next_allocation_vector      = NULL;
  char *next_finish_signal_vector = NULL;
  
  if( (root_state_ptr->type == ABSORBING) ||
      have_all_jobs_finished(root_finish_signal_vector) ) {
    
    /* terminal_state_action_pair->s = next_state; */

    root_state_action_pair = create_state_action_pair(submdp);
    root_state_action_pair->action = NO_ACTION;
    add_state_action_pair(submdp, root_state_ptr, root_state_action_pair);
    set_state_action_pair_dual_cost(submdp,
				    root_state_action_pair,
				    root_finish_signal_vector);
    _expand(submdp,
	    root_state_ptr,
	    root_state_action_pair,
	    root_allocation_vector,
	    root_finish_signal_vector,
	    NULL,
	    NULL);
    
    
  } else {
    
    for(int action_index = 0; action_index < NUM_JOBS_SUBMDP; ++action_index){
      
      if(!is_admissible_action(root_allocation_vector,
			       root_finish_signal_vector,
			       action_index)) {
	continue;
      }
      
      root_state_action_pair = create_state_action_pair(submdp);
      root_state_action_pair->action = action_index;
      add_state_action_pair(submdp, root_state_ptr, root_state_action_pair);
      set_state_action_pair_dual_cost(submdp,
				      root_state_action_pair,
				      root_finish_signal_vector);
      
      /* Then action "action_index" is a possible next action *
	 /* root_state_action_pair->s = root_state_ptr; */
      
      if(next_allocation_vector == NULL) {
	next_allocation_vector = clone_ul_array(root_allocation_vector,
						NUM_JOBS_SUBMDP);
      }
      
      next_allocation_vector[action_index]++;
      
      job_t *current_job = lookup_job_by_submdp_job_index(submdp,action_index); 
      
      int current_job_wcet = current_job->wcet;
      
      int num_possible_finish_signals_for_current_job;
      char arr_possible_finish_signals_for_current_job[2];
      
      if(current_job_wcet == next_allocation_vector[action_index]){
	num_possible_finish_signals_for_current_job = 1;
	arr_possible_finish_signals_for_current_job[0] = FINISHED;
      } else {
	assert(next_allocation_vector[action_index] < current_job_wcet);
	num_possible_finish_signals_for_current_job = 2;
	arr_possible_finish_signals_for_current_job[0] = NOT_FINISHED;
	arr_possible_finish_signals_for_current_job[1] = FINISHED;
      }
      
      
      for(int finish_signal_index = 0;
	  finish_signal_index < num_possible_finish_signals_for_current_job;
	  ++finish_signal_index) {
	
	int next_finish_signal =
	  arr_possible_finish_signals_for_current_job[finish_signal_index];
	
	if(next_finish_signal_vector == NULL) {
	  next_finish_signal_vector =
	    clone_char_array(root_finish_signal_vector, NUM_JOBS_SUBMDP);
	}
	
	next_finish_signal_vector[action_index] = next_finish_signal;
	
	_expand(submdp,
		root_state_ptr,
		root_state_action_pair,
		root_allocation_vector,
		root_finish_signal_vector,
		next_allocation_vector,
		next_finish_signal_vector);
	
	/* we are still inside the finish_signal loop */
	
	next_finish_signal_vector[action_index] =
	  root_finish_signal_vector[action_index];
      }
      next_allocation_vector[action_index] =
	root_allocation_vector[action_index];
    }
    
    
    if(next_finish_signal_vector != NULL){ free(next_finish_signal_vector); }
    if(next_allocation_vector    != NULL){ free(next_allocation_vector); }
  }
}

void build_state_tree(submdp_t *submdp) {
  
  state *initial_state_ptr = add_state(submdp,
				       NON_ABSORBING,
				       0,
				       initial_allocation_vector,
				       initial_finish_signal_vector);
  
  assert(initial_state_ptr != NULL);
  
  initial_state_ptr->is_initial = true;
  
  expand_state_tree_at(submdp,
		       initial_state_ptr,
		       initial_allocation_vector,
		       initial_finish_signal_vector);
  
}

double get_prob_job_demands_at_most(job_t *j, int demand_upper_bound) {
  
  double prob = 0.0;
  int max_demand = j->wcet;
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
/* bool compute_OCBP_priorities() { */
  
/*   int priority = 0; */
/*   int_list* job_indexes_list = NULL; */
/*   int_list* job_index = NULL; */
  
/*   for (int index = 0; index < NUM_JOBS; ++index) { */
/*     job_index = malloc(sizeof(int_list)); */
/*     job_index->index = index; */
/*     job_index->next = job_indexes_list; */
/*     job_index->prev = NULL; */
/*     if(job_indexes_list != NULL) { */
/*       job_indexes_list->prev = job_index; */
/*     }  */
/*     job_indexes_list = job_index; */
    
/*   } */
  
/*   char is_lowest_priority_job_found = false; */
/*   int num_remaining_jobs = NUM_JOBS; */
  
/*   while(num_remaining_jobs > 0) { */
    
/*     is_lowest_priority_job_found = false; */

/*     int_list *job_indexes_list_tail = job_indexes_list; */
    
/*     while(job_indexes_list_tail != NULL) { */
    
/*       int sum_exec_times = 0; */
/*       job *lowest_priority_job = &arr_jobs[job_indexes_list_tail->index]; */
      
/*       int_list *job_indexes_list_tail2 = job_indexes_list; */
      
/*       while(job_indexes_list_tail2 != NULL) { */
	
/* 	job *j = &arr_jobs[job_indexes_list_tail2->index]; */
/* 	if(lowest_priority_job->criticality < j->criticality) { */
/* 	  sum_exec_times += */
/* 	    j->arr_wcets_at_criticalities[lowest_priority_job->criticality]; */
/* 	} else { */
/* 	  sum_exec_times += j->arr_wcets_at_criticalities[j->criticality]; */
/* 	} */
/* 	job_indexes_list_tail2 = job_indexes_list_tail2->next; */
/*       } */
      
/*       if(sum_exec_times <= lowest_priority_job->deadline) { */
/* 	lowest_priority_job->priority = priority; */
/* 	++priority; */
/* 	// delete the job_indexes_list_tail */
/* 	if(job_indexes_list_tail->prev != NULL) { /\* middle or last *\/ */
/* 	  (job_indexes_list_tail->prev)->next = job_indexes_list_tail->next; */
/* 	} else { */
/* 	  /\* job_indexes_list_tail is the head  *\/ */
/* 	  assert(job_indexes_list == job_indexes_list_tail); */
/* 	  job_indexes_list = job_indexes_list_tail->next; */
/* 	} */
/* 	if(job_indexes_list_tail->next != NULL) { /\* not last element *\/ */
/* 	  (job_indexes_list_tail->next)->prev = job_indexes_list_tail->prev; */
/* 	} */
/* 	job_indexes_list_tail->next = NULL; */
/* 	job_indexes_list_tail->prev = NULL; */
/* 	free(job_indexes_list_tail); */
/* 	job_indexes_list_tail = NULL; */
/* 	is_lowest_priority_job_found = true; */
/* 	--num_remaining_jobs; */
/* 	break; */
/*       } else { */
/* 	is_lowest_priority_job_found = false; */
/* 	job_indexes_list_tail = job_indexes_list_tail->next; */
/*       } */
/*     } */
   
/*     if(!is_lowest_priority_job_found) { */
/*       break; */
/*     } */
    
/*   } */
  
/*   int_list *job_indexes_list_tail = job_indexes_list; */
/*   while(job_indexes_list_tail != NULL) { */
/*     job_indexes_list_tail = job_indexes_list_tail->next; */
/*     free(job_indexes_list); */
/*     job_indexes_list = job_indexes_list_tail; */
/*   } */
  
/*   return is_lowest_priority_job_found; */
/* } */

void print_exec_time_sample(int *arr_exec_time_sample) {

  printf("Job demands sample:     ");
  for(int index = 0; index < NUM_JOBS; ++index) {
    printf("J%d: %d", index, arr_exec_time_sample[index]);
  }
  printf("\n");
}

void print_job_set() {
  
  printf( "\nNumber of jobs: %d\n", NUM_JOBS);
  fflush(stdout);
  
  int index = 0;
  for(index = 0; index < NUM_JOBS; ++index) {
    job_t *j = &arr_jobs[index];
    printf( "J%d   deadline: %d, criticality: %d, WCET: %d\n",
	    index,
	    j->deadline,
	    j->criticality,
	    j->wcet);
    fflush(stdout);
  }
  
  
}

void print_submdp(submdp_t *submdp, int submdp_index) {

  printf("+===========================SubMDP-%d============================\n",
	 submdp_index);
  
  fflush(stdout);
  
  int index = 0;
  for(index = 0; index < NUM_JOBS_SUBMDP; ++index) {
    job_t *j = lookup_job_by_submdp_job_index(submdp, index);
    printf( "| J%d   deadline: %d, criticality: %d, WCET: %d\n",
	    submdp->arr_job_indexes[index],
	    j->deadline,
	    j->criticality,
	    j->wcet);
    fflush(stdout);
    
  }
  
  printf( "| Horizon: %d\n", submdp->horizon);
  fflush(stdout);
  printf("+================================================================\n\n");
  fflush(stdout);
}

/* void solve_unconstrained_primal_submdp_exact(submdp_t *submdp) { */
    
/*   /\**--------------------Begin Building and Solving LP----------------------**\/ */
  
/*   /\* TODO: add the variables by iterating through  */
/*      submdp->state_action_pair_ptrs_list and for each sap_ptr,  */
/*      consider the state in the sap with coeff. = 1, and the  */
/*      out_state_ptrs_list with coeff. that is -Q(s'|s,a) =  */
/*      out_state_ptr->prob_state_given_child_sap. The others are all */
/*      zeros and do not need to be added. */
/*   *\/ */
  
/*   GRBenv   *master_env   = NULL; */
/*   GRBmodel *model = NULL; */
/*   int       error = 0; */
  
/*   submdp->arr_primal_solution = calloc(get_total_submdp_state_count(submdp), */
/* 				       sizeof(double));  */
  
/*   // <NOTE> */
/*   // By BNA */
/*   // ====== */
/*   // state_action_pair->index corresponds to its index in arr_primal_solution */
/*   // and arr_objective_coefficients. */
/*   // arr_primal_variable_indexes contains the */
/*   // state_action_pair indexes that have nonzero coefficients in the */
/*   // constraint being constructed. */
/*   // </NOTE> */
/*   int *arr_primal_variable_indexes = */
/*     calloc(get_total_submdp_state_count(submdp), sizeof(int)); */
  
/*   double *arr_constraint_coefficients = */
/*     calloc(get_total_submdp_state_count(submdp), sizeof(double)); */
  
/*   double *arr_objective_coefficients = */
/*     calloc(get_total_submdp_state_count(submdp), sizeof(double)); */
  
/*   int       optimization_status; */
/*   double    objective_value; */
  
/*   /\* Create environment *\/ */
  
/*   error = GRBloadenv(&master_env, "lp.log"); */
/*   if (error) goto QUIT; */
      
/*   /\* Create an empty model *\/ */
  
/*   error = GRBnewmodel(master_env,&model, "lp", 0, NULL, NULL, NULL, NULL, NULL); */
/*   if (error) goto QUIT; */
  
/*   printf( "Creating objective coefficients array ... \n"); */
/*   fflush(stdout); */
  
  
/*   /\* Add variables *\/ */
/*   submdp->state_tree_tail = submdp->state_tree; */
/*   while(submdp->state_tree_tail != NULL) { */

/*     state *s = submdp->state_tree_tail; */
    
/*     arr_objective_coefficients[s->index] = s->is_initial ? 1 : 0; */
    
/*     submdp->state_tree_tail = submdp->state_tree_tail->next; */
/*   }  */
  
/*   printf( "Adding variables ... \n"); */
/*   fflush(stdout); */
  
/*   error = */
/*     GRBXaddvars(model, get_total_submdp_state_count(submdp), */
/* 		0, NULL, NULL, NULL, */
/* 		arr_objective_coefficients, NULL, NULL, NULL, NULL); */
/*   if (error) goto QUIT; */
  
/*   /\* Change objective sense to maximization *\/ */
  
/*   error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE); */
/*   if (error) goto QUIT; */
  
/*   GRBenv *model_env = GRBgetenv(model); */
  
/*   /\* error = GRBsetintparam(model_env, GRB_INT_PAR_THREADS, 10); *\/ */
/*   /\* if (error) goto QUIT; *\/ */
  
/*   error = GRBsetintparam(model_env, GRB_INT_PAR_CROSSOVER, 0); */
/*   if (error) goto QUIT; */

/*   /\* GRB_INT_PAR_METHOD = 0: Primal Simplex *\/ */
/*   /\*                      1: Primal Simplex *\/ */
/*   /\*                      2: Barrier *\/ */
/*   /\*                      3: Nondeterministic Concurrent *\/ */
/*   /\*                      4: Deterministic Concurrent *\/ */
/*   /\* error = GRBsetintparam(model_env, GRB_INT_PAR_METHOD, 0); *\/ */
			 
/*   if (error) goto QUIT; */
  
/*   /\* Integrate new variables *\/ */
  
/*   error = GRBupdatemodel(model); */
/*   if (error) goto QUIT; */

/*   printf( "Creating constraints ... \n"); */
/*   fflush(stdout); */
/*   /\* Create the constraints *\/ */

/*   ul num_nonzero_constraint_coefficients = 0; */
/*   double constraint_rhs = 0; */

/*   submdp->state_action_pair_ptrs_list_tail = */
/*     submdp->state_action_pair_ptrs_list; */
  
/*   while(submdp->state_action_pair_ptrs_list_tail != NULL) { */

/*     num_nonzero_constraint_coefficients = 0; */
    
/*     state_action_pair *sap = */
/*       submdp->state_action_pair_ptrs_list_tail->sap; */
    
/*     arr_primal_variable_indexes[num_nonzero_constraint_coefficients] = */
/*       sap->s->index; */

/*     arr_constraint_coefficients[num_nonzero_constraint_coefficients] = 1; */
    
/*     num_nonzero_constraint_coefficients++; */

/*     // Now loop through the out_states of the sap. */
/*     out_state_pointer *out_state_ptrs_list = sap->out_state_ptrs_list; */
/*     while(out_state_ptrs_list != NULL) { */
      
/*       out_state_pointer *out_state_ptr = out_state_ptrs_list; */
      
/*       arr_primal_variable_indexes[num_nonzero_constraint_coefficients] = */
/* 	out_state_ptr->out_state_index; */

/*       /\* we create an out_state for a sap only if  */
/* 	 out_state_ptr->prob_out_state_given_sap > 0 *\/ */
/*       arr_constraint_coefficients[num_nonzero_constraint_coefficients] = */
/* 	-1 * out_state_ptr->prob_out_state_given_sap; */
      
/*       ++num_nonzero_constraint_coefficients; */
      
/*     } */
    
    
/*     constraint_rhs = sap->dual_objective_cost; */
    
/*     error = GRBaddconstr(model, */
/*     			 num_nonzero_constraint_coefficients, */
/*     			 arr_primal_variable_indexes, */
/*     			 arr_constraint_coefficients, */
/*     			 GRB_LESS_EQUAL, */
/*     			 constraint_rhs, */
/*     			 NULL); */
/*     if (error) goto QUIT; */
    
    
/*     submdp->state_action_pair_ptrs_list_tail = */
/*       submdp->state_action_pair_ptrs_list_tail->next; */
    
    
    
/*   }  */
  
/*   /\* Optimize model *\/ */

/*   printf( "Solving LP ... \n"); */
/*   fflush(stdout); */
  
/*   error = GRBoptimize(model); */
/*   if (error) goto QUIT; */
  
/*   /\* Write model to 'lp.lp' *\/ */
  
/*   error = GRBwrite(model, "lp.lp"); */
/*   if (error) goto QUIT; */
  
/*   /\* Capture solution information *\/ */
  
/*   error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimization_status); */
/*   if (error) goto QUIT; */

/*   if (optimization_status == GRB_OPTIMAL) { */
/*     error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objective_value); */
/*     if (error) goto QUIT; */
/*     printf( "Optimal objective: %.4e\n", objective_value); */
/*     fflush(stdout); */
    
/*     error = GRBgetdblattrarray(model, */
/* 			       GRB_DBL_ATTR_X, */
/* 			       0, */
/* 			       get_total_submdp_state_count(submdp), */
/* 			       submdp->arr_primal_solution); */
/*     if (error) goto QUIT; */
    
/*     /\* printf( "  x=%.0f, y=%.0f, z=%.0f\n", sol[0], sol[1], sol[2]); *\/ */
/*   } else if (optimization_status == GRB_INF_OR_UNBD) { */
/*     printf( "Model is infeasible or unbounded\n"); */
/*     fflush(stdout); */
/*     goto QUIT; */
    
/*   } else { */
/*     printf( "Optimization was stopped early\n"); */
/*     fflush(stdout); */
/*     goto QUIT; */
/*   } */
  
  
/*  QUIT: */
  
/*   /\* Error reporting *\/ */
  
/*   if (error) { */
/*     printf( "ERROR: %s\n", GRBgeterrormsg(master_env)); */
/*     exit(1); */
/*   } */
  
/*   /\* Free model *\/ */
  
/*   GRBfreemodel(model); */
  
/*   /\* Free environment *\/ */
  
/*   GRBfreeenv(master_env); */
  
/* } */


double lookup_primal_variable_solution_by_index(submdp_t *submdp,
						ul state_index){
  double solution = 0.0;
  
  GRBgetdblattrelement(submdp->grb_model,
		       GRB_DBL_ATTR_PI,
		       state_index,
		       &solution);
  return solution;
  
}

double lookup_dual_variable_solution_by_index(submdp_t *submdp, ul sap_index){

  double solution = 0.0;

  GRBgetdblattrelement(submdp->grb_model,
		       GRB_DBL_ATTR_X,
		       sap_index,
		       &solution);
    
  return solution;
  
}

void solve_unconstrained_submdp_exact(submdp_t *submdp) {
  
  /**--------------------Begin Building and Solving LP----------------------**/
  
  GRBenv   *master_env   = NULL;
  submdp->grb_model      = NULL;
  int       error        = 0;
  
  /* submdp->arr_dual_solution = */
  /*   calloc(submdp->num_state_action_pairs, sizeof(double)); */
  
  // <NOTE>
  // By BNA
  // ======
  // state_action_pair->index corresponds to its index in arr_dual_solution
  // and arr_objective_coefficients.
  // arr_dual_variable_indexes contains the
  // state_action_pair indexes that have nonzero coefficients in the
  // constraint being constructed.
  // </NOTE>
  int *arr_variable_indexes = calloc(submdp->num_state_action_pairs,
					  sizeof(int));
  
  double *arr_constraint_coefficients = calloc(submdp->num_state_action_pairs,
					       sizeof(double));
  double *arr_objective_coefficients = calloc(submdp->num_state_action_pairs,
					      sizeof(double));
  
  int       optimization_status;
  double    objective_value;
  
  /* Create environment */
  
  error = GRBloadenv(&master_env, "lp.log");
  if (error) goto QUIT;
  
  /* Create an empty model */
  
  error = GRBnewmodel(master_env, &submdp->grb_model, "lp", 0,
		      NULL, NULL, NULL, NULL, NULL);
  if (error) goto QUIT;
  
  printf( "Creating objective coefficients array ... \n");
  fflush(stdout);
  
  
  /* Add variables */
  submdp->state_tree_tail = submdp->state_tree;
  while(submdp->state_tree_tail != NULL) {
    
    state_action_pair *sap_list =
      submdp->state_tree_tail->admissible_state_action_pairs_list;
    
    while(sap_list != NULL) {

      /* index i in arr_objective_coefficients corresponds to 
	 state_action_pair whose index is i. We did the tying here. */
      
      arr_objective_coefficients[sap_list->index] = sap_list->objective_cost;
      
      sap_list = sap_list->next;
    } 
    
    submdp->state_tree_tail = submdp->state_tree_tail->next;
  } 
  
  printf( "Adding variables ... \n");
  fflush(stdout);
  
  error =
    GRBXaddvars(submdp->grb_model, submdp->num_state_action_pairs,
		0, NULL, NULL, NULL,
		arr_objective_coefficients, NULL, NULL, NULL, NULL);
  if (error) goto QUIT;
  
  /* Change objective sense to minimization */
  
  error = GRBsetintattr(submdp->grb_model,
			GRB_INT_ATTR_MODELSENSE,
			GRB_MINIMIZE);
  if (error) goto QUIT;
  
  GRBenv *model_env = GRBgetenv(submdp->grb_model);
  
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_THREADS, 10); */
  /* if (error) goto QUIT; */
  
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_CROSSOVER, 0); */
  /* if (error) goto QUIT; */
  
  /* GRB_INT_PAR_METHOD = 0: Primal Simplex */
  /*                      1:  Simplex */
  /*                      2: Barrier */
  /*                      3: Nondeterministic Concurrent */
  /*                      4: Deterministic Concurrent */
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_METHOD, 0); */
  
  if (error) goto QUIT;
  
  /* Integrate new variables */
  
  error = GRBupdatemodel(submdp->grb_model);
  if (error) goto QUIT;

  printf("Creating constraints ... \n");
  fflush(stdout);
  /* Create the constraints */

  ul num_nonzero_constraint_coefficients = 0;
  double constraint_rhs = 0;
  
  submdp->state_tree_tail = submdp->state_tree;
  /* int num_constraints = 0; */


  /* Add the constraints per state, for all states */

  // <NOTE>
  // By BNA
  // ======
  // state_tree is actually a list, and states in this list are by
  // their indexes (increasing). Since the constraints below
  // are added per state in
  // the order of their appearance in state_tree and thus their index,
  // querying the grb model for the shadow price of the constraint whose index
  // is i corresponds to querying for the value of the  variable
  // (and hence the state) whose index is also i.
  // </NOTE>

  while(submdp->state_tree_tail != NULL) {
    
    state_action_pair *admissible_saps_list =
      submdp->state_tree_tail->admissible_state_action_pairs_list;
    
    num_nonzero_constraint_coefficients = 0;
    
    while(admissible_saps_list != NULL) {
      
      arr_variable_indexes[num_nonzero_constraint_coefficients] =
	admissible_saps_list->index;
      
      arr_constraint_coefficients[num_nonzero_constraint_coefficients] = 1.0;
      ++num_nonzero_constraint_coefficients;
      admissible_saps_list = admissible_saps_list->next;
      
      
    } 
    
    
    state_action_pair_pointer *in_sap_ptrs_list_tail =
      submdp->state_tree_tail->in_state_action_pairs_list;
    
    while(in_sap_ptrs_list_tail != NULL) {
      
      arr_variable_indexes[num_nonzero_constraint_coefficients] =
	in_sap_ptrs_list_tail->sap->index;
      
      arr_constraint_coefficients[num_nonzero_constraint_coefficients] =
	-1 * in_sap_ptrs_list_tail->prob_parent_state_given_sap;
      
      ++num_nonzero_constraint_coefficients;
      in_sap_ptrs_list_tail = in_sap_ptrs_list_tail->next;
      
    } 

    constraint_rhs = submdp->state_tree_tail->is_initial ? 1 : 0;
    
    error = GRBaddconstr(submdp->grb_model,
    			 num_nonzero_constraint_coefficients,
    			 arr_variable_indexes,
    			 arr_constraint_coefficients,
    			 GRB_EQUAL,
    			 constraint_rhs,
    			 NULL);
    if (error) goto QUIT;
    
    /* ++num_constraints; */
    /* printf( "Added constraint %d\n", num_constraints); */
    
    submdp->state_tree_tail = submdp->state_tree_tail->next;
  } 
  
  
  
  
  /* Optimize model */

  printf( "Solving LP ... \n");
  fflush(stdout);
  
  error = GRBoptimize(submdp->grb_model);
  if (error) goto QUIT;
  
  /* Write model to 'lp.lp' */
  
  error = GRBwrite(submdp->grb_model, "lp.lp");
  if (error) goto QUIT;
  
  /* Capture solution information */
  
  error = GRBgetintattr(submdp->grb_model, GRB_INT_ATTR_STATUS, &optimization_status);
  if (error) goto QUIT;

  if (optimization_status == GRB_OPTIMAL) {
    error = GRBgetdblattr(submdp->grb_model,
			  GRB_DBL_ATTR_OBJVAL,
			  &objective_value);
    
    if (error) goto QUIT;
    printf( "Optimal objective: %.4e\n", objective_value);
    fflush(stdout);
    
    /* error = GRBgetdblattrarray(submdp->grb_model, */
    /* 			       GRB_DBL_ATTR_X, */
    /* 			       0, */
    /* 			       submdp->num_state_action_pairs, */
    /* 			       submdp->arr_solution); */
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

  
 QUIT:
  
  /* Error reporting */
  
  if (error) {
    printf( "ERROR: %s\n", GRBgeterrormsg(master_env));
    exit(1);
  }
  
  /* Free model */
  
  /* GRBfreemodel(model); */
  
  /* Free environment */
  
  /* GRBfreeenv(master_env); */
  
  
}

/* Assumes: the required sap exists in submdp given the params. */
/* Here the vectors are for the 2 jobs as they appearin the submdp. 
   submdp_action is index of a job relative to the 
   submdp; i.e., 0 or 1, and the global job index in arr_jobs must be 
   determined */ 
state_action_pair *lookup_submdp_sap_by_params(submdp_t *submdp,
					       ul *submdp_allocation_vector,
					       char *submdp_finish_signal_vector,
					       int submdp_action) {
  state *submdp_state =
    lookup_state_ptr_by_params(submdp,
			       NON_ABSORBING,
			       NULL,
			       submdp_allocation_vector,
			       submdp_finish_signal_vector);
  
  assert(submdp_state != NULL && submdp_state != 0);
  state_action_pair *admissible_saps_list =
    submdp_state->admissible_state_action_pairs_list;

  state_action_pair *sap = NULL;
  
  while(admissible_saps_list != NULL) {
    if(admissible_saps_list->action == submdp_action) {
      sap = admissible_saps_list;
      break;
    }
    admissible_saps_list = admissible_saps_list->next;
  }
  
  assert(sap != NULL);
  admissible_saps_list = NULL;
  
  return sap;
  
}

/* here the vectors are for the n jobs as they appear in arr_jobs
   and not per submdp */ 
int get_num_legal_error_flags_for_crit(int time,
				       ul *allocation_vector,
				       char *finish_signal_vector,
				       int crit) {
  
  int num_legal_error_flags = 0;

  int num_jobs_of_crit = arr_num_jobs_of_criticality[crit];
  
  /* Loop over the deadlines of crit-criticality jobs only */
  int index; 
  for(index = 0; index < num_jobs_of_crit; ++index) {
    
    int job_index = arr_job_indexes_per_criticality[crit][index];

    if(finish_signal_vector[job_index] == NOT_FINISHED) {
      if(time >= arr_jobs[job_index].deadline) {
	/* either error or error', but one of them exclusively */
	num_legal_error_flags = 1;
	break;
      }
    }
  }

  if(num_legal_error_flags == 0) {/* time < d_i for every i for which y_i = 1*/
    
    for(index = 0; index < num_jobs_of_crit; ++index) {
      
      int job_index = arr_job_indexes_per_criticality[crit][index];
      
      if(finish_signal_vector[job_index] == FINISHED) {
	if(time > arr_jobs[job_index].deadline) {
	  num_legal_error_flags = 2;
	  break;
	}
      }
    }
    
    if(num_legal_error_flags == 0) { /* either no jobs finished at time or 
					time <= deadline of every job that 
					have finished */
      num_legal_error_flags = 1;
    }
    
  }
  
  return num_legal_error_flags;
  
  
}

/* m: variable (submdp) index
   l: criticality level */
double compute_var_coef_in_crit_constraint(int m, int crit) {
  
  double coef = 0.0;
  /* must loop through every t in {0,.., submdp m's horizon} */
  /* Consider only the states in S'; i.e., where there is an l-crit job 
     that has not finished (y^i = 1) */
  
  submdp_t *submdp = arr_submdps[m];
  int job1_index = submdp->arr_job_indexes[0];
  int job2_index = submdp->arr_job_indexes[1];
  

  job_t *job1 = lookup_job_by_submdp_job_index(submdp, 0);
  job_t *job2 = lookup_job_by_submdp_job_index(submdp, 1);
  
  int num_jobs = NUM_JOBS;
  
  int time;
  int alloc1;
  int alloc2;
  int num_jobs_of_crit = arr_num_jobs_of_criticality[crit];
  
  ul *allocation_vector       = calloc(num_jobs, sizeof(ul)); 
  char *finish_signal_vector  = calloc(num_jobs, sizeof(char));
  
  /* Loop over the deadlines of crit-criticality jobs only */
  int index;
  for(index = 0; index < num_jobs_of_crit; ++index) {
    
    int job_index = arr_job_indexes_per_criticality[crit][index];
    time = arr_jobs[job_index].deadline;
    
    if(time > submdp->horizon) { continue; }
    
    /* now generate pairs of execution times whose sum is t */
    if((job1->wcet + job2->wcet) < time) { continue; }
    
    
    /* here we are looping the submdp actions (0 and 1) only because 
       rho(..., a) = 0 for other actions. submdp_action is {0,1} and it 
       corresponds to the action's index from the submdp prespective; 
       the actual job corrresponding to submdp_action in arr_jobs is
       submdp->arr_job_indexes[submdp_action] */
    
    for(alloc1 = 0; alloc1 <= imin((int)time, (int)job1->wcet); ++alloc1) {
      /* printf("%d, %d\n", a, imin(sum - a, b_limit)); */
      alloc2 = imin((int)time - alloc1, (int)job2->wcet);
      
      
      if(alloc1 + alloc1 != time) { continue; }
      
      /* now we have two legal allocations alloc1 and alloc2 for jobs
	 job1 and job2 of the submdp, respectively. */
      
      allocation_vector[job1_index] = alloc1;
      allocation_vector[job2_index] = alloc2;
      
      int i;
      for(i = 0; i < num_jobs; ++i) {
	
	if(i == job1_index || i == job2_index) { continue; }
	
	finish_signal_vector[i] = NOT_FINISHED;
	
      }
      
      /* Now loop through all the possible values for the finish signals
	 of job1 and job2, which are (0, 0), (0, 1), (1,0), (1,1) */
      
      int num_possible_finish_signals_job1 = 0;
      int num_possible_finish_signals_job2 = 0;
      
      char *arr_possible_finish_signals_job1;
      char *arr_possible_finish_signals_job2;
      
      /* If allocn = 0 for any n, then jobn cannot have a finish signal
	 of FINISHED. Also if allocn = wcetn for any n, then job n cannot
	 have finish signal NOT_FINISHED */
      if(alloc1 == 0 || alloc1 == job1->wcet) {
	
	num_possible_finish_signals_job1 = 1;
	arr_possible_finish_signals_job1 =
	  calloc(num_possible_finish_signals_job1, sizeof(char));
	
	if(alloc1 == 0) {
	  arr_possible_finish_signals_job1[0] = NOT_FINISHED;
	} else {
	  arr_possible_finish_signals_job1[0] = FINISHED;
	}
	
      } else { /* 0 < alloc < wcet */
	
	num_possible_finish_signals_job1 = 2;
	arr_possible_finish_signals_job1 =
	  calloc(num_possible_finish_signals_job1, sizeof(char));
	arr_possible_finish_signals_job1[0] = NOT_FINISHED;
	arr_possible_finish_signals_job1[1] = FINISHED;
	
      }	
      
      if(alloc2 == 0 || alloc2 == job2->wcet) {
	
	num_possible_finish_signals_job2 = 1;
	arr_possible_finish_signals_job2 =
	  calloc(num_possible_finish_signals_job2, sizeof(char));
	
	if(alloc2 == 0) {
	  arr_possible_finish_signals_job2[0] = NOT_FINISHED;
	} else {
	  arr_possible_finish_signals_job2[0] = FINISHED;
	}
	
      } else { /* 0 < alloc < wcet */
	
	num_possible_finish_signals_job2 = 2;
	arr_possible_finish_signals_job2 =
	  calloc(num_possible_finish_signals_job2, sizeof(char));
	arr_possible_finish_signals_job2[0] = NOT_FINISHED;
	arr_possible_finish_signals_job2[1] = FINISHED;
	
      }	
      
      
      int j;
      for(i = 0; i < num_possible_finish_signals_job1; i++) {
	for(j = 0; j < num_possible_finish_signals_job2; j++) {
	  
	  finish_signal_vector[job1_index] =
	    arr_possible_finish_signals_job1[i];
	  
	  finish_signal_vector[job2_index] =
	    arr_possible_finish_signals_job2[j];
	  
	  ul *submdp_allocation_vector = calloc(NUM_JOBS_SUBMDP, sizeof(ul));
	  char *submdp_finish_signal_vector =
	    calloc(NUM_JOBS_SUBMDP, sizeof(char));
	  
	  submdp_allocation_vector[0] = allocation_vector[job1_index];
	  submdp_allocation_vector[1] = allocation_vector[job2_index];
	  
	  submdp_finish_signal_vector[0] = finish_signal_vector[job1_index];
	  submdp_finish_signal_vector[1] = finish_signal_vector[job2_index];
	  
	  
	  int submdp_action;
	  for(submdp_action = 0;
	      submdp_action < NUM_JOBS_SUBMDP;
	      ++submdp_action) {
	    
	    if(!is_admissible_action(submdp_allocation_vector,
				     submdp_finish_signal_vector,
				     submdp_action)) {
	      continue;
	    }
	    
	    state_action_pair* submdp_sap =
	      lookup_submdp_sap_by_params(submdp,
					  submdp_allocation_vector,
					  submdp_finish_signal_vector,
					  submdp_action);
	    
	    double submdp_variable_solution =
	      lookup_dual_variable_solution_by_index(submdp, submdp_sap->index);

	    
	    /* fflush(stdout); /\* debug *\/ */
	    
	    if(essentially_equal(submdp_variable_solution, 0)) { continue; }
	    
	    int num_legal_error_flags = 0;
	    
	    int crit_level;
	    for(crit_level = 0; crit_level < NUM_CRITICALITY_LEVELS;
		++crit_level){
	      
	      if (crit_level != crit) { /* because r^{crit} = error is fixed */
		
		num_legal_error_flags +=
		  get_num_legal_error_flags_for_crit(time,
						     allocation_vector,
						     finish_signal_vector,
						     crit_level);
	      }	
	      
	    }
	    
	    coef += (double) num_legal_error_flags * submdp_variable_solution;
	    
	  }

	  free(submdp_allocation_vector);
	  free(submdp_finish_signal_vector);
	  
	}
	
      }
    }
  }
  // TODO: free arr_possible_finish_signals_job1 and ..._job2
  free(allocation_vector);
  free(finish_signal_vector);
  
  return coef; 
  
}

double compute_var_coef_in_obj(int m) {

  double coef = 0.0;
  /* must loop through every t in {0,.., submdp m's horizon} */
  /* Consider only the states in S'; i.e., where there is a  job 
     that has not finished (y^i = 1) */
  
  submdp_t *submdp = arr_submdps[m];
  int job1_index = submdp->arr_job_indexes[0];
  int job2_index = submdp->arr_job_indexes[1];
  
  
  job_t *job1 = lookup_job_by_submdp_job_index(submdp, 0);
  job_t *job2 = lookup_job_by_submdp_job_index(submdp, 1);
  
  int num_jobs = NUM_JOBS;
  
  int time;
  int alloc1;
  int alloc2;
  
  ul *allocation_vector       = calloc(num_jobs, sizeof(ul)); 
  char *finish_signal_vector  = calloc(num_jobs, sizeof(char));
  
  /* Loop over all n jobs */
  
  int global_job_index;
  for(global_job_index = 0; global_job_index < NUM_JOBS; ++global_job_index){
    
    time = arr_jobs[global_job_index].deadline;
    
    if(time > submdp->horizon) { continue; }
    
    /* now generate pairs of execution times whose sum is t */
    if((job1->wcet + job2->wcet) < time) { continue; }
    
    
    /* here we are looping the submdp actions (0 and 1) only because 
       rho(..., a) = 0 for other actions. submdp_action is {0,1} and it 
       corresponds to the action's index from the submdp prespective; 
       the actual job corrresponding to submdp_action in arr_jobs is
       submdp->arr_job_indexes[submdp_action] */
    
    for(alloc1 = 0; alloc1 <= imin((int)time, (int)job1->wcet); ++alloc1) {
      /* printf("%d, %d\n", a, imin(sum - a, b_limit)); */
      alloc2 = imin((int)time - alloc1, (int)job2->wcet);
      
      
      if(alloc1 + alloc1 != time) { continue; }
      
      /* now we have two legal allocations alloc1 and alloc2 for jobs
	 job1 and job2 of the submdp, respectively. */
      
      allocation_vector[job1_index] = alloc1;
      allocation_vector[job2_index] = alloc2;
      
      int i;
      for(i = 0; i < num_jobs; ++i) {
	
	if(i == job1_index || i == job2_index) { continue; }
	
	finish_signal_vector[i] = NOT_FINISHED;
	
      }
      
      /* Now loop through all the possible values for the finish signals
	 of job1 and job2, which are (0, 0), (0, 1), (1,0), (1,1) */
      
      int num_possible_finish_signals_job1 = 0;
      int num_possible_finish_signals_job2 = 0;
      
      char *arr_possible_finish_signals_job1;
      char *arr_possible_finish_signals_job2;
      
      if(global_job_index == job1_index) {
	if(alloc1 == job1->wcet) {
	  continue;
	} else {
	  num_possible_finish_signals_job1 = 1;
	  arr_possible_finish_signals_job1 =
	    calloc(num_possible_finish_signals_job1, sizeof(char));
	  arr_possible_finish_signals_job1[0] = NOT_FINISHED;
	}
      }
      
      if(global_job_index == job2_index) {
	if(alloc2 == job2->wcet) {
	  continue;
	} else {
	  num_possible_finish_signals_job2 = 1;
	  arr_possible_finish_signals_job2 =
	    calloc(num_possible_finish_signals_job2, sizeof(char));
	  arr_possible_finish_signals_job2[0] = NOT_FINISHED;
	}
      }
      
      
      
      char are_possible_finish_signals_determined_job1 =
	num_possible_finish_signals_job1 > 0;

      if(!are_possible_finish_signals_determined_job1) {
	/* If allocn = 0 for any n, then jobn cannot have a finish signal
	   of FINISHED. Also if allocn = wcetn for any n, then job n cannot
	   have finish signal NOT_FINISHED */
	if(alloc1 == 0 || alloc1 == job1->wcet) {
	  
	  num_possible_finish_signals_job1 = 1;
	  arr_possible_finish_signals_job1 =
	    calloc(num_possible_finish_signals_job1, sizeof(char));
	  
	  if(alloc1 == 0) {
	    arr_possible_finish_signals_job1[0] = NOT_FINISHED;
	  } else {
	    arr_possible_finish_signals_job1[0] = FINISHED;
	  }
	  
	} else { /* 0 < alloc < wcet */
	  
	  num_possible_finish_signals_job1 = 2;
	  arr_possible_finish_signals_job1 =
	    calloc(num_possible_finish_signals_job1, sizeof(char));
	  arr_possible_finish_signals_job1[0] = NOT_FINISHED;
	  arr_possible_finish_signals_job1[1] = FINISHED;
	  
	}	
      }


      char are_possible_finish_signals_determined_job2 =
	num_possible_finish_signals_job2 > 0;

      if(!are_possible_finish_signals_determined_job2) {
	if(alloc2 == 0 || alloc2 == job2->wcet) {
	  
	  num_possible_finish_signals_job2 = 1;
	  arr_possible_finish_signals_job2 =
	    calloc(num_possible_finish_signals_job2, sizeof(char));
	  
	  if(alloc2 == 0) {
	    arr_possible_finish_signals_job2[0] = NOT_FINISHED;
	  } else {
	    arr_possible_finish_signals_job2[0] = FINISHED;
	  }
	  
	} else { /* 0 < alloc < wcet */
	  
	  num_possible_finish_signals_job2 = 2;
	  arr_possible_finish_signals_job2 =
	    calloc(num_possible_finish_signals_job2, sizeof(char));
	  arr_possible_finish_signals_job2[0] = NOT_FINISHED;
	  arr_possible_finish_signals_job2[1] = FINISHED;
	  
	}	
      }
	
      int j;
      for(i = 0; i < num_possible_finish_signals_job1; i++) {
	for(j = 0; j < num_possible_finish_signals_job2; j++) {
	  
	  finish_signal_vector[job1_index] =
	    arr_possible_finish_signals_job1[i];
	  
	  finish_signal_vector[job2_index] =
	    arr_possible_finish_signals_job2[j];
	  
	  ul *submdp_allocation_vector = calloc(NUM_JOBS_SUBMDP, sizeof(ul));
	  char *submdp_finish_signal_vector =
	    calloc(NUM_JOBS_SUBMDP, sizeof(char));
	  
	  submdp_allocation_vector[0] = allocation_vector[job1_index];
	  submdp_allocation_vector[1] = allocation_vector[job2_index];
	  
	  submdp_finish_signal_vector[0] = finish_signal_vector[job1_index];
	  submdp_finish_signal_vector[1] = finish_signal_vector[job2_index];
	  
	  
	  int submdp_action;
	  for(submdp_action = 0;
	      submdp_action < NUM_JOBS_SUBMDP;
	      ++submdp_action) {
	    
	    if(!is_admissible_action(submdp_allocation_vector,
				     submdp_finish_signal_vector,
				     submdp_action)) {
	      continue;
	    }
	    
	    state_action_pair* submdp_sap =
	      lookup_submdp_sap_by_params(submdp,
					  submdp_allocation_vector,
					  submdp_finish_signal_vector,
					  submdp_action);
	    
	    double submdp_variable_solution =
	      lookup_dual_variable_solution_by_index(submdp, submdp_sap->index);
	    
	    if(essentially_equal(submdp_variable_solution, 0)) { continue; }
	    
	    int num_legal_error_flags = 0;
	    
	    int crit_level;
	    for(crit_level = 0; crit_level < NUM_CRITICALITY_LEVELS;
		++crit_level){
	      
	      num_legal_error_flags +=
		get_num_legal_error_flags_for_crit(time,
						   allocation_vector,
						   finish_signal_vector,
						   crit_level);
	      
	    }
	    
	    coef += (double) num_legal_error_flags * submdp_variable_solution;
	    
	  }
	  
	  free(submdp_allocation_vector);
	  free(submdp_finish_signal_vector);
	  
	}
	
      }
    }
  }
  // TODO: free arr_possible_finish_signals_job1 and ..._job2
  free(allocation_vector);
  free(finish_signal_vector);
  
  return coef; 
  
  
}

double arr_error_upper_bounds[NUM_CRITICALITY_LEVELS] = {0.3, 0.4, 0.5};

void build_approximate_lp() {

  GRBenv   *master_env = NULL;
  GRBmodel *model      = NULL;
  int       error      = 0;
  
  int num_submpds = NUM_JOB_PAIRS;
  
  /* submdp->arr_solution = calloc(num_submpds, sizeof(double)); */
 
  int *arr_variable_indexes = calloc(num_submpds, sizeof(int));

  double *arr_objective_coefficients = calloc(num_submpds,  sizeof(double));
  
  double *arr_constraint_coefficients = calloc(num_submpds, sizeof(double));

  
  int       optimization_status;
  double    objective_value;
  
  /* Create environment */
  
  error = GRBloadenv(&master_env, "approx-lp.log");
  if (error) goto QUIT;
  
  /* Create an empty model */
  
  error = GRBnewmodel(master_env, &model, "lp", 0,
		      NULL, NULL, NULL, NULL, NULL);
  if (error) goto QUIT;

  int num_crit_levels = NUM_CRITICALITY_LEVELS; 

  printf( "Creating objective coefficients array for approximate LP ... \n");
  fflush(stdout);

  double coef_var_m_in_obj = 0.0;
  
  /* m: the index of variable v_m in the approximate LP CDAULP */
  int m;
  for(m = 0; m < num_submpds; ++m) {
    
    arr_objective_coefficients[m] = compute_var_coef_in_obj(m);
    
    fprintf(fdebug, "coef of variable m = %d in objective: %f\n",
	    m, arr_objective_coefficients[m]);
    
  }
  
  
  printf( "Adding variables ... \n");
  fflush(stdout);
  
  error =
    GRBXaddvars(model,
		num_submpds,
		0, NULL, NULL, NULL,
		arr_objective_coefficients, NULL, NULL, NULL, NULL);
  if (error) goto QUIT;

  
  double constraint_rhs = 0.0;
  
  /* Compute the coefficients of v_m in the L constraints kDv <= epsilon */
  int l;
  int num_nonzero_var_coefs = 0;
  double coef_var_m_in_l_constraint = 0.0;
  
  for(l = 0; l < num_crit_levels; ++l) {

    num_nonzero_var_coefs = 0;
    
    for(m = 0; m < num_submpds; ++m) {
      
      
      coef_var_m_in_l_constraint = compute_var_coef_in_crit_constraint(m,l);

      /* fprintf(fdebug, "coef of variable m = %d in constraint l = %d: %f\n", coef_var_m_in_l_constraint, m,l); */
      
      if(definitely_greater_than(coef_var_m_in_l_constraint, 0.0)) {
	
	arr_variable_indexes[num_nonzero_var_coefs++] = m;
	
	arr_constraint_coefficients[num_nonzero_var_coefs] =
	  coef_var_m_in_l_constraint;
      }
      
    }
    
    constraint_rhs = arr_error_upper_bounds[l];


    
    /* add the constraint */
    error = GRBaddconstr(model,
    			 num_nonzero_var_coefs,
    			 arr_variable_indexes,
    			 arr_constraint_coefficients,
    			 GRB_LESS_EQUAL,
    			 constraint_rhs,
    			 NULL);
    if (error) goto QUIT;
    
  }

  printf("\nAdded the %d criticality constraints for the approximate LP\n\n",
	 NUM_CRITICALITY_LEVELS); /* debug */
  
  /* Compute the coefficients of v_m for every m in the K_D constraints 
     B^{T}A^{T}Dv = B^{T}e_{s_0} */
  int k;
  for(k = 0; k < num_submpds; ++k) {
    
    num_nonzero_var_coefs = 0;
    
    for(m = 0; m < num_submpds; ++m) {
      
      double coef_var_m_in_k_constraint = 0.0;
      
      if(definitely_greater_than(coef_var_m_in_k_constraint, 0.0)) {
	
	arr_variable_indexes[num_nonzero_var_coefs++] = m;
	
	arr_constraint_coefficients[num_nonzero_var_coefs] =
	  coef_var_m_in_k_constraint;
      }
      
    }
    
    constraint_rhs = arr_error_upper_bounds[l]; //TODO: change
    
    /* add the constraint */
    error = GRBaddconstr(model,
    			 num_nonzero_var_coefs,
    			 arr_variable_indexes,
    			 arr_constraint_coefficients,
    			 GRB_EQUAL,
    			 constraint_rhs,
    			 NULL);
    if (error) goto QUIT;
    
  }
  
  
  
  
  
  /* Change objective sense to minimization */
  
  error = GRBsetintattr(model,
			GRB_INT_ATTR_MODELSENSE,
			GRB_MINIMIZE);
  if (error) goto QUIT;
  
  GRBenv *model_env = GRBgetenv(model);
  
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_THREADS, 10); */
  /* if (error) goto QUIT; */
  
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_CROSSOVER, 0); */
  /* if (error) goto QUIT; */
  
  /* GRB_INT_PAR_METHOD = 0: Primal Simplex */
  /*                      1:  Simplex */
  /*                      2: Barrier */
  /*                      3: Nondeterministic Concurrent */
  /*                      4: Deterministic Concurrent */
  /* error = GRBsetintparam(model_env, GRB_INT_PAR_METHOD, 0); */
  
  if (error) goto QUIT;
  
  /* Integrate new variables */
  
  error = GRBupdatemodel(model);
  if (error) goto QUIT;

  printf("Creating constraints ... \n");
  fflush(stdout);

  
  
  /* Optimize model */

  printf( "Solving LP ... \n");
  fflush(stdout);
  
  error = GRBoptimize(model);
  if (error) goto QUIT;
  
  /* Write model to 'approx-lp.lp' */
  
  error = GRBwrite(model, "approx-lp.lp");
  if (error) goto QUIT;
  
  /* Capture solution information */
  
  error = GRBgetintattr(model,
			GRB_INT_ATTR_STATUS,
			&optimization_status);
  if (error) goto QUIT;
  
  if (optimization_status == GRB_OPTIMAL) {
    error = GRBgetdblattr(model,
			  GRB_DBL_ATTR_OBJVAL,
			  &objective_value);
    
    if (error) goto QUIT;
    printf( "Optimal objective: %.4e\n", objective_value);
    fflush(stdout);
    
    /* error = GRBgetdblattrarray(submdp->grb_model, */
    /* 			       GRB_DBL_ATTR_X, */
    /* 			       0, */
    /* 			       submdp->num_state_action_pairs, */
    /* 			       submdp->arr_solution); */
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

  
 QUIT:
  
  /* Error reporting */
  
  if (error) {
    printf( "ERROR: %s\n", GRBgeterrormsg(master_env));
    exit(1);
  }
  
  /* Free model */
  
  /* GRBfreemodel(model); */
  
  /* Free environment */
  
  /* GRBfreeenv(master_env); */
  

  
  
}

submdp_t *create_submdp() {
  
  submdp_t *submdp = (submdp_t *)malloc(sizeof(submdp_t));
  submdp->horizon = 0;

  submdp->state_tree      = NULL;
  submdp->state_tree_tail = NULL;
  
  /* submdp->state_action_pair_ptrs_list      = NULL; */
  /* submdp->state_action_pair_ptrs_list_tail = NULL; */
  
  submdp->num_absorbing_states    = 0;
  submdp->num_nonabsorbing_states = 0;
    
  submdp->num_state_action_pairs = 0;

  return submdp;
  
}

void setup_jobs(ransampl_ws **arr_exec_times_ws) {
  
  
  for(int crit_index = 0; crit_index < NUM_CRITICALITY_LEVELS; ++crit_index){
    arr_num_jobs_of_criticality[crit_index] = 0;
  }
  
  int job_index;
  for(job_index = 0; job_index < NUM_JOBS; ++job_index) {
    
    job_t *job = &arr_jobs[job_index];
    ul wcet = (ul)job->wcet;
    arr_exec_times_ws[job_index] = ransampl_alloc( wcet );
    ransampl_set( arr_exec_times_ws[job_index], job->execution_time_pmf );
    
    arr_num_jobs_of_criticality[job->criticality] += 1;
    
  }
  
  for(int crit_index = 0; crit_index < NUM_CRITICALITY_LEVELS;++crit_index){
    
    arr_job_indexes_per_criticality[crit_index] =
      malloc(arr_num_jobs_of_criticality[crit_index] * sizeof(int));
  }

  char *arr_num_jobs_of_criticality_copy =
    clone_char_array(arr_num_jobs_of_criticality, NUM_CRITICALITY_LEVELS);
  
  for(job_index = 0; job_index < NUM_JOBS; ++job_index){
    
    job_t *j = &arr_jobs[job_index];
    int criticality = j->criticality;
    arr_job_indexes_per_criticality[criticality]
      [arr_num_jobs_of_criticality_copy[criticality] - 1] = job_index;
    --arr_num_jobs_of_criticality_copy[criticality];
    
  }
 
  free(arr_num_jobs_of_criticality_copy);
  
}

int *generate_int_pairs_for_sum(int a_limit, int b_limit, int sum) {

  int a;
  int b;
  
  if(a_limit + b_limit < sum) { return NULL; }

  
  
  for(a = 0; a <= imin(sum, a_limit); ++a) {
    printf("%d, %d\n", a, imin(sum - a, b_limit));
  }
}

int play() {
  
  int sum = 250;
  int a_limit = 100;
  int b_limit = 200;

  int a;
  int b;

  if(a_limit + b_limit < sum) {return;}
  
  for(a = 0; a <= imin(sum, a_limit); ++a) {
    b = imin(sum - a, b_limit);
    if(a + b == sum) {
      printf("%d, %d\n", a, b);
    }
  }
  //back
  /*       for(b = 0; b <= fmin(sum - a, arr_jobs[1].wcet); ++b) */
  
}
	     
int main(int argc, char *argv[]){
  /* play(); exit(0); */
  
  /* Input jobset */

  fdebug = fopen ("debug.txt", "w+");
  
  job_t *job = &arr_jobs[0];
  job->deadline    = 120;
  job->criticality = 1;
  job->wcet = 100;
  job->execution_time_pmf = generate_execution_time_pmf(job->wcet);
  
  job = &arr_jobs[1];
  job->deadline    = 250;
  job->criticality = 0;
  job->wcet = 200;
  job->execution_time_pmf = generate_execution_time_pmf(job->wcet);

  job = &arr_jobs[2];
  job->deadline    = 200;
  job->criticality = 2;
  job->wcet = 500;
  job->execution_time_pmf = generate_execution_time_pmf(job->wcet);

  
  job = &arr_jobs[3];
  job->deadline    = 200;
  job->criticality = 2;
  job->wcet = 150;
  job->execution_time_pmf = generate_execution_time_pmf(job->wcet);
  
  print_job_set(); 
  
  /* bool is_ocbp_schedulable = compute_OCBP_priorities(); */
  
  
  // Initialize random number generator:
  gsl_rng_env_setup();
  gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
  gsl_rng_set(rng, time(NULL));
  
  /* Create job execution time distribiutions */
  ransampl_ws *arr_exec_times_ws[NUM_JOBS];

  setup_jobs(arr_exec_times_ws);
  
  int submdp_index = 0;
  
  int i = 0;
  int j = 0;
  for (i = 0; i < NUM_JOBS; ++i) {
    for (j = i + 1; j < NUM_JOBS; ++j){
      submdp_t *submdp = create_submdp();
      submdp->arr_job_indexes[0] = i;
      submdp->arr_job_indexes[1] = j;
      
      
      initial_allocation_vector     = calloc(NUM_JOBS_SUBMDP, sizeof(ul));
      initial_finish_signal_vector  = calloc(NUM_JOBS_SUBMDP, sizeof(char));
      /* initial_error_flag            = NOT_ERROR; */
      
      submdp->horizon = 0;
           
      ul state_array_size = 1;
      
      int index;
      for(index = 0; index < NUM_JOBS_SUBMDP; ++index) {
	
	initial_finish_signal_vector[index] = NOT_FINISHED;
	
	job_t *job = lookup_job_by_submdp_job_index(submdp, index);
	ul wcet = (ul)job->wcet;
	submdp->horizon += wcet;
	submdp->dims[index] = wcet + 1;
	state_array_size *= submdp->dims[index];
	
      }
      
      for(index = NUM_JOBS_SUBMDP; index < NUM_DIMS; ++index) {
	submdp->dims[index] = MAX_FINISH_SIGNAL;
	state_array_size *= submdp->dims[index];
      }
      /* dims[ NUM_DIMS - 1 ] = MAX_ERROR_FLAG; */
      
      /* state_array_size *= dims[ NUM_DIMS - 1 ]; */
      
      submdp->arr_nonabsorbing_states =
	calloc(state_array_size, sizeof(state*));

      /* There is an absorbing state for each 
	 t = NUM_JOBS_SUBMDP + 1, ..., horizon + 1, so the total 
	 is (horizon + 1) - (NUM_JOBS_SUBMDP + 1) + 1 = 
	 horizon - NUM_JOBS_SUBMDP + 1 */
      submdp->arr_absorbing_states =
	calloc(submdp->horizon - NUM_JOBS_SUBMDP + 1, sizeof(state*));

      
      print_submdp(submdp, submdp_index + 1);
      
      printf( "Creating state tables ...\n");
      fflush(stdout);
      
      build_state_tree(submdp);
      
      setlocale(LC_NUMERIC, "");
      printf( "Building LP\nNumber of variables: %'lu\nNumber of constraints: %'lu\n", submdp->num_state_action_pairs, get_total_submdp_state_count(submdp));
      fflush(stdout);
      
      solve_unconstrained_submdp_exact(submdp);

      
      
      double primal_solution =
	lookup_primal_variable_solution_by_index(submdp, 0);
      printf("dual variable-0 value: %f\n", primal_solution);
      
      arr_submdps[submdp_index] = submdp;
      submdp_index++;
      
      
    }
  }

  build_approximate_lp();
  
  /* compute the coefficient of CDAULP */
  
  
  /**---------------------------Begin Simulation------------------------**/
  
  
  /* const int num_sim_samples = 100000; */
  
  /* sim_state_t *sim_state = malloc(sizeof(sim_state_t)); */
  /* sim_stats_t *sim_stats  = malloc(sizeof(sim_stats_t)); */
  
  /* sim_stats->arr_num_deadline_misses_per_job = */
  /*   calloc(NUM_JOBS, sizeof(int)); */
  
  /* /\* sim_stats->arr_num_errors_per_criticality = *\/ */
  /* /\*   calloc(NUM_CRITICALITY_LEVELS, sizeof(int)); *\/ */
  
  /* sim_stats->num_errors = 0; */
  
  /* printf( "\n\nSimulating execution over %d demand samples (vectors) \n", */
  /* 	 num_sim_samples); */
  /* fflush(stdout); */
  
  /* int sample_index; */
  
  /* for(sample_index = 0; sample_index < num_sim_samples; ++sample_index) { */
    
  /*   int arr_exec_time_samples[NUM_JOBS]; */
  /*   for(int index = 0; index < NUM_JOBS; ++index) { */
  /*     /\* ransampl_draw(...) generates an index in {0, ..., max_wcet - 1} *\/ */
  /*     arr_exec_time_samples[index] = */
  /* 	ransampl_draw( arr_exec_times_ws[index], */
  /* 		       gsl_rng_uniform(rng), */
  /* 		       gsl_rng_uniform(rng) ) + 1; */
      
  /*   } */

  /*   /\* print_exec_time_sample(arr_exec_time_samples); *\/ */
    
  /*   /\* printf( "sample ok ... \n"); *\/ */

  /*   // <TODO> */
  /*   // By BNA */
  /*   // Change vectors to be for the whole job system; i.e., of dim n */
  /*   // </TODO> */
  /*   sim_state->error_flag	    = initial_error_flag; */
  /*   sim_state->allocation_vector    = */
  /*     clone_ul_array(initial_allocation_vector, NUM_JOBS); */
  /*   sim_state->finish_signal_vector = */
  /*     clone_char_array(initial_finish_signal_vector, NUM_JOBS); */
    
  /*   state *s; */
  /*   ul sim_time = 0; */
    
  /*   bool *arr_job_missed_deadlines = calloc(NUM_JOBS, sizeof(bool)); */
    
  /*   /\* simulation loop *\/ */
  /*   while(!have_all_jobs_finished(sim_state->finish_signal_vector)) { */
      
      
      
  /*     /\* print_state(sim_time, *\/ */
  /*     /\* 		  sim_state->allocation_vector, *\/ */
  /*     /\* 		  sim_state->finish_signal_vector, *\/ */
  /*     /\* 		  sim_state->error_flag); *\/ */
      
  /*     /\* Lookup all the admissible state-action pairs for state s *\/ */

  /*     /\** TODO LATER: do this for each submdp *\/ */
  /*     s = lookup_state_ptr_by_params(sim_state->allocation_vector, */
  /* 				     sim_state->finish_signal_vector, */
  /* 				     sim_state->error_flag); */
      
  /*     assert(s != NULL); */
      
  /*     state_action_pair *admissible_saps_list_tail = */
  /* 	s->admissible_state_action_pairs_list; */
      
  /*     double *sap_pmf = malloc((s->num_admissible_saps) * sizeof(double)); */
  /*     state_action_pair **admissible_sap_pointers = */
  /* 	malloc((s->num_admissible_saps) * sizeof(state_action_pair*)); */
      
  /*     int index = 0; */
  /*     double sap_lp_solution = 0; */
  /*     double sum_admissible_saps_lp_solutions = 0; */
      
  /*     while(admissible_saps_list_tail != NULL) { */
  /* 	admissible_sap_pointers[index] = admissible_saps_list_tail; */
	
  /* 	sap_lp_solution = arr_dual_solution[admissible_saps_list_tail->index]; */
	
  /* 	/\* printf( "sap solution: %.*lf\n", 200, sap_lp_solution); *\/ */
	
  /* 	sap_pmf[index] = sap_lp_solution; */
  /* 	sum_admissible_saps_lp_solutions += sap_lp_solution; */
  /* 	++index; */
  /* 	admissible_saps_list_tail = admissible_saps_list_tail->next; */
  /*     } */

  /*     /\* sum_admissible_saps_lp_solutions can be zero. If so, choose the action  */
  /*     /\* arbitrarily. Choose the first action *\/ */

  /*     int action = 0; */
  /*     if(definitely_greater_than(sum_admissible_saps_lp_solutions, 0)) { */
	
  /* 	/\* Normalize sap_pmf *\/ */
  /* 	for(index = 0; index < s->num_admissible_saps; ++index) { */
  /* 	  /\* printf( "sap: %lf \n", sap_pmf[index]); *\/ */
  /* 	  sap_pmf[index] = (sap_pmf[index]/sum_admissible_saps_lp_solutions); */
  /* 	  /\* printf( "sap normalized: %lf \n", sap_pmf[index]); *\/ */
  /* 	  assert(definitely_less_than_or_equal(sap_pmf[index], 1)); */
  /* 	} */
  /* 	/\* printf( "ok\n"); *\/ */
  /* 	/\* Now sample a state-action pair (sap) pointer from sap_pmf *\/ */
  /* 	/\* action is the index of a job  *\/ */
  /* 	ransampl_ws *sap_index_ws = ransampl_alloc( s->num_admissible_saps ); */
  /* 	ransampl_set( sap_index_ws, sap_pmf ); */
  /* 	int sap_index_decision = ransampl_draw( sap_index_ws, */
  /* 						gsl_rng_uniform(rng), */
  /* 						gsl_rng_uniform(rng) ); */
	
  /* 	/\* printf( "number of admissible actions: %d\n", s->num_admissible_saps); *\/ */
	
  /* 	action = admissible_sap_pointers[sap_index_decision]->action; */

  /* 	ransampl_free( sap_index_ws ); */
  /* 	free(sap_pmf); */
  /* 	free(admissible_sap_pointers); */
	
  /*     } else { */
	
  /* 	action = admissible_sap_pointers[0]->action; */
	
  /*     } */
      
  /*     /\* printf( "Job (%d) chosen at time (%lu) with probability %lf\n\n", *\/ */
  /*     /\* 	     action, sim_time, sap_pmf[sap_index_decision]); *\/ */
      
  /*     /\* printf( "after acting\n"); *\/ */
      
  /*     assert(action != NO_ACTION); */
  /*     assert(sim_state->finish_signal_vector[action] != FINISHED); */
      
  /*     /\* Update the state *\/ */
  /*     sim_state->allocation_vector[action]++; */
      
  /*     if(sim_state->allocation_vector[action] == */
  /* 	 arr_exec_time_samples[action]) { */
  /* 	sim_state->finish_signal_vector[action] = FINISHED; */
  /*     } */
      
  /*     ++sim_time; */
      
      
  /* 	int e = */
  /* 	  get_next_error_flag(sim_state->error_flag, */
  /* 			      sim_time, */
  /* 			      sim_state->allocation_vector, */
  /* 			      sim_state->finish_signal_vector); */
	
  /* 	sim_state->error_flag = e; */

  /*     /\* Check if a job just missed its deadline *\/ */
  /*     for(int job_index = 0; job_index < NUM_JOBS; ++job_index) { */
  /* 	job *j = &arr_jobs[job_index]; */
  /* 	bool job_not_finished = */
  /* 	  (sim_state->finish_signal_vector[job_index] == NOT_FINISHED); */
	
  /* 	if((sim_time == j->deadline) && job_not_finished) { */
  /* 	  assert(arr_job_missed_deadlines[job_index] == false); */
  /* 	  sim_stats->arr_num_deadline_misses_per_job[job_index]++; */
  /* 	  arr_job_missed_deadlines[job_index] = true; */
  /* 	  if(j->criticality == HI) { */
  /* 	    assert(sim_state->error_flag == ERROR); */
  /* 	  } */
  /* 	} */
  /*     } */
      
  /*     if(sim_state->error_flag == ERROR) { */
  /* 	sim_stats->num_errors++; */
  /*     } */
      
      
      
  /*   } // end simulation per sample while loop */
    
  /*   free(sim_state->allocation_vector); */
  /*   free(sim_state->finish_signal_vector); */
  /*   free(arr_job_missed_deadlines); */
    
  /* } //end looping over samples for loop */
  
  
  /* /\* Display simulation statistics *\/ */
  /* printf( "+=======================Simulation Stats=======================\n"); */
  /* printf( "+ Number of samples: %d\n", num_sim_samples); */
  /* printf( "+ Number of deadline misses per job:\n     ", num_sim_samples); */
  /* for(int job_index = 0; job_index < NUM_JOBS; ++job_index) { */
    
  /*   printf( "J%d: (%d/%d)   ", */
  /* 	   job_index, */
  /* 	   sim_stats->arr_num_deadline_misses_per_job[job_index], */
  /* 	   num_sim_samples); */
    
  /* } */
  
  /* printf( "\n+ Total number of error executions: (%d/%d)\n", */
  /* 	 sim_stats->num_errors, */
  /* 	 num_sim_samples); */
  /* printf( "+==============================================================\n"); */
  /* fflush(stdout); */

  /*----------------------End Simulation---------------------*/

    
  /* char *sched_str = is_ocbp_schedulable ? "YES" : "NO"; */
  /* printf( "OCBP Schdulability: %s\n", sched_str); */
  /* fflush(stdout); */
  
  /* char str[20]; */
  /* printf( "Done. Press any key to continue ...\n"); */
  /* scanf("%s", str); */
  return 0;
  
}
