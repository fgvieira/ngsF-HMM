#pragma once

#include "read_data.hpp"
#include "threadpool.h"
#include "bfgs.h"

extern bool SIG_COND;
const uint64_t N_STATES = 2;


// Struct to store all input arguments
typedef struct {
  char *in_geno;
  char *in_pos;
  bool in_bin;
  bool in_lkl;
  bool in_loglkl;
  uint64_t n_ind;
  uint64_t n_sites;
  double ***geno_lkl; // n_ind * n_sites+1 * N_GENO
  double *pos_dist;   // n_sites+1
  bool call_geno;
  char *in_freq;
  bool freq_fixed;
  char *in_indF;
  bool indF_fixed;
  char *in_path;
  bool path_fixed;
  char *out_prefix;
  unsigned int log;
  bool log_bin;
  uint min_iters;
  uint max_iters;
  double min_epsilon;
  uint n_threads;
  bool version;
  uint verbose;
  uint seed;

  double *freq;          // n_sites+1
  double ***prior;       // n_sites+1 * N_STATES * N_GENO
  char **path;           // n_ind * n_sites+1
  double ***marg_prob;   // n_ind * n_sites+1
  double *indF;          // n_ind
  double *aa;            // n_ind
  double *ind_lkl;       // n_ind
  double tot_lkl;

  threadpool_t *thread_pool;
} params;


// parse_args.cpp
void init_pars(params* );
void parse_cmd_args(params*, int, char**);
int init_output(params*);

// EM.cpp
int EM(params*);
void iter_EM(params*);
void post_prob(double*, double*, double*, uint64_t);
void print_iter(char*, params*);
void dump_data(gzFile, params*, bool = false);
double calc_trans(char, char, double, double, double, bool = true);
void calc_prior(double *, double, uint64_t);
void update_priors(double ***, double *, uint64_t);
