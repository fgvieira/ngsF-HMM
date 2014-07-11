#pragma once

#include "shared.hpp"

using namespace std;

extern bool SIG_COND;

const uint64_t INF = 0xFFFFFFFFFFFFFFFF;
const uint64_t N_STATES = 2;
const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;

#define path_pos_get(seq, pos) ( seq[pos/8]  & (1 << pos%8) )
#define path_pos_set0(seq, pos) ( seq[pos/8] &= ~(1 << pos%8) )
#define path_pos_set1(seq, pos) ( seq[pos/8] |= (1 << pos%8) )

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
  uint64_t *pos_dist; // n_sites+1
  bool call_geno;
  char *in_freq;
  bool freq_fixed;
  char *in_trans;
  bool trans_fixed;
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
  threadpool_t *thread_pool;
} params;



// Output data
typedef struct {
  double ***a;           // n_ind * N_STATES * N_STATES
  double *freq;          // n_sites+1
  double ***prior;       // n_sites+1 * N_STATES * N_GENO
  char **path;           // n_ind * n_sites+1
  double ***marg_prob;   // n_ind * n_sites+1
  double *indF;          // n_ind
  double *lkl;           // n_ind
} out_data;



// parse_args.cpp
void init_pars(params* );
void parse_cmd_args(params*, int, char**);

// read_data.cpp
int init_output(params*, out_data*);

// EM.cpp
int EM(params*, out_data*);
void iter_EM(params*, out_data*);
void post_prob(double*, double*, double*, uint64_t);
void print_iter(char*, params*, out_data*);
int update_priors(out_data*, uint64_t);
