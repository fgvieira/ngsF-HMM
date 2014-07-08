#pragma once

#include "shared.hpp"

using namespace std;

extern bool SIG_COND;

const uint64_t N_STATES = 2;
const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;

#define path_pos_get(seq, pos) ( seq[pos/8]  & (1 << pos%8) )
#define path_pos_set0(seq, pos) ( seq[pos/8] &= ~(1 << pos%8) )
#define path_pos_set1(seq, pos) ( seq[pos/8] |= (1 << pos%8) )

// Struct to store all input arguments //GZIP
typedef struct {
  char *in_geno;
  bool in_bin;
  bool in_lkl;
  bool in_loglkl;
  double ***geno_lkl;
  uint64_t n_ind;
  uint64_t n_sites;
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
  double ***a;
  double *freq;
  double ***prior;
  char **path;
  double ***marg_prob;
  double *indF;
  double *lkl;
} out_data;


// Prepare structure
typedef struct {
  params* pars;
  uint64_t chunk_size;
  uint64_t chunk_abs_start_pos;
  double** chunk_data;
  int iter;
  unsigned int* n_run_threads;
  out_data* data;
} pth_params;


// parse_args.cpp
void init_pars(params* );
void parse_cmd_args(params*, int, char**);

// read_data.cpp
int init_output(params*, out_data*);
int read_geno(params*);
uint64_t read_chunk(double**, params*, uint64_t);

// EM.cpp
int EM(params*, out_data*);
void *run_chunk(void*);
void iter_EM(params*, out_data*);
void post_prob(double*, double*, double*, uint64_t);
void print_iter(char*, params*, out_data*);
int update_priors(out_data*, uint64_t);
