#pragma once

#include <semaphore.h>
#include <pthread.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>

#include "shared.hpp"

using namespace std;

extern bool SIG_COND;

const uint64_t N_STATES = 2;
const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;


// Struct to store all input arguments //GZIP
typedef struct {
  char* in_geno;
  bool in_bin;
  bool in_lkl;
  bool in_loglkl;
  double*** geno_lkl;
  uint64_t n_ind;
  uint64_t n_sites;
  int call_geno;
  char* in_freq;
  bool freq_fixed;
  char* in_trans;
  bool trans_fixed;
  char* in_path;
  bool path_fixed;
  char* out_prefix;
  bool log;
  bool log_bin;
  uint min_iters;
  uint max_iters;
  double min_epsilon;
  uint n_threads;
  uint n_chunks;
  uint max_chunk_size;
  bool version;
  uint verbose;
  uint seed;
  
  sem_t launch_thread_semaph;
  sem_t running_thread_semaph;
  pthread_mutex_t F_lock;
} params;


// Output data
typedef struct {
  double*** a;
  double* freq;
  double*** e;
  uint** path;
  double*** marg_prob;
  double* indF;
  double* lkl;
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
int parse_cmd_args(int, char**, params*);

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
int update_e(out_data*, uint64_t);
