#include <getopt.h>
#include "ngsF-HMM.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_bin = false;
  pars->in_lkl = false;
  pars->in_loglkl = false;
  pars->geno_lkl = NULL;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->call_geno = 0;
  pars->in_freq = NULL;
  pars->freq_fixed = false;
  pars->in_trans = NULL;
  pars->trans_fixed = false;
  pars->in_path = NULL;
  pars->path_fixed = false;
  pars->out_prefix = NULL;
  pars->log = false;
  pars->log_bin = false;
  pars->min_iters = 10;
  pars->max_iters = 1500;
  pars->min_epsilon = 1e-5;
  pars->n_threads = 1;
  pars->n_chunks = 0;
  pars->max_chunk_size = 10000;
  pars->version = false;
  pars->verbose = 1;
  pars->seed = time(NULL);
}


// Parses command line args and stores them into struct params
int parse_cmd_args(int argc, char** argv, params* pars) {
  
  static struct option long_options[] =
    {
      {"geno", required_argument, NULL, 'g'},
      {"lkl", no_argument, NULL, 'l'},
      {"loglkl", no_argument, NULL, 'L'},
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"call_geno", required_argument, NULL, 'G'},
      {"freq", required_argument, NULL, 'f'},      
      {"freq_fixed", no_argument, NULL, 'F'},
      {"trans", required_argument, NULL, 't'},
      {"trans_fixed", no_argument, NULL, 'T'},
      {"path", required_argument, NULL, 'p'},
      {"path_fixed", no_argument, NULL, 'P'},
      {"out", required_argument, NULL, 'o'},
      {"log", no_argument, NULL, 'X'},
      {"log_bin", no_argument, NULL, 'b'},
      {"min_iters", required_argument, NULL, 'm'},
      {"max_iters", required_argument, NULL, 'M'},
      {"min_epsilon", required_argument, NULL, 'e'},
      {"n_threads", required_argument, NULL, 'x'},
      {"chunk_size", required_argument, NULL, 'c'},
      {"version", no_argument, NULL, 'v'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'S'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:lLn:s:G:f:Ft:Tp:Po:Xbm:M:e:x:c:vV:S:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_geno = optarg;
      break;
    case 'l':
      pars->in_lkl = true;
      break;
    case 'L':
      pars->in_lkl = true;
      pars->in_loglkl = true;
      break;
    case 'n':
      pars->n_ind = atoi(optarg);
      break;
    case 's':
      pars->n_sites = atoi(optarg);
      break;
    case 'G':
      pars->call_geno = atoi(optarg);
      break;
    case 'f':
      pars->in_freq = optarg;
      break;
    case 'F':
      pars->freq_fixed = true;
      break;
    case 't':
      pars->in_trans = optarg;
      break;
    case 'T':
      pars->trans_fixed = true;
      break;
    case 'p':
      pars->in_path = optarg;
      break;
    case 'P':
      pars->path_fixed = true;
      break;
    case 'o':
      pars->out_prefix = optarg;
      break;
    case 'X':
      pars->log = true;
      break;
    case 'b':
      pars->log_bin = true;
      break;
    case 'm':
      pars->min_iters = atoi(optarg);
      break;
    case 'M':
      pars->max_iters = atoi(optarg);
      break;
    case 'e':
      pars->min_epsilon = atof(optarg);
      break;
    case 'x':
      pars->n_threads = atoi(optarg);
      break;
    case 'c':
      pars->max_chunk_size = atoi(optarg);
      break;
    case 'v':
      pars->version = true;
      break;
    case 'V':
      pars->verbose = atoi(optarg);
      break;
    case 'S':
      pars->seed = atoi(optarg);
      break;
    default:
      exit(-1);
    }

  // Default value for initial values
  if(pars->in_freq == NULL) {
    pars->in_freq = new char[2];
    pars->in_freq[0] = 'r';
    pars->in_freq[1] = '\0';
  }
  if(pars->in_trans == NULL) {
    pars->in_trans = new char[2];
    pars->in_trans[0] = 'r';
    pars->in_trans[1] = '\0';
  }
  if(pars->in_path == NULL) {
    pars->in_path = new char[2];
    pars->in_path[0] = 'r';
    pars->in_path[1] = '\0';
  }

  return 0;
}
