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
  pars->log = 0;
  pars->log_bin = false;
  pars->min_iters = 10;
  pars->max_iters = 1500;
  pars->min_epsilon = 1e-5;
  pars->n_threads = 1;
  pars->version = false;
  pars->verbose = 1;
  pars->seed = time(NULL);
  pars->thread_pool = NULL;
}


// Parses command line args and stores them into struct params
void parse_cmd_args(params* pars, int argc, char** argv){
  
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
      {"log", required_argument, NULL, 'X'},
      {"log_bin", required_argument, NULL, 'b'},
      {"min_iters", required_argument, NULL, 'm'},
      {"max_iters", required_argument, NULL, 'M'},
      {"min_epsilon", required_argument, NULL, 'e'},
      {"n_threads", required_argument, NULL, 'x'},
      {"version", no_argument, NULL, 'v'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'S'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:lLn:s:G:f:Ft:Tp:Po:X:b:m:M:e:x:vV:S:", long_options, NULL)) != -1 )
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
      pars->log = atoi(optarg);
      break;
    case 'b':
      pars->log = atoi(optarg);
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


  // Check Arguments
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "genotype input file (-geno) missing!");
  if(pars->out_prefix == NULL)
    error(__FUNCTION__, "output prefix (-out_prefix) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "number of individuals (-n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "number of sites (-n_sites) missing!");
  if(pars->log < 0)
    error(__FUNCTION__, "invalid LOG (-log) option!");


  // Print arguments to STDOUT
  if(pars->verbose >= 1){
    printf("==> Input Arguments:\n");
    printf("\tgeno file: %s\n\tgeno lkl: %s\n\tgeno loglkl: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tcall_geno: %d\n\tfreq: %s\n\tfreq_fixed: %s\n\ttrans: %s\n\ttrans_fixed: %s\n\tpath: %s\n\tpath_fixed: %s\n\tout prefix: %s\n\tlog: %u\n\tlog_bin: %s\n\tmin_iters: %d\n\tmax_iters: %d\n\tmin_epsilon: %.10f\n\tn_threads: %d\n\tversion: %s\n\tverbose: %d\n\tseed: %d\n\n",
           pars->in_geno, pars->in_lkl ? "true":"false", pars->in_loglkl ? "true":"false", pars->n_ind, pars->n_sites, pars->call_geno, pars->in_freq, pars->freq_fixed ? "true":"false", pars->in_trans, pars->trans_fixed ? "true":"false", pars->in_path, pars->path_fixed ? "true":"false", pars->out_prefix, pars->log, pars->log_bin ? "true":"false", pars->min_iters, pars->max_iters, pars->min_epsilon, pars->n_threads, pars->version ? "true":"false", pars->verbose, pars->seed);
  }
  if(pars->verbose >= 4)
    printf("==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");
}





int init_output(params* pars, out_data* data) {
  char* buf = new char[BUFF_LEN];
  double* t;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, pars->seed);


  ///////////////////////////////////
  // Set TRANSITION initial values //
  ///////////////////////////////////
  double trans_rng_max = 1;
  double trans_rng_min = 0;
  gzFile in_trans_fh;

  data->a = init_ptr(pars->n_ind, N_STATES, N_STATES, -INFINITY);

  if( strcmp("r", pars->in_trans) == 0 )
    for(uint64_t i = 0; i < pars->n_ind; i++){
      data->a[i][0][1] = trans_rng_min + gsl_rng_uniform(r) * (trans_rng_max - trans_rng_min);
      data->a[i][0][0] = 1 - data->a[i][0][1];
      data->a[i][1][0] = trans_rng_min + gsl_rng_uniform(r) * (trans_rng_max - trans_rng_min);
      data->a[i][1][1] = 1 - data->a[i][1][0];
    }
  
  else if( (in_trans_fh = gzopen(pars->in_trans, "r")) != NULL ){
    uint64_t i = 0;
    while( gzgets(in_trans_fh, buf, BUFF_LEN) != NULL ){
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;

      if( i > pars->n_ind || split(buf, (const char*) " ,-\t", &t) != 2)
	error(__FUNCTION__, "wrong TRANS file format!");

      data->a[i][0][1] = min(max(t[0], trans_rng_min), trans_rng_max);
      data->a[i][0][0] = 1 - data->a[i][0][1];
      data->a[i][1][0] = min(max(t[1], trans_rng_min), trans_rng_max);
      data->a[i][1][1] = 1 - data->a[i][1][0];
      i++;

      delete [] t;
    }
    gzclose(in_trans_fh);
  }

  else{
    if( split(pars->in_trans, (const char*) ",-", &t) != 2 )
      error(__FUNCTION__, "wrong TRANS parameters format!");

    for(uint64_t i = 0; i < pars->n_ind; i++){
      data->a[i][0][1] = min(max(t[0], trans_rng_min), trans_rng_max);
      data->a[i][0][0] = 1 - data->a[i][0][1];
      data->a[i][1][0] = min(max(t[1], trans_rng_min), trans_rng_max);
      data->a[i][1][1] = 1 - data->a[i][1][0];
    }
    delete [] t;
  }
  // Convert transition probs to log-scale
  for (uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t k = 0; k < N_STATES; k++)
      for(uint64_t l = 0; l < N_STATES; l++)
	data->a[i][k][l] = log(data->a[i][k][l]);



  /////////////////////////////
  // Set FREQ initial values //
  /////////////////////////////
  double freq_rng_max = 0.49;
  double freq_rng_min = 0.01;
  gzFile in_freq_fh;

  data->freq = init_ptr(pars->n_sites+1, freq_rng_min);
  // Initialize site 0 to invalid value
  data->freq[0] = -1;

  if( strcmp("r", pars->in_freq) == 0 )
    for(uint64_t s = 1; s <= pars->n_sites; s++)
      data->freq[s] = freq_rng_min + gsl_rng_uniform(r) * (freq_rng_max - freq_rng_min);

  else if( (in_freq_fh = gzopen(pars->in_freq, "r")) != NULL ){
    uint64_t s = 1;
    while( gzgets(in_freq_fh, buf, BUFF_LEN) != NULL ){
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;

      if( s > pars->n_sites || split(buf, (const char*) " ,-\t", &t) != 1)
        error(__FUNCTION__, "wrong FREQ file format!");

      data->freq[s] = min(max(t[0], freq_rng_min), freq_rng_max);
      s++;
      delete [] t;
    } 
    gzclose(in_freq_fh);
  }

  else
    for(uint64_t s = 1; s <= pars->n_sites; s++)
      data->freq[s] = min(max(atof(pars->in_freq), freq_rng_min), freq_rng_max);



  ///////////////////////////////////////////////
  // Set EMISSION probabilities initial values //
  ///////////////////////////////////////////////
  data->e = init_ptr(pars->n_sites+1, N_STATES, N_GENO, -INFINITY);
  // Update emission probs based on allele freqs
  update_e(data, pars->n_sites);
  


  /////////////////////////////
  // Set PATH initial values //
  /////////////////////////////
  gzFile in_path_fh;

  data->path = init_ptr(pars->n_ind, pars->n_sites+1, (const char*) '\0');

  if( strcmp("r", pars->in_path) == 0 )
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++)
	data->path[i][s] = gsl_rng_uniform(r) > 0.5 ? 1 : 0;

  else if( (in_path_fh = gzopen(pars->in_path, "r")) != NULL ){
    uint64_t i = 0;
    while(gzgets(in_path_fh, buf, BUFF_LEN) != NULL){
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;

      int* t = NULL;
      if(i >= pars->n_ind || split(buf, (const char*) "", &t) != pars->n_sites)
        error(__FUNCTION__, "wrong PATH file format!");

      for(uint64_t s = 0; s < pars->n_sites; s++)
	data->path[i][s] = t[s] > 0.5 ? 1 : 0;
      i++;
      delete [] t;
    }
    gzclose(in_path_fh);
  }

  else
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++)
	data->path[i][s] = atoi(pars->in_path) > 0.5 ? 1 : 0;



  ///////////////////////////////////////
  // Initialize Marginal Probabilities //
  ///////////////////////////////////////
  data->marg_prob = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, 0.0);
  // Initialize site 0 to invalid value
  for (uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t k = 0; k < N_STATES; k++)
      data->marg_prob[i][0][k] = -1;



  ///////////////////////////////////////
  // Initialize indF and lkl variables //
  ///////////////////////////////////////
  data->indF = init_ptr(pars->n_ind, 0.0);
  data->lkl = init_ptr(pars->n_ind, -INFINITY);



  delete [] buf;
  gsl_rng_free(r);
  return(0);
}
