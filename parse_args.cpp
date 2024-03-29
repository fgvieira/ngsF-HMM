#include <getopt.h>
#include "ngsF-HMM.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_pos = NULL;
  pars->in_bin = false;
  pars->in_lkl = false;
  pars->in_loglkl = false;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->geno_lkl = NULL;
  pars->pos_dist = NULL;
  pars->call_geno = false;
  pars->in_freq = NULL;
  pars->freq_est = 1;
  pars->e_prob_calc = 1;
  pars->in_indF = NULL;
  pars->indF_fixed = false;
  pars->alpha_fixed = false;
  pars->out_prefix = NULL;
  pars->log = 0;
  pars->log_bin = false;
  pars->min_iters = 10;
  pars->max_iters = 100;
  pars->min_epsilon = 1e-5;
  pars->n_threads = 1;
  pars->verbose = 1;
  pars->seed = rand() % 1000;
  pars->prev_tot_lkl = 0;
  pars->tot_lkl = 0;
  pars->thread_pool = NULL;
}



/////////////////////////////////////////////////////////////////
// Parses command line args and stores them into struct params //
/////////////////////////////////////////////////////////////////
void parse_cmd_args(params* pars, int argc, char** argv){
  
  static struct option long_options[] =
    {
      {"geno", required_argument, NULL, 'g'},
      {"pos", required_argument, NULL, 'Z'},
      {"lkl", no_argument, NULL, 'l'},
      {"loglkl", no_argument, NULL, 'L'},
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"call_geno", no_argument, NULL, 'G'},
      {"freq", required_argument, NULL, 'f'},
      {"freq_est", required_argument, NULL, 'F'},
      {"e_prob", required_argument, NULL, 'e'},
      {"indF", required_argument, NULL, 'i'},
      {"indF_fixed", no_argument, NULL, 'I'},
      {"alpha_fixed", no_argument, NULL, 'A'},
      {"out", required_argument, NULL, 'o'},
      {"log", required_argument, NULL, 'X'},
      {"log_bin", required_argument, NULL, 'b'},
      {"min_iters", required_argument, NULL, 'm'},
      {"max_iters", required_argument, NULL, 'M'},
      {"min_epsilon", required_argument, NULL, 'E'},
      {"n_threads", required_argument, NULL, 'x'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'S'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:Z:lLn:s:Gf:F:e:i:IAo:X:b:m:M:E:x:V:S:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_geno = optarg;
      break;
    case 'Z':
      pars->in_pos = optarg;
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
      pars->call_geno = true;
      break;
    case 'f':
      pars->in_freq = optarg;
      break;
    case 'F':
      pars->freq_est = atoi(optarg);
      break;
    case 'e':
      pars->e_prob_calc = atoi(optarg);
      break;
    case 'i':
      pars->in_indF = optarg;
      break;
    case 'I':
      pars->indF_fixed = true;
      break;
    case 'A':
      pars->alpha_fixed = true;
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
    case 'E':
      pars->min_epsilon = atof(optarg);
      break;
    case 'x':
      pars->n_threads = atoi(optarg);
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


  
  //////////////////////////////////////
  // Default value for initial values //
  //////////////////////////////////////
  if(pars->in_freq == NULL) {
    pars->in_freq = init_ptr(20, (const char*) '\0');
    strcat(pars->in_freq, "r");
  }
  if(pars->in_indF == NULL) {
    pars->in_indF = init_ptr(20, (const char*) '\0');
    strcat(pars->in_indF, "0.01-0.001");
  }



  ///////////////////////////////
  // Print arguments to STDOUT //
  ///////////////////////////////
  if(pars->verbose >= 1){
    printf("==> Input Arguments:\n");
    printf("\tgeno: %s\n\tpos: %s\n\tlkl: %s\n\tloglkl: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tcall_geno: %s\n\tfreq: %s\n\tfreq_est: %d\n\te_prob: %d\n\tindF: %s\n\tindF_fixed: %s\n\talpha_fixed: %s\n\tout: %s\n\tlog: %u\n\tlog_bin: %s\n\tmin_iters: %d\n\tmax_iters: %d\n\tmin_epsilon: %.10f\n\tn_threads: %d\n\tverbose: %d\n\tseed: %d\n\tversion: %s (%s @ %s)\n\n",
           pars->in_geno,
	   pars->in_pos,
	   pars->in_lkl ? "true":"false",
	   pars->in_loglkl ? "true":"false",
	   pars->n_ind,
	   pars->n_sites,
	   pars->call_geno ? "true":"false",
	   pars->in_freq,
	   pars->freq_est,
	   pars->e_prob_calc,
	   pars->in_indF,
	   pars->indF_fixed ? "true":"false",
	   pars->alpha_fixed ? "true":"false",
	   pars->out_prefix,
	   pars->log,
	   pars->log_bin ? "true":"false",
	   pars->min_iters,
	   pars->max_iters,
	   pars->min_epsilon,
	   pars->n_threads,
	   pars->verbose,
	   pars->seed,
	   version, __DATE__, __TIME__);
  }
  if(pars->verbose >= 4)
    printf("==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");



  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "genotype input file (--geno) missing!");
  if(pars->in_pos == NULL)
    error(__FUNCTION__, "positions input file (--pos) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "number of individuals (--n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "number of sites (--n_sites) missing!");
  if(pars->call_geno && !pars->in_lkl)
    error(__FUNCTION__, "can only call genotypes from likelihoods!");
  if(pars->freq_est < 0 || pars->freq_est > 2)
    error(__FUNCTION__, "invalid MAF estimation method!");
  if(pars->e_prob_calc < 0 || pars->e_prob_calc > 2)
    error(__FUNCTION__, "invalid emission probability calculation method!");
  if(pars->e_prob_calc > 1)
    warn(__FUNCTION__, "calculation of emission probabilities accounting for LD is still under development!");
  //if(pars->freq_est == 2 && pars->e_prob_calc == 1)
  //error(__FUNCTION__, "incompatible MAF and EMISSION algorithms!");
  if(pars->out_prefix == NULL)
    error(__FUNCTION__, "output prefix (--out) missing!");
  if(pars->log < 0)
    error(__FUNCTION__, "invalid LOG (--log) option!");
  if(pars->min_iters < 1 || pars->max_iters < 1 || pars->min_iters >= pars->max_iters)
    error(__FUNCTION__, "invalid number of iterations!");
  if(pars->n_threads < 1)
    error(__FUNCTION__, "invalid number of threads!");
}



int init_output(params* pars) {
  char* buf = new char[BUFF_LEN];
  double* t;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, pars->seed);


  ///////////////////////////////////
  // Set INBREEDING initial values //
  ///////////////////////////////////
  double indF_rng_min = 0.000001;
  double indF_rng_max = 1 - indF_rng_min;
  double alpha_rng_min = 0.000001;
  double alpha_rng_max = 1 - alpha_rng_min;
  gzFile in_indF_fh;

  pars->indF = init_ptr(pars->n_ind, (double) 0);
  pars->alpha = init_ptr(pars->n_ind, (double) 0);

  if( strcmp("r", pars->in_indF) == 0 ){
    if(pars->verbose >= 1)
      printf("==> Using random initial inbreeding values.\n");
    for(uint64_t i = 0; i < pars->n_ind; i++){
      pars->indF[i] = indF_rng_min + gsl_rng_uniform(r) * (indF_rng_max - indF_rng_min);
      pars->alpha[i] = alpha_rng_min + gsl_rng_uniform(r) * (alpha_rng_max - alpha_rng_min);
    }
  }else if( (in_indF_fh = gzopen(pars->in_indF, "r")) != NULL ){
    if(pars->verbose >= 1)
      printf("==> Reading initial inbreeding values from file \"%s\".\n", pars->in_indF);

    uint64_t i = 0;
    while( gzgets(in_indF_fh, buf, BUFF_LEN) != NULL ){
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;

      if( i > pars->n_ind || split(buf, (const char*) " ,-\t", &t) != 2)
	error(__FUNCTION__, "wrong INDF file format!");

      pars->indF[i] = min(max(t[0], indF_rng_min), indF_rng_max);
      pars->alpha[i] = min(max(t[1], alpha_rng_min), alpha_rng_max);
      i++;

      delete [] t;
    }
    gzclose(in_indF_fh);
  }else{
    if(pars->verbose >= 1)
      printf("==> Setting initial inbreeding values to: %s\n", pars->in_indF);

    if( split(pars->in_indF, (const char*) ",-", &t) != 2 )
      error(__FUNCTION__, "wrong INDF parameters format!");

    for(uint64_t i = 0; i < pars->n_ind; i++){
      pars->indF[i] = min(max(t[0], indF_rng_min), indF_rng_max);
      pars->alpha[i] = min(max(t[1], alpha_rng_min), alpha_rng_max);
    }
    delete [] t;
  }



  /////////////////////////////
  // Set FREQ initial values //
  /////////////////////////////
  double freq_rng_min = 0.01;
  double freq_rng_max = 0.5 - freq_rng_min;
  double hap_freq[4]; // P_AB, P_Ab, P_aB, P_ab
  gzFile in_freq_fh;

  pars->freq = init_ptr(pars->n_sites+1, freq_rng_min);
  // Initialize site 0 to invalid value
  pars->freq[0] = -1;

  if( strcmp("r", pars->in_freq) == 0 ) {
    if(pars->verbose >= 1)
      printf("==> Using random initial frequency values.\n");

    for(uint64_t s = 1; s <= pars->n_sites; s++)
      pars->freq[s] = freq_rng_min + gsl_rng_uniform(r) * (freq_rng_max - freq_rng_min);

  } else if( strcmp("e", pars->in_freq) == 0 ){
    if(pars->verbose >= 1)
      printf("==> Estimating initial frequency values assuming HWE.\n");

    for(uint64_t s = 1; s <= pars->n_sites; s++)
      if(pars->freq_est == 1 || s == 1){
	pars->freq[s] = est_maf(pars->n_ind, pars->geno_lkl_s[s], (double) 0);
      }else if(pars->freq_est == 2){
	double loglkl;
	uint64_t n_iter, n_ind_data;
	n_iter = haplo_freq(hap_freq, &loglkl, &n_ind_data, pars->geno_lkl_s[s-1], pars->geno_lkl_s[s], pars->freq[s-1], pars->freq[s], pars->n_ind, false);
	pars->freq[s] = hap_freq[1] + hap_freq[3];
      }

  } else if( (in_freq_fh = gzopen(pars->in_freq, "r")) != NULL ){
    if(pars->verbose >= 1)
      printf("==> Reading initial frequency values from file \"%s\".\n", pars->in_freq);

    uint64_t s = 1;
    uint64_t n_fields = 0;
    while( gzgets(in_freq_fh, buf, BUFF_LEN) != NULL ){
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;

      // Parse input line into array
      n_fields = split(buf, (const char*) " ,-\t", &t);

      // Check if header and skip
      if(!n_fields){
	printf("> Header found! Skipping line...\n");
	continue;
      }

      if( s > pars->n_sites || n_fields != 1)
        error(__FUNCTION__, "wrong FREQ file format!");

      pars->freq[s] = min(max(t[0], freq_rng_min), freq_rng_max);
      s++;
      delete [] t;
    } 
    gzclose(in_freq_fh);

  } else {
    if(pars->verbose >= 1)
      printf("==> Setting initial frequency values to: %s\n", pars->in_freq);

    for(uint64_t s = 1; s <= pars->n_sites; s++)
      pars->freq[s] = min(max(atof(pars->in_freq), freq_rng_min), freq_rng_max);
  }



  ///////////////////////////////////////
  // Initialize emission probabilities //
  ///////////////////////////////////////
  if(pars->verbose >= 1)
    printf("==> Calculating initial emission probabilities\n");
  pars->e_prob = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, (double) 0);

  for(uint64_t s = 1; s <= pars->n_sites; s++){
    if(pars->e_prob_calc == 2) {
      double loglkl;
      uint64_t n_iter, n_ind_data;
      n_iter = haplo_freq(hap_freq, &loglkl, &n_ind_data, pars->geno_lkl_s[s-1], pars->geno_lkl_s[s], pars->freq[s-1], pars->freq[s], pars->n_ind, false);
    }

    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t k = 0; k < N_STATES; k++)
	if(pars->e_prob_calc == 1 || s == 1)
	  pars->e_prob[i][s][k] = calc_emission(pars->geno_lkl[i][s], pars->freq[s], k);
	else if(pars->e_prob_calc == 2)
	  pars->e_prob[i][s][k] = calc_emissionLD(hap_freq, pars->geno_lkl[i][s-1], pars->geno_lkl[i][s], pars->freq[s-1], pars->freq[s], k);
  }



  //////////////////////////
  // Allocate PATH memory //
  //////////////////////////
  pars->path = init_ptr(pars->n_ind, pars->n_sites+1, (const char*) '\0');



  ///////////////////////////////////////
  // Initialize Marginal Probabilities //
  ///////////////////////////////////////
  pars->marg_prob = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, (double) 0);
  // Initialize site 0 to invalid value
  for (uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t k = 0; k < N_STATES; k++)
      pars->marg_prob[i][0][k] = -1;



  //////////////////////////
  // Initialize Lkl array //
  //////////////////////////
  pars->ind_lkl = init_ptr(pars->n_ind, (double) -INFINITY);



  delete [] buf;
  gsl_rng_free(r);
  return(0);
}
