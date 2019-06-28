
#include "ngsF-HMM.hpp"


// General structure for launching threads
struct pth_struct{
  int type;
  double **ptr;
  double *F;
  bool F_fixed;
  double *alpha;
  bool alpha_fixed;
  double **e_prob;
  char *path;
  double *pos_dist;
  uint64_t length;
};


// Function prototypes
void threadpool_add_task(threadpool_t *thread_pool, int type, double **ptr, double *F, bool F_fixed, double *alpha, bool alpha_fixed, double **e_prob, char *path, double *pos_dist, uint64_t length);
void thread_slave(void *ptr);
double lkl(const double*, const void*);



int EM (params *pars) {
  SIG_COND = true;
  catch_SIG();

  uint64_t iter = 0;
  double max_lkl_epsilon = -INFINITY;
  double *prev_ind_lkl = init_ptr(pars->n_ind, (double) -INFINITY);
  double *ind_lkl_epsilon = init_ptr(pars->n_ind, (double) -INFINITY);



  // Print out initial parameters
  if(pars->verbose >= 5){
    printf("==> Initial parameters:\n");
    // indF and alpha
    for(uint64_t i = 0; i < pars->n_ind; i++)
      printf("\t%.10f\t%f\n", pars->indF[i], pars->alpha[i]);

    // freq
    for(uint64_t s = 1; s <= pars->n_sites; s++)
      printf("\t%f", pars->freq[s]);
    printf("\n");
  }



  ////////////////////
  // Iteration loop //
  ////////////////////
  while((pars->prev_tot_lkl - pars->tot_lkl > pars->min_epsilon || max_lkl_epsilon > pars->min_epsilon || iter < pars->min_iters) && iter < pars->max_iters && SIG_COND) {

    // Dump previous iteration data
    if(pars->log && (iter == 1 || iter % pars->log == 0)){
      if(pars->verbose >= 1)    
	printf("==> Printing current iteration parameters\n");
      print_iter(pars->out_prefix, pars);
    }

    // Next Iteration...
    time_t iter_start = time(NULL);
    iter++;
    if(pars->verbose >= 1)
      printf("\nIteration %lu:\n", iter);

    // Run next EM iteration
    iter_EM(pars);

    // Check convergence criteria
    pars->prev_tot_lkl = pars->tot_lkl;
    pars->tot_lkl = 0;
    for (uint64_t i = 0; i < pars->n_ind; i++){
      // Get total Lkl
      pars->tot_lkl += pars->ind_lkl[i];
      // Get per-indiv lkl epsilon
      ind_lkl_epsilon[i] = (pars->ind_lkl[i] - prev_ind_lkl[i]) / fabs(prev_ind_lkl[i]);
    }
    uint64_t ind_max_lkl_epsilon  = array_max_pos(ind_lkl_epsilon, pars->n_ind);
    max_lkl_epsilon = ind_lkl_epsilon[ind_max_lkl_epsilon];
    // Save current LKLs
    cpy(prev_ind_lkl, pars->ind_lkl, pars->n_ind, sizeof(double));

    // Print iteration info..
    time_t iter_end = time(NULL);
    if(pars->verbose >= 1)
      printf("\tLogLkl: %.15f\t max lkl epsilon: %.15f\ttime: %.0f (s)\n", pars->tot_lkl, max_lkl_epsilon, difftime(iter_end, iter_start) );

    if(pars->verbose >= 3)
      for (uint64_t i = 0; i < pars->n_ind; i++)
	printf("\tInd %lu: %.15f\t lkl epsilon: %.15f%s\n", i+1, pars->ind_lkl[i], ind_lkl_epsilon[i], i == ind_max_lkl_epsilon ? " (max)" : "");

    fflush(stdout);
  }

  if(iter >= pars->max_iters)
    printf("WARN: Maximum number of iterations reached! Check if analysis converged... \n");



  /////////////
  // Viterbi //
  /////////////
  if(pars->verbose >= 1)
    printf("\n==> Decoding most probable path (Viterbi)\n");
  double ***Vi = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, (double) 0);

  for (uint64_t i = 0; i < pars->n_ind; i++)
    threadpool_add_task(pars->thread_pool, 3, Vi[i], &pars->indF[i], false, &pars->alpha[i], false, pars->e_prob[i], pars->path[i], pars->pos_dist, pars->n_sites);

  threadpool_wait(pars->thread_pool);
  free_ptr((void***) Vi, pars->n_ind, pars->n_sites+1);



  /////////////////////////
  // Print Final Results //
  /////////////////////////
  if(pars->verbose >= 1){
    printf("Final logLkl: %f\n", pars->tot_lkl);
    printf("Printing final results\n");
  }
  print_iter(pars->out_prefix, pars);



  // Free memory and return
  free_ptr((void*) ind_lkl_epsilon);
  free_ptr((void*) prev_ind_lkl);
  return 0;
}



void iter_EM(params *pars) {
  double ***marg_prob = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, -INF);
  cpy(marg_prob, pars->marg_prob, pars->n_ind, pars->n_sites+1, N_STATES, sizeof(double));
  double ***Fw = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, (double) 0);
  double ***Bw = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, (double) 0);



  // Forward recursion
  time_t fwd_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Forward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++)
    threadpool_add_task(pars->thread_pool, 1, Fw[i], &pars->indF[i], false, &pars->alpha[i], false, pars->e_prob[i], NULL, pars->pos_dist, pars->n_sites);
    
  // Backward recursion
  time_t bwd_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Backward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++)
    threadpool_add_task(pars->thread_pool, 2, Bw[i], &pars->indF[i], false, &pars->alpha[i], false, pars->e_prob[i], NULL, pars->pos_dist, pars->n_sites);

  threadpool_wait(pars->thread_pool);



  // Lkl check! - relaxed (to 0.001) due to precision issues on large datasets
  for (uint64_t i = 0; i < pars->n_ind; i++)
    if( abs(logsum(Fw[i][pars->n_sites],2) - logsum(Bw[i][0],2)) > 0.001 ){
      printf("Ind %lu: %.15f\t%.15f (%.15f)\n", i, logsum(Fw[i][pars->n_sites],2), logsum(Bw[i][0],2), abs(logsum(Fw[i][pars->n_sites],2) - logsum(Bw[i][0],2)) );
      error(__FUNCTION__, "Fw and Bw lkl do not match!");
    }



  // Marginal probabilities
  time_t mp_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Marginal probabilities\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    // Get per-indiv Lkl
    pars->ind_lkl[i] = logsum(Fw[i][pars->n_sites],2);
    // Get marg probs
    for (uint64_t s = 1; s <= pars->n_sites; s++)
      for(uint64_t k = 0; k < N_STATES; k++)
	pars->marg_prob[i][s][k] = check_interv(exp(Bw[i][s][k] + Fw[i][s][k] - pars->ind_lkl[i]), false);
  }


  
  // Estimate inbreeding and transition parameter
  time_t indF_t = time(NULL);
  if(pars->indF_fixed){
    if(pars->verbose >= 1)
      printf("==> Inbreeding and transition parameter not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update inbreeding and transition parameter\n");

    for(uint64_t i = 0; i < pars->n_ind; i++)
      threadpool_add_task(pars->thread_pool, 4, NULL, &pars->indF[i], pars->indF_fixed, &pars->alpha[i], pars->alpha_fixed, pars->e_prob[i], NULL, pars->pos_dist, pars->n_sites);

    threadpool_wait(pars->thread_pool);

    if(pars->verbose >= 4)
      for(uint64_t i = 0; i < pars->n_ind; i++)
	printf("\t%.10f\t%f\n", pars->indF[i], pars->alpha[i]);
  }



  // Estimate allele frequencies (EM)
  time_t freqs_t = time(NULL);
  if(pars->freq_est == 0){
    if(pars->verbose >= 1)
      printf("==> Alelle frequencies not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Estimating allele frequencies and calculating emission probabilities\n");

    double prior[N_GENO], hap_freq[4];
    double *indF = init_ptr(pars->n_ind, 0.0);
    double **prev_site = init_ptr(pars->n_ind, N_GENO, 0.0);
    double **curr_site = init_ptr(pars->n_ind, N_GENO, 0.0);

    for (uint64_t s = 1; s <= pars->n_sites; s++){
      for(uint64_t i = 0; i < pars->n_ind; i++){
	indF[i] = pars->marg_prob[i][s][1];

	calc_HWE(prior, pars->freq[s-1], pars->marg_prob[i][s-1][1]);
	post_prob(prev_site[i], pars->geno_lkl[i][s-1], prior, N_GENO);
	calc_HWE(prior, pars->freq[s], pars->marg_prob[i][s][1]);
	post_prob(curr_site[i], pars->geno_lkl[i][s], prior, N_GENO);
      }

      // Calculate haplotype frequency through an EM
      if(pars->freq_est == 2 || pars->e_prob_calc == 2)
	haplo_freq(hap_freq, prev_site, curr_site, pars->freq[s-1], pars->freq[s], pars->n_ind);

      // Estimated MAF
      if(pars->freq_est == 1 || s == 1){
	// Calculate MAF assuming independent sites through an EM
	pars->freq[s] = est_maf(pars->n_ind, pars->geno_lkl_s[s], indF);
      }else if(pars->freq_est == 2){
	// Calculate MAF through the haplotype frequency
	pars->freq[s] = hap_freq[1] + hap_freq[3];
      }else
	error(__FUNCTION__, "wrong MAF estimation method!");

      // Calculate emission probabilites
      if(pars->e_prob_calc == 1 || s == 1) {
	for(uint64_t i = 0; i < pars->n_ind; i++)
	  for(uint64_t k = 0; k < N_STATES; k++)
	    if(pars->e_prob_calc == 1 || s == 1)
	      // Calculate emission probability conditioned on MAF
	      pars->e_prob[i][s][k] = calc_emission(pars->geno_lkl[i][s], pars->freq[s], k);
	    else if(pars->e_prob_calc == 2)
	      // Calculate emission probability conditioned on previous site genotype
	      pars->e_prob[i][s][k] = calc_emissionLD(hap_freq, pars->geno_lkl[i][s-1], pars->geno_lkl[i][s], pars->freq[s-1], pars->freq[s], k);
	    else
	      error(__FUNCTION__, "wrong emission probability calculation method!");
      }

      if(pars->verbose >= 7){
	printf("Site %lu; freq: %f; emission: ", s, pars->freq[s]);
	for(uint64_t i = 0; i < pars->n_ind; i++)
	  printf("\t%f/%f", exp(pars->e_prob[i][s][0]), exp(pars->e_prob[i][s][1]));
	printf("\n");
      }
    }
  }



  time_t end_t = time(NULL);
  if(pars->verbose >= 3)
    printf("\nFw: %.1f\nBw: %.1f\nMP: %.1f\nindF: %.1f\nfreqs: %.1f\n", 
	   difftime(bwd_t,fwd_t), 
	   difftime(mp_t,bwd_t), 
	   difftime(indF_t,mp_t), 
	   difftime(freqs_t,indF_t),
	   difftime(end_t,freqs_t)
	   );

  free_ptr((void***) marg_prob, pars->n_ind, pars->n_sites+1);
  free_ptr((void***) Fw, pars->n_ind, pars->n_sites+1);
  free_ptr((void***) Bw, pars->n_ind, pars->n_sites+1);
}



void print_iter(char *out_prefix, params *pars){
  // Open filehandle to "indF" file
  char *tmp_out = strdcat(out_prefix, ".indF");
  gzFile out_fh = open_gzfile(tmp_out, "wT");
  delete [] tmp_out;

  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open INDF output file!");

  // Print total Lkl
  gzprintf(out_fh, "%.10f\n", pars->tot_lkl);

  // Print indF and transition pars
  for(uint16_t i = 0; i < pars->n_ind; i++){
    if(pars->indF[i] < EPSILON)
      gzprintf(out_fh, "%.5f\tNA\n", (double) 0);
    else if(pars->indF[i] > 1-EPSILON)
      gzprintf(out_fh, "%.5f\tNA\n", (double) 1);
    else
      gzprintf(out_fh, "%.5f\t%f\n", pars->indF[i], pars->alpha[i]);
  }

  // Print allele freqs
  for(uint64_t s = 1; s <= pars->n_sites; s++)
    gzprintf(out_fh, "%f\n", pars->freq[s]);

  // Close "indF" filehandle
  gzclose(out_fh);

  /////////////////////////////////////////////////
  // Open filehandle to "IBD" file
  tmp_out = strdcat(out_prefix, ".ibd");
  out_fh = open_gzfile(tmp_out, "wT", max(10000,pars->n_sites)+100);
  delete [] tmp_out;

  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open IBD output file!");

  // Print IBD info: most probable path (Viterbi) and IBD marg probs
  char *buf = join(pars->ind_lkl, pars->n_ind, "\t");
  // Print Lkl
  if(gzprintf(out_fh, "//\t%s\n", buf) <= 0)
    error(__FUNCTION__, "cannot write LKL info to file!");
  delete [] buf;

  // Print most probable path (Viterbi)
  for (uint64_t i = 0; i < pars->n_ind; i++){
    for (uint64_t s = 1; s <= pars->n_sites; s++)
      if(gzprintf(out_fh, "%c", pars->path[i][s]+48) <= 0)
	error(__FUNCTION__, "cannot write PATH info to file!");
    gzprintf(out_fh, "\n");
  }

  // Print marginal probs
  for (uint64_t i = 0; i < pars->n_ind; i++){
    // To avoid leading \t
    gzprintf(out_fh, "%f", pars->marg_prob[i][1][1]);
    for (uint64_t s = 2; s <= pars->n_sites; s++)
      gzprintf(out_fh, "\t%f", pars->marg_prob[i][s][1]);
    gzprintf(out_fh, "\n");
  }

  // Close "IBD" filehandle
  gzclose(out_fh);

  /////////////////////////////////////////////////
  // Print genotype posterior probabilities
  tmp_out = strdcat(out_prefix, ".geno");
  out_fh = open_gzfile(tmp_out, "wbT");
  delete [] tmp_out;

  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open GENO output file!");

  double pp[N_GENO], prior[N_GENO];

  for(uint64_t s = 1; s <= pars->n_sites; s++)
    for (uint64_t i = 0; i < pars->n_ind; i++){
      //calc_HWE(prior, pars->freq[s], pars->marg_prob[i][s][1]);
      calc_HWE(prior, pars->freq[s], (double) pars->path[i][s]);
      post_prob(pp, pars->geno_lkl[i][s], prior, N_GENO);
      conv_space(pp, N_GENO, exp);
      gzwrite(out_fh, pp, sizeof(double)*N_GENO);
    }

  // Close filehandle
  gzclose(out_fh);
}



// General thread function
void threadpool_add_task(threadpool_t *thread_pool, int type, double **ptr, double *F, bool F_fixed, double *alpha, bool alpha_fixed, double **e_prob, char *path, double *pos_dist, uint64_t length){
  pth_struct *p = new pth_struct;

  p->type = type;
  p->ptr = ptr;
  p->F = F;
  p->F_fixed = F_fixed;
  p->alpha = alpha;
  p->alpha_fixed = alpha_fixed;
  p->e_prob = e_prob;
  p->path = path;
  p->pos_dist = pos_dist;
  p->length = length;

  // Add task to thread pool
  int ret = threadpool_add(thread_pool, thread_slave, (void*) p, 0);
  if(ret == -1)
    error(__FUNCTION__, "invalid thread pool!");
  else if(ret == -2)
    error(__FUNCTION__, "thread pool lock failure!");
  else if(ret == -3)
    error(__FUNCTION__, "queue full!");
  else if(ret == -4)
    error(__FUNCTION__, "thread pool is shutting down!");
  else if(ret == -5)
    error(__FUNCTION__, "thread failure!");
}

void thread_slave(void *ptr){
  pth_struct* p = (pth_struct*) ptr;
  double F[N_STATES] = {1-*p->F, *p->F};

  if(p->type == 1)
    forward(p->ptr, F, *p->alpha, p->e_prob, p->pos_dist, p->length, N_STATES);
  else if(p->type == 2)
    backward(p->ptr, F, *p->alpha, p->e_prob, p->pos_dist, p->length, N_STATES);
  else if(p->type == 3)
    viterbi(p->ptr, F, *p->alpha, p->e_prob, p->path, p->pos_dist, p->length, N_STATES);
  else if(p->type == 4){
    double val[2] = {*p->F, *p->alpha};
    double l_bound[2] = {1/INF, 1/INF};
    double u_bound[2] = {1-l_bound[0], 10};
    int lims[2] = {2, 2};

    if(p->F_fixed){
      l_bound[0] = *p->F;
      u_bound[0] = *p->F;
    }
    if(p->alpha_fixed){
      l_bound[1] = *p->alpha;
      u_bound[1] = *p->alpha;
    }

    findmax_bfgs(2, val, (void*) p, &lkl, NULL, l_bound, u_bound, lims, -1);
    *p->F = val[0];
    *p->alpha = val[1];
  }else
    error(__FUNCTION__, "invalid thread task option!");

  delete p;
}



double lkl(const double *pars, const void *data){
  pth_struct* p = (pth_struct*) data;
  double **Fw = init_ptr(p->length+1, N_STATES, (double) 0);
  double lkl = 0;

  if(isnan(pars[0]) || isinf(pars[0]) ||
     isnan(pars[1]) || isinf(pars[1]) )
    lkl = INF; // Added due to a putative bug on the BFGS function
  else{
    double F[N_STATES] = {1-pars[0], pars[0]};
    lkl = forward(Fw, F, pars[1], p->e_prob, p->pos_dist, p->length, N_STATES);
  }

  free_ptr((void**) Fw, p->length+1);
  return -lkl;
}
