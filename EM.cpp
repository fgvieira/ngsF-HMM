
#include "ngsF-HMM.hpp"
#include "HMM.cpp"


int EM (params *pars) {
  SIG_COND = true;
  catch_SIG();

  uint64_t iter = 0;
  double max_lkl_epsilon = -INFINITY;
  double *prev_ind_lkl = init_ptr(pars->n_ind, (double) -INFINITY);
  double *ind_lkl_epsilon = init_ptr(pars->n_ind, (double) -INFINITY);
  gzFile log_fh = NULL;



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



  // Open filehandle for iteration log
  if(pars->log){
    char *tmp_out = strdcat(pars->out_prefix, ".log.gz");
    log_fh = open_gzfile(tmp_out, (pars->log_bin ? "wb9" : "w9"), pars->n_sites+100);
    delete [] tmp_out;

    if(log_fh == NULL)
      error(__FUNCTION__, "cannot open LOG file!");

    if(pars->verbose >= 1)
      printf("==> Dumping initial values to LOG file\n");
    dump_data(log_fh, pars, pars->log_bin);
  }



  ////////////////////
  // Iteration loop //
  ////////////////////
  while((max_lkl_epsilon > pars->min_epsilon || iter < pars->min_iters) && iter < pars->max_iters && SIG_COND) {
    time_t iter_start = time(NULL);

    // Next Iteration...
    iter++;
    if(pars->verbose >= 1)
      printf("\nIteration %lu:\n", iter);

    // Run next EM iteration
    iter_EM(pars);

    // Check convergence criteria
    pars->tot_lkl = 0;
    for (uint64_t i = 0; i < pars->n_ind; i++){
      // Get total Lkl
      pars->tot_lkl += pars->ind_lkl[i];
      // Get per-indiv lkl epsilon
      ind_lkl_epsilon[i] = (pars->ind_lkl[i] - prev_ind_lkl[i]) / fabs(prev_ind_lkl[i]);
    }
    max_lkl_epsilon = ind_lkl_epsilon[array_max_pos(ind_lkl_epsilon, pars->n_ind)];
    // Save current LKLs
    cpy(prev_ind_lkl, pars->ind_lkl, pars->n_ind, sizeof(double));

    if(pars->verbose >= 3)
      printf("Lkl epsilon: %s\n", join(ind_lkl_epsilon, pars->n_ind, "\t"));

    // Dump iteration data
    if(pars->log && (iter == 1 || iter % pars->log == 0)){
      if(pars->verbose >= 1)
	printf("==> Dumping iteration to log file\n");
      dump_data(log_fh, pars, pars->log_bin);
      
      if(pars->verbose >= 1)    
	printf("==> Printing current iteration parameters\n");
      print_iter(pars->out_prefix, pars);
    }

    // Print iteration info..
    if(pars->verbose >= 1){
      time_t iter_end = time(NULL);
      printf("\tLogLkl: %.15f\t lkl epsilon: %.15f\ttime: %.0f (s)\n", pars->tot_lkl, max_lkl_epsilon, difftime(iter_end, iter_start) );
    }

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
    threadpool_add_task(pars->thread_pool, 3, Vi[i], pars->geno_lkl[i], &pars->indF[i], &pars->alpha[i], pars->freq, pars->path[i], pars->pos_dist, pars->n_sites);

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
  if(pars->log)
    gzclose(log_fh);
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
    threadpool_add_task(pars->thread_pool, 1, Fw[i], pars->geno_lkl[i], &pars->indF[i], &pars->alpha[i], pars->freq, NULL, pars->pos_dist, pars->n_sites);
    
  // Backward recursion
  time_t bwd_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Backward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++)
    threadpool_add_task(pars->thread_pool, 2, Bw[i], pars->geno_lkl[i], &pars->indF[i], &pars->alpha[i], pars->freq, NULL, pars->pos_dist, pars->n_sites);

  threadpool_wait(pars->thread_pool);


  /*
  // Lkl check! - relaxed (to 0.001) due to precision issues on large datasets
  for (uint64_t i = 0; i < pars->n_ind; i++)
    if( abs(logsum(Fw[i][pars->n_sites],2) - logsum(Bw[i][0],2)) > 0.001 ){
      printf("Ind %lu: %.15f\t%.15f (%.15f)\n", i, logsum(Fw[i][pars->n_sites],2), logsum(Bw[i][0],2), abs(logsum(Fw[i][pars->n_sites],2) - logsum(Bw[i][0],2)) );
      error(__FUNCTION__, "Fw and Bw lkl do not match!");
    }
  */


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
      threadpool_add_task(pars->thread_pool, 4, NULL, pars->geno_lkl[i], &pars->indF[i], &pars->alpha[i], pars->freq, NULL, pars->pos_dist, pars->n_sites);

    threadpool_wait(pars->thread_pool);

    if(pars->verbose >= 4)
      for(uint64_t i = 0; i < pars->n_ind; i++)
	printf("\t%.10f\t%f\n", pars->indF[i], pars->alpha[i]);
  }



  // Estimate allele frequencies (EM)
  time_t freqs_t = time(NULL);
  if(pars->freq_fixed){
    if(pars->verbose >= 1)
      printf("==> Alelle frequencies not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Estimate allele frequencies\n");

    double *indF = init_ptr(pars->n_ind, 0.0);
    for (uint64_t s = 1; s <= pars->n_sites; s++){
      for(uint64_t i = 0; i < pars->n_ind; i++)
	indF[i] = pars->marg_prob[i][s][1];

      pars->freq[s] = est_maf(pars->n_ind, pars->geno_lkl_s[s], indF);

      if(pars->verbose >= 7)
	printf("%lu; freq: %f\n", s, pars->freq[s]);
    }
    free_ptr((void*) indF);
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
  char *tmp_out;
  gzFile out_fh;

  // Open filehandle to "indF" file
  tmp_out = strdcat(out_prefix, ".indF");
  out_fh = open_gzfile(tmp_out, "wT");
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
  out_fh = open_gzfile(tmp_out, "wT", pars->n_sites+100);
  delete [] tmp_out;

  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open IBD output file!");

  // Print IBD info: most probable path (Viterbi) and IBD marg probs
  dump_data(out_fh, pars, false);

  // Close "IBD" filehandle
  gzclose(out_fh);

  /////////////////////////////////////////////////
  // Print genotype posterior probabilities
  tmp_out = strdcat(out_prefix, ".geno");
  out_fh = open_gzfile(tmp_out, "wbT");
  delete [] tmp_out;

  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open GENO output file!");

  double pp[N_GENO];

  for(uint64_t s = 1; s <= pars->n_sites; s++)
    for (uint64_t i = 0; i < pars->n_ind; i++){
      double prior[3];
      //calc_prior(prior, pars->freq[s], pars->marg_prob[i][s][1]);
      calc_prior(prior, pars->freq[s], (double) pars->path[i][s]);
      post_prob(pp, pars->geno_lkl[i][s], prior, N_GENO);
      conv_space(pp, N_GENO, exp);
      gzwrite(out_fh, pp, sizeof(double)*N_GENO);
    }

  // Close filehandle
  gzclose(out_fh);
}



void dump_data(gzFile fh, params *pars, bool out_bin){
  char *buf;

  if(out_bin){
    // Print Lkl
    gzwrite(fh, pars->ind_lkl, sizeof(double)*pars->n_ind);

    // Print most probable path (Viterbi)
    for (uint64_t i = 0; i < pars->n_ind; i++)
      gzwrite(fh, pars->path[i]+1, sizeof(char)*pars->n_sites);

    // Print marginal probs
    for (uint64_t i = 0; i < pars->n_ind; i++)
      for (uint64_t s = 1; s <= pars->n_sites; s++){
	gzwrite(fh, &pars->marg_prob[i][s][1], sizeof(double));
      }
  }else{
    buf = join(pars->ind_lkl, pars->n_ind, "\t");
    // Print Lkl
    if(gzprintf(fh, "//\t%s\n", buf) <= 0)
      error(__FUNCTION__, "cannot write LKL info to file!");
    delete [] buf;

    // Print most probable path (Viterbi)
    for (uint64_t i = 0; i < pars->n_ind; i++){
      buf = join(pars->path[i]+1, pars->n_sites, "");
      if(gzprintf(fh, "%s\n", buf) <= 0)
	error(__FUNCTION__, "cannot write PATH info to file!");
      delete [] buf;
    }

    // Print marginal probs
    for (uint64_t i = 0; i < pars->n_ind; i++){
      // To avoid leading \t
      gzprintf(fh, "%f", pars->marg_prob[i][1][1]);
      for (uint64_t s = 2; s <= pars->n_sites; s++)
	gzprintf(fh, "\t%f", pars->marg_prob[i][s][1]);
      gzprintf(fh, "\n");
    }
  }
}
