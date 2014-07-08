
#include "ngsF-HMM.hpp"
#include "HMM.cpp"


int EM (params *pars, out_data *data) {
  SIG_COND = true;
  catch_SIG();

  uint64_t iter = 0;
  double  max_lkl_epsilon = -INFINITY;
  double *lkl_epsilon = init_ptr(pars->n_ind, -INFINITY);
  double *prev_lkl = init_ptr(pars->n_ind, -INFINITY);
  gzFile log_fh;

  // Open filehandle for iteration log
  if(pars->log){
    char *tmp_out = strdcat(pars->out_prefix, ".log.gz");
    if((log_fh = gzopen(tmp_out, pars->log_bin ? "wb" : "w")) == NULL)
      error(__FUNCTION__, "cannot open LOG file!");

    if(gzbuffer(log_fh, pars->n_sites + pars->n_ind*10) < 0)
      error(__FUNCTION__, "cannot increase ZLIB buffer size!");
    delete [] tmp_out;
  }

  // Loop through all iterations
  while((max_lkl_epsilon > pars->min_epsilon || iter < pars->min_iters) && iter < pars->max_iters && SIG_COND) {
	  
    time_t iter_start = time(NULL);
	  
    // Next Iteration...
    iter++;
    if(pars->verbose >= 1)
      printf("\nIteration %lu:\n", iter);

    // Run next EM iteration
    iter_EM(pars, data);

    // Check convergence criteria
    for (uint64_t i = 0; i < pars->n_ind; i++)
      lkl_epsilon[i] = (data->lkl[i] - prev_lkl[i]) / fabs(prev_lkl[i]);
    max_lkl_epsilon = lkl_epsilon[array_max_pos(lkl_epsilon, pars->n_ind)];
    // Save current LKLs
    cpy(prev_lkl, data->lkl, pars->n_ind, sizeof(double));

    if(pars->verbose >= 2)
      printf("Lkl epsilon: %s\n", join(lkl_epsilon, pars->n_ind, "\t"));

    if(pars->verbose >= 1){
      // Get total lkl
      double sum = 0;
      for(uint64_t i = 0; i < pars->n_ind; i++)
	sum += data->lkl[i];

      time_t iter_end = time(NULL);
      printf("\tLogLkl: %.15f\t lkl epsilon: %.15f\ttime: %.0f (s)\n", sum, max_lkl_epsilon, difftime(iter_end, iter_start) );
    }

    fflush(stdout);



    /////////////////////////
    // Dump iteration data //
    /////////////////////////
    char *buf;
    if(pars->log && (iter == 1 || iter % pars->log == 0)){
      if(pars->verbose >= 1)
	printf("==> Dumping iteration to log file\n");

      if(pars->log_bin){
	// Print Lkl
	gzwrite(log_fh, data->lkl, sizeof(double)*pars->n_ind);

	// Print most probable path (Viterbi)
	for (uint64_t i = 0; i < pars->n_ind; i++)
	  gzwrite(log_fh, data->path[i]+1, sizeof(char)*pars->n_sites);

	// Print marginal probs
	for (uint64_t i = 0; i < pars->n_ind; i++)
	  for (uint64_t s = 1; s <= pars->n_sites; s++){
	    double mp = exp(data->marg_prob[i][s][1]);
	    gzwrite(log_fh, &mp, sizeof(double));
	  }
      }else{
	buf = join(data->lkl, pars->n_ind, "\t");
	// Print Lkl
	if(gzprintf(log_fh, "//\t%s\n", buf) <= 0)
	  error(__FUNCTION__, "cannot write LKL info to LOG file!");
	delete [] buf;

	// Print most probable path (Viterbi)
	for (uint64_t i = 0; i < pars->n_ind; i++){
	  buf = join(data->path[i]+1, pars->n_sites, "");
	  /*
	  buf = init_ptr(pars->n_sites+1, '\0');
	  for(uint64_t s = 1; s <= pars->n_sites; s++)
	    sprintf(&buf[s-1], "%d", data->path[i][s]);
	  */
	  if(gzprintf(log_fh, "%s\n", buf) <= 0)
	    error(__FUNCTION__, "cannot write PATH info to LOG file!");
	  delete [] buf;
	}

	// Print marginal probs
	for (uint64_t i = 0; i < pars->n_ind; i++){
	  // To avoid leading \t
	  gzprintf(log_fh, "%f", exp(data->marg_prob[i][1][1]));
	  for (uint64_t s = 2; s <= pars->n_sites; s++)
	    gzprintf(log_fh, "\t%f", exp(data->marg_prob[i][s][1]));
	  gzprintf(log_fh, "\n");
	}
      }
    }

    // Disabled since printing data each iteration takes too long
    if(pars->verbose >= 1)    
      printf("==> Printing current iteration parameters\n");
    print_iter(pars->out_prefix, pars, data);
  }


  // Close filehandle for iteration log
  if(pars->log)
    gzclose(log_fh);

  
  if(iter >= pars->max_iters)
    printf("WARN: Maximum number of iterations reached! Check if analysis converged... \n");

  // Free memory and return
  free_ptr((void*) lkl_epsilon);
  free_ptr((void*) prev_lkl);
  return 0;
}



void iter_EM(params *pars, out_data *data) {
  double ***a = init_ptr(pars->n_ind, N_STATES, N_STATES, -INFINITY);
  cpy(a, data->a, pars->n_ind, N_STATES, N_STATES, sizeof(double));
  double ***prior = init_ptr(pars->n_sites+1, N_STATES, N_GENO, -INFINITY);
  cpy(prior, data->prior, pars->n_sites+1, N_STATES, N_GENO, sizeof(double));
  char **path = init_ptr(pars->n_ind, pars->n_sites+1, (const char*) '\0');
  cpy(path, data->path, pars->n_ind, pars->n_sites+1, sizeof(char));

  double ***Fw = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, 0.0);
  double ***Bw = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, 0.0);
  double ***Vi = init_ptr(pars->n_ind, pars->n_sites+1, N_STATES, 0.0);


  // Forward recursion
  time_t fwd_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Forward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++) {
    Fw[i][0][0] = log(1-data->indF[i]);
    Fw[i][0][1] = log(data->indF[i]);

    threadpool_add_task(pars->thread_pool, 1, NULL, pars->geno_lkl[i], Fw[i], NULL, NULL, prior, a[i], path[i], pars->n_sites);
  }
  threadpool_wait(pars->thread_pool);
    
  // Backward recursion
  time_t bwd_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Backward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    Bw[i][pars->n_sites][0] = log(1);
    Bw[i][pars->n_sites][1] = log(1);

    threadpool_add_task(pars->thread_pool, 2, NULL, pars->geno_lkl[i], NULL, Bw[i], NULL, prior, a[i], path[i], pars->n_sites);
  }
  threadpool_wait(pars->thread_pool);

  if(pars->verbose >= 1)
    printf("> Backward Recursion termination\n");
  for (uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t k = 0; k < N_STATES; k++)
      Bw[i][0][k] += Fw[i][0][k];

  // Marginal probabilities
  time_t mp_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Marginal probabilities\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    data->lkl[i] = logsum(Fw[i][pars->n_sites],2);
    for (uint64_t s = 1; s <= pars->n_sites; s++)
      for(uint64_t k = 0; k < N_STATES; k++)
	data->marg_prob[i][s][k] = Bw[i][s][k] + Fw[i][s][k] - data->lkl[i];
  }
  
  // Viterbi
  time_t path_t = time(NULL);
  if(pars->path_fixed){
    if(pars->verbose >= 1)
      printf("==> Most probable path not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update most probable path (Viterbi)\n");
    for (uint64_t i = 0; i < pars->n_ind; i++){
      Vi[i][0][0] = log(1-data->indF[i]);
      Vi[i][0][1] = log(data->indF[i]);

      threadpool_add_task(pars->thread_pool, 3, NULL, pars->geno_lkl[i], NULL, NULL, Vi[i], prior, a[i], path[i], pars->n_sites);
    }
    threadpool_wait(pars->thread_pool);

    if(pars->verbose >= 1)
      printf("> Back-tracking\n");
    for (uint64_t i = 0; i < pars->n_ind; i++)
      for (uint64_t s = 1; s <= pars->n_sites; s++)
	data->path[i][s] = (Vi[i][s][0] > Vi[i][s][1] ? 0 : 1);
  }
  
  // Estimate transition probabilities
  time_t trans_t = time(NULL);
  if(pars->trans_fixed){
    if(pars->verbose >= 1)
      printf("==> Transition probabilities not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update transition probabilities\n");

    for (uint64_t i = 0; i < pars->n_ind; i++)
      threadpool_add_task(pars->thread_pool, 4, data->a[i], pars->geno_lkl[i], Fw[i], Bw[i], NULL, prior, a[i], path[i], pars->n_sites);
  }
  threadpool_wait(pars->thread_pool);

  // Estimate allele frequencies
  time_t freqs_t = time(NULL);
  if(pars->freq_fixed){
    if(pars->verbose >= 1)
      printf("==> Alelle frequencies not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update allele frequencies\n");
    double pp[N_GENO];

    for (uint64_t s = 1; s <= pars->n_sites; s++){
      // Expected number minor alleles
      double num = 0;
      // Expected total number of alleles
      double den = 0;

      for (uint64_t i = 0; i < pars->n_ind; i++){
	//post_prob(pp, pars->geno_lkl[i][s], prior[s][(int) path[i][s]], N_GENO);
	post_prob(pp, pars->geno_lkl[i][s], NULL, N_GENO);

	num = num + exp(pp[1]) + exp(pp[2])*(2-exp(data->marg_prob[i][s][1]));
	den = den + 2*exp(pp[1]) + exp(logsum2(pp[0],pp[2]))*(2-exp(data->marg_prob[i][s][1]));
	if(pars->verbose >= 8)
	  printf("%lu %lu; num: %f; den; %f; pp: %f %f %f; marg: %f\n", s, i, num, den, exp(pp[0]), exp(pp[1]), exp(pp[2]), exp(data->marg_prob[i][s][1]));
      }
      data->freq[s] = num/den;
    }
  }

  // Update emission probabilities
  time_t emission_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Update emission probabilities\n");
  update_priors(data, pars->n_sites);

  // Calculate inbreeding coefficients
  time_t F_t = time(NULL);
  if(pars->verbose >= 1)
    printf("==> Update inbreeding coefficients\n");
  for (uint64_t i = 0; i < pars->n_ind; i++)
    data->indF[i] = exp(data->a[i][0][1] - logsum2(data->a[i][0][1], data->a[i][1][0]));



  time_t end_t = time(NULL);
  if(pars->verbose >= 3)
    printf("\nFw: %.1f\nBw: %.1f\nMP: %.1f\npath: %.1f\ntrans: %.1f\nfreqs: %.1f\nemission: %.1f\nF: %.1f\n", 
	   difftime(bwd_t,fwd_t), 
	   difftime(mp_t,bwd_t), 
	   difftime(path_t,mp_t), 
	   difftime(trans_t,path_t), 
	   difftime(freqs_t,trans_t), 
	   difftime(emission_t,freqs_t), 
	   difftime(F_t,emission_t), 
	   difftime(end_t,F_t)
	   );

  free_ptr((void***) a, pars->n_ind, N_STATES);
  free_ptr((void***) prior, pars->n_sites+1, N_STATES);
  free_ptr((void**) path, pars->n_ind);

  free_ptr((void***) Fw, pars->n_ind, pars->n_sites+1);
  free_ptr((void***) Bw, pars->n_ind, pars->n_sites+1);
  free_ptr((void***) Vi, pars->n_ind, pars->n_sites+1);
}




void post_prob(double *pp, double *lkl, double *prior, uint64_t n_geno){
  for(uint64_t cnt = 0; cnt < n_geno; cnt++){
    pp[cnt] = lkl[cnt];
    if(prior != NULL)
      pp[cnt] += prior[cnt];
  }

  double norm = logsum(pp, n_geno);

  for(uint64_t cnt = 0; cnt < n_geno; cnt++)
    pp[cnt] -= norm;
}




void print_iter(char *out_prefix, params *pars, out_data *data){
  char *tmp_out;
  FILE *out_fh;

  // Open filehandle to "indF" file
  tmp_out = strdcat(out_prefix, ".indF");
  out_fh = fopen(tmp_out, "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open \"indF\" output file!");

  // Print total Lkl
  double sum = 0;
  for(uint64_t i = 0; i < pars->n_ind; i++)
    sum += data->lkl[i];
  fprintf(out_fh,"%.10f\n", sum);

  // Print indF
  fprintf(out_fh,"%s\n", join(data->indF, pars->n_ind, "\n"));
  
  // Print transition prob
  for(uint16_t i = 0; i < pars->n_ind; i++)
    fprintf(out_fh,"%f\t%f\n", exp(data->a[i][0][1]), exp(data->a[i][1][0]));

  // Print allele freqs
  for(uint64_t s = 1; s <= pars->n_sites; s++)
    fprintf(out_fh,"%f\n", data->freq[s]);

  // Close "indF" filehandle
  fclose(out_fh);
  delete [] tmp_out;

  /////////////////////////////////////////////////
  // Open filehandle to "viterbi" file
  tmp_out = strdcat(out_prefix, ".viterbi");
  out_fh = fopen(tmp_out, "wb");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open \"viterbi\" output file!");
  
  // Print most probable path (Viterbi)
  for(uint64_t i = 0; i < pars->n_ind; i++)
    fwrite(data->path[i]+1, sizeof(char), pars->n_sites, out_fh);

  // Close "viterbi" filehandle
  fclose(out_fh);
  delete [] tmp_out;

  /////////////////////////////////////////////////
  // Print viterbi path in text mode
  /*  
  char *buf;
  // Open filehandle to "viterbi.txt" file
  tmp_out = strdcat(out_prefix, ".viterbi.txt");
  out_fh = fopen(tmp_out, "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open \"viterbi.txt\" output file!");

  // Print most probable path (Viterbi)
  for(uint64_t i = 0; i < pars->n_ind; i++){
    buf = join(data->path[i]+1, pars->n_sites, "");
    fprintf(out_fh, "%s\n", buf);
    delete [] buf;
  }
  
  // Close "viterbi.txt" filehandle
  fclose(out_fh);
  delete [] tmp_out;
  */
}




int update_priors(out_data *data, uint64_t n_sites){
  for(uint64_t s = 1; s <= n_sites; s++)
    for(uint64_t F = 0; F < N_STATES; F++){
      double f = data->freq[s];
      data->prior[s][F][0] = log(pow(1-f,2)+(1-f)*f*F);
      data->prior[s][F][1] = log(2*(1-f)*f-2*(1-f)*f*F);
      data->prior[s][F][2] = log(pow(f,2)+(1-f)*f*F);
    }

  return 0;
}
