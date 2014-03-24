
#include "ngsF-HMM.hpp"



int EM (params *pars, out_data *data) {
  SIG_COND = true;
  catch_SIG();

  char out_prefix[1000];

  uint64_t iter = 0;
  double  max_lkl_epsilon = -INFINITY;
  double *lkl_epsilon = init_double(pars->n_ind, -INFINITY);
  double *prev_lkl = init_double(pars->n_ind, -INFINITY);
  gzFile log_fh;

  // Open filehandle for iteration log
  if(pars->log){
    strcpy(out_prefix, pars->out_prefix);
    if((log_fh = gzopen( strcat(out_prefix,".log.gz"), pars->log_bin ? "wb" : "w")) == NULL)
      error(__FUNCTION__, "cannot open LOG file!");

    if(gzbuffer(log_fh, pars->n_sites + pars->n_ind*10) < 0)
      error(__FUNCTION__, "cannot increase ZLIB buffer size!");
  }

  while((max_lkl_epsilon > pars->min_epsilon || iter < pars->min_iters) && iter < pars->max_iters && SIG_COND) {
	  
    time_t iter_start = time(NULL);
	  
    // Next Iteration...
    iter++;
    if(pars->verbose >= 1)
      printf("\nIteration %lu:\n", iter);
	  
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
      double sum = 0;
      time_t iter_end = time(NULL);

      for(uint64_t i = 0; i < pars->n_ind; i++)
	sum += data->lkl[i];

      printf("\tLogLkl: %.15f\t lkl epsilon: %.15f\ttime: %.0f (s)\n", sum, max_lkl_epsilon, difftime(iter_end, iter_start) );
    }

    fflush(stdout);



    /////////////////////////
    // Dump iteration data //
    /////////////////////////
    char *buf;
    if(pars->log){
      if(pars->verbose >= 1)
	printf("==> Dumping iteration to log file\n");

      if(pars->log_bin){
	// Print Lkl
	gzwrite(log_fh, data->lkl, sizeof(double)*pars->n_ind);

	// Print most probable path (Viterbi)
	for (uint64_t i = 0; i < pars->n_ind; i++)
	  gzwrite(log_fh, data->path[i]+1, sizeof(uint)*pars->n_sites);

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

    /* Disabled since printing data each iteration takes too long
    printf("==> Printing current iteration parameters\n");
    strcpy(out_prefix, pars->out_prefix);
    print_iter(strcat(out_prefix,".out"), pars, data);
    */
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
  double inbreed = 0.1;
  double pp[N_GENO];

  double ***a = init_double(pars->n_ind, N_STATES, N_STATES, -INFINITY);
  cpy(a, data->a, pars->n_ind, N_STATES, N_STATES, sizeof(double));
  double ***e = init_double(pars->n_sites+1, N_STATES, N_GENO, -INFINITY);
  cpy(e, data->e, pars->n_sites+1, N_STATES, N_GENO, sizeof(double));
  uint **path = init_uint(pars->n_ind, pars->n_sites+1, 0);
  cpy(path, data->path, pars->n_ind, pars->n_sites+1, sizeof(uint));

  double ***Fw = init_double(pars->n_ind, pars->n_sites+1, N_STATES, 0);
  double ***Bw = init_double(pars->n_ind, pars->n_sites+1, N_STATES, 0);
  double ***Vi = init_double(pars->n_ind, pars->n_sites+1, N_STATES, 0);

  if(pars->verbose >= 1)
    printf("==> Forward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    // Initialise Forward table
    Fw[i][0][0] = log(1-inbreed);
    Fw[i][0][1] = log(inbreed);

    for (uint64_t s = 1; s <= pars->n_sites; s++){
      post_prob(pp, pars->geno_lkl[i][s], e[s][path[i][s]], N_GENO);
      if(pars->call_geno == 2) 
	call_geno(pp, N_GENO, true);

      for(uint64_t l = 0; l < N_STATES; l++){
	// logsum(k==0,k==1)
	Fw[i][s][l] = logsum2(Fw[i][s-1][0] + a[i][0][l],
			      Fw[i][s-1][1] + a[i][1][l]);
	Fw[i][s][l] += logsum3(e[s][l][0]+pp[0], e[s][l][1]+pp[1], e[s][l][2]+pp[2]);
      }
    }
  }
    

  if(pars->verbose >= 1)
    printf("==> Backward Recursion\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    // Initialise Backward table
    Bw[i][pars->n_sites][0] = log(1);
    Bw[i][pars->n_sites][1] = log(1);

    for (uint64_t s = pars->n_sites; s > 0; s--){ 
      post_prob(pp, pars->geno_lkl[i][s], e[s][path[i][s]], N_GENO);
      if(pars->call_geno == 2) 
	call_geno(pp, N_GENO, true);

      double LS_0 = logsum3(e[s][0][0]+pp[0], e[s][0][1]+pp[1], e[s][0][2]+pp[2]);
      double LS_1 = logsum3(e[s][1][0]+pp[0], e[s][1][1]+pp[1], e[s][1][2]+pp[2]);
	
      for(uint64_t k = 0; k < N_STATES; k++)
	// logsum(l==0,l==1)
	Bw[i][s-1][k] = logsum2(a[i][k][0] + LS_0 + Bw[i][s][0],
				a[i][k][1] + LS_1 + Bw[i][s][1]);
    }
  }
  if(pars->verbose >= 1)
    printf("==> Backward Recursion termination\n");
  for (uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t k = 0; k < N_STATES; k++)
      Bw[i][0][k] += Fw[i][0][k];



  if(pars->verbose >= 1)
    printf("==> Marginal probabilities\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    data->lkl[i] = logsum(Fw[i][pars->n_sites],2);
    for (uint64_t s = 1; s <= pars->n_sites; s++)
      for(uint64_t k = 0; k < N_STATES; k++)
	data->marg_prob[i][s][k] = Bw[i][s][k] + Fw[i][s][k] - data->lkl[i];
  }
    

    
  if(pars->path_fixed){
    if(pars->verbose >= 1)
      printf("==> Most probable path not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update most probable path (Viterbi)\n");
    for (uint64_t i = 0; i < pars->n_ind; i++){
      // Initialise Forward table
      Vi[i][0][0] = log(1-inbreed);
      Vi[i][0][1] = log(inbreed);
	
      for (uint64_t s = 1; s <= pars->n_sites; s++){
	post_prob(pp, pars->geno_lkl[i][s], e[s][path[i][s]], N_GENO);
	if(pars->call_geno == 2)
	  call_geno(pp, N_GENO, true);

	for(uint64_t l = 0; l < N_STATES; l++){
	  // max(k==0,k==1)
	  Vi[i][s][l] = max(Vi[i][s-1][0] + a[i][0][l], 
			    Vi[i][s-1][1] + a[i][1][l]);
	  Vi[i][s][l] += logsum3(e[s][l][0]+pp[0], e[s][l][1]+pp[1], e[s][l][2]+pp[2]);
	}
      }
    }

    if(pars->verbose >= 1)
      printf("> Back-tracking\n");
    for (uint64_t i = 0; i < pars->n_ind; i++)
      for (uint64_t s = 1; s <= pars->n_sites; s++)
	data->path[i][s] = (Vi[i][s][0] > Vi[i][s][1] ? 0 : 1);
  }
    
    
  if(pars->trans_fixed){
    if(pars->verbose >= 1)
      printf("==> Transition probabilities not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update transition probabilities\n");
    double *sPk = init_double(pars->n_sites+1, -INFINITY);

    for (uint64_t i = 0; i < pars->n_ind; i++){
      for(uint64_t k = 0; k < N_STATES; k++){
	// Get P(k)
	for (uint64_t s = 1; s <= pars->n_sites; s++)
	  sPk[s] = Fw[i][s][k] + Bw[i][s][k];
	// Sum all site_Pk, skipping site 0 and last
	double Pk = logsum(sPk+1, pars->n_sites-1);

	for(uint64_t l = 0; l < N_STATES; l++){
	  double tmp_a = -INFINITY;
	  for (uint64_t s = 1; s < pars->n_sites; s++){
	    post_prob(pp, pars->geno_lkl[i][s+1], e[s+1][path[i][s+1]], N_GENO);
	    if(pars->call_geno == 2)
	      call_geno(pp, N_GENO, true);
	    double LS = logsum3(e[s+1][l][0]+pp[0],
				e[s+1][l][1]+pp[1],
				e[s+1][l][2]+pp[2]);
	    tmp_a = logsum2(tmp_a, 
			    Fw[i][s][k] + a[i][k][l] + LS + Bw[i][s+1][l] - Pk);
	  }
	  data->a[i][k][l] = tmp_a;
	}
      }
    }

    delete [] sPk;
  }



  if(pars->freq_fixed){
    if(pars->verbose >= 1)
      printf("==> Alelle frequencies not estimated!\n");
  }else{
    if(pars->verbose >= 1)
      printf("==> Update allele frequencies\n");
    for (uint64_t s = 1; s <= pars->n_sites; s++){
      // Expected number minor alleles
      double num = 0;
      // Expected total number of alleles
      double den = 0;

      for (uint64_t i = 0; i < pars->n_ind; i++){
	post_prob(pp, pars->geno_lkl[i][s], e[s][path[i][s]], N_GENO);
	if(pars->call_geno == 2)
	  call_geno(pp, N_GENO, true);

	num = num + exp(pp[1]) + exp(pp[2])*(2-exp(data->marg_prob[i][s][1]));
	den = den + 2*exp(pp[1]) + exp(logsum2(pp[0],pp[2]))*(2-exp(data->marg_prob[i][s][1]));
	if(pars->verbose >= 8)
	  printf("%lu %lu; num: %f; den; %f; pp: %f %f %f; marg: %f\n", s, i, num, den, exp(pp[0]), exp(pp[1]), exp(pp[2]), exp(data->marg_prob[i][s][1]));
      }
      data->freq[s] = num/den;
    }
  }



  if(pars->verbose >= 1)
    printf("==> Update emission probabilities\n");
  update_e(data, pars->n_sites);



  if(pars->verbose >= 1)
    printf("==> Update inbreeding coefficients\n");
  for (uint64_t i = 0; i < pars->n_ind; i++){
    double indF = -INFINITY;
    for (uint64_t s = 1; s <= pars->n_sites; s++)
      indF = logsum2(indF, data->marg_prob[i][s][1]);
    data->indF[i] = exp(indF)/pars->n_sites;
  }


  free_ptr((void***) a, pars->n_ind, N_STATES);
  free_ptr((void***) e, pars->n_sites+1, N_STATES);
  free_ptr((void**) path, pars->n_ind);

  free_ptr((void***) Fw, pars->n_ind, pars->n_sites+1);
  free_ptr((void***) Bw, pars->n_ind, pars->n_sites+1);
  free_ptr((void***) Vi, pars->n_ind, pars->n_sites+1);
}




void post_prob(double *pp, double *lkl, double *prior, uint64_t n_geno){
  for(uint64_t cnt = 0; cnt < n_geno; cnt++)
    pp[cnt] = lkl[cnt];//+prior[cnt];

  double norm = logsum(pp, n_geno);

  for(uint64_t cnt = 0; cnt < n_geno; cnt++)
    pp[cnt] -= norm;
}




void print_iter(char *out_file, params *pars, out_data *data){
  char *buf;
  // Open filehandle to write
  FILE *out_fh = fopen(out_file, "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open iteration file!");

  // Print indF
  fprintf(out_fh,"%s\n", join(data->indF, pars->n_ind, "\t"));
  
  // Print transition prob
  for(uint16_t i = 0; i < pars->n_ind; i++)
    fprintf(out_fh,"%f\t%f\n", exp(data->a[i][0][1]), exp(data->a[i][0][1]));

  // Print most probable path (Viterbi)
  for (uint64_t i = 0; i < pars->n_ind; i++){
    buf = join(data->path[i]+1, pars->n_sites, "");
    fprintf(out_fh, "%s\n", buf);
    delete [] buf;
  }
  
  // Print allele freqs
  for(uint64_t s = 1; s <= pars->n_sites; s++)
    fprintf(out_fh,"%f\n", data->freq[s]);
  
  // Close filehandle
  fclose(out_fh);
}




int update_e(out_data *data, uint64_t n_sites){
  for(uint64_t s = 1; s <= n_sites; s++)
    for(uint64_t F = 0; F < N_STATES; F++){
      double f = data->freq[s];
      data->e[s][F][0] = log(pow(1-f,2)+(1-f)*f*F);
      data->e[s][F][1] = log(2*(1-f)*f-2*(1-f)*f*F);
      data->e[s][F][2] = log(pow(f,2)+(1-f)*f*F);
    }

  return 0;
}
