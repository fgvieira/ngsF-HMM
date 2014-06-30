
/*
 *
 * ngsF-HMM - NGS data individual inbreeding coefficients estimation.
 * Copyright (C) 2012  Filipe G. Vieira
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "ngsF-HMM.hpp"
#include "read_data.cpp"

char const* version = "0.0.1b";


int main (int argc, char** argv) {
  /////////////////////
  // Parse Arguments //
  /////////////////////
  params* pars = new params;
  init_pars(pars);
  parse_cmd_args(argc, argv, pars);
  if(pars->version) {
    printf("ngsF-HMM v%s\nCompiled on %s @ %s", version, __DATE__, __TIME__);
    exit(0);
  }

  if(pars->verbose >= 1){
    printf("==> Input Arguments:\n");
    printf("\tgeno file: %s\n\tgeno lkl: %s\n\tgeno loglkl: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tcall_geno: %d\n\tfreq: %s\n\tfreq_fixed: %s\n\ttrans: %s\n\ttrans_fixed: %s\n\tpath: %s\n\tpath_fixed: %s\n\tout prefix: %s\n\tlog: %s\n\tlog_bin: %s\n\tmin_iters: %d\n\tmax_iters: %d\n\tmin_epsilon: %.10f\n\tn_threads: %d\n\tversion: %s\n\tverbose: %d\n\tseed: %d\n\n",
	   pars->in_geno, pars->in_lkl ? "true":"false", pars->in_loglkl ? "true":"false", pars->n_ind, pars->n_sites, pars->call_geno, pars->in_freq, pars->freq_fixed ? "true":"false", pars->in_trans, pars->trans_fixed ? "true":"false", pars->in_path, pars->path_fixed ? "true":"false", pars->out_prefix, pars->log ? "true":"false", pars->log_bin ? "true":"false", pars->min_iters, pars->max_iters, pars->min_epsilon, pars->n_threads, pars->version ? "true":"false", pars->verbose, pars->seed);
  }
  if(pars->verbose >= 4)
    printf("==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");
  


  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "Genotype input file (-geno) missing!");
  if(pars->out_prefix == NULL)
    error(__FUNCTION__, "Output prefix (-out_prefix) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "Number of individuals (-n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "Number of sites (-n_sites) missing!");
  
  
  
  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  struct stat st;
  if(stat(pars->in_geno, &st) != 0)
    error(__FUNCTION__, "cannot check file size!");

  if(strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0){
    if(pars->verbose >= 1)
      printf("==> GZIP input file (never BINARY)\n");
    pars->in_bin = false;
  }else if(pars->n_sites == st.st_size/sizeof(double)/pars->n_ind/N_GENO){
    if(pars->verbose >= 1)
      printf("==> BINARY input file (always loglkl)\n");
    pars->in_bin = true;
    pars->in_lkl = true;
    pars->in_loglkl = true;
  }else
    error(__FUNCTION__, "invalid/corrupt genotype input file!");
  


  ////////////////////////////
  // Prepare initial values //
  ////////////////////////////
  // Create thread pool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, pars->n_ind, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");

  // Declare and initialize output variables
  out_data* data = new out_data;
  init_output(pars, data);

  // Read data from GENO file
  pars->geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_lkl, pars->n_ind, pars->n_sites);
  
  // If input not called genotypes, check whether to call genotypes and/or convert to log-space
  if(pars->in_lkl)
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++){
	// Call genotypes
	if(pars->call_geno == 1)
	  call_geno(pars->geno_lkl[i][s], N_GENO, pars->in_loglkl);
	// Convert space
	if(!pars->in_loglkl)
	  conv_space(pars->geno_lkl[i][s], N_GENO, log);
      }

  

  //////////////////
  // Analyze Data //
  //////////////////
  EM(pars, data);
  if(pars->verbose >= 1){
    double sum = 0;
    for(uint64_t i = 0; i < pars->n_ind; i++)
      sum += data->lkl[i];
    printf("\nFinal logLkl: %f\n", sum);
  }
  


  /////////////////////////
  // Print Final Results //
  /////////////////////////
  if(pars->verbose >= 1)
    printf("==> Printing final results\n");
  print_iter(pars->out_prefix, pars, data);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    printf("Freeing memory...\n");

  // thread pool
  threadpool_wait(pars->thread_pool);
  if(threadpool_destroy(pars->thread_pool, threadpool_graceful) != 0)
    error(__FUNCTION__, "cannot free thread pool!");

  // data struct
  free_ptr((void***) data->a, pars->n_ind, N_STATES);
  free_ptr((void*) data->freq);
  free_ptr((void***) data->e, pars->n_sites+1, N_STATES);
  free_ptr((void**) data->path, pars->n_ind);
  free_ptr((void***) data->marg_prob, pars->n_ind, pars->n_sites+1);
  free_ptr((void*) data->indF);
  free_ptr((void*) data->lkl);
  delete data;

  // pars struct
  //free_ptr((void*) pars->in_geno);
  free_ptr((void***) pars->geno_lkl, pars->n_ind, pars->n_sites+1);

  if(pars->verbose >= 1)
    printf("Done!\n");
  delete pars;

  return 0;
}
