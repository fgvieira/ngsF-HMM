
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

char const* version = "0.0.1b";


int main (int argc, char** argv) {
  /////////////////////
  // Parse Arguments //
  /////////////////////
  params* pars = new params;
  init_pars(pars);
  parse_cmd_args(pars, argc, argv);

  // Adjust number of threads
  if(pars->n_threads > pars->n_ind){
    warn(__FUNCTION__, "adjusting threads (--n_threads) to match number of individuals!");
    pars->n_threads = pars->n_ind;
  }



  ///////////////////////
  // Adjust parameters //
  ///////////////////////
  // Get file total size
  struct stat st;
  if(stat(pars->in_geno, &st) != 0)
    error(__FUNCTION__, "cannot check GENO file size!");

  if(strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0){
    if(pars->verbose >= 1)
      printf("==> GZIP input file (not BINARY)\n");
    pars->in_bin = false;
  }else if(pars->n_sites == st.st_size/sizeof(double)/pars->n_ind/N_GENO){
    if(pars->verbose >= 1)
      printf("==> BINARY input file (always lkl)\n");
    pars->in_bin = true;
    pars->in_lkl = true;
  }else
    error(__FUNCTION__, "invalid/corrupt genotype input file!");



  /////////////////////
  // Read input data //
  /////////////////////
  // Read positions from file
  if(pars->verbose >= 1)
    printf("==> Getting sites coordinates\n");
  if(pars->in_pos)
    pars->pos_dist = read_pos(pars->in_pos, pars->n_sites);
  else
    pars->pos_dist = init_ptr(pars->n_sites+1, (double) INFINITY);
  // Convert position distances to Mb
  for(uint64_t s = 1; s <= pars->n_sites; s++)
    pars->pos_dist[s] /= 1e6;

  // Read data from GENO file
  if(pars->verbose >= 1)
    printf("> Reading data from file...\n");
  pars->geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_lkl, pars->in_loglkl, pars->n_ind, pars->n_sites);
  pars->geno_lkl_s = transp_matrix(pars->geno_lkl, pars->n_ind, pars->n_sites+1);

  // Read_geno always returns genos in logscale
  pars->in_loglkl = true;

  // If input is not genotypes, check whether to call genotypes
  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
	call_geno(pars->geno_lkl[i][s], N_GENO);

      // Ensure minimum GL allowed
      if( pars->geno_lkl[i][s][array_min_pos(pars->geno_lkl[i][s], N_GENO)] < log(0.001) ){
	for(uint64_t g = 0; g < N_GENO; g++)
	  if(pars->geno_lkl[i][s][g] < log(0.001))
	    pars->geno_lkl[i][s][g] = log(0.001);
	// Re-normalize GL
	post_prob(pars->geno_lkl[i][s], pars->geno_lkl[i][s], NULL, N_GENO, true);
      }
    }



  /////////////////////////////////////////////
  // Declare and initialize output variables //
  /////////////////////////////////////////////
  init_output(pars);



  //////////////////
  // Analyze Data //
  //////////////////
  // Create thread pool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, 2*pars->n_ind, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");

  // Run EM!
  EM(pars);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    printf("Freeing memory...\n");

  // thread pool
  threadpool_wait(pars->thread_pool);
  if(threadpool_destroy(pars->thread_pool, threadpool_graceful) != 0)
    error(__FUNCTION__, "cannot free thread pool!");

  // pars struct
  //free_ptr((void*) pars->in_geno);
  free_ptr((void***) pars->geno_lkl, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) pars->geno_lkl_s, pars->n_sites+1);
  free_ptr((void*) pars->pos_dist);
  free_ptr((void*) pars->freq);
  free_ptr((void**) pars->path, pars->n_ind);
  free_ptr((void***) pars->marg_prob, pars->n_ind, pars->n_sites+1);
  free_ptr((void*) pars->indF);
  free_ptr((void*) pars->alpha);
  free_ptr((void*) pars->ind_lkl);

  if(pars->verbose >= 1)
    printf("Done!\n");
  delete pars;

  return 0;
}
