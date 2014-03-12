#include <gsl/gsl_rng.h>
#include "ngsF-HMM.hpp"



int init_output(params* pars, out_data* data) {
  char* buf = new char[BUFF_LEN];
  double* t;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, pars->seed);


  ////////////////////////////////////////
  // Set TRANSITION initial values ...
  ////////////////////////////////////////
  double trans_rng_max = 1;
  double trans_rng_min = 0;
  gzFile in_trans_fh;

  data->a = init_double(pars->n_ind, N_STATES, N_STATES, -INFINITY);

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
      // Check for empty lines
      if( strcmp(buf,"\r") == 0 ||
	  strcmp(buf,"\n") == 0 ||
	  strcmp(buf,"\r\n") == 0 ) break;

      if( i > pars->n_ind || split(buf, (const char*) " ,-\t\r\n", &t) != 2)
	error("wrong TRANS file format!");
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
      error("wrong TRANS parameters format!");

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



  ////////////////////////////////////////
  // Set FREQ initial values...
  ////////////////////////////////////////
  double freq_rng_max = 0.49;
  double freq_rng_min = 0.01;
  gzFile in_freq_fh;

  data->freq = init_double(pars->n_sites+1, freq_rng_min);
  // Initialize site 0 to invalid value
  data->freq[0] = -1;

  if( strcmp("r", pars->in_freq) == 0 )
    for(uint64_t s = 1; s <= pars->n_sites; s++)
      data->freq[s] = freq_rng_min + gsl_rng_uniform(r) * (freq_rng_max - freq_rng_min);

  else if( (in_freq_fh = gzopen(pars->in_freq, "r")) != NULL ){
    uint64_t s = 1;
    while( gzgets(in_freq_fh, buf, BUFF_LEN) != NULL ){
      // Check for empty lines
      if( strcmp(buf,"\r") == 0 ||
	  strcmp(buf,"\n") == 0 ||
	  strcmp(buf,"\r\n") == 0 ) break;

      if( s > pars->n_sites || split(buf, (const char*) " ,-\t\r\n", &t) != 1)
        error("wrong FREQ file format!");
      data->freq[s] = min(max(t[0], freq_rng_min), freq_rng_max);
      s++;
      delete [] t;
    } 
    gzclose(in_freq_fh);
  }

  else
    for(uint64_t s = 1; s <= pars->n_sites; s++)
      data->freq[s] = min(max(atof(pars->in_freq), freq_rng_min), freq_rng_max);



  ////////////////////////////////////////
  // Set EMISSION initial values ...
  ////////////////////////////////////////
  data->e = init_double(pars->n_sites+1, N_STATES, N_GENO, -INFINITY);
  // Update emission probs based on allele freqs
  update_e(data, pars->n_sites);
  


  ////////////////////////////////////////
  // Set PATH initial values ...
  ////////////////////////////////////////
  gzFile in_path_fh;

  data->path = init_uint64(pars->n_ind, pars->n_sites+1, 0);

  if( strcmp("r", pars->in_path) == 0 )
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++)
	data->path[i][s] = gsl_rng_uniform(r) > 0.5 ? 1 : 0;

  else if( (in_path_fh = gzopen(pars->in_path, "r")) != NULL ){
    uint64_t s = 1;
    while( gzgets(in_path_fh, buf, BUFF_LEN) != NULL ){
      // Check for empty lines
      if( strcmp(buf,"\r") == 0 ||
	  strcmp(buf,"\n") == 0 ||
	  strcmp(buf,"\r\n") == 0 ) break;

      int* t = NULL;
      if( s > pars->n_sites || split(buf, (const char*) " ,-\t\r\n", &t) != pars->n_ind )
        error("wrong PATH file format!");
      for(uint64_t i = 0; i < pars->n_ind; i++)
	data->path[i][s] = t[i] > 0.5 ? 1 : 0;
      s++;
      delete [] t;
    }
    gzclose(in_path_fh);
  }

  else
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++)
	data->path[i][s] = atoi(pars->in_path) > 0.5 ? 1 : 0;



  ////////////////////////////////////////
  // Initialize Marginal Probabilities
  ////////////////////////////////////////
  data->marg_prob = init_double(pars->n_ind, pars->n_sites+1, N_STATES, 0);
  // Initialize site 0 to invalid value
  for (uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t k = 0; k < N_STATES; k++)
      data->marg_prob[i][0][k] = -1;



  ////////////////////////////////////////
  // Initialize indF and lkl variables
  ////////////////////////////////////////
  data->indF = init_double(pars->n_ind, 0);
  data->lkl = init_double(pars->n_ind, -INFINITY);



  delete [] buf;
  gsl_rng_free(r);
  return(0);
}



int read_geno(params* pars){
  uint64_t n_fields;
  // Depending on input we will have either 1 or 3 genot
  uint64_t n_geno = (pars->in_lkl ? N_GENO : 1);
  double* t;
  double* ptr;
  char* buf = new char[BUFF_LEN];

  // Allocate memory
  pars->geno_lkl = init_double(pars->n_ind, pars->n_sites+1, N_GENO, -INFINITY);
  
  // Open GENO file
  gzFile in_geno_fh;
  if( (in_geno_fh = gzopen(pars->in_geno, pars->in_bin ? "rb" : "r")) == NULL )
    error("cannot open genotype file!");

  for(uint64_t s = 1; s <= pars->n_sites; s++){
    if(pars->in_bin){
      for(uint64_t i = 0; i < pars->n_ind; i++)
	if( gzread(in_geno_fh, pars->geno_lkl[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
	  error("cannot read GENO file!");
    }
    else{
      if( gzgets(in_geno_fh, buf, BUFF_LEN) == NULL)
	error("cannot read GENO file!");
      
      // Parse input line into array
      n_fields = split(buf, (const char*) " \t\r\n", &t);

      // Check if header and skip
      if(!n_fields){
	s--;
	continue;
      }

      if(n_fields < pars->n_ind * n_geno)
	error("wrong GENO file format!");
      
      // Use last "n_ind * n_geno" columns
      ptr = t + (n_fields - pars->n_ind * n_geno);
      
      if(pars->in_lkl)
	for(uint64_t i = 0; i < pars->n_ind; i++)
	  for(uint64_t g = 0; g < N_GENO; g++)
	    pars->geno_lkl[i][s][g] = pars->in_loglkl ? ptr[i*N_GENO+g] : log(ptr[i*N_GENO+g]);
      else
	for(uint64_t i = 0; i < pars->n_ind; i++){
	  int g = (int) ptr[i];
	  pars->geno_lkl[i][s][g] = log(1);
	}

      delete [] t;
    }
  }
  
  gzclose(in_geno_fh);
  delete [] buf;
  return 0;
}




uint64_t read_chunk(double** chunk_data, params* pars, uint64_t chunk) {
  uint64_t total_elems_read = 0;
  
  if(chunk >= pars->n_chunks)
    error("invalid chunk number!");
  
  // Define chunk start and end positions
  uint64_t start_pos = chunk * pars->max_chunk_size;
  uint64_t end_pos = start_pos + pars->max_chunk_size - 1;
  if(end_pos >= pars->n_sites)	end_pos = pars->n_sites - 1;
  uint64_t chunk_size = end_pos - start_pos + 1;
  if( pars->verbose >= 6 ) printf("\tReading chunk %lu from position %lu to %lu (%lu)\n", chunk+1, start_pos, end_pos, chunk_size);
  
  // Read data from file
  for(uint64_t c = 0; c < chunk_size; c++) {
    //    chunk_data[c] = pars->geno_lkl[start_pos+c];
    uint64_t elems_read = pars->n_ind * 3;

    if( elems_read != pars->n_ind * 3 )
      error("cannot read GLF file!");
    total_elems_read += elems_read;
  }

  return( total_elems_read/(pars->n_ind * 3) );
}
