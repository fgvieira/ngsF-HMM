#include "shared.hpp"

// Reads both called genotypes (1 field per site and indiv), genotype lkls or genotype post probs (3 fields per site and indiv)
double*** read_geno(char *in_geno, bool in_bin, bool in_probs, uint64_t n_ind, uint64_t n_sites){
  uint64_t n_fields;
  // Depending on input we will have either 1 or 3 genot
  uint64_t n_geno = (in_probs ? N_GENO : 1);
  double *t;
  double *ptr;
  char *buf = init_char(BUFF_LEN, '\0');

  // Allocate memory
  double ***geno = init_double(n_ind, n_sites+1, N_GENO, -INFINITY);
  
  // Open GENO file
  gzFile in_geno_fh = gzopen(in_geno, in_bin ? "rb" : "r");
  if(in_geno_fh == NULL)
    error(__FUNCTION__, "cannot open genotype file!");

  for(uint64_t s = 1; s <= n_sites; s++){
    if(in_bin){
      for(uint64_t i = 0; i < n_ind; i++)
	if( gzread(in_geno_fh, geno[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
	  error(__FUNCTION__, "cannot read GENO file!");
    }
    else{
      if( gzgets(in_geno_fh, buf, BUFF_LEN) == NULL)
	error(__FUNCTION__, "cannot read GENO file!");
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;

      // Parse input line into array
      n_fields = split(buf, (const char*) " \t", &t);

      // Check if header and skip
      if(!n_fields){
	s--;
	continue;
      }

      if(n_fields < n_ind * n_geno)
	error(__FUNCTION__, "wrong GENO file format!");
      
      // Use last "n_ind * n_geno" columns
      ptr = t + (n_fields - n_ind * n_geno);
      
      if(in_probs)
	for(uint64_t i = 0; i < n_ind; i++)
          for(uint64_t g = 0; g < N_GENO; g++)
            geno[i][s][g] = ptr[i*N_GENO+g];
      else
	for(uint64_t i = 0; i < n_ind; i++){
          int g = (int) ptr[i];
          if(g < 0 || g > 2)
            error(__FUNCTION__, "wrong GENO format!");
          geno[i][s][g] = log(1);
        }

      delete [] t;
    }
  }
  
  gzclose(in_geno_fh);
  delete [] buf;
  return geno;
}


// normalize GL!?!!?
