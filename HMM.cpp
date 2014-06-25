
#include "ngsF-HMM.hpp"

int forward(double **data, double **Fw, double start, double ***e, double **a, uint *path, uint64_t length){
  double pp[N_GENO];

  Fw[0][0] = log(1-start);
  Fw[0][1] = log(start);

  for (uint64_t s = 1; s <= length; s++){
    post_prob(pp, data[s], e[s][path[s]], N_GENO);

    for(uint64_t l = 0; l < N_STATES; l++){
      // logsum(k==0,k==1)
      Fw[s][l] = logsum2(Fw[s-1][0] + a[0][l],
			 Fw[s-1][1] + a[1][l]);
      Fw[s][l] += logsum3(e[s][l][0]+pp[0], e[s][l][1]+pp[1], e[s][l][2]+pp[2]);
    }
  }

  return 0;
}


int backward(double **data, double **Bw, double start, double ***e, double **a, uint *path, uint64_t length){
  double pp[N_GENO];

  // Initialise Backward table
  Bw[length][0] = log(1);
  Bw[length][1] = log(1);

  for (uint64_t s = length; s > 0; s--){
    post_prob(pp, data[s], e[s][path[s]], N_GENO);

    double LS_0 = logsum3(e[s][0][0]+pp[0], e[s][0][1]+pp[1], e[s][0][2]+pp[2]);
    double LS_1 = logsum3(e[s][1][0]+pp[0], e[s][1][1]+pp[1], e[s][1][2]+pp[2]);

    for(uint64_t k = 0; k < N_STATES; k++)
      // logsum(l==0,l==1)
      Bw[s-1][k] = logsum2(a[k][0] + LS_0 + Bw[s][0],
			   a[k][1] + LS_1 + Bw[s][1]);
  }

  return 0;
}



int viterbi(double **data, double **Vi, double start, double ***e, double **a, uint *path, uint64_t length){
  double pp[N_GENO];

  // Initialise Forward table
  Vi[0][0] = log(1-start);
  Vi[0][1] = log(start);

  for (uint64_t s = 1; s <= length; s++){
    post_prob(pp, data[s], e[s][path[s]], N_GENO);

    for(uint64_t l = 0; l < N_STATES; l++){
      // max(k==0,k==1)
      Vi[s][l] = max(Vi[s-1][0] + a[0][l],
		     Vi[s-1][1] + a[1][l]);
      Vi[s][l] += logsum3(e[s][l][0]+pp[0], e[s][l][1]+pp[1], e[s][l][2]+pp[2]);
    }
  }

  return 0;
}


int trans(double **new_a, double **data, double **Fw, double **Bw, double ***e, double **a, uint *path, uint64_t length){
  double pp[N_GENO];
  double sPk[length+1];
  
  for(uint64_t k = 0; k < N_STATES; k++){
    // Get P(k)                                                                                                                                                                            
    for (uint64_t s = 1; s <= length; s++)
      sPk[s] = Fw[s][k] + Bw[s][k];
    // Sum all site_Pk, skipping site 0 and last                                                                                                                                           
    double Pk = logsum(sPk+1, length-1);

    for(uint64_t l = 0; l < N_STATES; l++){
      double tmp_a = -INFINITY;
      for (uint64_t s = 1; s < length; s++){
	post_prob(pp, data[s+1], e[s+1][path[s+1]], N_GENO);

	double LS = logsum3(e[s+1][l][0]+pp[0],
			    e[s+1][l][1]+pp[1],
			    e[s+1][l][2]+pp[2]);
	tmp_a = logsum2(tmp_a,
			Fw[s][k] + a[k][l] + LS + Bw[s+1][l] - Pk);
      }
      new_a[k][l] = tmp_a;
    }
  }

  return 0;
}
