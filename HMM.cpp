
#include "ngsF-HMM.hpp"


// General structure for launching threads
struct pth_struct{
  int type;
  double **ptr;
  double *F;
  double *alpha;
  double **e_prob;
  char *path;
  double *pos_dist;
  uint64_t length;
};

// Function prototypes
void threadpool_add_task(threadpool_t *thread_pool, int type, double **ptr, double *F, double *alpha, double **e_prob, char *path, double *pos_dist, uint64_t length);
void thread_slave(void *ptr);


// General thread function
void threadpool_add_task(threadpool_t *thread_pool, int type, double **ptr, double *F, double *alpha, double **e_prob, char *path, double *pos_dist, uint64_t length){
  pth_struct *p = new pth_struct;

  p->type = type;
  p->ptr = ptr;
  p->F = F;
  p->alpha = alpha;
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

  if(p->type == 1)
    forward(p->ptr, *p->F, *p->alpha, p->e_prob, p->pos_dist, p->length);
  else if(p->type == 2)
    backward(p->ptr, *p->F, *p->alpha, p->e_prob, p->pos_dist, p->length);
  else if(p->type == 3)
    viterbi(p->ptr, *p->F, *p->alpha, p->e_prob, p->path, p->pos_dist, p->length);
  else if(p->type == 4){
    double val[2] = {*p->F, *p->alpha};
    double l_bound[2] = {1/INF, 1/INF};
    double u_bound[2] = {1-l_bound[0], 10};
    int lims[2] = {2, 2};

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
  else
    lkl = forward(Fw, pars[0], pars[1], p->e_prob, p->pos_dist, p->length);

  free_ptr((void**) Fw, p->length+1);
  return -lkl;
}


// Forward function
double forward(double **Fw, double F, double alpha, double **e_prob, double *pos_dist, uint64_t length){
  Fw[0][0] = log(1-F);
  Fw[0][1] = log(F);

  for (uint64_t s = 1; s <= length; s++)
    for(uint64_t l = 0; l < N_STATES; l++){
      // logsum(k==0,k==1)
      Fw[s][l] = logsum(Fw[s-1][0] + calc_trans(0,l,pos_dist[s],F,alpha),
			Fw[s-1][1] + calc_trans(1,l,pos_dist[s],F,alpha)) + e_prob[s][l];

      if(isnan(Fw[s][l])){
	printf("site: %lu\tdist: %f\tF: %.15f %.15f\tstate: %lu\tFw: %f %f %f\ttrans: %f %f\temission: %f\n", s, pos_dist[s], F, alpha, l, Fw[s-1][0], Fw[s-1][1], Fw[s][l], calc_trans(0,l,pos_dist[s],F,alpha), calc_trans(1,l,pos_dist[s],F,alpha), e_prob[s][l]);
	error(__FUNCTION__, "invalid Lkl found!");
      }
    }

  return logsum(Fw[length],2);
}



// Backward function
double backward(double **Bw, double F, double alpha, double **e_prob, double *pos_dist, uint64_t length){
  Bw[length][0] = log(1);
  Bw[length][1] = log(1);

  for (uint64_t s = length; s > 0; s--){
    for(uint64_t k = 0; k < N_STATES; k++){
      // logsum(l==0,l==1)
      Bw[s-1][k] = logsum(calc_trans(k,0,pos_dist[s],F,alpha) + e_prob[s][0] + Bw[s][0],
			  calc_trans(k,1,pos_dist[s],F,alpha) + e_prob[s][1] + Bw[s][1]);

      if(isnan(Bw[s-1][k])){
	printf("site: %lu\tdist: %f\tF: %.15f %.15f\tstate: %lu\tBw: %f %f %f\ttrans: %f %f\temission: %f %f\n", s, pos_dist[s], F, alpha, k, Bw[s][0], Bw[s][1], Bw[s-1][k], calc_trans(k,0,pos_dist[s],F,alpha), calc_trans(k,1,pos_dist[s],F,alpha), e_prob[s][0], e_prob[s][1]);
	error(__FUNCTION__, "invalid Lkl found!");
      }
    }
  }

  Bw[0][0] += log(1-F);
  Bw[0][1] += log(F);

  return logsum(Bw[0],2);
}



double viterbi(double **Vi, double F, double alpha, double **e_prob, char *path, double *pos_dist, uint64_t length){
  double Vi_prob[N_STATES] = {log(1-F), log(F)};

  for (uint64_t s = 1; s <= length; s++){
    for(uint64_t l = 0; l < N_STATES; l++){
      double vmax = -INF;
      uint64_t k_vmax = 0;

      for(uint64_t k = 0; k < N_STATES; k++){
	double pval = Vi_prob[k] + calc_trans(k,l,pos_dist[s],F,alpha);
	if( vmax < pval ) { vmax = pval; k_vmax = k; }
      }

      Vi[s][l] = k_vmax;
      Vi_prob[l] = vmax + e_prob[s][l];
    }
  }

  // Trace back the Viterbi path
  path[length] = array_max_pos(Vi_prob, N_STATES);
  for (int s = length-1; s >= 0; s--)
    path[s] = Vi[s+1][(int) path[s+1]];

  return max(Vi_prob[0], Vi_prob[1]);
}



// Calculates transition probabilities between states k and l, depending on distance, inbreeding and transition rate
double calc_trans(char k, char l, double pos_dist, double F, double alpha, bool cont){
  double trans = 0;

  if(cont){
    double coanc_change = exp(-alpha*pos_dist);
    // Continuous HMM
    if(k == 0 && l == 0)
      trans = (1-coanc_change) * (1-F) + coanc_change;
    else if(k == 0 && l == 1)
      trans = (1-coanc_change) * F;
    else if(k == 1 && l == 0)
      trans = (1-coanc_change) * (1-F);
    else if(k == 1 && l == 1)
      trans = (1-coanc_change) * F + coanc_change;
  } else {
    // HMM
    if(k == 0 && l == 0)
      trans = 1-alpha*F;
    else if(k == 0 && l == 1)
      trans = alpha*F;
    else if(k == 1 && l == 0)
      trans = alpha*(1-F);
    else if(k == 1 && l == 1)
      trans = 1-alpha*(1-F);
  }

  return log(trans);
}



// Calculates emission probabilities for "state" depending on the MAF. if state==0, assumes HWE, if not, assumes HWE with inbreeding.
double calc_emission(double gl[3], double maf, uint64_t state){
  if(maf < 0 || maf > 1)
    error(__FUNCTION__, "invalid MAF!");

  double emission[3];
  calc_HWE(emission, maf, state);

  return logsum(gl[0]+emission[0], gl[1]+emission[1], gl[2]+emission[2]);
}



// Calculates emission probabilities for "state" depending on the MAF. if state==0, assumes HWE, if not, assumes HWE with inbreeding.
// site _p is previous site
// site _c is current site
double calc_emissionLD(double hap_freq[4], double *gl_p, double *gl_c, double maf_p, double maf_c,  uint64_t state){
  if(maf_p < 0 || maf_p > 1 || maf_c < 0 || maf_c > 1)
    error(__FUNCTION__, "invalid MAF!");

  double sum, P_Gc, P_Gp[N_GENO];
  double s_p[N_GENO], s_c[N_GENO];
  for(uint64_t g = 0; g < N_GENO; g++){
    s_p[g] = exp(gl_p[g]);
    s_c[g] = exp(gl_c[g]);
  }

  sum = 0;
  if(0){
    for(uint64_t g_c = 0; g_c < N_GENO; g_c++){
      calc_HWE(P_Gp, maf_p, state);
      for(uint64_t g_p = 0; g_p < N_GENO; g_p++)
	P_Gc += joint_geno_prob(hap_freq, g_p, state, g_c, state) * P_Gp[g_p];

      sum += s_c[g_c] * P_Gc;
    }

    return log(sum);
  }else{
    // results are calculated with current algorithm
    // results _2 are calculated just taking the haplotype frequency into account (removing s_p[g_p] * s_c[g_c])

    for(uint64_t g_c = 0; g_c < N_GENO; g_c++)
      for(uint64_t g_p = 0; g_p < N_GENO; g_p++)
	sum += joint_geno_prob(hap_freq, g_p, state, g_c, state) * s_p[g_p] * s_c[g_c];

    return log(sum) - calc_emission(gl_p, maf_p, state);
  }
}



double joint_geno_prob(double hap_freq[4], uint64_t g_p, uint64_t state_p, uint64_t g_c, uint64_t state_c){
  if(state_p != state_c)
    error(__FUNCTION__, "so far only cases where prev and curr positions have same state!");

  if(g_p == 0 && g_c == 0){
    return (state_c == 0 ? pow(hap_freq[0], 2) : hap_freq[0]);
  }else if(g_p == 0 && g_c == 1){
    return (state_c == 0 ? 2*hap_freq[0]*hap_freq[1] : 0);
  }else if(g_p == 0 && g_c == 2){
    return (state_c == 0 ? pow(hap_freq[1], 2): hap_freq[1]);
  }else if(g_p == 1&& g_c == 0){
    return (state_c == 0 ? 2*hap_freq[0]*hap_freq[2] : 0);
  }else if(g_p == 1&& g_c == 1){
    return (state_c == 0 ? 2*(hap_freq[0]*hap_freq[3]+hap_freq[1]*hap_freq[2]) : 0);
  }else if(g_p == 1 && g_c == 2){
    return (state_c == 0 ? 2*hap_freq[1]*hap_freq[3] : 0);
  }else if(g_p == 2 && g_c == 0){
    return (state_c == 0 ? pow(hap_freq[2], 2) : hap_freq[2]);
  }else if(g_p == 2 && g_c == 1){
    return (state_c == 0 ? 2*hap_freq[2]*hap_freq[3] : 0);
  }else if(g_p == 2 && g_c == 2){
    return (state_c == 0 ? pow(hap_freq[3], 2) : hap_freq[3]);
  }

  return -1;
}
