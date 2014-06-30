
#include "ngsF-HMM.hpp"

// General structure for launching threads
struct pth_struct{
  int type;
  double **new_a;
  double **data;
  double **Fw;
  double **Bw;
  double **Vi;
  double ***e;
  double **a;
  unsigned short int *path;
  uint64_t length;
};



// Function prototypes
void threadpool_add_task(int type, double **new_a, double **data, double **Fw, double **Bw, double **Vi, double ***e, double **a, unsigned short int *path, uint64_t length);
void thread_slave(void *ptr);
void forward(double **data, double **Fw, double ***e, double **a, unsigned short int *path, uint64_t length);
void backward(double **data, double **Bw, double ***e, double **a, unsigned short int *path, uint64_t length);
void viterbi(double **data, double **Vi, double ***e, double **a, unsigned short int *path, uint64_t length);
void trans(double **new_a, double **data, double **Fw, double **Bw, double ***e, double **a, unsigned short int *path, uint64_t length);



// General thread functions
void threadpool_add_task(threadpool_t *thread_pool, int type, double **new_a, double **data, double **Fw, double **Bw, double **Vi, double ***e, double **a, unsigned short int *path, uint64_t length){
  pth_struct *p = new pth_struct;

  p->type = type;
  p->new_a = new_a;
  p->data = data;
  p->Fw = Fw;
  p->Bw = Bw;
  p->Vi = Vi;
  p->e = e;
  p->a = a;
  p->path = path;
  p->length = length;

  // Add task to thread pool
  //thread_slave((void*) p);

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
    forward(p->data, p->Fw, p->e, p->a, p->path, p->length);
  if(p->type == 2)
    backward(p->data, p->Bw, p->e, p->a, p->path, p->length);
  if(p->type == 3)
    viterbi(p->data, p->Vi, p->e, p->a, p->path, p->length);
  if(p->type == 4)
    trans(p->new_a, p->data, p->Fw, p->Bw, p->e, p->a, p->path, p->length);
}



// Forward functions
void forward(double **data, double **Fw, double ***e, double **a, unsigned short int *path, uint64_t length){
  double pp[N_GENO];

  for (uint64_t s = 1; s <= length; s++){
    post_prob(pp, data[s], e[s][path[s]], N_GENO);

    for(uint64_t l = 0; l < N_STATES; l++){
      // logsum(k==0,k==1)
      Fw[s][l] = logsum2(Fw[s-1][0] + a[0][l],
			 Fw[s-1][1] + a[1][l]);
      Fw[s][l] += logsum3(e[s][l][0]+pp[0], e[s][l][1]+pp[1], e[s][l][2]+pp[2]);
    }
  }
}


void backward(double **data, double **Bw, double ***e, double **a, unsigned short int *path, uint64_t length){
  double pp[N_GENO];

  for (uint64_t s = length; s > 0; s--){
    post_prob(pp, data[s], e[s][path[s]], N_GENO);

    double LS_0 = logsum3(e[s][0][0]+pp[0], e[s][0][1]+pp[1], e[s][0][2]+pp[2]);
    double LS_1 = logsum3(e[s][1][0]+pp[0], e[s][1][1]+pp[1], e[s][1][2]+pp[2]);

    for(uint64_t k = 0; k < N_STATES; k++)
      // logsum(l==0,l==1)
      Bw[s-1][k] = logsum2(a[k][0] + LS_0 + Bw[s][0],
			   a[k][1] + LS_1 + Bw[s][1]);
  }
}



void viterbi(double **data, double **Vi, double ***e, double **a, unsigned short int *path, uint64_t length){
  double pp[N_GENO];

  for (uint64_t s = 1; s <= length; s++){
    post_prob(pp, data[s], e[s][path[s]], N_GENO);

    for(uint64_t l = 0; l < N_STATES; l++){
      // max(k==0,k==1)
      Vi[s][l] = max(Vi[s-1][0] + a[0][l],
		     Vi[s-1][1] + a[1][l]);
      Vi[s][l] += logsum3(e[s][l][0]+pp[0], e[s][l][1]+pp[1], e[s][l][2]+pp[2]);
    }
  }
}


void trans(double **new_a, double **data, double **Fw, double **Bw, double ***e, double **a, unsigned short int *path, uint64_t length){
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
}


/*
// UNTESTED!!!!!
void emission(double ***new_e, double **data, double **Fw, double **Bw, double ***e, double **a, unsigned short int *path, uint64_t length){
  double pp[N_GENO];
  double sPk[length+1];

  for(uint64_t k = 0; k < N_STATES; k++){
    // Get P(k)
    for (uint64_t s = 1; s <= length; s++)
      sPk[s] = Fw[s][k] + Bw[s][k];
    // Sum all site_Pk, skipping site 0 and last
    double Pk = logsum(sPk+1, length-1);

    for (uint64_t s = 1; s <= length; s++){
      post_prob(pp, data[s], e[s][path[s]], N_GENO);

      for(uint64_t g = 0; g < N_GENO; g++)
	Fw[s][k]+Bw[s][k]+pp[g]-Pk;
	new_e[s][k][g] = logsum(new_e[s][k][g], );
    }
  }

  for (uint64_t s = 1; s <= length; s++)
    new_e[s] /= length;
}
*/
