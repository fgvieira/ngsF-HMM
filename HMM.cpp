
#include "ngsF-HMM.hpp"

// General structure for launching threads
struct pth_struct{
  int type;
  double **new_a;
  double **data;
  double **Fw;
  double **Bw;
  double **Vi;
  double ***prior;
  double **a;
  char *path;
  uint64_t length;
};



// Function prototypes
void threadpool_add_task(int type, double **new_a, double **data, double **Fw, double **Bw, double **Vi, double ***prior, double **a, char *path, uint64_t length);
void thread_slave(void *ptr);
void forward(double **data, double **Fw, double ***prior, double **a, char *path, uint64_t length);
void backward(double **data, double **Bw, double ***prior, double **a, char *path, uint64_t length);
void viterbi(double **data, double **Vi, double ***prior, double **a, char *path, uint64_t length);
void trans(double **new_a, double **data, double **Fw, double **Bw, double ***prior, double **a, char *path, uint64_t length);



// General thread functions
void threadpool_add_task(threadpool_t *thread_pool, int type, double **new_a, double **data, double **Fw, double **Bw, double **Vi, double ***prior, double **a, char *path, uint64_t length){
  pth_struct *p = new pth_struct;

  p->type = type;
  p->new_a = new_a;
  p->data = data;
  p->Fw = Fw;
  p->Bw = Bw;
  p->Vi = Vi;
  p->prior = prior;
  p->a = a;
  p->path = path;
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
    forward(p->data, p->Fw, p->prior, p->a, p->path, p->length);
  if(p->type == 2)
    backward(p->data, p->Bw, p->prior, p->a, p->path, p->length);
  if(p->type == 3)
    viterbi(p->data, p->Vi, p->prior, p->a, p->path, p->length);
  if(p->type == 4)
    trans(p->new_a, p->data, p->Fw, p->Bw, p->prior, p->a, p->path, p->length);
}



// Forward functions
void forward(double **data, double **Fw, double ***prior, double **a, char *path, uint64_t length){
  for (uint64_t s = 1; s <= length; s++)
    for(uint64_t l = 0; l < N_STATES; l++){
      double e_l = logsum3(data[s][0]+prior[s][l][0], data[s][1]+prior[s][l][1], data[s][2]+prior[s][l][2]);
      // logsum(k==0,k==1)
      Fw[s][l] = logsum2(Fw[s-1][0] + a[0][l], Fw[s-1][1] + a[1][l]) + e_l;
    }
}


void backward(double **data, double **Bw, double ***prior, double **a, char *path, uint64_t length){
  for (uint64_t s = length; s > 0; s--){
    double e_nIBD = logsum3(data[s][0]+prior[s][0][0], data[s][1]+prior[s][0][1], data[s][2]+prior[s][0][2]);
    double e_IBD  = logsum3(data[s][0]+prior[s][1][0], data[s][1]+prior[s][1][1], data[s][2]+prior[s][1][2]);

    for(uint64_t k = 0; k < N_STATES; k++)
      // logsum(l==0,l==1)
      Bw[s-1][k] = logsum2(a[k][0] + e_nIBD + Bw[s][0],
			   a[k][1] + e_IBD  + Bw[s][1]);
  }
}



void viterbi(double **data, double **Vi, double ***prior, double **a, char *path, uint64_t length){
  for (uint64_t s = 1; s <= length; s++)
    for(uint64_t l = 0; l < N_STATES; l++){
      double e_l = logsum3(data[s][0]+prior[s][l][0], data[s][1]+prior[s][l][1], data[s][2]+prior[s][l][2]);
      // max(k==0,k==1)
      Vi[s][l] = max(Vi[s-1][0] + a[0][l], Vi[s-1][1] + a[1][l]) + e_l;
    }
}


void trans(double **new_a, double **data, double **Fw, double **Bw, double ***prior, double **a, char *path, uint64_t length){
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
	double e_l = logsum3(data[s+1][0]+prior[s+1][l][0], data[s+1][1]+prior[s+1][l][1], data[s+1][2]+prior[s+1][l][2]);
	tmp_a = logsum2(tmp_a,
			Fw[s][k] + a[k][l] + e_l + Bw[s+1][l] - Pk);
      }
      new_a[k][l] = tmp_a;
    }
  }
}
