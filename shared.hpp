
#include <algorithm>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <signal.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>

#include "threadpool.h"

using namespace std;

void error(const char*, const char*);
void handler(int);
void catch_SIG();
double check_interv(double, bool);
int array_max_pos(double*, int);
double draw_rnd(gsl_rng*, uint64_t, uint64_t);
void call_geno(double*, int, bool);
void conv_space(double*, int, double (*func)(double));
double logsum(double*, uint64_t);
double logsum2(double, double);
double logsum3(double, double, double);
void chomp(char*);
int64_t read_file(char*, char***, uint64_t);
uint64_t split(char*, const char*, int**);
uint64_t split(char*, const char*, float**);
uint64_t split(char*, const char*, double**);
uint64_t split(char*, const char*, char***);
char *join(unsigned short int*, uint64_t, const char*);
char *join(uint64_t*, uint64_t, const char*);
char *join(double*, uint64_t, const char*);
unsigned short int *init_usint(uint64_t, unsigned short int);
unsigned short int **init_usint(uint64_t, uint64_t, unsigned short int);
uint64_t *init_uint64(uint64_t, uint64_t);
uint64_t **init_uint64(uint64_t, uint64_t, uint64_t);
double *init_double(uint64_t, double);
double **init_double(uint64_t, uint64_t, double);
double ***init_double(uint64_t, uint64_t, uint64_t, double);
char *strdcat(char*, const char*);
char *init_char(uint64_t, const char*);
char **init_char(uint64_t, uint64_t, const char*);
void free_ptr(void*);
void free_ptr(void**, uint64_t);
void free_ptr(void***, uint64_t, uint64_t);
void cpy(void*, void*, uint64_t, uint64_t);
void cpy(void*, void*, uint64_t, uint64_t, uint64_t);
void cpy(void*, void*, uint64_t, uint64_t, uint64_t, uint64_t);
