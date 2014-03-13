
#include <algorithm>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

void error(const char*);
void handler(int);
void catch_SIG();
double check_interv(double, bool);
int array_max_pos(double*, int);
void call_geno(double*, int);
double factorial(uint64_t);
double logsum(double*, uint64_t);
double logsum2(double, double);
double logsum3(double, double, double);
uint64_t split(char*, const char*, int**);
uint64_t split(char*, const char*, float**);
uint64_t split(char*, const char*, double**);
uint64_t split(char*, const char*, char***);
char* join(uint64_t*, uint64_t, const char*);
char* join(double*, uint64_t, const char*);
uint64_t* init_uint64(uint64_t, uint64_t);
uint64_t** init_uint64(uint64_t, uint64_t, uint64_t);
double* init_double(uint64_t, double);
double** init_double(uint64_t, uint64_t, double);
double*** init_double(uint64_t, uint64_t, uint64_t, double);
char* init_char(uint64_t A, const char* init);
char** init_char(uint64_t A, uint64_t B, const char* init);
void free_ptr(void*);
void free_ptr(void**, uint64_t);
void free_ptr(void***, uint64_t, uint64_t);
void cpy(void*, void*, uint64_t, uint64_t);
void cpy(void*, void*, uint64_t, uint64_t, uint64_t);
void cpy(void*, void*, uint64_t, uint64_t, uint64_t, uint64_t);
