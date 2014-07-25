#pragma once

#include "gen_func.hpp"


const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;


double*** read_geno(char*, bool, bool, uint64_t, uint64_t);
double* read_pos(char*, uint64_t, uint64_t);
