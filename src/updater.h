#include <ctime>
#include "timer.h"
#include "globals.h"
#include <cmath>
#include "mt.h"
#include "preparator.h"

void OutputAllEnergies(long long int step);
void OutputSumForce();
void OutputForces();

void update(long long int step, int* mt_len);
int change_conc(int* delta, int* mt_len);
void mt_length(long long int step, int* mt_len);