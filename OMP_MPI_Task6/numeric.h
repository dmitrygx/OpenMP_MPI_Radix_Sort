#ifndef NUMERIC_H
#define NUMERIC_H

#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <ctime>
#include "util.h"

using namespace std;

#define uint unsigned long long int

double* OmpNumRandomGenerate(double min, double max, unsigned long long int num);

double OmpRand(double min, double max);

#endif
