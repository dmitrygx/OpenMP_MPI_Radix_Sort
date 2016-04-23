#ifndef SORT_H
#define SORT_H

#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <stack>
#include <cmath>
#include "numeric.h"
#include "util.h"

#ifndef CHUNK 
#define CHUNK 10
#endif

#define NUM_VAL 2

using namespace std;

double* OmpRadixSortMSD(const double* array, const uint len, uint radix);

double* OmpMSDRadixSort(const double* array, const uint len, uint radix, uint full);

double* MPIGetTaskForSort(const double* array, uint len, uint &lenTask, int procRank, int procNum, uint degree);

double** MPIGetTasksForSort(const double* array, 
							uint len, 
							uint **lenTasks, 
							int &countOfTasks, 
							int procNum);

#endif