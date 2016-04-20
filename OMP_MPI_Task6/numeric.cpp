#include "numeric.h"

/**
* Allocate array of double numver and initialize them by
* random value
* @param min The minimum value of number
* @param max The maximum value of number
* @param num The numbers of double number
* @return The pointer to array of double number
*/
double* OmpNumRandomGenerate(double min, double max, uint num)
{
	if (0 >= num)
	{
		return NULL;
	}
	srand((unsigned int)time(NULL));
	double* result = OmpGetMemoryPool()->OmpAlloc(num);
	for (int i = 0; i < (int)num; ++i)
	{
		result[i] = OmpRand(min, max);
	}
	return result;
}

/**
* Generate random number from min value to max value
* @param min The minimum value of number
* @param max The maximum value of number
* @return The double number
*/
double OmpRand(double min, double max)
{
	double f = 0;
	{
		f = (double)rand() / RAND_MAX;
	}
	return min + f * (max - min);
}