#include "sort.h"

double* OmpRadixSortMSDStack(stack<double> st, uint radix)
{
	double* result;
	stack<double> stack[NUM_VAL];
	uint thr1 = 0;
	uint thr0 = 0;
	if (st.empty())
	{
		return NULL;
	}
	else
	{
		result = OmpGetMemoryPool()->OmpAlloc(st.size());
		if ((64 == radix) || (st.size() == 1))
		{
			uint counter = 0;
			while (!st.empty())
			{
				result[counter] = st.top();
				st.pop();
				counter++;
			}
			return result;
		}
	}

	while (!st.empty())
	{
		double value = st.top();
		st.pop();
		(OmpGetBit(value, radix)) == 0 ? thr0++ : thr1++;
		stack[OmpGetBit(value, radix)].push(value);
	}
	uint counter = 0;
	uint counter1 = 0;
	int i;
	double* res1;
	double* res2;
	//#pragma omp parallel shared(stack) num_threads(2)

		//#pragma omp sections
#pragma omp parallel shared(stack) num_threads(2)
	{
#pragma omp sections
	{
#pragma omp section
	{
		{
			i = 1;
			res1 = OmpRadixSortMSDStack(stack[i], radix + 1);
			if (NULL != res1)
			{
				for (int j = 0; j < (int)stack[i].size(); j++)
				{
					{
						result[counter] = res1[j];
						counter++;
					}
				}
			}
		}
	}
#pragma omp section
	{
		{
			int i1 = 0;
			res2 = OmpRadixSortMSDStack(stack[i1], radix + 1);
			if (NULL != res2)
			{
				for (int j = 0; j < (int)stack[i1].size(); j++)
				{
					{
						result[thr1 + counter1] = res2[j];
						counter1++;
					}
				}
			}
		}
	}
	}
	}
	OmpGetMemoryPool()->OmpFree(thr1, res1);
	OmpGetMemoryPool()->OmpFree(thr0, res2);
	return result;
}

double* OmpRadixSortMSD(const double* array, const uint len, uint radix)
{
	double* result;
	//#pragma omp critical
	{
		result = OmpGetMemoryPool()->OmpAlloc(len);
	}
	stack<double> stack[NUM_VAL];
	uint thr1 = 0;
	uint thr0 = 0;
	double val = 0;
	//create_lock(mutex);
	//init_lock(&mutex);
#pragma omp parallel shared(stack/*, mutex*/) private(val)
	{
#pragma omp for schedule(dynamic, CHUNK) nowait
		for (int i = 0; i < (int)len; ++i)
		{
			//lock(&mutex);
			{
#pragma omp critical
			{
				val = array[i];
				(OmpGetBit(val, radix)) == 0 ? thr0++ : thr1++;
				stack[OmpGetBit(val, radix)].push(val);
			}
			}
			//unlock(&mutex);
		}
	}
	//destroy_lock(&mutex);
	uint counter1 = 0;
	uint counter2 = 0;
	int j = 0;
	double* res1;
	double* res2;

	//#pragma omp parallel shared(stack) private(j) num_threads(2)
#pragma omp parallel shared(stack) private(j) num_threads(2)
	{
#pragma omp sections
	{
#pragma omp section
	{
		res1 = OmpRadixSortMSDStack(stack[0], radix + 1);
		if (NULL != res1)
		{
			//#pragma omp parallel shared(counter,  stack) private(j)
			{
				//#pragma omp for schedule(dynamic, CHUNK) nowait
				for (j = 0; j < (int)thr0; j++)
				{
					result[counter1] = res1[j];
					counter1++;
				}
			}
			//#pragma omp critical
			//OmpGetMemoryPool()->OmpFree(thr0, res1);
		}
	}
	//
#pragma omp section
	{
		res2 = OmpRadixSortMSDStack(stack[1], radix + 1);
		if (NULL != res2)
		{
			//#pragma omp parallel shared(counter,  stack) private(j)
			{
				//#pragma omp for schedule(dynamic, CHUNK) nowait
				{
					for (j = (int)thr1 - 1; j >= 0; --j)
					{
						result[thr0 + counter2] = res2[j];
						counter2++;
					}
				}
			}
			//#pragma omp critical
			//OmpGetMemoryPool()->OmpFree(thr1, res2);
		}
	}
	}
	}
	return result;
}

void OmpAddNewElemToAuxArr(double elem, double *arr, unsigned int bit, uint *left0, uint *right1)
{
	if (!bit)
	{
		arr[*left0] = elem;
		(*left0)++;
	}
	else
	{
		(*right1)--;
		arr[*right1] = elem;
	}
}

double* OmpMSDRadixSort(const double* array, const uint len, uint radix, uint full)
{
	double *result;
	double *aux;
	uint thr1 = 0;
	uint thr0 = 0;
	double val = 0;
	unsigned int bit = 0;
	uint left0 = 0;
	uint right1 = len;
	if (len == 0)
	{
		return NULL;
	}
	else
	{
		if ((64 == radix) || (len == 1))
		{
			result = OmpGetMemoryPool()->OmpAlloc(len);
			uint counter = 0;
			for (int i = 0; i < (int)len; ++i)
			{
				result[counter] = array[i];
				counter++;
			}
			return result;
		}
	}
	result = OmpGetMemoryPool()->OmpAlloc(len);
	aux = OmpGetMemoryPool()->OmpAlloc(len);
	for (int i = 0; i < (int)len; ++i)
	{
		val = array[i];
		bit = OmpGetBit(val, radix);
		OmpAddNewElemToAuxArr(val, aux, bit, &left0, &right1);
		(!bit) ? thr0++ : thr1++;
	}
	uint counter1 = 0;
	uint counter2 = 0;
	double* res1;
	double* res2;

	if (len > full / 2/*omp_get_num_threads()*/)
	{
#pragma omp parallel shared(aux) num_threads(2)
	{
#pragma omp sections
	{
#pragma omp section
	{
		res2 = OmpMSDRadixSort(aux + thr0, thr1, radix + 1, full);
		if (NULL != res2)
		{
			if (0 == radix)
			{
				for (int j = (int)thr1 - 1; j >= 0; --j)
				{
					result[counter2] = res2[j];
					counter2++;
				}
			}
			else
			{
				for (int j = 0; j < (int)thr1; j++)
				{
					result[thr0 + counter2] = res2[j];
					counter2++;
				}
			}
			OmpGetMemoryPool()->OmpFree(thr0, res2);
		}
	}
#pragma omp section
	{
		res1 = OmpMSDRadixSort(aux, thr0, radix + 1, full);
		if (NULL != res1)
		{
			if (0 == radix)
			{
				for (int j = 0; j < (int)thr0; j++)
				{
					result[thr1 + counter1] = res1[j];
					counter1++;
				}
			}
			else
			{
				for (int j = 0; j < (int)thr0; j++)
				{
					result[counter1] = res1[j];
					counter1++;
				}
			}
			OmpGetMemoryPool()->OmpFree(thr1, res1);
		}
	}
	}
	}
	}
	else
	{
		res2 = OmpMSDRadixSort(aux + thr0, thr1, radix + 1, full);
		if (NULL != res2)
		{
			if (0 == radix)
			{
				for (int j = (int)thr1 - 1; j >= 0; --j)
				{
					result[counter2] = res2[j];
					counter2++;
				}
			}
			else
			{
				for (int j = 0; j < (int)thr1; j++)
				{
					result[thr0 + counter2] = res2[j];
					counter2++;
				}
			}
			OmpGetMemoryPool()->OmpFree(thr0, res2);
		}
		res1 = OmpMSDRadixSort(aux, thr0, radix + 1, full);
		if (NULL != res1)
		{
			if (0 == radix)
			{
				for (int j = 0; j < (int)thr0; j++)
				{
					result[thr1 + counter1] = res1[j];
					counter1++;
				}
			}
			else
			{
				for (int j = 0; j < (int)thr0; j++)
				{
					result[counter1] = res1[j];
					counter1++;
				}
			}
			OmpGetMemoryPool()->OmpFree(thr1, res1);
		}
	}
	OmpGetMemoryPool()->OmpFree(len, aux);
	return result;
}

double* MPIGetTaskForSort(const double* array, uint len, uint &lenTask, int procRank, int procNum, uint degree)
{
	double val;
	double *Task = new double[len];
	uint *bit = new uint[degree];
	uint pos = 0;
	int j = degree;
	for (int i = 64 - degree; i < 64; ++i)
	{
		bit[degree - j] = OmpGetBit(procRank, i);
		//cout << bit[degree - j];
		j--;
	}
	//cout << endl;
	lenTask = 0;
	for (int i = 0; i < len; ++i)
	{
		j = 0;
		int b = OmpGetBit(array[i], j);
		while ((bit[j] == b) && (j != degree))
		{
			j++;
			b = OmpGetBit(array[i], j);
		}
		if (j == degree)
		{
			Task[lenTask++] = array[i];
		}
	}
	return Task;
}

double** MPIGetTasksForSort(const double* array, uint len, uint *lenTasks, int procRanks, int procNum, uint degree, uint count)
{
	double val;
	double **Tasks = new double*[count];
	uint **bit = new uint*[count];
	uint pos = 0;
	int j = degree;
	for (int k = 0; k < count; ++k)
	{
		Tasks[k] = new double[len];
		lenTasks[k] = 0;
	}
	for (int k = 0; k < procNum; k++)
	{
		procRanks = k;
		bit[k] = new uint[degree];
		j = degree;
		for (int i = 64 - degree; i < 64; ++i)
		{
			bit[k][degree - j] = OmpGetBit(procRanks, i);
			//cout << bit[degree - j];
			j--;
		}
	}
	//cout << endl;
	bool gotNewTask = false;
	for (int i = 0; i < len; ++i)
	{
		gotNewTask = false;
		for (int k = 0; ((k < count)); k++)
		{
			j = 0;
			int b = OmpGetBit(array[i], j);
			while ((bit[k][j] == b) && (j != degree))
			{
				j++;
				b = OmpGetBit(array[i], j);
			}
			if (j == degree)
			{
				gotNewTask = true;
				Tasks[k][lenTasks[k]++] = array[i];
			}
			if (gotNewTask)
				break;
		}
	}
	for (int i = 0; i < procNum; i++)
	{
		delete[] bit[i];
	}
	delete[] bit;
	return Tasks;
}