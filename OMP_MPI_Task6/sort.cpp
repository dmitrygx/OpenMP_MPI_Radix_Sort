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

double** MPIGetTasksForSort(const double* array, uint len, uint **lenTasks, int &countOfTasks, int procNum)
{
	uint radix = 0;
	double **Tasks = new double*[2 * (radix + 1)];
	double **NewTasks;
	uint *newLenTasks;
	(*lenTasks) = new uint[2 * (radix + 1)];
	double val;
	int bit;
	for (int i = 0; i < 2 * (radix + 1); ++i)
	{
		Tasks[i] = new double[len];
		(*lenTasks)[i] = 0;
	}
	for (int i = 0; i < len; ++i)
	{
		val = array[i];
		bit = OmpGetBit(val, radix);
		Tasks[bit][(*lenTasks)[bit]++] = val;
	}
	for (int i = 0; i < 2 * (radix + 1); ++i)
	{
		if (((*lenTasks)[i]) > 0)
		{
			countOfTasks++;
		}
	}
	radix++;
	
	while ((countOfTasks < procNum))
	{
		cout << "Radix = " << radix << " CT" << countOfTasks << endl;
		NewTasks = new double*[(int)pow(2, (int)radix + 1)];
		newLenTasks = new uint[(int)pow(2, (int)radix + 1)];
		for (int i = 0; i < pow(2, (int)radix + 1); ++i)
		{
			NewTasks[i] = new double[len];
			newLenTasks[i] = 0;
		}
		for (int i = 0; i < (int)pow(2, (int)radix); ++i)
		{
			for (int j = 0; j < (*lenTasks)[i]; ++j)
			{
				val = Tasks[i][j];
				bit = OmpGetBit(val, radix);
				NewTasks[i * 2 + bit][newLenTasks[i * 2 + bit]++] = val;
			}
			delete[] Tasks[i];
		}
		delete[] Tasks;
		delete[](*lenTasks);
		countOfTasks = 0;
		Tasks = new double*[(int)pow(2, (int)radix + 1)];
		(*lenTasks) = new uint[(int)pow(2, (int)radix + 1)];
		for (int i = 0; i < (int)pow(2, (int)radix + 1); ++i)
		{
			if (newLenTasks[i] > 0)
			{
				countOfTasks++;
				cout << "TASK+1 => COT = " << countOfTasks << endl;
			}
			Tasks[i] = new double[newLenTasks[i]];
			(*lenTasks)[i] = newLenTasks[i];
			//cout << i << ": NEWLEN = " << newLenTasks[i] << "; LEN = " << (*lenTasks)[i] << endl;
			for (int j = 0; j < (*lenTasks)[i]; ++j)
			{
				Tasks[i][j] = NewTasks[i][j];
			}
			delete[] NewTasks[i];
		}
		delete[] newLenTasks;
		delete[] NewTasks;
		radix++;
		if (countOfTasks >= len)
		{
			break;
		}
	}
	cout << "COT = " << countOfTasks << endl;
	return Tasks;
}