#include "main.h"
#include <ctime>
#include "mpi.h"

using namespace std;

int main(int argc, char **argv)
{

	double min, max;
	uint num, precision;
	int procNum, procRank;
	MPI_Status status;

	(void)OmpParseArgs(argc, argv, min, max, num, precision);

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	//OmpInitMemoryPool(num * 64);

	omp_set_nested(true);

	uint degree = 0;
	int counter = 1;
	while (counter < procNum)
	{
		degree++;
		counter *= 2;
	}

	if (0 == procRank)
	{
		double* array = OmpNumRandomGenerate(min, max, num);
		double* result;
		if (NULL == array)
		{
			cout << "Random array hasn't been obtained" << endl;
			return -1;
		}

		OmpOutput(true, "omp_array.txt", "Source array:", array, num, precision);
		time_t start = clock();
		int radix = 0;
		double **ResultTasks;
		uint *ResultLen;
		uint *lenTasks;

		ResultTasks = new double*[counter];
		lenTasks = new uint[counter];
		ResultLen = new uint[counter];

		uint numOfSend = 0;
		uint first = 0;
		int procRanks = 0;
		double **Tasks = MPIGetTasksForSort(array, num, lenTasks, procRanks, procNum, degree, counter);
		for (int i = 0; i < counter; ++i)
		{
			if (i % procNum != procRank)
			{
				MPI_Send(&lenTasks[i], 1, MPI_INT, i % procNum, (counter - 1) - i, MPI_COMM_WORLD);
				MPI_Send(Tasks[i], (int)lenTasks[i], MPI_DOUBLE, i % procNum, (counter - 1) - i, MPI_COMM_WORLD);
				numOfSend++;
				delete[] Tasks[i];
			}
			else
			{
				ResultTasks[(counter - 1) - i] = OmpMSDRadixSort(Tasks[i], lenTasks[i], degree, num);
				ResultLen[(counter - 1) - i] = lenTasks[i];
				delete[] Tasks[i];
				/*cout << procRank << " - " << (counter - 1) - i << endl;*/
				//first++;
			}
		}
		delete[] Tasks;
		delete[] lenTasks;
		int ResultTasksLen = 0;
		int numOfRecv = (int)first;
		while (numOfSend != 0)
		{
			MPI_Recv(&ResultTasksLen, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			ResultLen[status.MPI_TAG] = ResultTasksLen;
			ResultTasks[status.MPI_TAG] = new double[ResultTasksLen];
			MPI_Recv(ResultTasks[status.MPI_TAG], ResultTasksLen, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			numOfRecv++;
			numOfSend--;
		}
		result = new double[num];
		int count = 0;
		for (int i = 0; i < counter / 2; ++i)
		{
			for (int j = (int)ResultLen[i] - 1; j >= 0; --j)
			{
				result[count++] = ResultTasks[i][j];
			}
			delete[] ResultTasks[i];
		}
		for (int i = counter - 1; i >= counter / 2; --i)
		{
			for (int j = 0; j < ResultLen[i]; ++j)
			{
				result[count++] = ResultTasks[i][j];
			}
			delete[] ResultTasks[i];
		}
		time_t end = clock();
		time_t diff = end - start;
		cout << "Time: " << ((double)diff) / CLOCKS_PER_SEC << " sec." << endl;
		delete[] ResultLen;
		delete[] ResultTasks;
		if (NULL == result)
		{
			cout << "Array hasn't been sorted" << endl;
			return -1;
		}
		OmpOutput(true, "omp_result.txt", "Result of array sorting:", result, num, precision);
		OmpGetMemoryPool()->OmpFree(num, result);
	}
	else
	{
		uint numOfTasks = counter / procNum;
		if (counter > numOfTasks * procNum + procRank)
		{
			numOfTasks++;
		}
		//cout << "Proc " << procRank << " has numOfTasks = " << numOfTasks << endl;
		int *lenTask = new int[numOfTasks];
		double **Task = new double*[numOfTasks];
		double **ResultTask = new double*[numOfTasks];
		uint *tags = new uint[numOfTasks];
		uint currentTask = 0;
		while (currentTask != numOfTasks)
		{
			MPI_Recv(&lenTask[currentTask], 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			Task[currentTask] = new double[lenTask[currentTask]];
			MPI_Recv(Task[currentTask], lenTask[currentTask], MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			/*cout << "Proc " << procRank << " got lenTask = " << lenTask[currentTask] << endl;
			for (int i = 0; i < lenTask[currentTask]; ++i)
			{
				cout << procRank << ": " << Task[currentTask][i] << "; ";
			}*/
			tags[currentTask] = status.MPI_TAG;
			//cout << procRank << " - " << status.MPI_TAG << endl;
			currentTask++;
		}
		for (int i = 0; i < numOfTasks; ++i)
		{
			ResultTask[i] = OmpMSDRadixSort(Task[i], lenTask[i], degree, num);
			MPI_Send(&lenTask[i], 1, MPI_INT, 0, (int)tags[i], MPI_COMM_WORLD);
			MPI_Send(ResultTask[i], (int)lenTask[i], MPI_DOUBLE, 0, (int)tags[i], MPI_COMM_WORLD);
			delete[] ResultTask[i];
			delete[] Task[i];
			//char buf[100];
			//sprintf_s(buf, "omp_result_for_%d.txt", procRank);
			//OmpOutput(true, buf, "Result of array sorting:", ResultTask[i], lenTask[i], precision);
		}
		delete[] lenTask;
		delete[] Task;
		delete[] ResultTask;
		delete[] tags;
	}

	//OmpTerminateMemoryPool();
	MPI_Finalize();

	return 0;
}
