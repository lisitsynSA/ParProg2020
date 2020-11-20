#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <ctime>
#include <unistd.h>

#define TABLE_SIZE 20


int main(int argc, char** argv)
{
	omp_set_nested(1);
	
	clock_t start_time, end_time;
	clock_t arr_time[TABLE_SIZE][TABLE_SIZE];

	std::ofstream fout("heatTable.txt");
	fout << TABLE_SIZE << '\n';

	for(int i = 1; i <= TABLE_SIZE; ++i)
		for(int j = 1; j <= TABLE_SIZE; ++j)
		{
			start_time = clock();
			#pragma omp parallel num_threads(i)
			{
				#pragma omp parallel num_threads(j)
				{
					usleep(10);
				}
			}
			end_time = clock();
      		arr_time[i - 1][j - 1] = end_time - start_time;
		}
	
	for(int i = 0; i < TABLE_SIZE; ++i)
	{
		for(int j = 0; j < TABLE_SIZE; ++j)
			fout << arr_time[i][j] << ' ';
		fout << '\n';
	}
	fout.close();
	return 0;
}
