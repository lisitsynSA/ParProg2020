#include <stdio.h>
#include <omp.h>

int create_thread(int n) 
{
	int res = 0;
	#pragma omp parallel num_threads(n) 
	{
		res++;
	}
	return res;
}

int main() 
{
	FILE * fp;
	fp = fopen("heatmap.csv","w");
	for (int x = 0; x <= 30; x++) {
		if (x != 30) {
			fprintf(fp, "%d,", x);
		} else {
			fprintf(fp, "%d\n", x);
		}
	}

	omp_set_nested(1);
	for (int y = 1; y <= 30; y++) {
		fprintf(fp, "%d,", y);
		for (int x = 1; x <= 30; x++) {
			double start = omp_get_wtime(); 
			int res = 0;
			#pragma omp parallel for num_threads(x) 
			for (int i = 0; i < 1000; i++) {
				res += create_thread(y);
			}
			double end = omp_get_wtime();
			if (x != 30) {
				fprintf(fp, "%lf,", end - start);
			} else {
				fprintf(fp, "%lf\n", end - start);
			}
		}
	}
	
	int res = 0;

	

	return 0;
}