#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

double calc(uint32_t x_last, uint32_t num_threads)
{
  double* loc_res = (double*)calloc(num_threads, sizeof(double));
  double* loc_fact = (double*)malloc(num_threads * sizeof(double));

  #pragma omp parallel num_threads(num_threads) 
  {
    int tid = omp_get_thread_num();
    loc_fact[tid] = 1.0;

    #pragma omp for 
    for (uint32_t i = 1; i < x_last; i++) {
      loc_fact[tid] /= double(i);
      loc_res[tid] += loc_fact[tid];
    } 
  }

  double res = 1.0;
  double fact = 1.0;

  for (uint32_t tid = 0; tid < num_threads; tid++) {
    res += loc_res[tid] * fact;
    fact *= loc_fact[tid];
  }
  
  free(loc_fact);
  free(loc_res);
  return res;
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 3)
  {
    std::cout << "[Error] Usage <inputfile> <output file>\n";
    return 1;
  }

  // Prepare input file
  std::ifstream input(argv[1]);
  if (!input.is_open())
  {
    std::cout << "[Error] Can't open " << argv[1] << " for write\n";
    return 1;
  }

  // Prepare output file
  std::ofstream output(argv[2]);
  if (!output.is_open())
  {
    std::cout << "[Error] Can't open " << argv[2] << " for read\n";
    input.close();
    return 1;
  }

// Read arguments from input
  uint32_t x_last = 0, num_threads = 0;
  input >> x_last >> num_threads;

  // Calculation
  double res = calc(x_last, num_threads);

  // Write result
  output << std::setprecision(16) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
