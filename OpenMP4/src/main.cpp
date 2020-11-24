#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

double calc(uint32_t x_last, uint32_t num_threads) //Failed at Test05, the last number is 5, but in test it is 6
{
  double sum = 0.0;
  double pre_sum = 0.0;
  double arr_sum[num_threads];
  double k = 1.0;
  double arr_k[num_threads];
  #pragma omp parallel num_threads(num_threads) firstprivate(k, pre_sum) shared(arr_sum, arr_k)
  {
    #pragma omp for
    for(uint32_t i = 0; i < x_last; ++i)
    {
      pre_sum += k;
      k /= i + 1.0;
    }
    int thread_id = omp_get_thread_num();
    arr_k[thread_id] = k;
    arr_sum[thread_id] = pre_sum; 
  }
  for(uint32_t i = 1; i < num_threads - 1; ++i) //the last of arr_k isn't necessary
    arr_k[i] *= arr_k[i - 1];
  #pragma omp parallel for num_threads(num_threads) reduction(+:sum)
  for(uint32_t i = num_threads - 1; i > 0; --i)
    sum += arr_sum[i] * arr_k[i - 1];
  sum += arr_sum[0]; 
  return sum;
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
