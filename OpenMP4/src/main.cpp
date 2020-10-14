#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

enum BORDERS
{
 BORDER_DEN = 10,
 BORDER_NUM = 2
};


double inverse_fact(int n)
{
  double fact(1.0);
  for(int i(1); i <= n; ++i)
  {
    fact /= double(i);
  }
  // std::cout << "inverse fact[" << n << "] " << fact << std::endl;
  return fact;
}

double inverse_fact_from_border(int n, double border_value, int BORDER)
{
  double fact(1.0);
  for(int i(BORDER + 1); i <= n; ++i)
  {
    fact /= double(i);
    // std::cout << "inverse fact from border: " << fact << std::endl;
  }
  // std::cout << "border value: " << border_value << std::endl;
  // std::cout << "inverse fact from border: " << fact * border_value << std::endl;
  return fact * border_value;
}

double calc(int x_last, int num_threads)
{
  std::cout << "STEPS: " << x_last << std::endl;
	num_threads++;
	num_threads--;
 /* int BORDER = x_last / BORDER_DEN;
  int BORDER1 = x_last / 10;*/
  double sum(0.0);
  // double border_value(0.0);

  int* borders = new int[BORDER_NUM];
  double* border_value = new double[BORDER_NUM];
  int border_size = x_last / BORDER_NUM;
  for(int i(0); i < BORDER_NUM; ++i)
  {
    borders[i] = i * border_size;
    // std::cout << "borders[" << i << "] = " << borders[i] << std::endl;
  }

  double* res = new double[num_threads];
  #pragma omp parallel num_threads(num_threads)
  {
    #pragma omp for
    for(int i = 0; i < x_last; ++i)
    {
      double value(0.0);
      if(i <= borders[1])
        value = inverse_fact(i);

      for(int j = 1; j < BORDER_NUM; ++j)
      {
        /*if(i > borders[j + 1])
          value = inverse_fact_from_border(i, border_value[j + 1], borders[j + 1]);
        if((i <= borders[j + 1]) && (i > borders[j]))
          value = inverse_fact_from_border(i, border_value[j], borders[j]);*/
        if(i == borders[j])
        {
          border_value[j] = value;
          // std::cout << " i = " << i << " border value[" << j << "] = " << border_value[j] << std::endl;
        } else if(i > borders[j])
          value = inverse_fact_from_border(i, border_value[j], borders[j]);
      }
      int tid = omp_get_thread_num();
      res[tid] += value;
      // sum += value;
    }
  }
  sum = 0.0;
  for(int i(0); i < num_threads; ++i)
    sum += res[i];

  delete[] borders;
  delete[] border_value;
  delete[] res;
  
  // std::cout << "SUM = " << sum << std::endl;
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
  int x_last = 0, num_threads = 0;
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
