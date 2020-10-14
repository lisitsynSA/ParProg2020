#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

double func(double x)
{
  return sin(x);
}

enum { TOKEN_NUM = 100 };

double calc(double x0, double x1, double dx, uint32_t num_threads)
{
  double sum(0.0);

  // allocate store buffer
  // buffer size is significantly less than x1-x0 to avoid frequents RAM access

  double* token = new double[TOKEN_NUM];

  // int i(0);
  int num_steps = (x1 - x0) / dx; // 11
  // std::cout << "num steps " << num_steps << std::endl;
  int token_size = num_steps / TOKEN_NUM; // 5
  // std::cout << "token size " << token_size << std::endl;

  #pragma omp parallel num_threads(num_threads)
  {
    #pragma omp for
    for(int i = (1); i < token_size * TOKEN_NUM; ++i)
    {
      int token_idx = i / token_size; // 0 0 0 0 0 1 1 1 1 1 1 2
      token[token_idx] += (func(x0 + i * dx));
    }
  }
  for(int i(token_size * TOKEN_NUM); i < num_steps; ++i)
    sum += (func(x0 + i * dx));// + func(x0 + (i + 1) * dx));

  for(int i(0); i < TOKEN_NUM; ++i)
  {
    sum += token[i];
    // std::cout << token[i] << std::endl;
  }
  sum += (func(x0) + func(x1)) / 2.0;
  sum *= dx;

  delete token;
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
  double x0 = 0.0, x1 =0.0, dx = 0.0;
  uint32_t num_threads = 0;
  input >> x0 >> x1 >> dx >> num_threads;

  // Calculation
  double res = calc(x0, x1, dx, num_threads);

  // Write result
  output << std::setprecision(13) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
