#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

double func(double x)
{
  return sin(x);
}

double calc(double x0, double x1, double dx, uint32_t num_threads)
{
  int steps_num = round((x1 - x0) / dx);
  double sum = (- func(x0) / 2.0 + func(x0 + steps_num * dx) / 2.0) * dx;  //<-*dx
  #pragma omp parallel for num_threads(num_threads) reduction(+:sum)
  for(int i = 0; i < steps_num; ++i)
    sum += func(x0 + i * dx) * dx;  //<--*dx
  //sum *= dx; //At first, i used this, and without "* dx" in marked "<--*dx" strings. But tests 06 and 07 were failed (1.99999999 != 2)
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
