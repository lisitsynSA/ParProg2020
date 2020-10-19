#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

double func(double x)
{
  return sin(x);
}

double trapezoid(double x0, double x1, double dx) {
  return dx * (func(x0) + func(x1)) * 0.5;
}

double calc(double x0, double x1, double dx, uint32_t num_threads)
{
  int num_of_dots = (x1 - x0) / dx + 1;
  double* dots = (double*)malloc(num_of_dots * sizeof(double));
  double res = 0.0;
  double trap = 0.0;
  double c = 0.0;
  double y = 0.0;
  double t = 0.0;
  #pragma omp parallel num_threads(num_threads) 
  {
    #pragma omp for 
    for (int i = 0; i < num_of_dots; i++) {
        dots[i] = x0 + i * dx;
    }
    #pragma omp for reduction(+:res) private(y, t, c, trap)
    for (int i = 0; i < num_of_dots - 1; i++) { 
      trap = trapezoid(dots[i], dots[i+1], dx);
      y = trap - c;
      t = res + y;
      c = (t - res) - y;
      res = t;
    }
  }
  free(dots);
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
