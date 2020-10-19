#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

enum Const
{
  NUM_PLACES = 100
};
double func(double x)
{
  return sin(x);
}

double calc(double x0, double x1, double dx, uint32_t num_threads)
{
  double res = 0.0;
  double * part_sum = new double[NUM_PLACES];
  int n_of_steps = (x1 - x0) / dx;
  int n_of_elems = n_of_steps / NUM_PLACES;
  #pragma omp parallel num_threads(num_threads)
  {
    #pragma omp for
    for(int i = 0; i < n_of_elems * NUM_PLACES; i++)
    {
      int id = i / n_of_elems;
      part_sum[id] += func(x0 + i * dx);
    }
  }
  for(int i = n_of_elems * NUM_PLACES; i < n_of_steps; i++)
    res += func(x0 + i * dx);
  for(int i = 0; i < NUM_PLACES; i++)
    res += part_sum[i];
  res *= dx;
  res += func(x1) * (x1 - x0 - n_of_steps * dx);
  delete [] part_sum;
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
