#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

enum { BORDER_DEN = 10 };


double inverse_fact(int n)
{
  double fact(1.0);
  for(int i(1); i <= n; ++i)
  {
    fact /= double(i);
  }
  return fact;
}

double inverse_fact_from_border(int n, double border_value, int BORDER)
{
  double fact(1.0);
  for(int i(BORDER + 1); i <= n; ++i)
  {
    fact /= double(i);
  }
  return fact * border_value;
}

double calc(uint32_t x_last, uint32_t num_threads)
{
  int BORDER = x_last / BORDER_DEN;
  double sum(0.0);
  double border_value(0.0);
 
  for(int i(0); i < x_last; ++i)
  {
    double value(0.0);
    if(i > BORDER)
      value = inverse_fact_from_border(i, border_value, BORDER);
    else
      value = inverse_fact(i);
    sum += value;
    if(i == BORDER)
    {
      border_value = value;
      std::cout << "border value: " << 1.0 / border_value << std::endl;
    }
  }
  
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
