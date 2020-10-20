#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

enum BORDERS
{
  N_of_borders = 100,
  Border_gate = 4
};


double inverted_fact(int n)
{
  double fact = 1.0;
  for(int i = 1; i <= n; ++i)
  {
    fact /= double(i);
  }
  return fact;
}

double inverted_fact_from_border(int n, double border_value, int border)
{
  double fact = 1.0;
  for(int i = border + 1; i <= n; ++i)
  {
    fact /= double(i);
  }
  return fact * border_value;
}

double calc(uint32_t x_last, uint32_t num_threads)
{
  double sum = 0.0;
  double* res = new double[num_threads];
  double* borders = new double[N_of_borders];
  double* border_value = new double[N_of_borders];

  borders[0] = Border_gate;
  border_value[0] = inverted_fact(Border_gate);
  //std::cout << borders[0] << " " << border_value[0] << std::endl;
  for(int i = 1; i < N_of_borders; i++)
  {
    borders[i] = Border_gate * (i + 1);
    border_value[i] = inverted_fact_from_border(Border_gate * (i + 1), border_value[i - 1], Border_gate * i);
    //std::cout << borders[i] << " " << border_value[i] << std::endl;
  }

  int new_x_last = std::min(int(x_last), Border_gate * N_of_borders);
  #pragma omp parallel num_threads(num_threads)
  {
    #pragma omp for
    for(int i = 0; i < new_x_last; ++i)
    {
      double temp_value = 0.0;
      if(i < borders[1])
        temp_value = inverted_fact(i);

      for(int j = 1; j < N_of_borders - 1; ++j)
      {
        if(i == borders[j])
          temp_value = border_value[j]; 
        if(i > borders[j])
          temp_value = inverted_fact_from_border(i, border_value[j], borders[j]);
      }
      int thread_id = omp_get_thread_num();
      res[thread_id] += temp_value;
    }
  }
  for(int i = 0; i < int(num_threads); ++i)
    sum += res[i];

  delete [] res;
  delete [] borders;
  delete [] border_value;
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
