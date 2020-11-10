#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

double calc3(uint32_t x_last, uint32_t num_threads)
{
    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    double res = 0.0;
    uint32_t i;
             
    double* res_buffer = (double*)calloc(x_last, sizeof(double));

    #pragma omp parallel for num_threads( num_threads) 
    for( i = x_last; i > 0; i--)
    {
        res_buffer[ i - 1] = 1.0 / i;
    }

    for ( i = x_last; i > 0; i--)
        res += res_buffer[ i - 1];
    
    return res;
}


double calc2(uint32_t x_last, uint32_t num_threads)
{
    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    double res = 0.0;
    uint32_t i;
             
    #pragma omp parallel for reduction(+:res) num_threads( num_threads) 
    for( i = x_last; i > 0; i--)
    {
        res += 1.0 / i;
    }
    
    return res;
}

       
double calc(uint32_t x_last, uint32_t num_threads)
{
    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    double res = 0.0;
    uint32_t i;
    
    double* res_buffer = (double*)calloc(num_threads, sizeof(double));
        
    #pragma omp parallel for num_threads( num_threads) 
    for( i = x_last; i > 0; i--)
    {
        res_buffer[ omp_get_thread_num()] += 1.0 / i;
    }
    
    for (int i = num_threads; i >= 0; i--)
        res += res_buffer[ i];
    
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
  double res = calc3(x_last, num_threads);

  // Write result
  output << std::setprecision(15) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
