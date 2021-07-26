#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

double func(double x)
{
  return sin(x);
}

double calc_reduce(double x0, double x1, double dx, uint32_t num_threads)
{
    
    if ( x1 == x0)
        return 0;

    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    int     count      = ( int)floor( ( x1 - x0) / dx);
    double  res = 0.0;


    /* Подсчет 2 "хвостиков" */
    res += ( func( x0) + func( count)) * dx / 2;

    /* sum = [f(0) + f(x1)] * dx / 2 + (f(1) + f(2) + ... + f(x1-1)) * dx */
    /* Основной ход метода трапеций */
    #pragma for reduce(+:res) num_threads( num_threads) 
    for ( int i = 1; i <= count; i++)
    {
        res += func( i * dx);
    }

    res *= dx;
    res += + (func( dx * ( count - 1)) + func( x1))/2 * ( (x1 - x0) - (dx * ( count - 1)));
    
  
    return res;
}


double calc(double x0, double x1, double dx, uint32_t num_threads)
{
    
    if ( x1 == x0)
        return 0;

    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    int     count      = ( int)floor( ( x1 - x0) / dx);
    double* res_buffer = ( double*)calloc( ( count + 1), sizeof( double));
    double  res = 0.0;

    /* sum = [f(0) + f(x1)] * dx / 2 + (f(1) + f(2) + ... + f(x1-1)) * dx */
    /* Основной ход метода трапеций */
    #pragma omp parallel for num_threads( num_threads) 
    for ( int i = 0; i <= count; i++)
    {
        res_buffer[i] = func( i * dx);
    }
    

    /* Сложение результатов */
    for( int i = 1; i < count; i++)
        res += res_buffer[i];

    res *= dx;
    
    /* Подсчет 2 "хвостиков" */
    res += (res_buffer[0] + res_buffer[count]) * dx / 2 + 
            ( func( dx * ( count - 1)) + func( x1))/2 * ( (x1 - x0) - (dx * ( count - 1)));

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
