#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <gmpxx.h>


void fact(int n){
  int i;
  mpz_t p;

  mpz_init_set_ui(p,1); /* p = 1 */

  for (i=1; i <= n ; ++i){
    mpz_mul_ui(p,p,i); /* p = p * i */
  }
  printf ("%d!  =  ", n);

}

#if 0
/**
 * Осуществляет принцип "разделяй и властвуй", где превращает список множителей в бинарное дерево
 * и перемножает их по 2.
 *
 * Рекурсивный алгоритм.
 */
double pod_tree_rec( int l, int r)
{
    if( l > r)
        return 1;

    if( l == r)
        return l;

    if( ( r - l) == 1)
        return (double)(l * r);

    int m = (l + r) / 2;
    
    return (pod_tree_rec( l, m) * pod_tree_rec( m + 1, r));
}

double tree_factorial( int n)
{
    if ( ( n == 0))
        return 1;

    if ( ( n == 2) || ( n == 1))
        return n;

    return ( pod_tree_rec( 2, n));
}

#endif
double calc(uint32_t x_last, uint32_t num_threads)
{
    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);
    
    //double res = 1.0;
    mpz_t res;
    double *factorial_array = (double*)calloc( x_last + 1, sizeof(double));
    factorial_array[0] = 1;

    #pragma omp parallel for num_threads( num_threads)
    for( uint32_t i = 1; i < x_last; i++)
        factorial_array[i] = tree_factorial( i);

    #pragma omp parallel for reduction(+:res) num_threads( num_threads)
    for( uint32_t i = 1; i < x_last - 1; i++)
    {
        res += 1.0 / factorial_array[i];
    }

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
  double res = calc(x_last, num_threads);

  // Write result
  output << std::setprecision(16) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
