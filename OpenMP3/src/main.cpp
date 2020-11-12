#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

<<<<<<< HEAD
double func(double x)
{
  return sin(x);
}

double calc(double x0, double x1, double dx, uint32_t num_threads)
{
  return 0;
=======
double calc() {
    size_t n_thread = 0;
    double x1, x2, dx, res = 0.0;
    printf("Enter first and last x, and dx\n");
    scanf("%lf %lf %lf", &x1, &x2, &dx);
    printf("Enter num of thread\n");
    scanf("%ld", &n_thread);
    int n = ceil((x2 - x1) / dx); // #pragma omp for don't work with double type
    
    #pragma omp parallel num_threads(n_thread)
    {
        #pragma omp for reduction(+:res) 
        for(int _x = 0; _x < n; ++_x) {
            res += (f(x1 + _x * dx) + f(x1 + (_x + 1)*dx)) / 2 * dx; 
        }
    }
    return res;
>>>>>>> for rebase
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
