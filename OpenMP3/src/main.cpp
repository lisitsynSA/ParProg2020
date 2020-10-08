#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <math.h>

#define PI (2 * asin(1))
#define DELTA(dx) (PI / dx)

double (*f)(double) = sin;

double calc(double x0, double x1, double dx, uint32_t num_threads)
{  
    if(x1 - x0 == 0) {
        return 0.0;
    }
    double res = 0.0;
    size_t n = size_t((x1 - x0) / dx) + 1; // #pragma omp for don't work with double type
    double* resbuf = (double*)calloc(n, sizeof(double));
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for 
        for(size_t _x = 0; _x < n; ++_x) {
            if(_x == 0 || _x == n - 1) {
                resbuf[_x] = f(x1 + _x * dx) / 2 * dx; 
            } else {
                resbuf[_x ]= f(x1 + _x * dx) * dx; 
            }
        }
    }
    size_t j = 0;
    size_t s = 0;
    for(size_t i = 0; i < n; ++i) {
        res += resbuf[j];
        j += DELTA(dx);
        if(j * dx > x1 - x0) {
            s++;
            j = s;
        }
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
  double x0 = 0.0, x1 =0.0, dx = 0.0;
  uint32_t num_threads = 0;
  input >> x0 >> x1 >> dx >> num_threads;

  // Calculation
  double res = calc(x0, x1, dx, num_threads);

  // Write result
<<<<<<< HEAD
  output << std::setprecision(15) << res << std::endl;
=======
  output << std::setprecision(13) << res << std::endl;
>>>>>>> 0f2d8930dcb55854e681a8927a619bc0a698c550
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
