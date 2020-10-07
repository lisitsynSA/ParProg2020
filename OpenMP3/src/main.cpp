#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

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

  // Calculation
  double res = calc();

  // Write result
  output << std::setprecision(15) << res;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
