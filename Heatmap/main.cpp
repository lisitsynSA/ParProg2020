#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <sys/time.h>
#include <climits>
#include <string>


int create_threads(int n_of_threads);
long int ReadArg(char * str);

int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 3)
  {
    std::cout << "[Error] Usage <n_of_external_processes> <n_of_internal_processes>\n";
    return 1;
  }

  int n_of_external_threads = ReadArg(argv[1]);
  int n_of_internal_threads = ReadArg(argv[2]);

  // enables to run nested parallelism
  omp_set_nested(true);
  int res = 0;
  double start = omp_get_wtime();

  // start of parallel code
  #pragma omp parallel num_threads(n_of_external_threads)
  {
    for(int i = 0; i < 100; ++i)
    {
      res += create_threads(n_of_internal_threads);
    }
  }

  // end of parallel code

  double end = omp_get_wtime();

  double delta = end - start;
  std::cout << delta << std::endl;

  return 0;
}

int create_threads(int n_of_threads)
{
  int res = 0;
  #pragma omp parallel num_threads(n_of_threads)
  {
    res++;
  }
  return res;
}

long int ReadArg(char * str)
{
  char* endptr;
  errno = 0;

  long int number = strtol(str, &endptr, 10);

  
  if ((errno == ERANGE && (number == LONG_MAX || number == LONG_MIN)) || (errno != 0 && number == 0)) 
  {
          perror("strtol");
          exit(EXIT_FAILURE);
    }

  if (endptr == str)
  {
          fprintf(stderr, "Error!\n");
          exit(EXIT_FAILURE);
    }
  if (*endptr != '\0')
  {
          fprintf(stderr, "Error!\n");
          exit(EXIT_FAILURE);
    }

  return number;
}

