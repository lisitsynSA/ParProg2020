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

