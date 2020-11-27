#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>
#include <cassert>

void st_calc(double* arr, uint32_t y_size, uint32_t x_size, int rank, int size)
{
  for (uint32_t y = 4; y < y_size; y++)
    for (uint32_t x = 0; x < x_size; x++)
      arr[y*x_size + x] = sin(arr[(y - 4)*x_size + x]);
}

void calc(double* arr, uint32_t y_size, uint32_t x_size, int rank, int size)
{
  if (size == 0)
    return;
  if (size == 1) {
    st_calc(arr, y_size, x_size, rank, size);
    return;
  }
  MPI_Bcast(&x_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  // In this implementation cannot parallel further than x_size cells
  // (theoretical limit (4*x_size))
  MPI_Comm my_comm;
  int color = (rank >= (int)x_size) ? MPI_UNDEFINED : 0;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &my_comm);
  if (rank >= (int)x_size)
    return;
  if (size > (int)x_size)
    size = x_size;

  double* transponded;
  if (rank == 0) {
    transponded = (double*)calloc(x_size * y_size, sizeof(double));
    assert(transponded);

    for (uint32_t i = 0; i < y_size; i++)
      for (uint32_t j = 0; j < x_size; j++)
        transponded[i*x_size + j] = arr[j*y_size + i];
  }

  uint32_t lines_per_process = x_size / size;
  uint32_t cells_per_process = lines_per_process * y_size;

  double* mycells = (double*)calloc(cells_per_process, sizeof(double));
  assert(mycells);

  MPI_Scatter(transponded, cells_per_process, MPI_DOUBLE,
              mycells,     cells_per_process, MPI_DOUBLE,
              0, my_comm);

  for (uint32_t x = 0; x < lines_per_process; x++)
    for (uint32_t y = 4; y < y_size; y++)
      mycells[x*y_size + y] = sin(mycells[x*y_size + (y - 4)]);

  MPI_Gather(mycells,     cells_per_process, MPI_DOUBLE,
             transponded, cells_per_process, MPI_DOUBLE,
             0, my_comm);

  free(mycells);

  // Compute remainder
  if (rank == 0) {
    for (uint32_t x = lines_per_process * size; x < x_size; x++)
      for (uint32_t y = 4; y < y_size; y++)
        transponded[x*y_size + y] = sin(transponded[x*y_size + (y - 4)]);

    for (uint32_t i = 0; i < y_size; i++)
      for (uint32_t j = 0; j < x_size; j++)
        arr[i*x_size + j] = transponded[j*y_size + i];

    free(transponded);
  }

  MPI_Comm_free(&my_comm);
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t y_size = 0, x_size = 0;
  double* arr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> y_size >> x_size;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[y_size * x_size];

    for (uint32_t y = 0; y < y_size; y++)
    {
     for (uint32_t x = 0; x < x_size; x++)
      {
        input >> arr[y*x_size + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (buf != 0)
    {
      return 1;
    }
  }

  calc(arr, y_size, x_size, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete arr;
      return 1;
    }
    for (uint32_t y = 0; y < y_size; y++)
    {
      for (uint32_t x = 0; x < x_size; x++)
      {
        output << " " << arr[y*x_size + x];
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
