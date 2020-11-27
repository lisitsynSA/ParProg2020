#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>
#include <cassert>

void st_calc(double* frame, uint32_t ySize, uint32_t xSize, double delta, int rank, int size)
{
  if (rank == 0 && size > 0)
  {
    double diff = 0;
    double* tmpFrame = new double[ySize * xSize];
    // Prepare tmpFrame
    for (uint32_t y = 0; y < ySize; y++)
    {
      tmpFrame[y*xSize] = frame[y*xSize];
      tmpFrame[y*xSize + xSize - 1] = frame[y*xSize + xSize - 1];
    }
    for (uint32_t x = 1; x < xSize - 1; x++)
    {
      tmpFrame[x] = frame[x];
      tmpFrame[(ySize - 1)*xSize + x] = frame[(ySize - 1)*xSize + x];
    }
    // Calculate first iteration
    for (uint32_t y = 1; y < ySize - 1; y++)
    {
      for (uint32_t x = 1; x < xSize - 1; x++)
      {
        tmpFrame[y*xSize + x] = (frame[(y + 1)*xSize + x] + frame[(y - 1)*xSize + x] +\
                                frame[y*xSize + x + 1] + frame[y*xSize + x - 1])/4.0;
        diff += std::abs(tmpFrame[y*xSize + x] - frame[y*xSize + x]);
      }
    }

    double* currFrame = tmpFrame;
    double* nextFrame = frame;
    uint32_t iteration = 1;
    // Calculate frames
    while (diff > delta)
    {
      diff = 0;
      for (uint32_t y = 1; y < ySize - 1; y++)
      {
        for (uint32_t x = 1; x < xSize - 1; x++)
        {
          nextFrame[y*xSize + x] = (currFrame[(y + 1)*xSize + x] + currFrame[(y - 1)*xSize + x] +\
                                  currFrame[y*xSize + x + 1] + currFrame[y*xSize + x - 1])/4.0;
          diff += std::abs(nextFrame[y*xSize + x] - currFrame[y*xSize + x]);
        }
      }
      std::swap(currFrame, nextFrame);
      iteration++;
    }

    // Copy result from tmp
    if (iteration % 2 == 1)
    {
      for (uint32_t i = 0; i < xSize*ySize; i++)
      {
        frame[i] = tmpFrame[i];
      }
    }
    delete[] tmpFrame;
  }
}

void calc(double* frame, uint32_t y_size, uint32_t x_size, double delta, int rank, int size)
{
  if (size == 0)
    return;
  if (size == 1) {
    st_calc(frame, y_size, x_size, delta, rank, size);
    return;
  }
  MPI_Bcast(&x_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta,  1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  // In this implementation cannot parallel further than y_size/2,
  // because each process whould own at least two lines
  MPI_Comm my_comm;
  int color = (rank >= (int)y_size/2) ? MPI_UNDEFINED : 0;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &my_comm);
  if (rank >= (int)y_size/2)
    return;
  if (size > (int)y_size/2)
    size = y_size/2;


  // y_size-2 because 0 always owns first and last lines
  // size-1 because 0 process is exception
  uint32_t lines_per_process = (y_size-2)/(size-1);

  // INITIALIZATION

#define PROCESS_OWNED_LINE_NUM(RANK) \
  ({ \
    uint32_t tmp; \
    if ((RANK) == 0) \
      tmp = -1; \
    else if ((RANK) != size - 1) \
      tmp = lines_per_process; \
    else \
      tmp = (y_size-2) - (size-2) * lines_per_process; \
    tmp; \
  })

#define PROCESS_LINE_NUM(RANK) (((RANK) == 0) ? -1 : (PROCESS_OWNED_LINE_NUM(RANK)+2))

  // For all ranks except 0
  uint32_t my_owned_lines_num = PROCESS_OWNED_LINE_NUM(rank);
  uint32_t my_lines_num = PROCESS_LINE_NUM(rank);
  uint32_t my_cells_num = PROCESS_LINE_NUM(rank) * x_size;
  double* my_cells_1 = NULL, *my_cells_2 = NULL;
  double* input = NULL, *output = NULL;

  // Only for rank 0
  double* trash = (double*)calloc(x_size, sizeof(double)); // For receiving what we don't need
  assert(trash);

  if (rank == 0) {
    for (int prank = 1; prank < size; prank++) {
      printf("prank = %d, off = %d, size = %d\n", prank, x_size * (prank-1) * lines_per_process, PROCESS_LINE_NUM(prank) * x_size);
      MPI_Send(&frame[x_size * (prank-1) * lines_per_process], PROCESS_LINE_NUM(prank) * x_size,
               MPI_DOUBLE, prank, 0, my_comm);
    }
  } else {
    input = my_cells_1 = (double*)calloc(my_cells_num, sizeof(double));
    assert(my_cells_1);
    output = my_cells_2 = (double*)calloc(my_cells_num, sizeof(double));
    assert(my_cells_2);

    MPI_Recv(my_cells_1, my_cells_num, MPI_DOUBLE,
        0, MPI_ANY_TAG, my_comm, NULL);
    memcpy(my_cells_2, my_cells_1, sizeof(double) * my_cells_num);
  }

  // CALCULATION

  double diff = delta + 1;
  int i = 0;
  while (diff > delta) {
    diff = 0;
    if (rank) {
      for (uint32_t y = 1; y < my_lines_num - 1; y++) {
        for (uint32_t x = 1; x < x_size - 1; x++) {
          output[y*x_size + x] = (input[(y + 1)*x_size + x] + input[(y - 1)*x_size + x] +
                                  input[y*x_size + x + 1]   + input[y*x_size + x - 1])/4.0;
          diff += fabs(output[y*x_size + x] - input[y*x_size + x]);
        }
      }
      MPI_Recv(&output[x_size*0],                x_size, MPI_DOUBLE, rank - 1,          MPI_ANY_TAG, my_comm, NULL);
      MPI_Send(&output[x_size*1],                x_size, MPI_DOUBLE, rank - 1,          0,           my_comm);
      MPI_Send(&output[x_size*(my_lines_num-2)], x_size, MPI_DOUBLE, (rank + 1) % size, 0,           my_comm);
      MPI_Recv(&output[x_size*(my_lines_num-1)], x_size, MPI_DOUBLE, (rank + 1) % size, MPI_ANY_TAG, my_comm, NULL);

      // Swap, C-way
      double* tmp = output; output = input; input = tmp;
    }
    if (!rank) {
      MPI_Send(&frame[x_size*0],          x_size, MPI_DOUBLE, 1,        0,           my_comm);
      MPI_Recv(trash,                     x_size, MPI_DOUBLE, 1,        MPI_ANY_TAG, my_comm, NULL);
      MPI_Recv(trash,                     x_size, MPI_DOUBLE, size - 1, MPI_ANY_TAG, my_comm, NULL);
      MPI_Send(&frame[x_size*(y_size-1)], x_size, MPI_DOUBLE, size - 1, 0,           my_comm);
    }

    if (rank == 0)
      MPI_Reduce(MPI_IN_PLACE, &diff, 1, MPI_DOUBLE, MPI_SUM, 0, my_comm);
    else
      MPI_Reduce(&diff,        &diff, 1, MPI_DOUBLE, MPI_SUM, 0, my_comm);
    MPI_Bcast(&diff, 1, MPI_DOUBLE, 0, my_comm);
    if (rank == 0) {
      i++;
      if (i % 100 == 0)
        printf("diff = %lg\n", diff);
    }
  }

  // COLLECTION

  if (rank == 0) {
    for (int prank = 1; prank < size; prank++)
      MPI_Recv(&frame[x_size * ((prank-1) * lines_per_process + 1)], PROCESS_OWNED_LINE_NUM(prank) * x_size,
               MPI_DOUBLE, prank, MPI_ANY_TAG, my_comm, NULL);
  } else {
    MPI_Send(&input[x_size], my_owned_lines_num * x_size,
             MPI_DOUBLE, 0, 0, my_comm);
  }

  if (rank == 0)
    free(trash);
  else {
    free(my_cells_1);
    free(my_cells_2);
  }
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  double delta = 0;
  uint32_t ySize = 0, xSize = 0;
  double* frame = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> ySize >> xSize >> delta;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    frame = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> frame[y*xSize + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(frame, ySize, xSize, delta, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete[] frame;
      return 1;
    }
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << frame[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete[] frame;
  }

  MPI_Finalize();
  return 0;
}
