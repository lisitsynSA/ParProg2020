#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* arr, uint32_t zSize, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  MPI_Status status;
  MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&zSize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  uint32_t num_of_tasks, arr_size = 0, arr_addr = 0;

  if (rank == 0) {
    double* new_arr = (double *) malloc(zSize * ySize * xSize * sizeof(double));
    bool* added = (bool *) calloc(zSize * ySize * xSize, sizeof(bool));
    int z, y, x;
    uint32_t* tasks = (uint32_t*) malloc(zSize * ySize * xSize * sizeof(uint32_t));
    uint32_t counter = 0, task_counter = 0, local_counter = 0, send_num_of_tasks, send_task_addr, send_arr_size[size], send_arr_addr;
    for (int i = 0; i < (int)zSize; i++) {
      for (int j = 0; j < (int)ySize; j++) {
        for (int k = 0; k < (int)xSize; k++) {
          if (added[i * ySize * xSize + j * xSize + k] == false) {
            z = i; 
            y = j;
            x = k;
            local_counter = 0;
            while (z < (int)zSize && y >= 0 && x >= 0) {
              added[z * ySize * xSize + y * xSize + x] = true;
              new_arr[counter] = arr[z * ySize * xSize + y * xSize + x];
              counter++;
              local_counter++;
              z++;
              y--;
              x--;
            }
            tasks[task_counter] = local_counter;
            task_counter++;
          }
        }
      }
    }
    num_of_tasks = task_counter / size;
    if (task_counter % size > 0) {
      num_of_tasks++;
    }
    send_arr_addr = 0;
    for (int i = 0; i < (int)num_of_tasks; i++) {
      arr_size += tasks[i];
    }
    send_arr_addr = arr_size;
    send_num_of_tasks = num_of_tasks;
    send_task_addr = num_of_tasks;
    for (int i = 1; i < size; i++) {
      send_arr_size[i] = 0;
      if (i == (int)task_counter % size) {
        send_num_of_tasks--;
      }
      for (int j = 0; j < (int)send_num_of_tasks; j++) {
        send_arr_size[i] += tasks[send_task_addr + j];
      }
      MPI_Send(&send_num_of_tasks, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&send_arr_size[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&tasks[send_task_addr], send_num_of_tasks, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&new_arr[send_arr_addr], send_arr_size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      send_arr_addr += send_arr_size[i];
      send_task_addr += send_num_of_tasks;
    }
    for (int i = 0; i < (int)num_of_tasks; i++) {
      for (int j = 0; j < (int)(tasks[i] - 1); j++) {
        new_arr[j + arr_addr + 1] = sin(new_arr[j + arr_addr]); 
      }
      arr_addr += tasks[i];
    }
    send_arr_addr = arr_size;
    for (int i = 1; i < size; i++) {
      MPI_Recv(&new_arr[send_arr_addr], send_arr_size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      send_arr_addr += send_arr_size[i];
    }
    counter = 0;
    for (int i = 0; i < (int)zSize; i++) {
      for (int j = 0; j < (int)ySize; j++) {
        for (int k = 0; k < (int)xSize; k++) {
          if (added[i * ySize * xSize + j * xSize + k] == true) {
            z = i; 
            y = j;
            x = k;
            while (z < (int)zSize && y >= 0 && x >= 0) {
              added[z * ySize * xSize + y * xSize + x] = false;
              arr[z * ySize * xSize + y * xSize + x] = new_arr[counter];
              counter++;
              z++;
              y--;
              x--;
            }
          }
        }
      }
    }
  } else {
    MPI_Recv(&num_of_tasks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&arr_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    uint32_t* tasks = (uint32_t *) malloc(num_of_tasks * sizeof(uint32_t));
    arr = (double *) malloc(arr_size * sizeof(double));
    MPI_Recv(tasks, num_of_tasks, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(arr, arr_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

    std::ofstream output("output" + std::to_string(rank) + ".txt");
    for (uint32_t z = 0; z < arr_size; z++) {
      output << " " << arr[z];
    }
    output.close();
    for (int i = 0; i < (int)num_of_tasks; i++) {
      for (int j = 0; j < (int)(tasks[i] - 1); j++) {
        arr[j + arr_addr + 1] = sin(arr[j + arr_addr]); 
      }
      arr_addr += tasks[i];
    }
    MPI_Send(arr, arr_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t zSize = 0, ySize = 0, xSize = 0;
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
    input >> zSize >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[zSize * ySize * xSize];
    for (uint32_t z = 0; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize; y++) {
        for (uint32_t x = 0; x < xSize; x++) {
          input >> arr[z*ySize*xSize + y*xSize + x];
        }
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

  calc(arr, zSize, ySize, xSize, rank, size);

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

    for (uint32_t z = 0; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize; y++) {
        for (uint32_t x = 0; x < xSize; x++) {
          output << " " << arr[z*ySize*xSize + y*xSize + x];
        }
        output << std::endl;
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
