#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#define ROOT 0

double acceleration(double t) {
    return sin(t);
}
/*
  // Sighting shot
  double v0 = 0;
  if (rank == 0 && size > 0)
  {
    trace[0] = y0;
    trace[1] = y0 + dt*v0;
    for (uint32_t i = 2; i < traceSize; i++)
    {
      trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*trace[i - 1] - trace[i - 2];
    }
  }

  // The final shot
  if (rank == 0 && size > 0)
  {
    v0 = (y1 - trace[traceSize - 1])/(dt*traceSize);
    trace[0] = y0;
    trace[1] = y0 + dt*v0;
    for (uint32_t i = 2; i < traceSize; i++)
    {
      trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*trace[i - 1] - trace[i - 2];
    }
  }
 

 */

void calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
{
    int* recv_dist = NULL;
    int* recv_range = NULL;

    int real_rank_size = 0;

    if(rank == ROOT) {
        size_t my_calc_size = std::max(ceil((traceSize) * 1.0 / size), 3.0);
        // if  proc_rank > ceil(traceSize / my_calc_size) it will do nothing
        real_rank_size = ceil((traceSize) * 1.0 / my_calc_size);
        size_t now_size = 0; 
    
        recv_dist = (int*)calloc(size, sizeof(int));
        recv_range = (int*)calloc(size, sizeof(int));
    
        for(int i = 0; i < real_rank_size; ++i) {
            recv_dist[i] = i * my_calc_size;
            recv_range[i] = std::min(my_calc_size, traceSize - now_size);
            now_size += std::min(my_calc_size, traceSize - now_size);
        }
        
        for(int i = 1; i < real_rank_size; ++i) {
            recv_dist[i] -= 1;
            recv_range[i] += 1;
        }
        for(int i = 0; i < real_rank_size - 1; ++i) {
            recv_range[i] += 1;
        }
    }    
    int recv_step = 0;
    int my_delta_t = 0;

    MPI_Scatter(
            recv_range, 1, MPI_INT,
            &recv_step, size, MPI_INT, 
            ROOT, MPI_COMM_WORLD);

    MPI_Scatter(
            recv_dist, 1, MPI_INT,
            &my_delta_t, size, MPI_INT,
            ROOT, MPI_COMM_WORLD);

    double* my_trace = (double*)calloc(recv_step, sizeof(double));
    // [0](math)(start_cond) -> [1](math)(start_cond) -> [2](math(start_cond)) -> ...
    
    // Sighting shot
    double my_y0 = 0;
    double my_y1 = 0;

    printf("[%d] step = %d\n", rank, recv_step);

    if(rank == ROOT) {
        MPI_Send(&y0, 1, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
        MPI_Send(&y0, 1, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
    }

    MPI_Recv(&my_y0, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&my_y1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, NULL);
    my_trace[0] = my_y0;
    my_trace[1] = my_y1;
    for (uint32_t i = 2; i < (uint32_t)recv_step; i++) {
        my_trace[i] = dt*dt*acceleration(t0 + (my_delta_t + i - 1)*dt) + 2*my_trace[i - 1] - my_trace[i - 2];
    }
    if(rank != size - 1) {
        MPI_Send(&my_trace[recv_step - 2], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(&my_trace[recv_step - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    } else {
        // send smth in y1 place!
        MPI_Send(&my_trace[recv_step - 1], 1, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
    }
    
    // The final shot
    double last_y = 0;
    if(rank == ROOT) {
        MPI_Recv(&last_y, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, NULL);
        
        MPI_Send(&y0, 1, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
        next_y0 = y0 + (y1 - last_y)/(1.0 * traceSize);
        MPI_Send(&next_y0, 1, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
    }
    
    MPI_Recv(&my_y0, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&my_y1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, NULL);
    my_trace[0] = my_y0;
    my_trace[1] = my_y1;
    
    for (uint32_t i = 2; i < (uint32_t)recv_step; i++) {
        my_trace[i] = dt*dt*acceleration(t0 + (my_delta_t + i - 1)*dt) + 2*my_trace[i - 1] - my_trace[i - 2];
    }
    if(rank != size - 1) {
        MPI_Send(&my_trace[recv_step - 2], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(&my_trace[recv_step - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }

    MPI_Gatherv(
            my_trace, recv_step, MPI_DOUBLE,
            trace, recv_range, recv_dist, MPI_DOUBLE, 
            ROOT, MPI_COMM_WORLD); 

    free(my_trace);
    if(rank == ROOT) {
        free(recv_range);
        free(recv_dist);
    }
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  uint32_t traceSize = 0;
  double t0 = 0, t1 = 0, dt = 0, y0 = 0, y1 = 0;
  double* trace = 0;

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
    input >> t0 >> t1 >> dt >> y0 >> y1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    traceSize = (t1 - t0)/dt;
    trace = new double[traceSize];

    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(trace, traceSize, t0, dt, y0, y1, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete trace;
      return 1;
    }

    for (uint32_t i = 0; i < traceSize; i++)
    {
      output << " " << trace[i];
    }
    output << std::endl;
    output.close();
    delete trace;
  }

  MPI_Finalize();
  return 0;
}
