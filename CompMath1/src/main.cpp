#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

double acceleration(double t)
{
  return sin(t);
}

void calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
{
  int count, pre_count;
  int mod;
  double v, u;
  double start_t;

  // Sighting shot
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0 && size > 0)
  {
    count = traceSize / size;
    mod = traceSize % size;
    
    ++count;
    for(int i = 0; i < mod; ++i)
    {
      start_t = t0 + dt * ((count - 1) + i * count);
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&start_t, 1 , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }

    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      start_t = t0 + dt * (count + mod + i * count);
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&start_t, 1 , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }
    start_t = t0;
    
  }
  if (rank != 0 && size > 0)
  {
    MPI_Recv(&count, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&start_t, 1 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  v = 0;
  double* local_array = new double[count + 1];
  if (rank != 0 && size > 0)
    local_array[0] = 0;
  if (rank == 0 && size > 0)
    local_array[0] = y0;
  local_array[1] = local_array[0] + dt * v;
  for (int i = 2; i < (count + 1); i++)
    local_array[i] = dt*dt*acceleration(start_t + (i - 1)*dt) + 2*local_array[i - 1] - local_array[i - 2];
  u = (local_array[count] - local_array[count - 1]) / dt; // b = local_array[count]

  double tmp_a, tmp_b, tmp_u, tmp_v;

  if (rank != 0 && size > 0)
  {
    MPI_Recv(&pre_count, 1 , MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tmp_a, 1 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tmp_b, 1 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tmp_u, 1 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tmp_v, 1 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array[0] = tmp_a + tmp_b + dt * pre_count * tmp_v;
    v = tmp_u + tmp_v;   
  }
  if ((rank != size - 1) && size > 1)
  {
    MPI_Send(&count, 1 , MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(&local_array[0], 1 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(&local_array[count], 1 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(&u, 1 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(&v, 1 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  }
  if (rank == size - 1)
  {
    if (size > 1)
    {
      double tmp_y1_miss = local_array[0] + local_array[count - 1] + dt * (count - 1) * v; //because we dont need local_array[count]
      MPI_Send(&tmp_y1_miss, 1 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  }


  // The final shot
  if (rank == 0 && size > 0)
  {
    double y1_miss;
    if (size > 1)
      MPI_Recv(&y1_miss, 1 , MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else
      y1_miss = local_array[0] + local_array[count - 1] + dt * (count - 1) * v;
    printf("y1_miss: %f\n", y1_miss);
    v = (y1 - y1_miss) / (dt * traceSize);
    local_array[0] = y0;
  }
  if (rank != 0 && size > 0)
  {
    double tmp_a_2, tmp_v_2;
    MPI_Recv(&tmp_a_2, 1 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tmp_v_2, 1 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array[0] = tmp_b + tmp_a_2 + dt * pre_count * tmp_v_2;
    v = tmp_u + tmp_v_2;
  }
  if ((rank != size - 1) && size > 0)
  { 
    MPI_Send(&local_array[0], 1 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(&v, 1 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  }
  local_array[1] = local_array[0] + dt * v;
  for(int i = 2; i < count; ++i)
    local_array[i] = dt*dt*acceleration(start_t + (i - 1)*dt) + 2*local_array[i - 1] - local_array[i - 2];
  
  if (rank != 0 && size > 0)
    MPI_Send(local_array, count , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  if (rank == 0 && size > 0)
  {
    for(int i = 0; i < count; ++i)
      trace[i] = local_array[i];
    ++count;
    for(int i = 0; i < mod; ++i)
      MPI_Recv(&trace[(count - 1) + i * count], count , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    --count;
    for(int i = mod; i < size - 1; ++i)
      MPI_Recv(&trace[count + mod + i * count], count , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  } 
  delete local_array;
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
