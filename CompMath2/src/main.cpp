#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

/**
 * Функция, которая высчитывает границы 
 */
void calcIterations( int* leftBorder, int* rightBorder, int rank, double sizeOfBlock, int length)
{
    *leftBorder = rank * sizeOfBlock + 1;
    
    /* Если справа хватает элементов, то проблем нет, просто выделяем их */
    if ( (length - *leftBorder) >= sizeOfBlock)
    {
	*rightBorder = *leftBorder + sizeOfBlock - 1;
    } else
	/* Если справа не хватает элементов (например последний блок), то присваиваем правую границу */
    {
	*rightBorder = length - 1;
    }
}

/**
 * Функция для отладки. Печатаем фрейм
 */
void printFrame(double* frame, uint32_t ySize, uint32_t xSize)
{
    for (uint32_t y = 0; y < ySize; y++)
    {
	for (uint32_t x = 0; x < xSize; x++)
	{
	    printf("%f ", frame[y*xSize + x]);
	}
	printf("\n");
    }
    printf("\n");
}

void calc(double* frame, uint32_t ySize, uint32_t xSize, double delta, int rank, int size)
{
    double* tmpFrame = new double[ySize * xSize]();

    /* Приходится создать ещё один, так как
     * MPI_Reduce не поддерживает send = recv
     */
    double* recvFrame = new double[ySize * xSize]();  

    double diff = 100;
    double diff_part = 0;

    if (rank != 0)
    {
	frame = new double[ySize * xSize]();
    }

    MPI_Bcast(frame, xSize*ySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

    /* Для каждого процесса расчитаем его "область действия" */
    /* Найдем границы для каждого процесса */
    int	xLeftBorder = 0,
	xRightBorder = 0;
    
    int xSizeOfBlock = (int)ceil( (float)(xSize - 2) / size);

    calcIterations( &xLeftBorder, &xRightBorder, rank, xSizeOfBlock, xSize - 1);

    double* currFrame = tmpFrame;
    double* nextFrame = frame;

    uint32_t iteration = 1;

    // Calculate frames
    while (diff > delta)
    {
	diff = 0;
	diff_part = 0;
	 
	std::fill(nextFrame, nextFrame + ySize*xSize, 0);
	std::fill(recvFrame, recvFrame + ySize*xSize, 0);

	for (uint32_t y = 0; y < ySize; y++)
	{
	    nextFrame[y*xSize] = currFrame[y*xSize];
	    nextFrame[y*xSize + xSize - 1] = currFrame[y*xSize + xSize - 1];
	}
	
	for (uint32_t x = 1; x < xSize - 1; x++)
	{
	    nextFrame[x] = currFrame[x];
	    nextFrame[(ySize - 1)*xSize + x] = currFrame[(ySize - 1)*xSize + x];
	}
	
	for (uint32_t y = 1; y < ySize - 1; y++)
	{
	    for (int x = xLeftBorder; x <= xRightBorder; x++)
	    {
		nextFrame[y*xSize + x] = (currFrame[(y + 1)*xSize + x] + currFrame[(y - 1)*xSize + x] +\
                                  currFrame[y*xSize + x + 1] + currFrame[y*xSize + x - 1])/4.0;
		diff_part += std::abs(nextFrame[y*xSize + x] - currFrame[y*xSize + x]);
	    }
	}

	MPI_Allreduce( nextFrame, recvFrame, xSize * ySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
	MPI_Allreduce( &diff_part, &diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
	for (uint32_t y = 0; y < ySize; y++)
	{
	    recvFrame[y*xSize] = currFrame[y*xSize];
	    recvFrame[y*xSize + xSize - 1] = currFrame[y*xSize + xSize - 1];
	}
	
	for (uint32_t x = 1; x < xSize - 1; x++)
	{
	    recvFrame[x] = currFrame[x];
	    recvFrame[(ySize - 1)*xSize + x] = currFrame[(ySize - 1)*xSize + x];
	}

	memcpy(currFrame, recvFrame, xSize*ySize*sizeof(double));
	iteration++;
    }

    if( rank == 0)
    {
	memcpy(frame, recvFrame, xSize*ySize*sizeof(double));
    }

    delete tmpFrame;

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

  /* Разошлём всем процессам xSize, ySize и zSize*/
  MPI_Bcast( &xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  calc(frame, ySize, xSize, delta, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete frame;
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
    delete frame;
  }

  MPI_Finalize();
  return 0;
}
