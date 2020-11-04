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
    *leftBorder = rank * sizeOfBlock;
    
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

/* Здесь вектор зависимости (+4, 0) (>, =) - истинная зависимость, значения сначала вычисляются, 
 * а затем используются. Можно распараллелить только по внутреннему циклу */

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    double* bufArr = new double[ySize * xSize]();
    double* reduceArr = new double[ySize * xSize]();

    if ( rank != 0)
    {
	arr = new double[ySize * xSize];
    }

    /* Передадим всем процессам изначальный массив arr */
    MPI_Bcast( arr, xSize * ySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

     /* Найдем границы для каждого процесса */
    int	xLeftBorder = 0,
	xRightBorder = 0;
    
    int xSizeOfBlock = (int)ceil( (float)xSize / size);

    calcIterations( &xLeftBorder, &xRightBorder, rank, xSizeOfBlock, (int)xSize);

    for (uint32_t y = 4; y < ySize; y++)
    {
      for (int x = xLeftBorder; x <= xRightBorder; x++)
      {
	bufArr[y*xSize + x] = sin(arr[(y - 4)*xSize + x]);
	arr[y*xSize + x] = bufArr[y*xSize + x];
      }
    }

    /* В каждой области имеем посчитанное значение, теперь надо их просто сложить в 0 процессе */
    MPI_Reduce( bufArr, reduceArr, (int)(xSize * ySize), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    for (int i = 0; i < 4; i++)
	for ( int j = 0; j < (int)xSize; j++)
	    reduceArr[i*xSize + j] = arr[i*xSize + j];
    
    if( rank == 0)
    {
	memcpy( arr, reduceArr, sizeof(double) * xSize * ySize);
    }

    delete bufArr;
    delete reduceArr;

}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t ySize = 0, xSize = 0;
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
    input >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> arr[y*xSize + x];
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

    /* Разошлём всем процессам xSize и ySize */
  MPI_Bcast( &xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  calc(arr, ySize, xSize, rank, size);

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
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << arr[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
