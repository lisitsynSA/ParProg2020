#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>
#include <string.h>

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

/**
 * Функция, которая сдвигает массив обратно по оси x
 */
void BackShiftArray3D( double* arr, double* shifted_arr, int xSize, int ySize, int zSize)
{
    int lenX = 2 * xSize - 1;
    for ( int z = 0; z < zSize; z++)
    {
	for ( int y = 0; y < ySize; y++)
	{
	    for ( int x = 0; x < xSize; x++)
	    {
		arr[z*xSize*ySize + y*xSize + x] = shifted_arr[z*lenX*ySize + y*lenX + (ySize - y + x - 1)];
	    }
	}
    }

}

/**
 * Функция, которая выполняет сдвиг массив "послойно" в направлении оси x
 *
 * WARNING: Массив удаляется из вызывающей функции
 */
double* ShiftArray3D( double* arr, int xSize, int ySize, int zSize)
{
    int lenX = 2*xSize - 1;
    /* Массив для сдвига должен быть в 2 раза больше по оси x*/
    double* shifted_arr = new double[lenX * ySize * zSize];

    /* Заполним массив значениями -2, так как изначальная заполненность положительными числами, а значения
     * синуса лежат на отрезке [-1; 1] */
    for ( int z = 0; z < zSize; z++)
    {
	for ( int y = 0; y < ySize; y++)
	{
	    for (int x = 0; x < lenX; x++)
	    {
		shifted_arr[z*lenX*ySize + y*lenX + x] = -2;
	    }
	}
    }

    /* Сдвигаем массив по оси x */
    for ( int z = 0; z < zSize; z++)
    {
	for ( int y = 0; y < ySize; y++)
	{
	    for (int x = 0; x < lenX; x++)
	    {
		if( x < xSize)
		{
		    shifted_arr[z*lenX*ySize + y*lenX + (ySize - y + x - 1)] = arr[z*xSize*ySize + y*xSize + x];
		}
	    }
	}
    }

    return ( shifted_arr);
}

/* Здесь вектор расстояний (1, -1, -1) (<, >, >). Т.к. первый символ < - то это истинная зависимость */
void calc(double* arr, uint32_t zSize, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    int lenX = 2 * xSize - 1;
    double* bufArr = new double[ySize * lenX * zSize]();
    double* reduceArr = new double[ySize * lenX * zSize]();

    /* Сначала рассылаем изначальный массив, а потом сдвигаем его, иначе нам потребуется в 3 раза больше операций */
    if ( rank != 0)
    {
	arr = new double[ySize * xSize * zSize];
    }

    /* Передадим всем процессам изначальный массив arr */
    MPI_Bcast( arr, xSize * ySize * zSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Сдвинем массив по y и x */
    double* shifted_arr = ShiftArray3D( arr, xSize, ySize, zSize);
    
    /* Так как мы сдвинули массив, что у него нет зависимости по x, то мы можем 
     * разрезать на слои наш параллелепипед, со сдвинутыми слоями на отдельные части 
     * и посчитать в отдельных потоках */

    /* Найдем границы для каждого процесса */
    int	xLeftBorder = 0,
	xRightBorder = 0;
    
    int xSizeOfBlock = (int)ceil( (float)lenX / size);

    calcIterations( &xLeftBorder, &xRightBorder, rank, xSizeOfBlock, lenX);

    /* Распараллелим по x */
    for (uint32_t z = 1; z < zSize; z++) {
	for (uint32_t y = 0; y < ySize - 1; y++) {
	    for (int x = xLeftBorder; x <= xRightBorder; x++) 
	    {
		if ( (shifted_arr[(z - 1)*ySize*lenX + (y + 1)*lenX + x] != -2) && (shifted_arr[z*ySize*lenX + y*lenX + x] != -2))
		{
		    bufArr[z*ySize*lenX + y*lenX + x] = sin(shifted_arr[(z - 1)*ySize*lenX + (y + 1)*lenX + x]);
		    shifted_arr[z*ySize*lenX + y*lenX+ x] = bufArr[z*ySize*lenX + y*lenX + x];
		}
	    }
	}
    }

    /* В каждой области имеем посчитанное значение, теперь надо их просто сложить в 0 процессе */
    MPI_Reduce( bufArr, reduceArr, (int)(lenX * ySize * zSize), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    /* Возвращаем обратно к виду правильного параллелепипеда */
    delete bufArr;
    bufArr = new double[ySize * xSize * zSize]();

    BackShiftArray3D( bufArr, reduceArr, xSize, ySize, zSize);

    /* Там где не было вычислений записываем предыдущие невычесленные значения */
    for (uint32_t z = 0; z < zSize; z++)
	/* Первый срез z = 0 */
	if ( z == 0)
	    for ( int y = 0; y < (int)ySize; y++)
		for ( int x = 0; x < (int)xSize; x++)
		    bufArr[z*ySize*xSize + y*xSize + x] = arr[z*ySize*xSize + y*xSize + x];
	else
	    for ( int y = 0; y < (int)ySize; y++)
		/* Правый столбец в срезе по z */
		if ( y < (int)(ySize - 1))
		    for ( int x = xSize - 1; x < (int)xSize; x++)
			bufArr[z*ySize*xSize + y*xSize + x] = arr[z*ySize*xSize + y*xSize + x];
		else
		    /* Последняя строка в срезе по z */ 
		    for ( int x = 0; x < (int)xSize ; x++)
			bufArr[z*ySize*xSize + y*xSize + x] = arr[z*ySize*xSize + y*xSize + x];

    if( rank == 0)
    {
	memcpy( arr, bufArr, sizeof(double) * xSize * ySize * zSize);
    }

    delete bufArr;
    delete reduceArr;
    delete shifted_arr;
    
    if ( rank != 0)
	delete arr;
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
  /* Разошлём всем процессам xSize, ySize и zSize*/
  MPI_Bcast( &xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &zSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

  
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
