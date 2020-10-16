#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

#define DEAD '0'
#define ALIVE '1'

/**
 * Конвертирует массив в 2D для более удобной работы
 *
 * WARNING Массив удаляется из вызывающей функции
 */
uint8_t** ConvertFrom1Dto2D( uint32_t xSize, uint32_t ySize, uint8_t* input_frame)
{
    uint8_t** temp_frame = (uint8_t**)calloc( ySize, sizeof(uint8_t*));
    for ( uint32_t i = 0; i < xSize; i++)
    {
        temp_frame[i] = (uint8_t*)calloc( ySize, sizeof(uint8_t));
    }
    
    /* x - столбцы, y - строки */
    for ( uint32_t i = 0; i < ySize; i++)
        for ( uint32_t j = 0; j < xSize; j++)
            temp_frame[i][j] = input_frame[i * ySize + j];
    return ( temp_frame);
}

/**
 * Переводит матрицу в линейный массив
 */
uint8_t* ConvertFrom2Dto1D( uint32_t xSize, uint32_t ySize, uint8_t** input_frame)
{
    
    uint8_t* temp_frame = (uint8_t*)calloc( ySize * xSize, sizeof(uint8_t));

    for ( uint32_t i = 0; i < ySize; i++)
	for( uint32_t j = 0; j < xSize; j++)
	    temp_frame[i * xSize + j] = input_frame[i][j];

    return ( temp_frame);
}

/**
 * Вычисляет количество соседей возле указанной точки (x, y)
 */
int NumOfNeighbours( uint32_t xSize, uint32_t ySize, uint32_t x, uint32_t y, uint8_t** input_frame)
{
    int neighbours_num = 0;
    xSize--;
    ySize--;
    /* Необходимо вокруг посмотреть 8 клеток, если там, кто-то есть, то количество соседей увеличиваем на 1 */

    /* Левый верхний угол */
    if ( ( x >= 1) && ( y >= 1))
    	if ( input_frame[x - 1][y - 1] == ALIVE)
	    neighbours_num++;
    
    /* Верхний средний */
    if ( x >= 1)
	if ( input_frame[x - 1][y] == ALIVE)
	    neighbours_num++;

    /* Правый сверху */
    if ( ( x >= 1) && ( y <= ( ySize - 1)))
	if ( input_frame[x - 1][y + 1] == ALIVE)
	     neighbours_num++;

    /* Слева средний */
    if ( y >= 1)
	if ( input_frame[x][y - 1] == ALIVE)
	    neighbours_num++;

    /* Справа средний */
    if ( y <= ( ySize - 1))
	if ( input_frame[x][y + 1] == ALIVE)
	    neighbours_num++;

    /* Слева нижний */
    if ( ( x <= ( xSize - 1)) && ( y >= 1))
	if ( input_frame[x + 1][y - 1] == ALIVE) 
	    neighbours_num++;

    /* Нижний средний */
    if ( x <= ( xSize - 1))
	if ( input_frame[x + 1][y] == ALIVE)
	    neighbours_num++;

    /* Справа нижний */      
    if ( ( y <= ( ySize - 1)) && ( x <= (xSize - 1)))
	if ( input_frame[x + 1][y + 1] == ALIVE)
	    neighbours_num++;

    return ( neighbours_num);
}


uint8_t* calc(uint32_t xSize, uint32_t ySize, uint32_t iterations, uint32_t num_threads, uint8_t* inputFrame)
{
    if (iterations == 0)
    {
	return inputFrame;
    }


    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    /* Для каждой ячейки необходимо проверить соседние и посчитать количество соседей */
    uint8_t** temp_frame = ConvertFrom1Dto2D( xSize, ySize, inputFrame);
   
    uint8_t** result_frame = (uint8_t**)calloc( xSize, sizeof(uint8_t*));
    for ( uint32_t i = 0; i < xSize; i++)
    {
        result_frame[i] = (uint8_t*)calloc( ySize, sizeof(uint8_t));
    }
    
    //#pragma omp parallel for num_threads( num_threads) 
    for ( uint32_t k = 0; k < iterations; k++)
    {
	for ( uint32_t i = 0; i < xSize; i++)
	{
	    for ( uint32_t j = 0; j < ySize; j++)
	    {
		/* В пустой клетке, рядом с которой есть ровно 3 соседа, зарождается жизнь */
		if ( ( NumOfNeighbours( xSize, ySize, i, j, temp_frame) == 3) && (temp_frame[i][j] == DEAD))
		    result_frame[i][j] = ALIVE;
		else
		/* Если у живой клетки 2 или 3 соседа, то она продолжает жить */
		if ( ( ( NumOfNeighbours( xSize, ySize, i, j, temp_frame) == 2) || ( NumOfNeighbours( xSize, ySize, i, j, temp_frame) == 3)) && (temp_frame[i][j] == ALIVE))
		    result_frame[i][j] = ALIVE;
		else
		    result_frame[i][j] = DEAD;
	    }
	}

	temp_frame = result_frame;
    }


    return ( ConvertFrom2Dto1D( xSize, ySize, result_frame));
    
/*
    for ( uint32_t i = 0; i < xSize; i++)
    {
	free( temp_frame[i]);
	free( result_frame[i]);
    }

    free( temp_frame);
    free( result_frame);
    */
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 3)
  {
    std::cout << "[Error] Usage <inputfile> <output file>\n";
    return 1;
  }

  // Prepare input file
  std::ifstream input(argv[1]);
  if (!input.is_open())
  {
    std::cout << "[Error] Can't open " << argv[1] << " for write\n";
    return 1;
  }

  // Prepare output file
  std::ofstream output(argv[2]);
  if (!output.is_open())
  {
    std::cout << "[Error] Can't open " << argv[2] << " for read\n";
    input.close();
    return 1;
  }

  // Read arguments from input
  uint32_t xSize = 0, ySize = 0, iterations = 0, num_threads = 0;
  input >> xSize >> ySize >> iterations >> num_threads;
  uint8_t* inputFrame = new uint8_t[xSize*ySize];
  uint8_t* outputFrame = new uint8_t[xSize*ySize];
  for (uint32_t y = 0; y < ySize; y++)
  {
    for (uint32_t x = 0; x < xSize; x++)
    {
      input >> inputFrame[y*xSize + x];
      inputFrame[y*xSize + x] -= '0';
    }
  }


  // Calculation
 outputFrame = calc(xSize, ySize, iterations, num_threads, inputFrame);

  // Write result
  for (uint32_t y = 0; y < ySize; y++)
  {
    for (uint32_t x = 0; x < xSize; x++)
    {
      outputFrame[y*xSize + x] += '0';
      output << " " << outputFrame[y*xSize + x];
    }
    output << "\n";
  }
#if 0
  // Prepare to exit
  delete outputFrame;
  delete inputFrame;
#endif  
  output.close();
  input.close();
  return 0;
}
