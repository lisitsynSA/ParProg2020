#include <stdio.h>
#include <omp.h>
#include <unistd.h>

#define TIME(code) {\
    double start = omp_get_wtime();\
    {code}\
    printf("[TIME %s] %lf\n\n", #code, omp_get_wtime() - start);\
}

// + barrier
void timeTest(int numThreads)
{
    #pragma omp parallel num_threads(numThreads)
    {
        printf("[%d] SLEEP\n", omp_get_thread_num());
        sleep(omp_get_thread_num());
        #pragma omp barrier
        printf("[%d] WAKE UP\n", omp_get_thread_num());
    }
}

void privateTest(int numThreads)
{
    int val = 13;
    printf("[%d] \t\t val = %d \t\t%d\n", omp_get_thread_num(), val, &val);
    #pragma omp parallel num_threads(numThreads) private(val)
    {
        printf("[%d] SLEEP \t val = %d \t%d\n", omp_get_thread_num(), val, &val);
        val = omp_get_thread_num();
        sleep(omp_get_thread_num());
        printf("[%d] WAKE UP \t val = %d \t%d\n", omp_get_thread_num(), val, &val);
    }
    printf("[%d] \t\t val = %d \t\t%d\n", omp_get_thread_num(), val, &val);
}

int create_threads(int numThreads)
{
    int res = 0;
    #pragma omp parallel num_threads(numThreads)
    {
        res++;
        //printf("[%d]\n", omp_get_thread_num());
    }
    return res;
}

void creationTest(int numThreads, int numSubThreads)
{
    int res = 0;
    #pragma omp parallel num_threads(numThreads)
    {
        for (int i = 0; i < 100; i++)
            res += create_threads(numSubThreads);
    }

}

void raceConditionTest(int numThreads, int iterations)
{
    int res = 0;
    #pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < iterations; i++)
        res++;
    printf("Iterations %d result %d\n", iterations, res);
}

void raceConditionTest_critical(int numThreads, int iterations)
{
    int res = 0;
    #pragma omp parallel num_threads(numThreads)
    {
        //....
        #pragma omp for
        for (int i = 0; i < iterations; i++)
        {
            //value = func();
            #pragma omp critical
            {
                res++;
            }
        }
        //....
    }

    printf("Iterations %d result %d\n", iterations, res);
}

void raceConditionTest_atomic(int numThreads, int iterations)
{
    int res = 0;
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp for
        for (int i = 0; i < iterations; i++)
        #pragma omp atomic
            res++;
    }
    printf("Iterations %d result %d\n", iterations, res);
}

// + firstprivate
void raceConditionTest_atomicSplit(int numThreads, int iterations)
{
    int res = 0;
    int * res_p = &res;
    #pragma omp parallel num_threads(numThreads) firstprivate(res)
    {
        //res = 0; use first private
        #pragma omp for
        for (int i = 0; i < iterations; i++)
            res++;
        #pragma omp atomic
        *res_p += res;
        //printf("[%d] %d_(%d)\n", omp_get_thread_num(), res, &res);
    }
    printf("Iterations %d result %d_(%d)\n", iterations, res, &res);
}

void raceConditionTest_reduction(int numThreads, int iterations)
{
    int res = 0;
    #pragma omp parallel for num_threads(numThreads) reduction(+:res)
    for (int i = 0; i < iterations; i++)
        res++;
    printf("Iterations %d result %d\n", iterations, res);
}

void lastPrivateTest(int numThreads, int iterations)
{
    int res = 0;
    #pragma omp parallel for num_threads(numThreads) lastprivate(res)
    for (int i = 0; i < iterations; i++)
    {
        sleep(numThreads - omp_get_thread_num());
        res = i*i;
    }
    printf("Iterations %d result %d\n", iterations, res);
}

void sectionsTest(int numThreads)
{
    #pragma omp parallel sections num_threads(numThreads)
    {
        //#pragma omp section
        printf("[%d] section0 \n", omp_get_thread_num());
        #pragma omp section
        {
            printf("[%d] section1 \n", omp_get_thread_num());
        }
        #pragma omp section
        {
            printf("[%d] section2 \n", omp_get_thread_num());
        }
    }
}

//singleTest(0);
void singleTest(int numThreads)
{
    #pragma omp parallel num_threads(numThreads)
    {
        printf("[%d] all\n", omp_get_thread_num());
        #pragma omp single
        {
            printf("[%d] single\n", omp_get_thread_num());
            if (numThreads - 1 != 0)
            {
                singleTest(numThreads - 1);
            }
        }
    }
}

void zeroThreadsTest(int numThreads)
{
    omp_set_num_threads(8);
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp single
        printf("NUM %d threads %d\n", numThreads, omp_get_num_threads());
    }
}

void ifTest(int numThreads)
{
    #pragma omp parallel num_threads(numThreads) if (numThreads<4)
    {
        printf("[%d]\n", omp_get_thread_num());
    }
}

void orderedTest(int numThreads)
{
    #pragma omp parallel for ordered num_threads(numThreads)
    for (int i = 0; i < 10; i++)
    {
        printf("(%d) %d\n", omp_get_thread_num(), i);
        #pragma omp ordered
        printf("[%d] %d\n", omp_get_thread_num(), i);
    }
}

void sheduleTest(int numThreads)
{
    //(static, 4)
    //(dynamic, 1)
    //(guided, 1)
    #pragma omp parallel for num_threads(numThreads) schedule(guided, 1)
    for (int i = 0; i < 10; i++)
    {
        printf("[%d] %d\n", omp_get_thread_num(), i);
    }
}

int main()
{
    omp_set_nested(1);
    //TIME(timeTest(4);)
    //TIME(privateTest(4);)
    //TIME(creationTest(4,10);)
    //TIME(raceConditionTest(40, 1000000);)
    //TIME(raceConditionTest_critical(40, 1000000);)
    //TIME(raceConditionTest_atomic(40, 1000000);)
    //TIME(raceConditionTest_atomicSplit(40, 1000000);)
    //TIME(raceConditionTest_reduction(40, 1000000);)
    //lastPrivateTest(4, 4);
    //sectionsTest(4);
    //singleTest(4);
    //singleTest(0);
    //zeroThreadsTest(0);
    //ifTest(5);
    //TIME(orderedTest(4);)
    sheduleTest(4);
    return 0;
}
