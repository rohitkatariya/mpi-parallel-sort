#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "psort.h"
/* includes MPI library code specs */
#define MAXSIZE 1000
// #define NUMDATA 1294967296
#define NUMDATA 20000000
// #define MAXSIZE RAND_MAX
using namespace std;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 
    
    
    long num_data;
    srand(time(0));
    int this_incr = 2000+rand()%100;
    #ifdef DEBUG
    if(myRank==0)
        printf("this_incr:%d",this_incr);
    #endif
    // this_incr=2028;
    // this_incr=2056;
    // this_incr=2034;
    // this_incr=2088;
    // this_incr = 2079;
    // this_incr = 2025;
    
    
    srand (myRank*10+this_incr);//+time(0));
    // srand (myRank*10+time(0));
    num_data= NUMDATA/nProcs;
    printf("\nnum_data %ld",num_data);
    // num_data=30;
    pSort::dataType *data = new pSort::dataType[num_data];
    // cout<<"\nrank"<<myRank<<"\t";
    for(int i=0;i<num_data;i++){
        data[i].a=97+rand()%26;
        data[i].b=97+rand()%26;
        data[i].x=rand()%MAXSIZE;
        // Set payload here 
        for (int v_i=0 ;v_i<58;v_i++){
            data[i].payload[v_i]=rand()%26+97;
        }
    }
    
    pSort this_p;
    this_p.sort(data,int32_t(num_data));
    delete[] data;
    #ifdef DEBUG
    printf("\nMPI_Finalize");
    #endif
    MPI_Finalize();
}
