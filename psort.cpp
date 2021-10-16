#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "psort.h"
#include <mpi.h>
using namespace std;

// #define createRandomData

void pSort::close(){}
void pSort::init(){}

int compareData(pSort::dataType d1,pSort::dataType d2){
    if (d1.a<d2.a){
        return 0;
    }
    else if(d2.a<d1.a){
        return 1;
    }
    if (d1.b<d2.b){
        return 0;
    }
    else if(d2.b<d1.b){
        return 1;
    }
    if(d1.x<=d2.x){
        return 0;
    }
    return 1;
}

void pSort::sort(dataType *data, int32_t n){
    int nProcs; 
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank

    int *all_counts = new int(nProcs)
    cout<<"\nsort called with rank="<<myRank<<"\tn="<<n; 
    MPI_Allgather( &n, 1, MPI_INT, all_counts, 1, MPI_INT, MPI_COMM_WORLD);
    printf("\n%d",nProcs);
    for(int i=0;i<nProcs;i++)
        cout<<"\t"<<all_counts[i];
}

