#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "psort.h"
#include <mpi.h>
#include <stddef.h>
#include <fstream>
using namespace std;

// #define createRandomData

#define DEBUG
typedef struct {
    uint32_t x;
    char a, b;
  } pivotStruct;
void pSort::close(){}
void pSort::init(){}

int compareData(pSort::dataType d1,pivotStruct d2){
    if (d1.a<d2.a){
        return 1;
    }
    else if(d2.a<d1.a){
        return 0;
    }
    if (d1.b<d2.b){
        return 1;
    }
    else if(d2.b<d1.b){
        return 0;
    }
    if(d1.x<=d2.x){
        return 1;
    }
    return 0;
}

int get_proc_from_ele_idx(int ele_idx, int *all_counts,int nProcs){
    int proc_num=0;
    for(int i=0;i<nProcs;i++){
        if(ele_idx<all_counts[i])
            return i;
        ele_idx-=all_counts[i];
    }
    return -1;
}

void sQuickSort(pSort::dataType *data, int *all_counts, int start, int end){
    
}

int internalPivoting(pSort::dataType *data,int start, int end,pivotStruct pivot_this){
    int i=start-1;
    int j=i;
    int count_less=0;
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    printf("\ninternal pivoting (rank,start,end), (%d,%d,%d)",myRank,start,end);
    for( int j=i+1;j<=end;j++){
        if( compareData(data[j],pivot_this)==1 ){
            i++;
            count_less++;
            // swap data[i] data[j]
            pSort::dataType temp = data[j];
            data[j]=data[i];
            data[i]=temp;
        }
    }
    return count_less;
}

void pquickSort(pSort::dataType *data, int *all_counts, int *all_offsets, int start, int end,int tree_level){
    int nProcs; 
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank
    
    int startproc=get_proc_from_ele_idx(start,all_counts,nProcs);
    int endproc=get_proc_from_ele_idx(end,all_counts,nProcs);
    if(startproc==endproc){
        // run sequential quicksort
        // change indexes
        sQuickSort(data, all_counts, start, end);
        cout<<"start proc= endproc";
        return;
    }

    MPI_Datatype PivotDataType_MPI;
    MPI_Datatype PivotDataType_T[3] = {MPI_CHAR, MPI_CHAR, MPI_UNSIGNED};
    int PivotDataType_B[3]  = {1,1,1};//block lengths
    // MPI_Aint PivotDataType_D[3]  = {0,1,2};//offsets
    MPI_Aint PivotDataType_D[3]  = {offsetof(pivotStruct, b), offsetof(pivotStruct, a), offsetof(pivotStruct, x)};//offsets
    MPI_Type_create_struct(3, PivotDataType_B, PivotDataType_D, PivotDataType_T, &PivotDataType_MPI);
    MPI_Type_commit(&PivotDataType_MPI);
    
    
    pivotStruct pivot_this;  
    MPI_Request *all_send_requests= new MPI_Request[nProcs];
    if (myRank==startproc){
        /*
         MPI_Comm com;
         printf("\ncom:%d\t%d",MPI_COMM_WORLD,com);
        */
        // Send pivot element to all other processors
        int curr_starting_idx = max(0,start-all_offsets[myRank]);
        pivot_this.a= data[curr_starting_idx].a;
        pivot_this.b= data[curr_starting_idx].b;
        pivot_this.x= data[curr_starting_idx].x;
        for( int i=startproc+1; i<=endproc;i++){
            MPI_Isend(&pivot_this, 1, PivotDataType_MPI, i, tree_level*10, MPI_COMM_WORLD,all_send_requests+i);
        }
    }
    else{
        MPI_Status recvStatus;
        MPI_Recv(&pivot_this, 1, PivotDataType_MPI, startproc, tree_level*10, MPI_COMM_WORLD,&recvStatus);
    }
    
    printf("\n%d\tafter send/recv, pivot struct val(a,b,x)=(%c,%c,%d)",myRank,pivot_this.a,pivot_this.b,pivot_this.x);

    // do internal pivoting
    // call internalPivoting after offsetting start end indexes
    int num_lt_pivot = 0;
    if(myRank!=startproc)
        num_lt_pivot = internalPivoting(data,0,min(end-all_offsets[myRank],all_counts[myRank]-1),pivot_this);
    else
        num_lt_pivot = internalPivoting(data,start-all_offsets[myRank],min(end-all_offsets[myRank],all_counts[myRank]-1),pivot_this);
    
    #ifdef DEBUG
        ofstream fout;
        fout.open("output_dir/ipvt_"+ to_string(myRank)+".txt");
        // fprintf(fout,"pivot (%c,%c,%d)\n",pivot_this.a,pivot_this.b,pivot_this.x);
        fout<<"\npivot:("<<pivot_this.a-97<<","<<pivot_this.b<<","<<pivot_this.x<<",nltp:"<<num_lt_pivot<<")\n";
        for(int i=max(0,start-all_offsets[myRank]); i<=min(end-all_offsets[myRank],all_counts[myRank]-1) ;i++){
            fout<<"\t("<<data[i].a-97<<","<<data[i].b<<","<<data[i].x<<","<<i <<")";
            if(i%10==9)
                fout<<"\n"<<myRank;
        }
        fout.close();
    #endif
    // send receive number of elements less than pivot
    int *num_lt_pivot_arr = new int(nProcs) ;
    num_lt_pivot_arr[myRank] = num_lt_pivot;
    MPI_Status status_this;
    MPI_Status recvStatus;
    for(int i=startproc;i<=endproc;i++){
        if(i==myRank)
            continue;
        if(myRank == startproc)
            MPI_Wait(all_send_requests+i,&status_this);
        MPI_Isend(&num_lt_pivot, 1, MPI_INT, i, tree_level*10, MPI_COMM_WORLD,all_send_requests+i-startproc);
    }
    
    for(int i=startproc;i<=endproc;i++){
        if(i==myRank)
            continue;
        MPI_Recv(num_lt_pivot_arr+i, 1, MPI_INT, i, tree_level*10, MPI_COMM_WORLD,&recvStatus);
    }
    

    // calculate pivot location
    int pivotLocation = 0;
    for( int i=startproc;i<=endproc;i++){
        pivotLocation+=num_lt_pivot_arr[i];
    }
    
    int pivotProc=get_proc_from_ele_idx(pivotLocation,all_counts,nProcs);
    printf("\npivot Location Global, pivot proc: %d, %d",pivotLocation,pivotProc);
    int start_idx_local=max(0,start-all_offsets[myRank]);
    int end_idx_local=max(all_offsets[myRank]-1,end-all_offsets[myRank]);

    MPI_Datatype myDataType_MPI;
    MPI_Datatype myDataType_T[4] = {MPI_CHAR, MPI_CHAR, MPI_UNSIGNED, MPI_INT8_T};
    int myDataType_B[4]  = {1,1,1,58};//block lengths
    // MPI_Aint PivotDataType_D[3]  = {0,1,2};//offsets
    MPI_Aint myDataType_D[4]  = {offsetof(pSort::dataType, b), offsetof(pSort::dataType, a), offsetof(pSort::dataType, x),offsetof(pSort::dataType, payload)};//offsets
    MPI_Type_create_struct(4, myDataType_B, myDataType_D, myDataType_T, &myDataType_MPI);
    MPI_Type_commit(&myDataType_MPI);
    
    MPI_Request pivot_sendRequest;
    if(myRank==startproc){
        MPI_Isend(data+start-all_offsets[myRank], 1, myDataType_MPI, pivotProc, tree_level*10, MPI_COMM_WORLD,&pivot_sendRequest);
        pSort::dataType pivotData_sent= data[start-all_offsets[myRank]];
        // printf("\npivot data(%d) sent(%c,%c,%d,(%d,%d)) from %d",pivotLocation,pivotData_sent.a,pivotData_sent.b,pivotData_sent.x,pivotData_sent.payload[0],pivotData_sent.payload[1],myRank);
    }
    pSort::dataType pivotData_received;
    if(myRank==pivotProc){
        MPI_Recv(&pivotData_received, 1, myDataType_MPI, startproc, tree_level*10, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // printf("\npivot data received(%c,%c,%d,(%d,%d)) at %d",pivotData_received.a,pivotData_received.b,pivotData_received.x,pivotData_received.payload[0],pivotData_received.payload[1],myRank);
    }

    int right_procNo = endproc;
    int left_procNo = startproc;
    int num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo]; 
    int num_swap_ele_left_procNo = all_counts[left_procNo]-(start-all_offsets[left_procNo])-num_lt_pivot_arr[left_procNo];
    int send_reqs_sent = 0;
    // MPI_Status *recvStatus_list = new MPI_Status(nProcs);
    while(right_procNo>pivotProc && left_procNo<pivotProc){
        if(num_swap_ele_right_procNo<0 or num_swap_ele_left_procNo<0){
            printf("\n##################count less than 0 :(%d,%d)",num_swap_ele_right_procNo,num_swap_ele_left_procNo);
        }
        if(num_swap_ele_right_procNo>num_swap_ele_left_procNo){
            //swap num_swap_ele_left_procNo elements
            
            if(myRank==right_procNo){//this is right processor
                // send and receive
                MPI_Sendrecv_replace(data+num_swap_ele_right_procNo-num_swap_ele_left_procNo, num_swap_ele_left_procNo, myDataType_MPI,
                         left_procNo, tree_level*10, left_procNo, tree_level*10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Sendrecv(data+num_swap_ele_right_procNo-num_swap_ele_left_procNo, num_swap_ele_left_procNo, myDataType_MPI,
                //  left_procNo, tree_level*10,
                //  data+num_swap_ele_right_procNo-num_swap_ele_left_procNo, num_swap_ele_left_procNo, myDataType_MPI,
                //  left_procNo, tree_level*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send_reqs_sent++;
                
            }
            if(myRank==left_procNo){ //this is the left processor
                // send and receive

                MPI_Sendrecv_replace(data+num_lt_pivot_arr[myRank], num_swap_ele_left_procNo, myDataType_MPI,
                         right_procNo, tree_level*10, right_procNo, tree_level*10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // MPI_Sendrecv(data+num_lt_pivot_arr[myRank], num_swap_ele_left_procNo, myDataType_MPI,
                //  right_procNo, tree_level*10,
                //  data+num_lt_pivot_arr[myRank], num_swap_ele_left_procNo, myDataType_MPI,
                //  right_procNo, tree_level*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send_reqs_sent++;
            }
            num_swap_ele_right_procNo=num_swap_ele_right_procNo-num_swap_ele_left_procNo;
            printf("\n###%d completed processor %d\n",myRank,left_procNo);
            left_procNo = left_procNo+1;
            num_swap_ele_left_procNo = all_counts[left_procNo]-num_lt_pivot_arr[left_procNo];
        }else if(num_swap_ele_right_procNo==num_swap_ele_left_procNo){
            if(myRank==right_procNo){//this is right processor
                // send and receive
                MPI_Sendrecv_replace(data, num_swap_ele_left_procNo, myDataType_MPI,
                         left_procNo, tree_level*10, left_procNo, tree_level*10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Sendrecv(data, num_swap_ele_left_procNo, myDataType_MPI,
                //  left_procNo, tree_level*10,
                //  data, num_swap_ele_left_procNo, myDataType_MPI,
                //  left_procNo, tree_level*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send_reqs_sent++;
            }
            if(myRank==left_procNo){ //this is the left processor
                // send and receive
                MPI_Sendrecv_replace(data+num_lt_pivot_arr[myRank], num_swap_ele_left_procNo, myDataType_MPI,
                         right_procNo, tree_level*10, right_procNo, tree_level*10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Sendrecv(data+num_lt_pivot_arr[myRank], num_swap_ele_left_procNo, myDataType_MPI,
                //  right_procNo, tree_level*10,
                //  data+num_lt_pivot_arr[myRank], num_swap_ele_left_procNo, myDataType_MPI,
                //  right_procNo, tree_level*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send_reqs_sent++;
            }
            printf("\n###%d completed processor %d,%d\n",myRank,left_procNo,right_procNo);
            left_procNo = left_procNo+1;
            right_procNo = right_procNo-1;
            
            num_swap_ele_left_procNo = all_counts[left_procNo]-num_lt_pivot_arr[left_procNo];
            num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo];
        }else if(num_swap_ele_right_procNo<num_swap_ele_left_procNo){
            if(myRank==right_procNo){//this is right processor
                // send and receive
                MPI_Sendrecv_replace(data, num_swap_ele_right_procNo, myDataType_MPI,
                         left_procNo, tree_level*10, left_procNo, tree_level*10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Sendrecv(data, num_swap_ele_right_procNo, myDataType_MPI,
                //  left_procNo, tree_level*10,
                //  data, num_swap_ele_right_procNo, myDataType_MPI,
                //  left_procNo, tree_level*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send_reqs_sent++;
            }
            if(myRank==left_procNo){ //this is the left processor
                // send and receive
                MPI_Sendrecv_replace(data+num_lt_pivot_arr[myRank]+num_swap_ele_left_procNo-num_swap_ele_right_procNo, num_swap_ele_right_procNo, myDataType_MPI,
                         right_procNo, tree_level*10, right_procNo, tree_level*10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // MPI_Sendrecv(data+num_lt_pivot_arr[myRank]+num_swap_ele_left_procNo-num_swap_ele_right_procNo, num_swap_ele_right_procNo, myDataType_MPI,
                //  right_procNo, tree_level*10,
                //  data+num_lt_pivot_arr[myRank]+num_swap_ele_left_procNo-num_swap_ele_right_procNo, num_swap_ele_right_procNo, myDataType_MPI,
                //  right_procNo, tree_level*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send_reqs_sent++;
            }
            // printf("\n###%d completed processor %d\n",myRank,right_procNo);
            right_procNo = right_procNo-1;
            // printf("$$$$%d-%d",num_swap_ele_left_procNo,num_swap_ele_right_procNo);
            num_swap_ele_left_procNo=num_swap_ele_left_procNo-num_swap_ele_right_procNo;
            num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo];
        }
    }
    if(right_procNo>pivotProc){
        //right not complete

    }elseif(left_procNo<pivotProc){
        //left not complete
    }
    #ifdef DEBUG
        // ofstream fout;
        fout.open("output_dir/pivoted_"+ to_string(myRank)+".txt");
        // fprintf(fout,"pivot (%c,%c,%d)\n",pivot_this.a,pivot_this.b,pivot_this.x);
        fout<<"\npivot:("<<(pivot_this.a-97)<<","<<pivot_this.b<<","<<pivot_this.x<<",nltp:"<<num_lt_pivot<<")\n";
        for(int i=max(0,start-all_offsets[myRank]); i<=min(end-all_offsets[myRank],all_counts[myRank]-1) ;i++){
            fout<<"\t("<<(data[i].a-97)<<","<<data[i].b<<","<<data[i].x<<","<<i <<")";
            if(i%10==9)
                fout<<"\n"<<myRank;
        }
        fout.close();
    #endif
    // for(int i=0;i<send_reqs_sent;i++){
    //     MPI_Wait(recvStatus+send_reqs_sent)
    // }
    // if(myRank==pivotProc){
    //     int pivot_idx_local = pivotLocation-all_offsets[myRank];

    // }
    // else if(myRank<pivotProc){
    //     // send all the right elements
    //     numberOfElementsToTransmit = (end_idx_local+1)-start_idx_local-num_lt_pivot
    //     numberOfElementsToReceive =  

    //     for(int i=max(0,start-all_offsets[myRank])+num_lt_pivot;i<=min(all_counts[myRank]-1 ,end-all_offsets[myRank]);i++){
    //         // send this element to a right processor
            
    //     }
    //     // receive all right elements

    // }
    // else{
    //     // send all the left elements
        
    //     // receive all left elements

    // }
    // #ifdef DEBUG
    //     fout.open("output_dir/lpvt_"+ to_string(myRank)+".txt");
    //     fout<<"\nat rank:"<<myRank;
    //     for(int i=startproc;i<=endproc;i++){
    //         fout<<"\t"<<num_lt_pivot_arr[i];
    //     }
    //     fout<<"\npivotLocation:"<<pivotLocation;
    //     fout.close();
    // #endif
    // pass elements between processors
    
    // Check this processor is left/right to pivot processor or this is the pivot processor
    
    // sort left and right arrays

    

}

void pSort::sort(dataType *data, int32_t n){
    int nProcs; 
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank

    // Gathering information of number of elements at each processor 
    int *all_counts = new int(nProcs);
    MPI_Allgather( &n, 1, MPI_INT, all_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    // Counting total number of elements to be sorted and creating an offset array useful for converting global index to processor's local index
    int num_ele_sort = 0;
    int *all_offsets=new int(nProcs);
    #ifdef DEBUG
    if(myRank==0)
        cout<<"\nPrinting sizes of arrays at each processor:";
    #endif
    for( int i =0;i<nProcs;i++){
        all_offsets[i]=num_ele_sort;
        num_ele_sort+=all_counts[i];
        #ifdef DEBUG
        if(myRank==0)
            cout<<all_counts[i]<<"\t";
        #endif
    }


    MPI_Barrier(MPI_COMM_WORLD);
    pquickSort(data,all_counts,all_offsets,0,num_ele_sort-1,0);

    // #ifdef DEBUG
    //     cout<<"\nsort called with rank="<<myRank<<"\tn="<<n; 
    //     printf("\n%d",nProcs);
    //     for(int i=0;i<nProcs;i++)
    //         cout<<"\t"<<all_counts[i];
    // #endif
}

