#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "psort.h"
#include <mpi.h>
#include <stddef.h>
#include <sstream>
#include <fstream>
using namespace std;

// #define createRandomData

// #define DEBUG
//#define DEBUGOUT
typedef struct {
    uint32_t x;
    char a, b;
  } pivotStruct;
void pSort::close(){}
void pSort::init(){}
template<class T>
int compareData(pSort::dataType d1,T d2){
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
    if(d1.x<d2.x){ 
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


string print_arr(int *arr,int len_arr){
    
    std::ostringstream os;
    for (int i=0;i<len_arr;i++) {
        os << i<<",";
    }
 
    std::string str(os.str());
    return str;
}

string print_arr2(pSort::dataType *arr,int len_arr){
    
    std::ostringstream os;
    os<<"[";
    for (int i=0;i<len_arr;i++) {
        os << arr[i].a-97<<",";
    }
    os<<"]";
    std::string str(os.str());
    return str;
}

template<class L>
int internalPivoting(pSort::dataType *data,int start, int end,L pivot_this){
    int i=start-1;
    int count_less=0;
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // printf("\n %d \t%s",myRank, print_arr2(data+start,end-start+1).c_str());
    // printf("\ninternal pivoting (rank,start,end), (%d,%d,%d)",myRank,start,end);
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

void sQuickSort(pSort::dataType *data, int start, int end){
    // return;
    if(start>=end)
        return;
    // printf("\nstart pivot(%d, %d)",start,end);
    int pivot_loc = start+internalPivoting(data,start+1,end,data[start]);
    pSort::dataType temp = data[pivot_loc];
    data[pivot_loc] = data[start];
    data[start]=temp;
    // printf("\ndone pivot:%d, start,end (%d,%d)",pivot_loc,start,end);
    // printf("\array after p:%s",print_arr2(data,end-start+1).c_str());
    // return;
    if(pivot_loc>start){
        sQuickSort(data,start,pivot_loc-1);
    }
    if(pivot_loc<end){
        sQuickSort(data,pivot_loc+1,end);
    }
}


void pquickSort(pSort::dataType *data, int *all_counts, int *all_offsets, int start, int end,int tree_level){
    int nProcs; 
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank
    
    int startproc=get_proc_from_ele_idx(start,all_counts,nProcs);
    int endproc=get_proc_from_ele_idx(end,all_counts,nProcs);
    
    if(myRank<startproc || myRank>endproc){
        return;
    }
    #ifdef DEBUG
    if(myRank==startproc)
        printf("\nstarting:(start,end):(%d,%d) (%d, %d)",start,end,startproc,endproc );
    #endif
    if(startproc==endproc){
        // change indexes
        if(end==start)
            return;
        // printf("\nqsStart: \t%s",print_arr2(data+start-all_offsets[myRank],end-start+1).c_str());
        sQuickSort(data, start-all_offsets[myRank], end-all_offsets[myRank]);
        // printf("\nqsEnd: \t%s",print_arr2(data+start-all_offsets[myRank],end-start+1).c_str());
        return;
    }
    
    MPI_Datatype PivotDataType_MPI;
    MPI_Datatype PivotDataType_T[3] = {MPI_CHAR, MPI_CHAR, MPI_UNSIGNED};
    int PivotDataType_B[3]  = {1,1,1};//block lengths
    MPI_Aint PivotDataType_D[3]  = {offsetof(pivotStruct, b), offsetof(pivotStruct, a), offsetof(pivotStruct, x)};//offsets
    MPI_Type_create_struct(3, PivotDataType_B, PivotDataType_D, PivotDataType_T, &PivotDataType_MPI);
    MPI_Type_commit(&PivotDataType_MPI);
    
    
    pivotStruct pivot_this;  
    MPI_Request *all_send_requests= new MPI_Request[nProcs];
    if (myRank==startproc){
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
    
    // do internal pivoting
    // call internalPivoting after offsetting start end indexes
    
    int num_lt_pivot = 0;
    if(myRank!=startproc)
        num_lt_pivot = internalPivoting(data,0,min(end-all_offsets[myRank],all_counts[myRank]-1),pivot_this);
    else
        num_lt_pivot = internalPivoting(data,start-all_offsets[myRank]+1,min(end-all_offsets[myRank],all_counts[myRank]-1),pivot_this);
    
    #ifdef DEBUG
        ofstream fout;
        fout.open("output_dir/ipvt_"+ to_string(myRank)+".txt");
        // fprintf(fout,"pivot (%c,%c,%d)\n",pivot_this.a,pivot_this.b,pivot_this.x);
        if(myRank==0)
        fout<<"\npivot:("<<pivot_this.a-97<<","<<pivot_this.b<<","<<pivot_this.x<<",nltp:"<<num_lt_pivot<<")\n";
        fout<<"\n";
        for(int i=max(0,start-all_offsets[myRank]); i<=min(end-all_offsets[myRank],all_counts[myRank]-1) ;i++){
            fout<<"\t("<<data[i].a-97<<","<<data[i].b-97;//<<","<<data[i].x<<","<<i <<")";
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
    #ifdef DEBUG
        if(myRank==startproc)
            printf("\nnum_lt_pivot_arr %s",print_arr(num_lt_pivot_arr+startproc,endproc-startproc+1).c_str());
    #endif
    // calculate pivot location
    int pivotLocation = start;
    for( int i=startproc;i<=endproc;i++){
        pivotLocation+=num_lt_pivot_arr[i];
    }
    #ifdef DEBUG
    if(myRank==startproc)
        printf("\npivot (%d,%d), pivotLocation %d",pivot_this.a-97,pivot_this.b-97,pivotLocation);
    #endif
    int pivotProc=get_proc_from_ele_idx(pivotLocation,all_counts,nProcs);
    
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
    }
    pSort::dataType pivotData_received;
    if(myRank==pivotProc){
        MPI_Recv(&pivotData_received, 1, myDataType_MPI, startproc, tree_level*10, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    int right_procNo = endproc;
    int left_procNo = startproc;
    int send_reqs_sent = 0;

    int num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo]; 
    int num_swap_ele_left_procNo = all_counts[left_procNo]-(start-all_offsets[left_procNo])-num_lt_pivot_arr[left_procNo]-1;//subtracting 1 because pivot is in startproc
    int right_procNo_data_start = num_lt_pivot_arr[right_procNo];
    int left_procNo_data_start = (start-all_offsets[left_procNo]) + num_lt_pivot_arr[left_procNo]+1; //+1 because of pivot in startproc
    // printf("\nprocno %d,%d,%d",left_procNo,pivotProc,right_procNo);
    while((right_procNo>pivotProc) || (left_procNo<pivotProc)){
        // if(myRank==right_procNo){
        //     printf("\n%d(%d),%d,%d(%d)",left_procNo,num_swap_ele_left_procNo,pivotProc,right_procNo,num_swap_ele_right_procNo);
        // }
        if(num_swap_ele_right_procNo>num_swap_ele_left_procNo){
            //swap num_swap_ele_left_procNo elements
            if(num_swap_ele_left_procNo!=0){
                if(myRank==right_procNo){//this is right processor
                    // send and receive
                    MPI_Sendrecv_replace(data+right_procNo_data_start-num_swap_ele_left_procNo, num_swap_ele_left_procNo, myDataType_MPI,
                            left_procNo, tree_level*10, left_procNo, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                }
                if(myRank==left_procNo){ //this is the left processor
                    // send and receive
                    MPI_Sendrecv_replace(data+left_procNo_data_start, num_swap_ele_left_procNo, myDataType_MPI, 
                            right_procNo, tree_level*10, right_procNo, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                    // printf("\n###%d completed processor %d\n",myRank,left_procNo);
                }
                right_procNo_data_start-=num_swap_ele_left_procNo;
                left_procNo_data_start+=num_swap_ele_left_procNo;
            }
            num_swap_ele_right_procNo=num_swap_ele_right_procNo-num_swap_ele_left_procNo;
            left_procNo = left_procNo+1;
            left_procNo_data_start = num_lt_pivot_arr[left_procNo];
            num_swap_ele_left_procNo = all_counts[left_procNo]-num_lt_pivot_arr[left_procNo];
        }else if(num_swap_ele_right_procNo==num_swap_ele_left_procNo){
            if(num_swap_ele_left_procNo>0){
                if(myRank==right_procNo){//this is right processor
                    // send and receive
                    MPI_Sendrecv_replace(data+right_procNo_data_start-num_swap_ele_left_procNo, num_swap_ele_left_procNo, myDataType_MPI,
                            left_procNo, tree_level*10, left_procNo, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                }
                if(myRank==left_procNo){ //this is the left processor
                    // send and receive
                    MPI_Sendrecv_replace(data+left_procNo_data_start, num_swap_ele_left_procNo, myDataType_MPI,
                            right_procNo, tree_level*10, right_procNo, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                // printf("\n###%d completed processor %d,%d\n",myRank,left_procNo,right_procNo);
                }
            }
            left_procNo = left_procNo+1;
            right_procNo = right_procNo-1;
            left_procNo_data_start = num_lt_pivot_arr[left_procNo];
            right_procNo_data_start = num_lt_pivot_arr[right_procNo];
            num_swap_ele_left_procNo = all_counts[left_procNo]-num_lt_pivot_arr[left_procNo];
            num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo];
        }else if(num_swap_ele_right_procNo<num_swap_ele_left_procNo){
            if(num_swap_ele_right_procNo>0){
                if(myRank==right_procNo){//this is right processor
                    // send and receive
                    MPI_Sendrecv_replace(data+right_procNo_data_start-num_swap_ele_right_procNo, num_swap_ele_right_procNo, myDataType_MPI,
                            left_procNo, tree_level*10, left_procNo, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                // printf("\n###%d completed processor %d\n",myRank,right_procNo);
                }
                if(myRank==left_procNo){ //this is the left processor
                    // send and receive
                    MPI_Sendrecv_replace(data+left_procNo_data_start, num_swap_ele_right_procNo, myDataType_MPI,
                            right_procNo, tree_level*10, right_procNo, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                }
                right_procNo_data_start-=num_swap_ele_right_procNo;
                left_procNo_data_start+=num_swap_ele_right_procNo;
            }
            num_swap_ele_left_procNo=num_swap_ele_left_procNo-num_swap_ele_right_procNo;
            right_procNo = right_procNo-1;
            num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo];
            right_procNo_data_start = num_lt_pivot_arr[right_procNo];
            
        }
    }
    
    delete[] num_lt_pivot_arr;
    //replace pivot 
    
    if(true){
        if(start==pivotLocation){
            ;
        }else if(pivotProc==startproc ){
            if( myRank==pivotProc){
                pSort::dataType temp = data[start-all_offsets[pivotProc]];
                data[start-all_offsets[pivotProc]] = data[pivotLocation-all_offsets[pivotProc]];
                data[pivotLocation-all_offsets[pivotProc]]=temp;
                // printf("\nmanual swap done %d",myRank);
            }
        }else{
            if(myRank==startproc){
                MPI_Sendrecv_replace(data+ (start-all_offsets[startproc]) , 1, myDataType_MPI,
                            pivotProc, tree_level*10, pivotProc, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if(myRank==pivotProc){
                MPI_Sendrecv_replace(data+(pivotLocation-all_offsets[pivotProc]), 1, myDataType_MPI,
                            startproc, tree_level*10, startproc, tree_level*10,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    
    if(pivotLocation>start){
        // printf("\ncalling on left: (%d,%d) ",start,pivotLocation-1);
        // if(tree_level<=1)
        pquickSort(data,all_counts,all_offsets,start,pivotLocation-1,tree_level+1);
    }
    if(pivotLocation<end){
        // printf("\ncalling on right: (%d,%d) ",pivotLocation+1,end);
        // if(true or tree_level<=1)
            pquickSort(data,all_counts,all_offsets,pivotLocation+1,end,tree_level+1);
        
    }
    
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
    
    for( int i =0;i<nProcs;i++){
        all_offsets[i]=num_ele_sort;
        num_ele_sort+=all_counts[i];
    }

    #ifdef DEBUG
        ofstream fout;
        fout.open("output_dir/inp_"+ to_string(myRank)+".txt");
        fout<<"\n"<<myRank;
        for(int i=0; i<=all_counts[myRank]-1 ;i++){
            // if(i%10==0)
            //     fout<<"\n"<<myRank;
            fout<<"\t"<<data[i].a-97<<","<<data[i].b-97;//<<","<<data[i].x<<","<<i <<")";
        }
        fout.close();
    #endif
    MPI_Barrier(MPI_COMM_WORLD);
    pquickSort(data,all_counts,all_offsets,0,num_ele_sort-1,0);
    //MPI_Barrier(MPI_COMM_WORLD);
    #ifdef DEBUGOUT
        printf("writing output to output_dir/out_%d.txt",myRank);
        ofstream fout;
        fout.open("output_dir/out_"+ to_string(myRank)+".txt");
        // fprintf(fout,"pivot (%c,%c,%d)\n",pivot_this.a,pivot_this.b,pivot_this.x);
        fout<<"\n";
        for(int i=0; i<=all_counts[myRank]-1 ;i++){
            // fout<<"\t"<<(data[i].a-97);//<<","<<data[i].b<<","<<data[i].x<<","<<i <<")";
            fout<<"\t"<<(data[i].a-97)<<","<<data[i].b-97; //<<","<<data[i].x<<","<<")";
            if(i%10==9)
                fout<<"\n";
        }
        fout.close();
    #endif
    delete[] all_counts;
    delete[] all_offsets;
    
}

