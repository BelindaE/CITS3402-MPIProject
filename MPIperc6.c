#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include "omp.h"
#include "mpi.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"

//  compile:  mpicc -fopenmp -o myMPI threadedMPI.c
//  run: syncCluster		/* needs to be run with every change to input file
//  run: mpirun --hostfile host myMPI seedType probability percType
// depending on host file (every node at uwa has 12 processers per node - 1 for master)
//


#define L 48            /* Linear edge dimension of map */
#define MAPSIZE L*L
#define NTHREADS 1
#define NPROCESSORS 12
#define NODES 1
#define threads 12
#define P L/threads    /* side length of single cluster map */
#define N L*P          // total number of nodes in a single map
#define MASTER 0

typedef struct Node NODE;

struct Node {
    int flag;                   // occupied=1, unoccupied=0, already in a cluster=2+
    int up;                    // (0=not visited,1=visited,2=all directions visited,-1=when parent node)
    int down;
    int left;
    int right;
    int x;                          //co-ords to be easily accessible by cluster
    int y;
};

NODE map[L][L];

typedef struct Stack STACK;

struct Stack{
    int x;
    int y;
    char dir;
} cluster[MAPSIZE];

int clusterMap[L][L];
int clusterArray[MAPSIZE] = {0};
int bDown[L] = {0};
int bRight[L*threads] = {0};

int top = 0;
int count = 0;
int clusterID = 2;
bool popped = false;

int idCounter[N*NTHREADS] = {0};

/*      this functions seeds the map with occupied or unoccupied nodes depending on probability entered
 by the user and initialises all node variables with 0.
 */

void initialise_stack(STACK* stackPointer){
	//#pragma omp critical
	{
		stackPointer->x = 0;
		stackPointer->y = 0;
		stackPointer->dir = 'z';
	}
}

int **allocate_2d(int rows, int cols) {
    int j;
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (j=0; j<rows; j++){
        array[j] = &(data[cols*j]);
    }
    
    return array;
}

int *allocate_array(int length) {
    int j;
    int *array = (int *)malloc(length*sizeof(int));
    for (j=0; j<length; j++){
        array[j] = 0;
    }
    
    return array;
}

// modifies 1d array sent over mpi into a 2d array
int **transform_array_toMap(int * array_1d) {
    int i,j;
    int **array= (int **)malloc(P*sizeof(int*));
    for (i=0; i< P; i++){
		for (j = 0; j < L; j++){
        array[j] = &(array_1d[j]);
		}
    }
    
    return array;
}

NODE **allocate_map(int rows, int cols){
	int i;
    struct Node** map = (struct Node**) malloc(rows*sizeof(struct Node*));
	map[0] = (struct Node*)malloc(sizeof(struct Node)*rows*cols);
	for (i = 0; i < rows; ++i){
		map[i] = (*map + cols * i);
	}

	return map;
}

// Free every row of a 2D array
void free_array(int** array, int size) 
{
    int i; 
    for (i=0; i<size; i++) 
    {
        free(array[i]); // Free the memory allocated to each cell
    }
    free(array); // Free the matrix allocation
    //array = NULL; // Nullify the pointer
}

// Free every row of a 2D array
void free_map(NODE** map, int size) 
{
    int i; 
    for (i=0; i<size; i++) 
    {
        free(map[i]); // Free the memory allocated to each cell
    }
    free(map); // Free the matrix allocation
    //array = NULL; // Nullify the pointer
}



void seedSite(double probability, int offset)
{
    double r;
    int seeded;
    int i, j;

    //map = allocate_map(P,L);
    
    for (i = offset; i < (offset + P); i++){  	                    //rows
        for (j = 0; j < L; j++){                //columns
            
            map[i][j].flag = 0;				//reset flag to 0
        }
    }
    
    for (i = offset; i < (P+offset); i++){  	                    //rows
        for (j = 0; j < L; j++){                //columns
            
            r = (double)rand()/RAND_MAX*1.0;
            //printf("pr: %f\t", r);                                //check print of probability values
            if (r<=probability) seeded = 1;
            else seeded = 0;
            
            map[i][j].flag = seeded;
            map[i][j].up = 0;
            map[i][j].down = 0;
            map[i][j].left = 0;
            map[i][j].right = 0;
            map[i][j].x = j;
            map[i][j].y = i;
            /*the following indicates that the top row of the group need not search up, and the bottom not search down, likewise left column need not search left or right column search right.*/
            if(i==0 || i%P==0){
                map[i][j].up = 2;
            }
            else if(i == (P-1) || i == (offset + P -1)){
                map[i][j].down = 2;
            }
            else if(j == 0){
                map[i][j].left = 2;
            }
            else if(j == L-1){
                map[i][j].right = 2;
            }
        }
    }
}

void seedBond(double probability, int offset, int rightOffset, int downOffset)
{
    double r1, r2;
    int seeded1, seeded2;
    int i, j, k;

    //map = allocate_map(P,L);

    for (i = offset; i < (offset + P); i++){  	                    //rows
        for (j = 0; j < L; j++){                //columns
            
            map[i][j].flag = 0;				//reset flag to 0
        }
    }
    
    for (i = offset; i < (offset + P); i++){  	                    //rows
        for (j = 0; j < L; j++){                //columns
            
            r1 = (double)rand()/RAND_MAX*1.0;
            r2 = (double)rand()/RAND_MAX*1.0;
            if (r1>probability) seeded1 = 2;
            else seeded1 = 0;
            if (r2>probability) seeded2 = 2;
            else seeded2 = 0;
            
            map[i][j].right = seeded1;
            if (j+1 < L) map[i][j+1].left = seeded1;
            if (map[i][j].right == 0) {
                map[i][j].flag = 1;
                if (j+1 < L) map[i][j+1].flag = 1;
            }
            
            map[i][j].down = seeded2;
            if (i+1 < P) map[i+1][j].up = seeded2;
            if (map[i][j].down == 0) {
                map[i][j].flag = 1;
                if (i+1 < P) map[i+1][j].flag = 1;
            }
            
            map[i][j].x = j;
            map[i][j].y = i;
            /*the following indicates that the top row of the group need not search up, and the bottom not search down, likewise left column need not search left or right column search right.*/
            
            if(i==0 || i%P==0){
                map[i][j].up = 2;
            }
            if(j == 0){
                map[i][j].left = 2;
            }
            if(i == (P-1) || i == (offset + P -1)){
                bDown[downOffset] =  map[i][j].down;
				
            }
            if(j == L -1){
                bRight[rightOffset] =  map[i][j].right;
				++k;
            }
        }
    }
    
}

// helper functions for depth first search

void push(int i, int j){
    cluster[top].x = j;
    cluster[top].y = i;
    popped = false;
}

STACK pop(){
    popped = true;
    return cluster[top--];
}

STACK peek(){
    return cluster[top];
}

bool isStackEmpty(){
    return top <= -1;
}

void deadEnd(int i, int j){                    //function to change paths to level 2 . no more options
    
    map[i][j].right = 2;                   // change entire node to 2 if it's a dead end
    map[i][j].left = 2;
    map[i][j].down = 2;
    map[i][j].up = 2;
}

void displayNode(){
    printf("(%i, %i)\t", peek().y, peek().x);               //prints out map co-ordinate
    
}

// determines next node in the cluster
// if no more paths from the current node returns previous node

NODE getNextNode(int taskid, int offset) {
    popped = false;
    int j = peek().x ;
    int i = peek().y ;
    //printf("check1 (%i, %i)\t",i,j);
    
    /*the following first checks if flag on, then path is 0 and if 1 it also checks previous direction to make sure not cycling back on itself.  If no options available then 2 is applied to the node and stack drops the top node
     note that I have not assumed the lattice can cycle around itself as I think we want to retain edge nodes assuming it's a small part of a larger area - can change this easily with %L if required*/
    
    
    
    if ((i+1 < (offset+P)) && map[i+1][j].flag == 1 && map[i+1][j].up == 0) {
        cluster[top].dir = 'd';
        map[i][j].down = 1;
        map[i+1][j].up = 1;
        i = i+1;
    }
    else if ((j+1 < L) && map[i][j+1].flag == 1 && map[i][j+1].left == 0){
        cluster[top].dir = 'r';
        map[i][j].right = 1;
        map[i][j+1].left = 1;
        j = j+1;
    }
    else if ((j-1 >= 0) && (map[i][j-1].flag == 1) && (map[i][j-1].right == 0)) {
        cluster[top].dir = 'l';
        map[i][j].left = 1;
        map[i][j-1].right = 1;
        j = j-1;
    }
    else if ((i-1 >= offset) && (map[i-1][j].flag == 1) && (map[i-1][j].down == 0)) {
        cluster[top].dir = 'u';
        map[i][j].up = 1;
        map[i-1][j].down = 1;
        i = i-1;
    }
    
    else if ((i+1 < (offset+P)) && map[i+1][j].flag == 1 && map[i+1][j].up == 1
             && cluster[top-1].dir != 'u' && top >= 1){
        cluster[top].dir = 'd';
        i = i+1;
        
    }
    else if ((j+1 < L) && map[i][j+1].flag == 1 && map[i][j+1].left ==1
             && cluster[top-1].dir != 'l' && top >= 1){
        cluster[top].dir = 'r';
        j = j+1;
    }
    else if ((j-1 >= 0) && (map[i][j-1].flag == 1) && (map[i][j-1].right == 1)
             && cluster[top-1].dir != 'r' && top >= 1)                       {
        cluster[top].dir = 'l';
        j = j-1;
    }
    else if ((i-1 >= offset) && (map[i-1][j].flag == 1) && (map[i-1][j].down == 1)
             && cluster[top-1].dir != 'd' && top >= 1) {
        cluster[top].dir = 'u';
        i = i-1;
    }
    
    else {
        //printf("CurrentID = %d\n",currentID);
        clusterMap[peek().y][peek().x] = clusterID;
        deadEnd(peek().y, peek().x);                    // changes path status to 2
        pop();
        // reduces top index by 1
        if (top == -1) return map[i][j];
        {
            return map[cluster[top].y][cluster[top].x];
        }
    }
    return map[i][j];
}


/*
 This searches on cluster of nodes
 */
int depthFirstSearch(int i, int j, int taskid, int offset){
    //the cluster ID will be the reference point for the side-arrays, and cluster sizes.
    top = 0;
    int tempX = 0;
    int tempY = 0;
    int tempC = 0;
    //STACK* cluster = (STACK*) malloc(sizeof(STACK));
	//initialise_stack(cluster);

    if (map[i][j].flag == 0 || map[i][j].flag >1) return 0;
    // check origin flag before commencing DFS,
    //aborting if either empty or already searched.
    count = 0;
    
    
    // reset cluster percolation arrays to zero
    
    while(!isStackEmpty()){
        
        //initialise temporary holding stack
        NODE temp[1][1];
        temp[0][0].x = 0;
        temp[0][0].y = 0;
        
        cluster[top].x = map[i][j].x;      //store co-ordinates in the cluster stack
        cluster[top].y = map[i][j].y;
        
        if(map[i][j].flag == 1){
            map[i][j].flag = clusterID + taskid*L*P/2;
            //marking the node as both explored, and to which cluster of nodes it belongs-
            count++;         //will only update if there's a new filled node attached.
        }
        //displayNode();                              // print check of top node in the cluster
        push(cluster[top].y, cluster[top].x);                   // push node into the stack
        tempX = cluster[top].x;
        tempY = cluster[top].y;
        tempC = clusterID;
        
        temp[0][0] = getNextNode(taskid, offset);                  // allocate new node coords
        i = temp[0][0].y;
        j = temp[0][0].x;
        int f = temp[0][0].flag;
        printf("f = %d, taskid = %d\n",f, taskid);
        
        if(tempC != f)
        {
            clusterMap[tempY][tempX] = f + taskid*L*P/2;
            //currentID = f + taskid*L*P/2;
        }
        if(f != 1)
        {
            clusterMap[i][j] = f + taskid*L*P/2;
            //currentID = f + taskid*L*P/2;
           
        }
        
        if (!popped) ++top;            //want to go back to previous node if no connection made
        
    }

	//free(cluster);
    idCounter[clusterID-2] = count;
    return count;
    
}



void searchControlBond(int taskid, double p, char seedType, int offset, int rightOffset, int downOffset){

	int i, j, k, dfs;
	//int** clusterMap = transform_array_toMap(clusterArray);

    seedBond(p, offset, rightOffset, downOffset);

    for(i =offset; i<(P+offset); i++){
    	for(j =0; j< L; j++){
        	dfs = depthFirstSearch(i, j, taskid, offset);
				if (dfs>0) clusterID++;
        }
    }

    printf("Cluster Map %d\n",taskid);

	// update clusterArray to transfer back to master
	k = offset;
	for (i=offset; i < (P+offset); i++){
		for (j=0; j < L; j++){
			clusterArray[k] = map[i][j].flag;
			++k;
		}
	}
	//free_map(map, P);
}

void searchControlSite(int taskid, double p, char seedType, int offset){

	int i, j, k, dfs;
	int** clusterMap = transform_array_toMap(clusterArray);

    seedSite(p, offset);

    for(i = offset; i < (P+offset); j++){
    	for(j =0; j<L; j++){
        	dfs = depthFirstSearch(i, j, taskid, offset);
				if (dfs>0) clusterID++;
        }
    }

	// update clusterArray to transfer back to master
	k = offset;
	for (i=offset; i < (P+offset); i++){
		for (j=0; j < L; j++){
			clusterArray[k] = map[i][j].flag;
			++k;
		}
	}
	//free_map(map, P);
}            


int main(int argc, char *argv[]){
    
        int sum;
        int numtasks, i;
        int taskid, dest, source, tag1, tag2, tag3, tag4, tag5, tag6, broadcast, rc;
        int offset, downOffset, rightOffset;
        double delta;
        srand(time(NULL));                              //should only be called once
        
        struct timeval start, end;
        gettimeofday(&start, NULL);

		// **note that argument position is different for MPI commands**

        char seedType = argv[1][0];            // enter 's' for site percolation, 'b' for bond percolation
        double p = atof(argv[2]);            // enter probability between 0 - 1, eg. 0.55
        int perc = atoi(argv[3]);            // enter 0 = column percolation
        //         1 = row percolation
        //          2 = both row and column percolation
        //       3 = see results for all three percolation types
            if (seedType == 'b') {

             void searchControlBond(int taskid, double p, char seedType, int offset, int rightOffset, int downOffset);
            } else {
             void searchControlSite(int taskid, double p, char seedType, int offset);            
            }
        
        //Start MPI
        MPI_Status status;
        
        // Initialisations //
        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

        tag1 = 1;
        tag2 = 2;
        tag3 = 3;
        tag4 = 4;
		tag5 = 5;
		tag6 = 6;
        // not sure we need this //
        /*if (L % (numtasks*NTHREADS -1) != 0) {
            printf("Quitting. Number of MPI tasks must be divisible by number of maps.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(0);
        }
        */
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        printf("MPI task %d has started...\n", taskid);
        
        int part = L * L / numtasks;
		printf("numtasks = %d\n", numtasks);

            // initialise the map array  MPIperc5.c changes clustermap to a 1d array
            //int *clusterArray;
	        /*clusterArray = allocate_array(MAPSIZE);
            int *bRight;
            bRight = allocate_array(L*numtasks);
            int *bDown;
            bDown = allocate_array(L);*/
        
        if (taskid == MASTER) {
		
		offset = part;   
        downOffset = P;  
		rightOffset = L;
            // send each task its portion of the map
            // master keeps 1st part    

            for (dest=1; dest < numtasks; dest++){
                 MPI_Send(&offset, 1 , MPI_INT, dest, tag1, MPI_COMM_WORLD);
                 //MPI_Send(&clusterArray[offset], part , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                if(seedType == 'b'){
	                 MPI_Send(&downOffset, 1 , MPI_INT, dest, tag3, MPI_COMM_WORLD);
                     //MPI_Send(&bDown[downOffset], P , MPI_INT, dest, tag4, MPI_COMM_WORLD);
	                 MPI_Send(&rightOffset, 1 , MPI_INT, dest, tag5, MPI_COMM_WORLD);
                     //MPI_Send(&bRight[rightOffset], L , MPI_INT, dest, tag6, MPI_COMM_WORLD);
                }
                printf("Sent elements to task %d offset, P = %d\n", dest,P);
				offset = offset + part;
				downOffset = downOffset + P;
				rightOffset = rightOffset + L;
            }

            //master does its part of the work
			offset = 0;
            downOffset = 0;
			rightOffset = 0;

            if (seedType == 'b') {

             searchControlBond(taskid, p, seedType, offset, rightOffset, downOffset);
            } else {
             searchControlSite(taskid, p, seedType, offset);            
            }

            /* Wait to receive results from each task */

            for (i=1; i<numtasks; i++) {
                source = i;
                MPI_Recv(&offset, 1 , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
                //MPI_Recv(&clusterArray[offset], part , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                if(seedType == 'b'){
                    MPI_Recv(&downOffset, 1 , MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
                    //MPI_Recv(&bDown[downOffset], P , MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
                    MPI_Recv(&rightOffset, 1 , MPI_INT, source, tag5, MPI_COMM_WORLD, &status);                     
					//MPI_Recv(&bRight[rightOffset], L , MPI_INT, source, tag6, MPI_COMM_WORLD, &status);
                }
	            printf("Master received tasks from %d\n", source);

            }
            
            
            // this has not been updated yet - we will need for threading file//
            int p , g;
            printf("MAP:\n");
                for(g = 0; g < MAPSIZE; g++){
                    printf("%d",clusterArray[g]);
					++p;
					if (p == L) {
						printf("\n");
						p = 0;
					}
                }
                printf("\n");
            printf("\n");
        } //end of master section

        //***** NON-MASTER tasks ******///

        if (taskid > MASTER) {
            /* Receive my portion of array from the master task */
            source = MASTER;
            MPI_Recv(&offset, 1 , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            //MPI_Recv(&clusterArray[offset], part , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
            if(seedType == 'b'){
                MPI_Recv(&downOffset, 1 , MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
                //MPI_Recv(&bDown[downOffset], P , MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
                MPI_Recv(&rightOffset, 1 , MPI_INT, source, tag5, MPI_COMM_WORLD, &status);
                //MPI_Recv(&bRight[rightOffset], L , MPI_INT, source, tag6, MPI_COMM_WORLD, &status);
            }
            printf("Slave %d received tasks from Master\n", taskid);

            if (seedType == 'b') {

             searchControlBond(taskid, p, seedType, offset, rightOffset, downOffset);
            } else {
             searchControlSite(taskid, p, seedType, offset);            
            }

    /*#pragma omp parallel num_threads(NTHREADS)
            {
                
                int thread_num = o
				omp_get_thread_num();
                printf("threaded num = %d\n",thread_num);
                searchControl(thread_num, p, seedType, map, clusterMap);
            }
            // Send my results back to the master task */
            dest = MASTER;
            MPI_Send(&offset, 1 , MPI_INT, dest, tag1, MPI_COMM_WORLD);
            //MPI_Send(&clusterArray[offset], part , MPI_INT, dest, tag2, MPI_COMM_WORLD);
            if(seedType == 'b'){
                MPI_Send(&downOffset, 1 , MPI_INT, dest, tag3, MPI_COMM_WORLD);
                //MPI_Send(&bDown[downOffset], P , MPI_INT, dest, tag4, MPI_COMM_WORLD);
                MPI_Send(&rightOffset, 1 , MPI_INT, dest, tag5, MPI_COMM_WORLD);
                //MPI_Send(&bRight[rightOffset], L , MPI_INT, dest, tag6, MPI_COMM_WORLD);
            }
            printf("Slave %d sent tasks to Master\n", taskid);
            
        } /* end of non-master */
        
        
        
        printf("perc = %i\n", perc);
        gettimeofday(&end, NULL);
        delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                 end.tv_usec - start.tv_usec) / 1.e6;

		if (taskid == MASTER) printf("time=%12.10f\n",delta);

        MPI_Finalize();

        return 0;
        exit(EXIT_SUCCESS);
}
