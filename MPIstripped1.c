#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include "omp.h"
#include "mpi.h"



#define L 24            /* Linear edge dimension of map */
#define MAPSIZE L*L
#define NTHREADS 1
#define NPROCESSORS 12
#define NODES 1
#define threads NTHREADS*NPROCESSORS*NODES
#define P L/threads	/* side length of single cluster map */
#define N L*P          // total number of nodes in a single map
#define MASTER 0

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

void searchControlBond(int taskid, double p, char seedType, int clusterMap[threads], int bRight[], int offset, int downOffset[]){
	// do something
}

void searchControl(int taskid, double p, char seedType, int clusterMap[threads], int offset){

	// do something
}			


int main(int argc, char *argv[]){
    
        
        int numtasks;
        int taskid, dest, source, offset, tag1, tag2, tag3, tag4, broadcast, rc;
		int downOffset, rightOffset;
        double delta;
        srand(time(NULL));                              //should only be called once
        
        struct timeval start, end;
        gettimeofday(&start, NULL);
        
        char seedType = argv[1][0];			// enter 's' for site percolation, 'b' for bond percolation
        double p = atof(argv[2]);			// enter probability between 0 - 1, eg. 0.55
        int perc = atoi(argv[3]);			// enter 0 = column percolation
        //		 1 = row percolation
        // 		 2 = both row and column percolation
        //       3 = see results for all three percolation types
        
        //Start MPI
        MPI_Status status;
        
        // Initialisations //
        
        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
        tag1 = 1;
        tag2 = 2;
        tag3 = 3;
        tag4 = 4;
        // not sure we need this //
        /*if (L % (numtasks*NTHREADS -1) != 0) {
            printf("Quitting. Number of MPI tasks must be divisible by number of maps.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(0);
        }
        */
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        printf("MPI task %d has started...\n", taskid);
        
        int chunksize = L * L / numtasks;
		offset = chunksize;   
		downOffset = P;     

        if (taskid == MASTER) {

			// initialise the map array
			int **clusterMap;
			clusterMap = allocate_2d(L,L);
			int **bRight;
			bRight = allocate_2d(threads,P);
			int *bDown;
			bDown = allocate_array(L);
        
			// send each task its portion of the map
			// master keeps 1st part	

            for (dest=1; dest < numtasks; dest++){
                MPI_Send(&offset, 1, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                MPI_Send(&clusterMap[dest][offset], chunksize , MPI_INT, dest, tag1, MPI_COMM_WORLD);
                if(seedType == 'b'){
                    MPI_Send(&bDown[downOffset], P , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                    MPI_Send(&bRight[dest][0], L , MPI_INT, dest, tag3, MPI_COMM_WORLD);
                }
                printf("Sent elements to task %d offset\n", dest);
				offset = offset + chunksize;
				downOffset = downOffset + P;
            }

			//master does its part of the work

			offset = 0;
			downOffset = 0;
			if (seedType == 'b') {
			//***note need to create new function***//
	 		searchControlBond(taskid, p, seedType, clusterMap[0], bRight[0], offset, downOffset);
			} else {
	 		searchControl(taskid, p, seedType, clusterMap[0], offset);			
			}

            /* Wait to receive results from each task */
            int i;
            for (i=1; i<numtasks; i++) {
                source = i;
                MPI_Recv(&offset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
                MPI_Recv(&clusterMap[i][offset], chunksize , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
                if(seedType == 'b'){
                    MPI_Recv(&bDown[downOffset], P , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&bRight[dest][0], L , MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
                }
            }
            
            
            // this has not been updated yet - we will need for threading file//
            int p , g;
            printf("MAP:\n");
            for(p = 0; p<L; p++)
            {
                for(g = 0; g < L; g++){
                    printf("%d",clusterMap[g][p]);
                }
                printf("\n");
            }
            printf("\n");
        } //end of master section

        //***** NON-MASTER tasks ******///

        if (taskid > MASTER) {
            /* Receive my portion of array from the master task */
            source = MASTER;
            MPI_Recv(&offset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
            MPI_Recv(&clusterMap[taskid][offset], chunksize , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            if(seedType == 'b'){
                MPI_Recv(&bDown[downOffset], P , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                MPI_Recv(&bRight[taskid][0], L , MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
            }

			if (seedType == 'b') {
			//***note need to create new function***//
	 		searchControlBond(taskid, p, seedType, clusterMap[0], bRight[0], offset, downOffset);
			} else {
	 		searchControl(taskid, p, seedType, clusterMap[0], offset);			
			}

	/*#pragma omp parallel num_threads(NTHREADS)
            {
                
                int thread_num = omp_get_thread_num();
                printf("threaded num = %d\n",thread_num);
                searchControl(thread_num, p, seedType, map, clusterMap);
            }
            // Send my results back to the master task */
            dest = MASTER;
            MPI_Send(&offset, 1, MPI_INT, dest, tag4, MPI_COMM_WORLD);
            MPI_Send(&clusterMap[taskid][offset], chunksize , MPI_INT, dest, tag1, MPI_COMM_WORLD);
            if(seedType == 'b'){
                MPI_Send(&bDown[downOffset][0], P , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                MPI_Send(&bRight[taskid][0], L , MPI_INT, dest, tag3, MPI_COMM_WORLD);
            }
            
        } /* end of non-master */
        
        
        
        MPI_Finalize();
    
        printf("perc = %i\n", perc);
        gettimeofday(&end, NULL);
        delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                 end.tv_usec - start.tv_usec) / 1.e6;
        printf("time=%12.10f\n",delta);

		free(clusterMap);
		free(bRight);
		free(bDown);
        
        return 0;
        exit(EXIT_SUCCESS);
}
