#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include "omp.h"
#include "mpi.h"



#define L 24            /* Linear edge dimension of map */
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


int main(int argc, char *argv[]){
    
        int **clusterMap;
        clusterMap = allocate_2d(L,L);
        int **bRight;
        bRight = allocate_2d(P,threads);
        int **bDown;
        bDown = allocate_2d(L,threads);
        
        
        int numtasks;
        int taskid, dest, source, offset, tag1, tag2, tag3, tag4, broadcast, rc;
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
        
        int chunksize = threads * P * P;
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
        
        
        
        //struct Node map[NTHREADS][P][L]; //dont get this !!
        //struct Clusters clusterMapLocal[NTHREADS];
        
        if (taskid == MASTER) {
            for (dest=1; dest < numtasks; dest++){
                int offset = (dest -1) * P;
                MPI_Send(&offset, 1, MPI_INT, dest,tag4, MPI_COMM_WORLD);
                MPI_Send(&clusterMap[offset][0], chunksize , MPI_INT, dest, tag1, MPI_COMM_WORLD);
                if(seedType == 'b'){
                    MPI_Send(&bDown[offset][0], P , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                    MPI_Send(&bDown[offset][0], L , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                }
                printf("Sent elements to task %d offset\n",dest);
            }
            /* Wait to receive results from each task */
            int i;
            for (i=1; i<numtasks; i++) {
                int source = i;
                int offset = (i -1) *P;
                MPI_Recv(&offset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
                MPI_Recv(&clusterMap[offset][0], chunksize , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
                if(seedType == 'b'){
                    MPI_Recv(&bDown[offset][0], P , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&bDown[offset][0], L , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                }
            }
            
            
            // this has not been updated yet - we will need for threading file//
            int p , g;
            printf("MAP:\n");
            for(p = 0; p<L; p++)
            {
                for(g = 0; g < L; g++){
                    printf("%d",clusterMap[p][g]);
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
            MPI_Recv(&clusterMap[offset][0], chunksize , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            if(seedType == 'b'){
                MPI_Recv(&bDown[offset][0], P , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                MPI_Recv(&bDown[offset][0], L , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
            }
/*#pragma omp parallel num_threads(NTHREADS)
            {
                
                int thread_num = omp_get_thread_num();
                printf("threaded num = %d\n",thread_num);
                searchControl(thread_num, p, seedType, map, clusterMap);
            }
            // Send my results back to the master task */
            int p, g;
            for(p = 0; p<P; p++)
            {
                for(g = 0; g < L; g++){
                    clusterMap[offset+p][g] = g*2;
                }
            }
            dest = MASTER;
            MPI_Send(&offset, 1, MPI_INT, dest,tag4, MPI_COMM_WORLD);
            MPI_Send(&clusterMap[offset][0], chunksize , MPI_INT, dest, tag1, MPI_COMM_WORLD);
            if(seedType == 'b'){
                MPI_Send(&bDown[offset][0], P , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                MPI_Send(&bDown[offset][0], L , MPI_INT, dest, tag2, MPI_COMM_WORLD);
            }
            
        } /* end of non-master */
        
        
        
        MPI_Finalize();
    
        printf("perc = %i\n", perc);
        gettimeofday(&end, NULL);
        delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                 end.tv_usec - start.tv_usec) / 1.e6;
        printf("time=%12.10f\n",delta);
        
        return 0;
        exit(EXIT_SUCCESS);
}

