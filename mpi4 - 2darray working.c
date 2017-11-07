

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


#define L 24            /* Linear edge dimension of map */
#define MAPSIZE L*L
#define NTHREADS 1
#define NPROCESSORS 12
#define NODES 1
#define threads 12
#define P L/threads    /* side length of single cluster map */
#define N L*P          // total number of nodes in a single map
#define MASTER 0

struct Node {
    int flag;                   // occupied=1, unoccupied=0, already in a cluster=2+
    int up;                    // (0=not visited,1=visited,2=all directions visited,-1=when parent node)
    int down;
    int left;
    int right;
    int x;                          //co-ords to be easily accessible by cluster
    int y;
};

struct Node** map;

struct Stack{
    int x;
    int y;
    char dir;
} cluster[N];

int top = 0;
int count = 0;
bool popped = false;

int idCounter[N*NTHREADS] = {0};

/*      this functions seeds the map with occupied or unoccupied nodes depending on probability entered
 by the user and initialises all node variables with 0.
 */

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

struct Node **allocate_map(int length){
    int i;
    struct Node** map = (struct Node**) malloc(length*sizeof(struct Node*));
    map[0] = (struct Node*)malloc(sizeof(struct Node)*length*length);
    for (i = 0; i < length; ++i){
        map[i] = (*map + length * i);
    }
    
    return map;
}

void free_array(int** array, int size)
{
    int i; // Declare a matrix index
    for (i=0; i<size; i++) // Iterate though every cell
    {
        free(array[i]); // Free the memory allocated to each cell
    }
    free(array); // Free the matrix allocation
    array = NULL; // Nullify the pointer
}





int main(int argc, char *argv[]){
    
    
    int numtasks;
    int taskid, dest, source, offset, tag1, tag2, tag3, tag4, broadcast, rc;
    int downOffset, rightOffset;
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
    
    //Start MPI
    MPI_Status status;
    
    // Initialisations //
    MPI_Init( &argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    tag2 = 2;
    tag4 = 3;
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
    
    int bDown[L][L];
    int i, j;
    for(i = 0 ; i<L;i++){
        for(j = 0 ; j<L; j++){
            bDown[i][j] = 0;
        }
    }
    
    if (taskid == MASTER) {
        
        // send each task its portion of the map
        // master keeps 1st part
        
        for (dest=1; dest < numtasks; dest++){
            MPI_Send(&downOffset, 1, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                MPI_Send(&bDown[downOffset][0], P*L , MPI_INT, dest, tag2, MPI_COMM_WORLD);
            printf("Sent %d elements to task %d offset= %d\n",P,dest,downOffset);
            printf("Sent elements to task %d offset\n", dest);
            downOffset = downOffset + P;
        }
        
        /* Wait to receive results from each task */
        int i;
        for (i=1; i<numtasks; i++) {
            source = i;
            MPI_Recv(&downOffset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
            MPI_Recv(&bDown[downOffset][0], P*L , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
            printf("Master received tasks from %d\n", source);
            
        }
        
        
        // this has not been updated yet - we will need for threading file//
        int p , g;
        for(p= 0; p<L; p++){
            for(g = 0; g < L; g++){
                printf("%d",bDown[p][g]);
            }
            printf("\n");
        }
        printf("\n");
        printf("perc = %i\n", perc);
        gettimeofday(&end, NULL);
        delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                 end.tv_usec - start.tv_usec) / 1.e6;
        printf("time=%12.10f\n",delta);
    } //end of master section
    
    
    //***** NON-MASTER tasks ******///
    
    if (taskid > MASTER) {
        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&downOffset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
        MPI_Recv(&bDown[downOffset][0], P*L , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
        printf("Slave %d received tasks from Master\n", taskid);
        int y, x;
        for(x = 0; x< L ; x++)
        {
            for(y=downOffset; y<downOffset +P; y++){
                bDown[y][x] = x;
            }
        }
         // Send my results back to the master task */
        dest = MASTER;
        MPI_Send(&downOffset, 1, MPI_INT, dest, tag4, MPI_COMM_WORLD);
        MPI_Send(&bDown[downOffset][0], P*L , MPI_INT, dest, tag2, MPI_COMM_WORLD);
        printf("Slave %d sent tasks to Master\n", taskid);
        
    } /* end of non-master */
    
    
    
    MPI_Finalize();
    
    return 0;
    exit(EXIT_SUCCESS);
}
