
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

//  compile:  gcc -fopenmp -o s sequential1.c
//  run:      ./s seedType probability
// this version changes the Node map back to a non-pointer to enable new seeded maps

#define L 24            /* Linear edge dimension of map */
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

void seedSite(double probability, int thread_num, struct Node map[P][L])
{
    double r;
    int seeded;
    int i, j;
    
    for (i = 0; i < P; i++){                        //rows
        for (j = 0; j < L; j++){                //columns
            
            map[i][j].flag = 0;                         //reset flag to 0
        }
    }
    
    for (i = 0; i < P; i++){                        //rows
        for (j = 0; j < L; j++){                //columns
            
            r = (double)rand()/RAND_MAX*1.0;
            //printf("pr: %f\t", r);                                //check print of probability values
            if (r<=probability) seeded = 1;
            else seeded = 0;
            
            //struct Node node = (struct Node) malloc(sizeof(struct Node));
            
            map[i][j].flag = seeded;
            map[i][j].up = 0;
            map[i][j].down = 0;
            map[i][j].left = 0;
            map[i][j].right = 0;
            map[i][j].x = j;
            map[i][j].y = i;
            /*the following indicates that the top row of the group need not search up, and the bottom not search down, likewise left column need not search left or right column search right.*/
            if(i==0){
                map[i][j].up = 2;
            }
            else if(i == L-1){
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

void seedBond(double probability, int thread_num, struct Node map[P][L], int bRight[][P],int  bDown[][L],int  taskid)
{
    double r1, r2;
    int seeded1, seeded2;
    int i, j;
    
    for (i = 0; i < P; i++){                        //rows
        for (j = 0; j < L; j++){                //columns
            
            map[i][j].flag = 0;                         //reset flag to 0
        }
    }
    
    for (i = 0; i < P; i++){                        //rows
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
            
            if(i==0){
                map[i][j].up = 2;
            }
            if(j == 0){
                map[i][j].left = 2;
            }
            if(i == P -1){
                bDown[taskid][j] =  map[i][j].down;
            }
            if(j == L -1){
                bRight[taskid][j] =  map[i][j].right;
            }
        }
    }
    
}


void printSites(struct Node map[P][L]){
    int i,j;
    for (i = 0; i < (L); i++){
        for (j = 0; j < L; j++){
            if (map[i][j].flag>0){
                if (map[i][(j+1)%L].flag > 0) printf("%i ----- ", map[i][j].flag);
                
                else printf("%i\t", map[i][j].flag);
            }
            else printf("%i\t", map[i][j].flag);
            
        }
        printf("\n");
        for (j = 0; j < L; j++){
            if (map[i][j].flag>0){
                
                if (map[(i+1)%L][j].flag >0) printf("|\t");
                else printf("\t");
            }
            else printf("\t");
        }
        printf("\n");
    }
    printf("\n\n");
}
void printBonds(struct Node map[P][L]){
    int i,j;
    for (i = 0; i < (P); i++){
        for (j = 0; j < L; j++){
            if (map[i][j].right==0){
                printf("%i ----- ", map[i][j].flag);
            }
            else printf("%i\t", map[i][j].flag);
            
        }
        printf("\n");
        for (j = 0; j < L; j++){
            if (map[i][j].down==0) printf("|\t");
            else printf("\t");
        }
        printf("\n");
    }
}


// helper functions for depth first search

void push(int i, int j){
    cluster[top].x = j;
    cluster[top].y = i;
    popped = false;
}

struct Stack pop(){
    popped = true;
    return cluster[top--];
}

struct Stack peek(){
    return cluster[top];
}

bool isStackEmpty(){
    return top <= -1;
}

void deadEnd(int i, int j, struct Node map[P][L]){                    //function to change paths to level 2 . no more options
    
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



struct Node getNextNode(int thread_num, struct Node map[P][L], int clusterMap[][L] ,int currentID, int offset, int taskid, int bRight[][P], int bDown[][L]) {
    popped = false;
    int j = peek().x ;
    int i = peek().y ;
    //printf("check1 (%i, %i)\t",i,j);
    
    /*the following first checks if flag on, then path is 0 and if 1 it also checks previous direction to make sure not cycling back on itself.  If no options available then 2 is applied to the node and stack drops the top node
     note that I have not assumed the lattice can cycle around itself as I think we want to retain edge nodes assuming it's a small part of a larger area - can change this easily with %L if required*/
    
    
    
    if ((i+1 < P) && map[i+1][j].flag == 1 && map[i+1][j].up == 0) {
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
    else if ((i-1 >= 0) && (map[i-1][j].flag == 1) && (map[i-1][j].down == 0)) {
        cluster[top].dir = 'u';
        map[i][j].up = 1;
        map[i-1][j].down = 1;
        i = i-1;
    }
    
    else if ((i+1 < P) && map[i+1][j].flag == 1 && map[i+1][j].up == 1
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
    else if ((i-1 >= 0) && (map[i-1][j].flag == 1) && (map[i-1][j].down == 1)
             && cluster[top-1].dir != 'd' && top >= 1) {
        cluster[top].dir = 'u';
        i = i-1;
    }
    
    else {
        //printf("CurrentID = %d\n",currentID);
        clusterMap[offset + peek().y][peek().x] = currentID;
        deadEnd(peek().y, peek().x, map);                    // changes path status to 2
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
int depthFirstSearch(int i, int j, int clusterID, int thread_num, struct Node map[P][L], int clusterMap[][L], int bRight[][P], int bDown[][L] , int offset, int taskid){
    //the cluster ID will be the reference point for the side-arrays, and cluster sizes.
    top = 0;
    int tempX = 0;
    int tempY = 0;
    int tempC = 0;
    int currentID = 0;
    
    if (map[i][j].flag == 0 || map[i][j].flag >1) return 0;
    // check origin flag before commencing DFS,
    //aborting if either empty or already searched.
    count = 0;
    
    
    // reset cluster percolation arrays to zero
    
    while(!isStackEmpty()){
        
        //initialise temporary holding stack
        struct Node temp[1][1];
        temp[0][0].x = 0;
        temp[0][0].y = 0;
        
        cluster[top].x = map[i][j].x;      //store co-ordinates in the cluster stack
        cluster[top].y = map[i][j].y;
        
        if(map[i][j].flag == 1){
            map[i][j].flag = clusterID + thread_num*L*P/2;
            //marking the node as both explored, and to which cluster of nodes it belongs-
            count++;         //will only update if there's a new filled node attached.
        }
        //displayNode();                              // print check of top node in the cluster
        push(cluster[top].y, cluster[top].x);                   // push node into the stack
        tempX = cluster[top].x;
        tempY = cluster[top].y;
        tempC = currentID;
        
        temp[0][0] = getNextNode(thread_num, map, clusterMap, currentID, offset , taskid, bRight, bDown);                  // allocate new node coords
        i = temp[0][0].y;
        j = temp[0][0].x;
        int f = temp[0][0].flag;
        
        if(tempC != f)
        {
            clusterMap[offset + tempY][tempX] = f + taskid*L*P/2;
            currentID = f + thread_num*L*P/2;
        }
        if(f != 1)
        {
            clusterMap[offset + tempY][tempX] = f + taskid*L*P/2;
            currentID = f + thread_num*L*P/2;
            
        }
        
        if (!popped) ++top;            //want to go back to previous node if no connection made
        
    }
    idCounter[currentID] = count;
    return count;
    
}

void matchClusters(float probability, char seedType, int perc, int clusterMap[][L], int bRight[][P], int bDown[][L]){
    int n, i, j, l , m , g;
    int max = 0;
    int rows = 0;
    int cols = 0;
    int both = 0;
    int connected[2][N*2];
    int biggestCluster[N*2] = {0};
    
    int tempCluster[N*P/2] =  {0};
    int k = 0;
    int tempRows[L] = {0};
    int tempCols[L] = {0};
    int tempMax = 0;
    int colP = 0;
    int rowP = 0;
    
    /*for (int i = 0; i < clusterID; i++){
     printf("clusterID %i, count %i\n", i, idCounter[i]);
     }
     for (n = 0; n < NTHREADS; n++){
     for(int i = 0; i < L; i++){
     printf("thread_num %i, clusterMap ID %i\n", n, clusterMap[n].m[i][L-1]);
     }
     }
     */
    for(i = 0; i<N*2;i++)
    {
        connected[0][i] = 0;
        connected[1][i] = 0;
    }
    
    int c = 0;
    
    for (n=0; n<threads; n++){
        for(i = 0; i < P; i++){
            if (clusterMap[n*P + i][L-1] > 0 && clusterMap[n*P+i][0] >0){
                if(seedType == 'b' && bRight[n][i] == 2) continue; // checking seeded bond between
                connected[0][c] = clusterMap[n*P+i][L-1];
                connected[1][c] = clusterMap[n*P+i][0];
                c++;
            }//R0
        }
        for(i = 0; i < L; i++){
            if(n == NTHREADS-1 && clusterMap[0][i] > 0 && clusterMap[n*P + P-1][i] > 0){
                if(seedType == 'b' && bDown[n][i] == 2) continue; // checking seeded bond between
                connected[0][c] = clusterMap[0][i] ;
                connected[1][c] = clusterMap[n*P + P-1][i];
                c++;
            }
            else if (clusterMap[n+1][i] > 0 && clusterMap[n*P + P-1][i] > 0){
                if(seedType == 'b' && bDown[n][i] == 2) continue; // checking seeded bond between
                connected[0][c] = clusterMap[n*P+1][i] ;
                connected[1][c] = clusterMap[n*P + P-1][i];
                c++;
            } //C2
        }
    }
    /*for (int i = 0; i < clusterID; i++){
     printf("connected ( %i, %i)\n", connected[0][i], connected[1][i]);
     }
     */
    for(i = 0; i < c ; i++)
    {
        if(connected[0][0] == 0 && i ==0)
        {
            break;
        }
        for(j = 0 ; j < c; j++)
        {
            if(tempCluster[0] == 0 && connected[0][j] !=0)
            {
                tempCluster[0] = connected[0][j];
                tempCluster[1] = connected[1][j];
                tempMax = tempMax+  idCounter[connected[0][j]] + idCounter[connected[1][j]];
                idCounter[connected[1][j]] = 0;
                idCounter[connected[0][j]] = 0;
                connected[0][j] = 0;
                connected[1][j] = 0;
                k =2;
            }
            else{
                for(l = 0; l < k; l++)
                {
                    if(tempCluster[l] == 0) continue;
                    else if(tempCluster[l] == connected[0][j]) {
                        tempCluster[k] = connected[1][j];
                        tempMax = tempMax + idCounter[connected[1][j]];
                        idCounter[connected[1][j]] = 0;
                        connected[0][j] =0;
                        connected[1][j] = 0;
                        l = -1;
                        j = -1;
                        k++;
                    }
                    else if(tempCluster[l] == connected[1][j])
                    {
                        tempCluster[k] = connected[0][j];
                        tempMax = tempMax + idCounter[connected[0][j]];
                        idCounter[connected[0][j]] = 0;
                        connected[0][j] = 0;
                        connected[1][j] = 0;
                        l = -1;
                        j = -1;
                        k++;
                    }
                }
                
            }
        }
        for(n = 0; n < NTHREADS ; n++)
        {
            for(j = 0; j < P; j++)
            {
                for(m = 0; m < L ;m++)
                {
                    for(l = 0; l < k ; l++)
                    {
                        if(clusterMap[n*P+ j][m] == tempCluster[l] )
                        {
                            tempRows[j +P*n] = 1;
                            tempCols[m] = 1;
                        }
                    }
                }
            }
        }
        for(j = 0; j < L;j++)
        {
            if(tempCols[j] == 0)
            {
                colP++;
            }
            if(tempRows[j] == 0)
            {
                rowP++;
            }
        }
        if(tempMax > max)
        {
            max = tempMax;
            for(j = 0 ; j < c ; j++) biggestCluster[j] = tempCluster[j];
        }
        if(colP == 0) cols = 1;
        if(rowP == 0) rows = 1;
        if(cols ==1 && rows == 1) both =1;
    }
    for(g = 0; g< NTHREADS; g++)
    {
        for(i = 0; i < P; i++)
        {
            for( j = 0; j  < L ; j++)
            {
                if(j <L  && i < L && clusterMap[g*P+ i][j] > 0){
                    int q = 0;
                    for(l = 0 ; l<c ; l++)
                    {
                        if(biggestCluster[l] ==clusterMap[g*P+ i][j])
                        {
                            printf(ANSI_COLOR_RED     "1"    ANSI_COLOR_RESET);
                            q = 1;
                            break;
                        }
                    }
                    if(q == 0)
                    {
                        printf("1");
                    }
                }
                else printf("0");
            }
            printf("\n");
        }
    }
    
    printf("Lattice size = %i x %i \n", L, L);
    printf("Number of threads = %i \n", NTHREADS);
    printf("Seeding type = %c, Probability = %f \n\n", seedType, probability);
    printf("perc in matchClusters = %i\n", perc);
    if (perc == 0) printf("columns percolated = %d \n", cols);
    else if (perc == 1) printf("rows percolated = %d \n", rows);
    else if (perc == 2) printf("both percolated = %d \n", both);
    else {
        printf("columns percolated = %d \n", cols);
        printf("rows percolated = %d \n", rows);
        printf("both percolated = %d \n", both);
    }
    
    printf("Biggest Cluster = %d \n",max);
    
    
    printf("\n");
}

void searchControl(int thread_num, double probability, char seedType, struct Node map[P][L], int clusterMap[][L], int bRight[][P], int bDown[][L], int offset, int taskid){
    int i, j, dfs;
    int clusterID = 2;
    if (seedType == 's') {
        seedSite(probability, thread_num, map);
        //printSites();
    }
    else {
        seedBond(probability, thread_num, map, bRight, bDown, taskid);
        //printBonds();
    }
    for(j = 0; j<P; j++)
    {
        for(i = 0; i<L; i++)
        {
            printf("%d",map[j][i].flag);
        }
        printf("\n");
    }
    for(j = 0; j<P; j++)
    {
        for(i = 0; i<L; i++)
        {
            clusterMap[j][i]= map[j][i].flag
        }
    }
    
    /*for(j =0; j<L; j++){
     for(i =0; i<P; i++){
     dfs = depthFirstSearch(i,j, clusterID, thread_num, map, clusterMap, bRight, bDown, offset, taskid);
     if (dfs>0) clusterID++;
     }
     }
    
    printf("Cluster Map %d\n",thread_num);
    
    /*printf("Cluster Map %d\n",thread_num);
     for(int r = 0 ; r < P ; r++)
     {
     for(int q = 0; q< L ; q++)
     {
     printf("%d\t",clusterMap[thread_num].m[r][q]);
     }
     printf("\n");
     }
     printf("\n");*/
    
}

int main(int argc, char *argv[]){
    
    if(argc != 4)
    {
        fprintf(stderr, " ERROR\n Usage: ./perc seedType probability\n");
        // example: "./perc s 0.55 3"
        exit(EXIT_FAILURE);             // Exit program indicating failure
    }
    else if(L % NTHREADS != 0)
    {
        fprintf(stderr, " ERROR\n L %% NTHREADS != 0\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        char seedType = argv[1][0];                     // enter 's' for site percolation, 'b' for bond percolation
        double p = atof(argv[2]);                       // enter probability between 0 - 1, eg. 0.55
        int perc = atoi(argv[3]);                       // enter 0 = column percolation
        //               1 = row percolation
        //               2 = both row and column percolation
        //       3 = see results for all three percolation types
        
        
        int numtasks;
        int taskid, dest, source, offset, tag1, tag2, tag3, tag4, broadcast, rc;
        int j, i, g, m;
        
        double delta;
        srand(time(NULL));                              //should only be called once
        
        struct timeval start, end;
        gettimeofday(&start, NULL);
        
        //Start MPI
        MPI_Status status;
        
        // Initialisations //
        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
        // initialise the map array
        tag1 = 1;
        tag2 = 2;
        tag3 = 3;
        tag4 = 4;
        // not sure we need this //
        /*if (L % (numtasks*NTHREADS -1) != 0)
         printf("Quitting. Number of MPI tasks must be divisible by number of maps.\n");
         MPI_Abort(MPI_COMM_WORLD, rc);
         exit(0);
         }
         */
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        
        
        int clusterMap[L][L];
        for(i = 0 ; i<L;i++){
            for(j = 0 ; j<L; j++){
                clusterMap[i][j] = 0;
            }
        }
        int bRight[threads][P];
        for(i = 0 ; i<threads;i++){
            for(j = 0 ; j<P; j++){
                bRight[threads][j] = 0;
            }
        }
        int bDown[threads][L];
        for(i = 0 ; i<threads;i++){
            for(j = 0 ; j<L; j++){
                bDown[i][j] = 0;
            }
        }
        
        int ntp = L/numtasks;
        offset = ntp;
        int partition = L * L / numtasks;
        if(taskid == MASTER) {
            
            
            // send each task its portion of the map
            // master keeps 1st part
            
            for(dest=1; dest < numtasks; dest++){
                MPI_Send(&offset, 1, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                MPI_Send(&clusterMap[offset][0], partition , MPI_INT, dest, tag1, MPI_COMM_WORLD);
                if(seedType == 'b'){
                    MPI_Send(&bDown[dest][0], L , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                    MPI_Send(&bRight[dest][0], P , MPI_INT, dest, tag3, MPI_COMM_WORLD);
                }
                offset = offset + ntp;
            }
            
            struct Node map[taskid][P][L];
            omp_set_num_threads(NTHREADS);
#pragma omp parallel
            {
                
                int thread_num = omp_get_thread_num();
                printf("0. threaded num = %d\n",thread_num);
                searchControl(thread_num, p, seedType, map[taskid], clusterMap, bRight, bDown, 0, taskid);
            }
            /* Wait to receive results from each task */
            
            for(i=1; i<numtasks; i++) {
                source = i;
                MPI_Recv(&offset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
                MPI_Recv(&clusterMap[offset][0], partition , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
                if(seedType == 'b'){
                    MPI_Recv(&bDown[dest][0], P , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&bRight[dest][0], L , MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
                }
                
                printf("Master received tasks from %d\n", source);
                
            }
            
            
            
            // this has not been updated yet - we will need for threading file//
            int p , g;
            printf("MAP:\n");
            for(m = 0; m<L; m++)
            {
                for(g = 0; g < L; g++){
                    printf("%d",clusterMap[m][g]);
                }
                printf("\n");
            }
            printf("\n");
            
            //matchClusters(p, seedType, perc, clusterMap, bRight, bDown);
            printf("perc = %i\n", perc);
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                     end.tv_usec - start.tv_usec) / 1.e6;
            printf("time=%12.10f\n",delta);
            
        } //end of master section
        
        //***** NON-MASTER tasks ******///
        
        if(taskid > MASTER) {
            /* Receive my portion of array from the master task */
            source = MASTER;
            MPI_Recv(&offset, 1, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
            MPI_Recv(&clusterMap[offset][0], partition , MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            if(seedType == 'b'){
                MPI_Recv(&bDown[taskid][0], P , MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
                MPI_Recv(&bRight[taskid][0], L , MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
            }
            struct Node map[taskid][P][L];
            omp_set_num_threads(NTHREADS);
#pragma omp parallel
            {
                
                int thread_num = omp_get_thread_num();
                searchControl(thread_num, p, seedType, map[taskid], clusterMap, bRight, bDown , offset, taskid);
                printf("1. thread num = %d\n", thread_num);
            }
            // Send my results back to the master task */
            dest = MASTER;
            MPI_Send(&offset, 1, MPI_INT, dest, tag4, MPI_COMM_WORLD);
            MPI_Send(&clusterMap[offset][0], partition , MPI_INT, dest, tag1, MPI_COMM_WORLD);
            if(seedType == 'b'){
                MPI_Send(&bDown[taskid][0], P , MPI_INT, dest, tag2, MPI_COMM_WORLD);
                MPI_Send(&bRight[taskid][0], L , MPI_INT, dest, tag3, MPI_COMM_WORLD);
            }
            
        } /* end of non-master */
        MPI_Finalize();
        return 0;
        exit(EXIT_SUCCESS);
    }
}
