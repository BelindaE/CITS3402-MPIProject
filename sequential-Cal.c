#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"

//  compile:  gcc -fopenmp -o s sequential1.c
//  run:      ./s seedType probability
// this version changes the Node map back to a non-pointer to enable new seeded maps

#define L 8             /* Linear edge dimension of map */
#define N L*L           // total number of nodes in a single map
#define NTHREADS 4

struct Node {
    int flag;                       // occupied=1, unoccupied=0, already in a cluster=2+
    int up;                         // (0=not visited,1=visited,2=all directions visited,-1=when parent node)
    int down;
    int left;
    int right;
    int x;                          //co-ords to be easily accessible by cluster
    int y;
};

struct Node map[L][L];

struct Stack{
    int x;
    int y;
    char dir;
} cluster[N];

struct Percolate{
    int ID;
    int size;
    char rows;
    char cols;
} perc[N];

struct Interface{
    char rows;
    char cols;
    int size;
} interface[4];

struct Combination{
    int a[L*4];
    int b[L*4];
    int c[L*4];
    int d[L*4];
    int aSize;
    int bSize;
    int cSize;
    int dSize;
    int maxSize;
}combination;
struct Clusters{
    int m[L][L];
}clusterMap[4];
struct bonds{
    int flag;               // occupied=1, unoccupied=0, already in a cluster=2+
    int up;                 // (0=not visited,1=visited,2=all directions visited,-1=when parent node)
    int down;
    int left;
    int right;
    
}bondMap[NTHREADS][L][L];

// Initialising the side arrays to 0-empty. Will be filled with clusterID numbers as they are found during the search/navigation.
int upSide[L]={0};
int downSide[L]={0};
int leftSide[L]={0};
int rightSide[L]={0};
int upSideP[L*2]={0};
int downSideP[L*2]={0};
int leftSideP[L*2]={0};
int rightSideP[L*2]={0};

// initialise arrays to store whether cluster has spanned all rows and/or all columns
int rowPerc[L]={0};
int columnPerc[L]={0};

int top = 0;
int count = 0;
bool popped = false;

int thread_num;
int currentID;
int idCounter[N*2] = {0};
int idSpace[N*2] = {0};

/*      this functions seeds the map with occupied or unoccupied nodes depending on probability entered
 by the user and initialises all node variables with 0.
 */
void seedSite(double probability)
{
    double r;
    int seeded;
    int i, j;
    
    for (i = 0; i < L; i++){  	                    //rows
        for (j = 0; j < L; j++){                //columns
            
            map[i][j].flag = 0;				//reset flag to 0
        }
    }
    
    for (i = 0; i < L; i++){  	                    //rows
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
            printf("%d",map[i][j].flag);
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
        printf("\n");
    }
    printf("\n");
}

void seedBond(double probability)
{
    double r1, r2;
    int seeded1, seeded2;
    int i, j;
    
    for (i = 0; i < L; i++){  	                    //rows
        for (j = 0; j < L; j++){                //columns
            
            map[i][j].flag = 0;				//reset flag to 0
        }
    }
    
    for (i = 0; i < L; i++){  	                    //rows
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
            if (i+1 < L) map[i+1][j].up = seeded2;
            if (map[i][j].down == 0) {
                map[i][j].flag = 1;
                if (i+1 < L) map[i+1][j].flag = 1;
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
            printf("%d= (%d,%d) = l%d, r%d, U%d, D%D\n",map[i][j].flag,i,j,map[i][j].left,map[i][j].right,map[i][j].up,map[i][j].down);
            bondMap[thread_num][i][i].flag = map[i][j].flag;
            bondMap[thread_num][i][i].up = map[i][j].up;
            bondMap[thread_num][i][i].down = map[i][j].down;
            bondMap[thread_num][i][i].left = map[i][j].left;
            bondMap[thread_num][i][i].right = map[i][j].right;
        }
    }
    
}

void printSites(){
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
void printBonds(){
    
    int i,j;
    for (i = 0; i < (L); i++){
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

struct Node getNextNode() {
    popped = false;
    
    int j = peek().x ;
    int i = peek().y ;
    //printf("check1 (%i, %i)\t",i,j);
    
    /*the following first checks if flag on, then path is 0 and if 1 it also checks previous direction to make sure not cycling back on itself.  If no options available then 2 is applied to the node and stack drops the top node
     note that I have not assumed the lattice can cycle around itself as I think we want to retain edge nodes assuming it's a small part of a larger area - can change this easily with %L if required*/
    
    
    
    if ((i+1 < L) && map[i+1][j].flag == 1 && map[i+1][j].up == 0) {
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
    
    else if ((i+1 < L) && map[i+1][j].flag == 1 && map[i+1][j].up == 1
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
        clusterMap[thread_num].m[peek().y][peek().x] = currentID;
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
int depthFirstSearch(int i, int j, int clusterID){
    
    //the cluster ID will be the reference point for the side-arrays, and cluster sizes.
    
    top = 0;
    int tempX = 0;
    int tempY = 0;
    int tempC = 0;
    
    if (map[i][j].flag == 0 || map[i][j].flag >1) return 0;
    // check origin flag before commencing DFS,
    //aborting if either empty or already searched.
    count = 0;
    
    
    // reset cluster percolation arrays to zero
    int q;
    for (q=0; q<L; q++){
        rowPerc[q]=0;
        columnPerc[q]=0;
    }
    
    
    while(!isStackEmpty()){
        
        //initialise temporary holding stack
        struct Node temp[1][1];
        temp[0][0].x = 0;
        temp[0][0].y = 0;
        
        cluster[top].x = map[i][j].x;      //store co-ordinates in the cluster stack
        cluster[top].y = map[i][j].y;
        
        rowPerc[i] = clusterID;            //store record of which row and column cluster reaches
        columnPerc[j] = clusterID;
        
        if(map[i][j].flag == 1){
            map[i][j].flag = clusterID;
            //marking the node as both explored, and to which cluster of nodes it belongs-
            count++;         //will only update if there's a new filled node attached.
        }
        displayNode();                              // print check of top node in the cluster
        push(cluster[top].y, cluster[top].x);                   // push node into the stack
        tempX = cluster[top].x;
        tempY = cluster[top].y;
        tempC = currentID;
        
        temp[0][0] = getNextNode();                  // allocate new node coords
        i = temp[0][0].y;
        j = temp[0][0].x;
        int f = temp[0][0].flag;
        
        if(tempC != f)
        {
            clusterMap[thread_num].m[tempY][tempX] = f;
            currentID = f;
        }
        if(f != 1)
        {
            clusterMap[thread_num].m[i][j] = f;
            currentID = f;
        }
        
        if (!popped) ++top;            //want to go back to previous node if no connection made
        
    }
    idCounter[currentID] = count;
    idSpace[currentID] = thread_num;
    return count;
    
}


/*	layout of 4 x map array combinations
	0	1	2	3	4	5	6	7						0	1	2	3	4	5	6	7
 31									8				31									8
 30									9				30									9
 29									10				29									10
 28			combination.a	  		11				28			combination.b			11
 27									12		R0		27									12		R1
 26									13				26									13
 25									14				25									14
 24									15				24									15
	23	22	21	20	19	18	17	16						23	22	21	20	19	18	17	16
 
 C0														C1
 
	0	1	2	3	4	5	6	7						0	1	2	3	4	5	6	7
 31									8				31									8
 30									9				30									9
 29									10				29									10
 28									11		R2		28									11		R3
 27			combination.c			12				27			combination.d			12
 26									13				26									13
 25									14				25									14
 24									15				24									15
	23	22	21	20	19	18	17	16						23	22	21	20	19	18	17	16
 
 
 C2														C3
 */
void matchClusters(float probability, char seedType, struct Combination combination){
    int n;
    int max = 0;
    int rows = 0;
    int cols = 0;
    int both = 0;
    int connected[2][N*2];
    int biggestCluster[N*2] = {0};
    for(int i = 0; i<N*2;i++)
    {
        connected[0][i] = 0;
        connected[1][i] = 0;
    }
    int c = 0;
    for (n=0; n<L; n++){
        if (clusterMap[0].m[n][L-1] > 0 && clusterMap[1].m[n][0] > 0){
            if(seedType == 'b' && bondMap[0][n][L-1].right == 2) continue; // checking seeded bond between
            connected[0][c] = clusterMap[0].m[n][L-1];
            connected[1][c] = clusterMap[1].m[n][0];
            c++;
        }//R0
        if (clusterMap[0].m[0][n] > 0 && clusterMap[2].m[L-1][n] > 0){
            if(seedType == 'b' && bondMap[2][L-1][n].down == 2) continue;
            connected[0][c] = clusterMap[0].m[0][n] ;
            connected[1][c] = clusterMap[2].m[L-1][n];
            c++;
        } //C2
        if (clusterMap[0].m[L-1][n] > 0 && clusterMap[2].m[0][n] > 0){
            if(seedType == 'b' && bondMap[0][L-1][n].down == 2) continue;
            connected[0][c] = clusterMap[0].m[L-1][n];
            connected[1][c] = clusterMap[2].m[0][n];
            c++;
            
        }//C0
        if (clusterMap[0].m[n][0] > 0 && clusterMap[1].m[n][L-1] > 0){
            if(seedType == 'b' && bondMap[1][n][L-1].right == 2) continue;
            connected[0][c] = clusterMap[0].m[n][0];
            connected[1][c] = clusterMap[1].m[n][L-1];
            c++;
        }//R1
        if (clusterMap[2].m[n][L-1] > 0 && clusterMap[3].m[n][0] > 0){
            if(seedType == 'b' && bondMap[2][n][L-1].right == 2) continue;
            connected[0][c] = clusterMap[2].m[n][L-1];
            connected[1][c] = clusterMap[3].m[n][0];
            c++;
        }//R2
        if (clusterMap[2].m[n][0] > 0 && clusterMap[3].m[n][L-1] > 0){
            if(seedType == 'b' && bondMap[3][n][L-1].right == 2) continue;
            connected[0][c] = clusterMap[2].m[n][0];
            connected[1][c] = clusterMap[3].m[n][L-1];
            c++;
        } //R3
        if (clusterMap[1].m[L-1][n] > 0 && clusterMap[3].m[0][n] > 0){
            if(seedType == 'b' && bondMap[1][L-1][n].down == 2) continue;
            connected[0][c] = clusterMap[1].m[L-1][n];
            connected[1][c] = clusterMap[3].m[0][n];
            c++;
        }//C1
        if (clusterMap[1].m[0][n] > 0 && clusterMap[3].m[L-1][n] > 0){
            if(seedType == 'b' && bondMap[1][L-1][n].down == 2) continue;
            connected[0][c] = clusterMap[1].m[0][n];
            connected[1][c] = clusterMap[3].m[L-1][n];
            c++;//C3
        }
    }
    /*printf("\n\nc = %d\n\n",c);
    for(int i = 0 ; i<c ; i++)
    {
        printf("%d,%d\n",connected[0][i],connected[1][i]);
    }*/
    for(int i = 0; i < 4 ; i++)
    {
        if(connected[0][0] == 0 && i ==0)
        {
            break;
        }
        int tempCluster[N*2] =  {0};
        int k = 0;
        int tempRows[N*2] = {0};
        int tempCols[N*2] = {0};
        int tempMax = 0;
        int colP = 0;
        int rowP = 0;
        for(int j = 0 ; j < c; j++)
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
                for(int l = 0; l < k; l++)
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
        for(int j = 0; j < L; j++)
        {
            for(int m = 0; m < L ;m++)
            {
                for(int l = 0; l < k ; l++)
                {
                    if(clusterMap[0].m[j][m] == tempCluster[l] )
                    {
                        tempRows[j] = 1;
                        tempCols[m] = 1;
                    }
                    if(clusterMap[2].m[j][m] == tempCluster[l])
                    {
                        tempRows[j+L] = 1;
                        tempCols[m] = 1;
                    }
                    if(clusterMap[1].m[j][m] == tempCluster[l])
                    {
                        tempRows[j] = 1;
                        tempCols[m+L] = 1;
                    }
                    if (clusterMap[3].m[j][m] == tempCluster[l]) {
                        tempRows[j +L] = 1;
                        tempCols[m+L] = 1;
                    }
                }
            }
        }
        for(int j = 0; j < L*2;j++)
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
            for(int j = 0 ; j < c ; j++) biggestCluster[j] = tempCluster[j];
        }
        if(colP == 0) cols = 1;
        if(rowP == 0) rows = 1;
        if(cols ==1 && rows == 1) both =1;
    }
    for(int j = 0 ; j < c ; j++) printf("%d,",biggestCluster[j]);
    printf("\n");
    for(int i = 0; i < L*2; i++)
    {
        for( int j = 0; j  < L *2; j++)
        {
            if(j <L  && i < L && clusterMap[0].m[i][j] > 0){
                int q = 0;
                for(int l = 0 ; l<c ; l++)
                {
                    if(biggestCluster[l] ==clusterMap[0].m[i][j])
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
            else if(j >=L  && i < L && clusterMap[1].m[i][j - L] > 0){
                int q = 0;
                for(int l = 0 ; l<c ; l++)
                {
                    if(biggestCluster[l] ==clusterMap[1].m[i][j- L])
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
            else if(j <L  && i >= L && clusterMap[2].m[i -L][j] > 0) {
                int q = 0;
                for(int l = 0 ; l<c ; l++)
                {
                    if(biggestCluster[l] ==clusterMap[2].m[i-L][j])
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
            else if(j >= L  && i >= L && clusterMap[3].m[i-L][j-L] > 0) {
                int q = 0;
                for(int l = 0 ; l<c ; l++)
                {
                    if(biggestCluster[l] ==clusterMap[3].m[i-L][j-L])
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
    
    printf("Lattice size = %i x %i \n", L*2, L*2);
    printf("Seeding type = %c, Probability = %f \n\n", seedType, probability);
    
    printf("columns percolated = %d \n", cols);
    printf("rows percolated = %d \n", rows);
    printf("both percolated = %d \n", both);
    printf("Biggest Cluster = %d \n",max);

    
    printf("\n");
}

void searchControl(double probability, char seedType){
    
    int dfs;
    
    struct Combination combination;
    
    for (int g = 0; g < L*4; g++){
        combination.a[g] = 0;
        combination.b[g] = 0;
        combination.c[g] = 0;
        combination.d[g] = 0;
    }
    
    int i, j, n, k;
    
    int clusterID = 2;
    char row = 'y';
    char column = 'y';
    
    //omp_set_num_threads(4);
    
    //#pragma omp parallel
    {
        //#pragma omp for private(clusterID, i, j, n, m, k, t, row, column, combination, size)
        for (k = 0; k<NTHREADS; k++)
        {
            thread_num=k;
            for(int r = 0 ; r < L ; r++)
            {
                for(int q = 0; q< L ; q++)
                {
                    clusterMap[k].m[r][q] = 0;
                }
            }
            if (seedType == 's') {
                seedSite(probability);
                //printSites();
            }
            else {
                seedBond(probability);
                printBonds();
                printf("\n\n");
            }
            
            {
                for(j =0; j<L; j++){
                    for(i =0; i<L; i++){
                        
                        dfs = depthFirstSearch(i,j, clusterID);
                        row = 'y';
                        column = 'y';
                        if(dfs>0){
                            // the following determines if the cluster has percolated over
                            //rows and/or columns and records details for each cluster
                            for (n=0; n<L; n++){
                                if (rowPerc[n]==0) row = 'n';
                                if (columnPerc[n]==0) column = 'n';
				                        }
                            perc[clusterID-2].ID = clusterID;
                            perc[clusterID-2].size = dfs;
                            perc[clusterID-2].rows = row;
                            perc[clusterID-2].cols = column;
                            clusterID++;
                        }
                    }
                }
            }
            printf("Cluster Map %d\n",thread_num);
            for(int r = 0 ; r < L ; r++)
            {
                for(int q = 0; q< L ; q++)
                {
                    printf("%d\t",clusterMap[thread_num].m[r][q]);
                }
                printf("\n");
            }
        }
    }
    
    for(i=0;i<(clusterID-2);i++){
        //printf("Cluster No. %i, size: %i, rows percolating: %c, columns percolating: %c \n",
               //perc[i].ID, perc[i].size, perc[i].rows, perc[i].cols);
    }
    //printf("\n");
    
}

int main(int argc, char *argv[]){
    
    if(argc != 3)
    {
        fprintf(stderr, " ERROR\n Usage: ./perc seedType probability\n");
        // example: "./perc s 0.55"
        exit(EXIT_FAILURE);             // Exit program indicating failure
    }
    else
    {
        char seedType = argv[1][0];			// enter 's' for site percolation, 'b' for bond percolation
        double p = atof(argv[2]);			// enter probability between 0 - 1, eg. 0.55
        
        double delta;
        
        srand(time(NULL));                              //should only be called once
        
        struct timeval start, end;
        gettimeofday(&start, NULL);
        
        searchControl(p, seedType);
        
        matchClusters(p, seedType, combination);
        gettimeofday(&end, NULL);
        delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                 end.tv_usec - start.tv_usec) / 1.e6;
        printf("time=%12.10f\n",delta);
        
        return 0;
    }
}

