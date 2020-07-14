#include <search.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>
#include <assert.h>



//Definition of the structures
struct edge{
  unsigned long long int node_i;
  unsigned long long int node_j;
  unsigned long long int key;
  unsigned int total_overlap;
  unsigned int overlap_aggregate;
};

struct layer{
  struct edge *links_table;
  unsigned long long int  max_node;
  unsigned long long int nlinks;
  unsigned int l1;
  unsigned int l2;
};


void usage(char *argv[]){
  printf(
  "\n\n\n"
  "******************************************************************************\n"
  "**                                                                          **\n"
  "**       Rewire multiplex in order to reduce / increase the structural      **\n"
  "** overlap from CURRENT to MIN/MAX saving the multiplex every 0.1 intervals **\n"
  "**                                                                          **\n"
  "**                                                                          **\n"
  "**           <filein> should be a list of the layers composing the          **\n"
  "**                            multiplex network                             **\n"
  "**                                                                          **\n"
  "**                                                                          **\n"
  "**             <unique_nodes> list of the nodes of the multiplex            **\n"
  "**                                                                          **\n"
  "**                                                                          **\n"
  "**  <I/D> I -> for increasing the overlap; D -> for decreasing the overlap  **\n"
  "**                                                                          **\n"
  "**                     ----   Optional Variables  ----                      **\n"
  "**                                                                          **\n"
  "**     <-t #thresh> represents the threshold value when computing the       **\n"
  "**    values of the structural overlap (i.e. |o_current -o* | < TOLL )      **\n"
  "**                                                                          **\n"
  "**   <-m #max_overlap> if inserted, this represents the maximum value of    **\n"
  "**             overlap desired.  After reaching the condition:              **\n"
  "**      |o_current -o_{maxoverlap} | < TOLL, the algorithm will stop        **\n"
  "**                                                                          **\n"
  "******************************************************************************\n");
  printf("Usage: %s <filein> <unique_nodes>  I/D [-t #thesh=0.001] [-m max_overlap=none] \n\n" , argv[0]);
}


unsigned long long int *load_unique_nodes(char *, unsigned long long int *, unsigned long long int *);

struct layer *read_slap(char *, unsigned long long int *, unsigned long long int , unsigned int, unsigned int *,
  unsigned long long int ***, unsigned long long int ***, unsigned long long int **, unsigned long long int **);

unsigned int read_ij(FILE *, unsigned long long int *,unsigned long long int, unsigned long long int **, unsigned long long int **,struct edge **);


int convert_ij2slap(unsigned long long int *, unsigned long long int ,
  unsigned long long int *, unsigned long long int *, unsigned  long long int, unsigned  long long int,
  unsigned long long int **, unsigned long long int **);

unsigned int find_max(unsigned long long int *, unsigned long long int );


void list_prime(unsigned long long int **,unsigned int );

int compare_keys(const void *, const void *);

double map_in_unitary(unsigned int,double);

void free_structure(unsigned long long int *,struct layer *, unsigned int,unsigned long long int **,
  unsigned long long int **,unsigned long long int *,unsigned long long int *);

unsigned int binary_search(unsigned long long int *,unsigned long long int, unsigned long long int );

unsigned long long int check_overlap(unsigned long long int, unsigned long long int, unsigned long long int,
  unsigned long long int **, unsigned long long int **,unsigned long long int *, unsigned int);

int is_neigh(unsigned long long int *, unsigned long long int *, unsigned long long int,
  unsigned long long int, unsigned long long int);

unsigned int check_existence(unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int,
  unsigned long long int  **,unsigned long long int **, unsigned long long int);

int evaluate_Delta_denom(unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int,unsigned long long int **,
  unsigned long long int **,unsigned long long int *, unsigned int, unsigned long long int);

void change_link(unsigned long long int, unsigned long long int *, unsigned long long int *, unsigned long long int, unsigned long long int);

int cmpfunc (const void *, const void *);

int delete_root(const void *, const void *);

void tdestroy(void *);

void write_edges_onfiles( unsigned long long int  **, unsigned long long int **, unsigned int, unsigned long long int, double);

void write_edges_onfiles_top_to_bottom( unsigned long long int  **, unsigned long long int **, unsigned int, unsigned long long int, double);

void write_edges_onfiles_bottom_to_top( unsigned long long int  **, unsigned long long int **, unsigned int, unsigned long long int, double);

unsigned long long int count_number_edge_aggregate(struct layer *,unsigned int, unsigned long long int);

void decrease_overlap(unsigned long long int  **J_slap_layers, unsigned long long int **r_slap_layers, unsigned long long int *v_slap,
                unsigned long long int *r_slap_total, unsigned long long int N, unsigned int M, unsigned long long int *nodes_unique,
                struct layer *multiplex,unsigned long long int length_unique, unsigned int flag_decrease_increase,double TOLL);



void increase_overlap(unsigned long long int  **J_slap_layers, unsigned long long int **r_slap_layers, unsigned long long int *v_slap,
                unsigned long long int *r_slap_total, unsigned long long int N, unsigned int M, unsigned long long int *nodes_unique,
                struct layer *multiplex,unsigned long long int length_unique, unsigned int flag_decrease_increase,double TOLL,double max_overlap);




int mpz_cmp_flag_to_01(int mpz_cmp_flag){
  if(mpz_cmp_flag <= 0)
    return 0;
  else
    return 1;
}



int main(int argc, char *argv[])
{
  //List of variables used in the main
  unsigned int M,i,j,k,flag_decrease_increase=-1;
  unsigned long long int *r_slap_total,**J_slap_layers,**r_slap_layers,*v_slap,*nodes_unique,length_unique, MAX_NODE,N;
  struct layer *multiplex;
  double TOLL=0.001,max_overlap=-1;


  if (argc <= 3){
    usage(argv);
    exit(1);
  }

  //Check if the procedure will have to increase or decrease the overlap
  if (isalpha(argv[3][0]))
  {
    if (argv[3][0] == 'D' || argv[3][0] == 'd')
    {
      flag_decrease_increase=0;
      fprintf(stderr,"Decreasing the overlap!\n");
    }
    else
    if (argv[3][0] == 'I' || argv[3][0] == 'i')
    {
      flag_decrease_increase=1;
      fprintf(stderr,"Increasing the overlap!\n");
    }
    if (flag_decrease_increase == -1)
    {
      fprintf(stderr,"Char inserted not valid! (Use instead I or D)\n EXIT NOW!\n");
      exit(1);
    }
  }

  for (i = 1; i < argc; i++) 
  {
    if (argv[i][0] == '-'  && argv[i][1] == 't' || argv[i][1] == 'T') 
    {
      if (argc > i+1)
        TOLL = (double)atof(argv[i+1]);
      fprintf(stderr,"Tollerance inserted: %lf\n", TOLL);
    }
    if (argv[i][0] == '-'  && argv[i][1] == 'm' || argv[i][1] == 'M') 
    {
      if (argc > i+1)
        max_overlap = (double)atof(argv[i+1]);
      fprintf(stderr,"Max overlap inserted: %lf\n", max_overlap);
    }
  }
  // //Check if a tollerance value and/or a max overlap value has been inserted
  // if (argc > 4)
  // {
  //   if (argv[4][0] !='-')
  //   {
  //     TOLL = (double)atof(argv[4]);
  //     fprintf(stderr,"Tollerance inserted: %lf\n", TOLL);
  //   }
  // }

  //Randomize the seed
  srand(time(NULL));

  //Load in nodes_unique the list of nodes in the multiplex network
  nodes_unique=load_unique_nodes(argv[2],&length_unique, &N);

  //Create the multiplex structure using a compress-row format on each layer (J_slap_layers,r_slap_layers).
  //It also determines the number of layers and map each node of the multiplex into the index of nodes_unique using
  //a binary search 
  multiplex=read_slap(argv[1],nodes_unique, length_unique, N, &M, &J_slap_layers, &r_slap_layers,&r_slap_total,&v_slap);
  //Biased random rewiring to increase (resp. decrease) the structural overlap of a multiplex networl
  //with M layers and same number of nodes N, yet maintaining the same degree sequence on each layer.
  if(flag_decrease_increase == 0)
    decrease_overlap(J_slap_layers, r_slap_layers, v_slap, r_slap_total,N, M,nodes_unique,multiplex,length_unique,flag_decrease_increase,TOLL);
  else
    increase_overlap(J_slap_layers, r_slap_layers, v_slap, r_slap_total,N, M,nodes_unique,multiplex,length_unique,flag_decrease_increase,TOLL,max_overlap);
  free_structure(nodes_unique,multiplex,M,J_slap_layers,r_slap_layers,v_slap,r_slap_total);
}


unsigned long long int *load_unique_nodes(char *filename, unsigned long long int *length_unique, unsigned long long int *MAX_NODE)
{
  FILE *fp;
  unsigned long long int *unique_nodes, size=10000, k=0;
  int scan=0;
  fp = fopen(filename,"r");
  unique_nodes=(unsigned long long int *)malloc(size*sizeof(unsigned long long int));
  while(scan!=EOF)
  {
    if(k==size)
    {
      size+=10000;
      unique_nodes=(unsigned long long int *)realloc(unique_nodes,size*sizeof(unsigned long long int));
    }
    scan=fscanf(fp,"%llu",&(unique_nodes[k]));
    k++;
  }
  k=k-1;
  unique_nodes=(unsigned long long int *)realloc(unique_nodes,k*sizeof(unsigned long long int));
  fprintf(stderr,"Number of unique nodes: %llu\n\n",k);
  fclose(fp);
  *MAX_NODE = unique_nodes[k-1];
  *length_unique=k;
  return(unique_nodes);
}



struct layer *read_slap(char *multiplex_filename, unsigned long long int *nodes_unique, unsigned long long int length_unique, unsigned int N, unsigned int *M_total_layer,
  unsigned long long int ***J_slap_layers, unsigned long long int ***r_slap_layers, unsigned long long int **r_slap_total, unsigned long long int **v_slap)
{
  unsigned long long int *I,*J;
  unsigned long long int i, k,counter_layer,K,numberoflinks_layer, cumulative,M;
  M=10;
  FILE *fp_multiplex,*fp_layer;
  char buff[256], *ptr, *token, line[100], stringa[200];
  struct  layer *array;


  array = (struct layer *)malloc(M*sizeof(struct layer ));
  fp_multiplex= fopen(multiplex_filename,"r");
  counter_layer=0;
  cumulative=0;
  *v_slap = (unsigned long long int *) calloc( (M+1), sizeof(unsigned long long int));
  *r_slap_total=(unsigned long long int *) malloc(M * sizeof(unsigned long long int));
  *r_slap_layers = (unsigned long long int **) malloc( M * sizeof(unsigned long long int *));
  *J_slap_layers = (unsigned long long int **) malloc( M * sizeof(unsigned long long int *));
  (*v_slap)[0]=0;
    while (1)
    {
      if (counter_layer == M-1)
      {
        M+=10;
        array = (struct layer *)realloc(array,M*sizeof(struct layer ));
        *v_slap = (unsigned long long int *)realloc(*v_slap, (M+1)* sizeof(unsigned long long int));
        *r_slap_total=(unsigned long long int *) realloc(*r_slap_total,M * sizeof(unsigned long long int));
        *r_slap_layers = (unsigned long long int **) realloc(*r_slap_layers, M * sizeof(unsigned long long int *));
        *J_slap_layers = (unsigned long long int **) realloc(*J_slap_layers, M * sizeof(unsigned long long int *));
      }
      I=NULL;
      J=NULL;
      //If it reaches the end of the file then break!
      if (fgets(line,100, fp_multiplex) == NULL) break;
      //token = strtok(NULL,"\n");
      //Read the name of the layer composing the multiplex list and print on the terminal and on the string "stringa"
      token = strtok(line,"\n");
      fprintf(stderr,"%s  ",token);
      sprintf(stringa,"%s",token);
      fp_layer= fopen(stringa,"r");
      k = read_ij(fp_layer, nodes_unique, length_unique, &I, &J,&(array[counter_layer].links_table));
      array[counter_layer].nlinks=k;
      K = 2 * k;
      cumulative+=K;
      I = realloc(I, K * sizeof(unsigned long long  int));
      J = realloc(J, K * sizeof(unsigned long long  int));
      for (i=k; i<K; i++)
      {
        I[i] = J[i-k];
        J[i] = I[i-k];
      }
      numberoflinks_layer =convert_ij2slap(nodes_unique, length_unique, I, J, K, M,&(*r_slap_layers)[counter_layer], &(*J_slap_layers)[counter_layer]);
      //printf("N_max: %u -- %u\n",numberoflinks_layer,counter_layer);
      fclose(fp_layer);
      fprintf(stderr,"--> Number of links: %llu\n",K);
      (*r_slap_total)[counter_layer]=numberoflinks_layer;
      array[counter_layer].max_node=numberoflinks_layer;
      counter_layer++;
      free(I);
      free(J);
      (*v_slap)[counter_layer]=cumulative;
      if(counter_layer==M) break;
    }
    fclose(fp_multiplex);
    M = counter_layer;
    *M_total_layer = M;
    array = (struct layer *)realloc(array,M*sizeof(struct layer ));
    *v_slap = (unsigned long long int *)realloc(*v_slap, (M+1)*sizeof(unsigned long long int));
    *r_slap_total=(unsigned long long  int *) realloc(*r_slap_total,M * sizeof(unsigned long long int));
    *r_slap_layers = (unsigned long long int **) realloc(*r_slap_layers, M * sizeof(unsigned long long int *));
    *J_slap_layers = (unsigned long long int **) realloc(*J_slap_layers, M * sizeof(unsigned long long int *));

    fprintf(stderr,"Double of the number of total link: %llu\n",(*v_slap)[counter_layer]);
    for (i = 0; i < M+1; i++)
    {
      fprintf(stderr,"%llu  ",(*v_slap)[i] );
    }
    fprintf(stderr,"\n");
    fprintf(stderr, "Number of layers: %llu\n",M);
    return(array);
}


/*
 * Read a file in ij format
 */
unsigned int read_ij(FILE *filein, unsigned long long int *unique_nodes,unsigned long long int length_unique, unsigned long long int **I, unsigned  long long int **J,
          struct edge **current){

  unsigned long long int size, K,flag=0,temp1,temp2,auxiliary,*item;
  char buff[256];
  char *ptr;
  unsigned long long int *p_a=&unique_nodes[0];
  unsigned long long int *p_b;
  unsigned long long int differenceInBytes;

  size = 10000;
  K = 0;
  *(current)=(struct edge *) malloc(size * sizeof(struct edge ));
  *I = malloc(size * sizeof(unsigned long long int));
  *J = malloc(size * sizeof(unsigned long long int));
      while(fgets(buff, 256, filein)){
        if (buff[0] == '#')
          continue;
        if (K == size){
          size += 10000;
          *(current)=(struct edge *) realloc(*(current),size * sizeof(struct edge ));
          *I = realloc(*I, size*sizeof(unsigned long long int));
          *J = realloc(*J, size*sizeof(unsigned long long int));
        }
        ptr = strtok(buff, " "); /* read the first node */
        temp1 = atoi(ptr);
        ptr = strtok(NULL, " "); /* read the second node */
        temp2 = atoi(ptr);
        // PARSING DATA WITH NODES_UNIQUE -> Take the corresponding index -1 of that element
        item = (unsigned long long int *) bsearch (&temp1, unique_nodes, length_unique, sizeof (unsigned long long int), cmpfunc);
        p_b=item;
        differenceInBytes = (p_b - p_a);
        temp1=differenceInBytes;
        (*I)[K] =temp1;

        item = (unsigned long long int *) bsearch (&temp2, unique_nodes, length_unique, sizeof (unsigned long long int), cmpfunc);
        p_b=item;
        differenceInBytes = (p_b - p_a);
        temp2=differenceInBytes;
        (*J)[K] = temp2;


        if(temp1>temp2)
        {
          auxiliary=temp2;
          temp2=temp1;
          temp1=auxiliary;
        }
        (*current)[K].node_i = temp1;
        (*current)[K].node_j = temp2;
        /////////////////////////////////////////
        //printf("nodei: %llu    nodej: %llu\n",(*I)[K],(*J)[K] );
        K += 1;
   // if(K==1 && flag ==0)
   //     {K=0; flag=1;}
  }

  *I = realloc(*I, K * sizeof(unsigned long long int));
  *J = realloc(*J, K * sizeof(unsigned long long int));
  *(current)=(struct edge *) realloc(*(current),K * sizeof(struct edge ));
  return K;
}



int convert_ij2slap(unsigned long long int *unique_nodes, unsigned long long int length_unique,
                    unsigned long long int *I, unsigned long long int *J, unsigned long long int K, unsigned  long long int  M,
                    unsigned long long int **r_slap, unsigned long long int **J_slap){

  unsigned long long int tmp, max;
  unsigned long long int N;
  unsigned long long int i, pos;
  unsigned long long int *p;

  max = find_max(I, K) + 1;
  tmp = find_max(J, K) + 1;
  if (tmp > max){
    max = tmp ;
  }
  //fprintf(stderr,"Max inside the layer: %u\n\n",max);
  (*r_slap)=(unsigned long long int *) malloc( (max+1) * sizeof(unsigned long long int));
  p = (unsigned long long int *)malloc(max * sizeof(unsigned long long int));

  *J_slap =(unsigned long long int *) malloc(K * sizeof(unsigned long long int));
  memset(*r_slap, 0, (max+1) * sizeof(unsigned long long int));
  for(i=0; i<max + 1; i++)
    (*r_slap)[i] = 0;
  memset(p, 0, max * sizeof(unsigned long long int));
  (*r_slap)[0] = 0;
  //fprintf(stderr, "WARNING!!!! R_SLAP[0] NOW IS SET TO ZERO!!!!!\n");
  for(i=0; i<K; i++){
    (*r_slap)[ I[i] + 1] += 1;
  }
  for(i=1; i<=max; i++){
    (*r_slap)[i] += (*r_slap)[i-1];
  }
  for(i=0; i<K; i++){
    pos = (*r_slap) [ I[i] ] + p[ I[i] ];
    (*J_slap)[pos] = J[i];
    p[ I[i] ] += 1;
  }
  free(p);
  return max;
}


unsigned int find_max(unsigned long long int *v, unsigned long long int N){

  unsigned long long int i, max;

  max = v[0];
  i = 0;
  while(++i < N){
    if (v[i] > max)
      max = v[i];
  }
  return max;
}


int cmpfunc (const void * a, const void * b)
{
   return ( *(long long int*)a - *(long long int*)b );
}





//Count the number of distinct edge in the multiplex structure having M number of layers. The counting procedure is done using a tree structure (tsearch,tfind)
unsigned long long int count_number_edge_aggregate(struct layer *multiplex,unsigned int M, unsigned long long int MAX_NODE)
{

  unsigned long long int i,j,k, nodei,nodej,temp,*shuffle, *distr1,distinct_edges=0;
  void *root=NULL;
  struct edge *update;

  //Creating a binary tree where the key of each leaf is given by key = MAX_NODE * nodei + nodej, where nodei < nodej
  for (k=0;k<M;k++)
  {
    for (i = 0; i < multiplex[k].nlinks; i++)
    {
      nodei=multiplex[k].links_table[i].node_i;
      nodej=multiplex[k].links_table[i].node_j;
      if(nodei > nodej)
      {
        multiplex[k].links_table[i].node_i = multiplex[k].links_table[i].node_j;
        multiplex[k].links_table[i].node_j = nodei;
        nodei=multiplex[k].links_table[i].node_i;
        nodej=multiplex[k].links_table[i].node_j;
      }

      multiplex[k].links_table[i].key= MAX_NODE * nodei + nodej;
      update=tfind (& multiplex[k].links_table[i], &root, compare_keys);

      // If the element does not belong to the tree, then update is egual to NULL and the element will be
      // inserted with the function tsearch
      if( update == NULL )
      {
        (void) tsearch(& multiplex[k].links_table[i], &root, compare_keys);
        distinct_edges++;
      }
    }
  }

  //Destroying the binary tree
  tdestroy(root);
  return(distinct_edges);
}


void tdestroy(void *root)
{
   struct edge * elementptr;
    while (root != NULL) {
        elementptr = *(struct edge **)root;
        tdelete((void *)elementptr, &(root), delete_root);
    }
}

int delete_root(const void *node1, const void *node2)
{
    return 0;
}




int compare_keys(const void *e1p, const void *e2p)
{
  const struct edge *e1, *e2;
  int last, first;

  e1 = (const struct edge *) e1p;
  e2 = (const struct edge *) e2p;
  /* check key values */
  if (e1->key < e2->key)
    return -1;
  else if (e1->key == e2->key)
    return 0;
  else
    return 1;
}


double map_in_unitary(unsigned int M,double x)
{
  double new_point;
  new_point=((1.0*M)/(1.0*M-1.0))*x - (1.0)/(1.0*M-1.0);
  return(new_point);
}

void free_structure(unsigned long long int *nodes_unique,struct layer *multiplex, unsigned int M,
  unsigned long long int **J_slap_layers,unsigned long long int **r_slap_layers,unsigned long long int *v_slap,unsigned long long int *r_slap_total)
{
  unsigned long long int i,k;
  for (k=0;k<M;k++)
  {
    free(J_slap_layers[k]);
    free(r_slap_layers[k]);
    free(multiplex[k].links_table);
  }
  free(nodes_unique);
  free(multiplex);
  free(v_slap);
  free(r_slap_total);
  free(r_slap_layers);
  free(J_slap_layers);
}



unsigned int binary_search(unsigned long long int *array,unsigned long long int size_array,unsigned long long int target)
{
  unsigned long long int middle;
  unsigned long long int left = 0, right = size_array,counter=0;
  while (left != right)
  {
      middle = (left + right)/ 2;
      //printf("%u \n",middle);

      if (array[middle] > target)
       {
        //printf("first -> %u   --- %u\n",left,right);
        if (middle == 0)
          right = left;
        else
        right =middle -1;

    }
      else
        {
          left = middle;
          //printf("second -> %u   --- %u\n",left,right);
          if (left == right -1)
          {
            if(array[right] <= target)
            {
              //printf("%u   --- %u\n",left,right);
              left = left+1;
            }
            else
            {
              //printf("%u   --- %u\n",left,right);
              right=right -1;
            }
            //printf("%u   --- %u\n",left,right);
          }

        }
  }
/* Now, left and right should be the index before the target selected */
  //printf("I found %u between %u - %u - %u -> the target was: %u\n",left,array[left-1],array[left],array[left+1],target);
  return(left);
}

unsigned long long int check_overlap(unsigned long long int j1, unsigned long long int j2, unsigned long long int alpha,
      unsigned long long int **r_slap_layers, unsigned long long int **J_slap_layers,unsigned long long int *r_slap_total, unsigned int M)
{
  unsigned long long int counter_overlap=0,i,k,id_layers,j1pos;
  for (id_layers=0;id_layers<M;id_layers++)
  {
    //j1pos=binary_search(r_slap_layers[id_layers],r_slap_total[id_layers],j1);
    if (j2 < r_slap_total[id_layers])
    {
    //printf("Nel layer %u sto cercando il link (%u,%u)\n", id_layers,j1,j2);
    for(k = r_slap_layers[id_layers][j1]; k< r_slap_layers[id_layers][j1+1]; k++)
    {
      //printf("Il k iniziale e'%u corrispondente alla coppia (%u, %u)\n", k, j1, J_slap_layers[id_layers][k]);
      if (J_slap_layers[id_layers][k]==j2)
      {
        counter_overlap++;
        break;
      }
    }
    }

  }
  return(counter_overlap);

}

/* Check if j is a neighbour of i */
int is_neigh(unsigned long long int *J_slap, unsigned long long int *r_slap, unsigned long long int N,
             unsigned long long int i, unsigned long long int j){

  unsigned long long int l;
  unsigned long long int count;
  count = 0;
  if (i >=N || j >=N)
    return 0;
  for(l=r_slap[i]; l<r_slap[i+1]; l++){
    if (J_slap[l] == j)
      count ++;
  }
  return count;
}

unsigned int check_existence(unsigned long long int j1, unsigned long long int k1, unsigned long long int j2, unsigned long long int k2,
       unsigned long long int  **r_slap_layers,unsigned long long int **J_slap_layers, unsigned long long int alpha)
{
  unsigned long long int i,k,existence=0;
  for (i = r_slap_layers[alpha][j1]; i < r_slap_layers[alpha][j1+1]; i++)
  {
    if (J_slap_layers[alpha][i]==k1)
      existence++;
  }

  for (i = r_slap_layers[alpha][j2]; i < r_slap_layers[alpha][j2+1]; i++)
  {
    if (J_slap_layers[alpha][i]==k2)
      existence++;
  }
  return(existence);

}

//Evaluating the the variations in the denominator for the computation of the total edge overlap. In particular, if we consider the link j1-j2 and k1-k2, we then have these
//possibilities for the variation of edge overlap
//                                       - the aggregate network will remain the same  <-> 0 <-> if  o_j1j2 > 1
//  remove j1 - j2 from the multiplex  -
//                                       - the aggregate network will have 1 less link <-> -1 <-> if  o_j1j2  = 1
//
//                                       - the aggregate network will remain the same  <-> 0 <-> if  o_k1k2 > 1
//  remove k1 - k2 from the multiplex  -
//                                       - the aggregate network will have 1 less link <-> -1 <-> if  o_k1k2  = 1
//
//                                       - the aggregate network will remain the same  <-> 0 <-> if  o_j1k1 > 0
//     add j1 - k1 in the multiplex  -
//                                       - the aggregate network will have  1 more link <-> +1 <-> if  o_j1k1  = 1
//
//                                       - the aggregate network will remain the same  <-> 0 <-> if  o_j2k2 > 0
//     add j2 - k2 in the multiplex  -
//                                       - the aggregate network will have  1 more link <-> +1 <-> if  o_j2k2  = 1
//
//All these steps mean that the variation of the overlap \delta o ={-2,-1,0,1,2} will assume only these 5 values
int evaluate_Delta_denom(unsigned long long int k1, unsigned long long int k2, unsigned long long int j1, unsigned long long int j2,
  unsigned long long int **r_slap_layers, unsigned long long int **J_slap_layers,unsigned long long int *r_slap_total,
  unsigned int M,unsigned long long int overlap_j1j2)
{
  int i,k,id_layers,j1pos, flag1=0,flag2=0,swap,temp;
  int counter_overlapk1k2=0,counter_overlapk1j1=0,counter_overlapk2j2=0,diff,counter_overlapj1j2=0;


  //printf("Sto facendo il check per il link %d %d e %d %d -> %d %d e %d %d\n",j1,j2,k1,k2, nodes_unique[j1],nodes_unique[j2],nodes_unique[k1], nodes_unique[k2] );
  for (id_layers=0;id_layers<M;id_layers++)
  {
    flag1=0;
    flag2=0;

    if (k2 < r_slap_total[id_layers] && k1 < r_slap_total[id_layers] )
    {
      for(k = r_slap_layers[id_layers][k1]; k< r_slap_layers[id_layers][k1+1]; k++)
      {
        //printf("J_slap_layers: %d\n",J_slap_layers[id_layers][k]);
        if (J_slap_layers[id_layers][k]==k2)
        {
          //printf("Ho trovato il link %d %d\n",k1,J_slap_layers[id_layers][k]);
          counter_overlapk1k2++;
          flag1=1;
        }
        //if (J_slap_layers[id_layers][k]==j1 && nodes_unique[k1] <= nodes_unique[r_slap_total[id_layers]] && nodes_unique[j1] <= nodes_unique[r_slap_total[id_layers]] )
        if (J_slap_layers[id_layers][k]==j1  && j1 < r_slap_total[id_layers])
        {
          //printf("Sto entrando qui per il link:  %d-%d  ma in realtà non dovrei\n",nodes_unique[k1],nodes_unique[j1]);
          counter_overlapk1j1++;
          flag2=1;
        }
        if(flag1 == 1 && flag2==1)
          break;

      }
    }
  }

  if (j2 > k2)
  {
    swap = j2;
    j2=k2;
    k2=swap;
  }
  for (id_layers=0;id_layers<M;id_layers++)
  {
    if (k2 < r_slap_total[id_layers]  && j2 < r_slap_total[id_layers])
    {
      for(k = r_slap_layers[id_layers][j2]; k< r_slap_layers[id_layers][j2+1]; k++)
      {
        //if (J_slap_layers[id_layers][k]==k2  && nodes_unique[k2] <= nodes_unique[r_slap_total[id_layers]] && nodes_unique[j2] <= nodes_unique[r_slap_total[id_layers]])
        if (J_slap_layers[id_layers][k]==k2)
        {
          //printf("Sto entrando qui per il link:  %d-%d  ma in realtà non dovrei\n",nodes_unique[k2],nodes_unique[j2]);
          counter_overlapk2j2++;
          break;
        }

      }
    }

  }

  //Final count for the denominator
  if(counter_overlapk1k2 >1)  counter_overlapk1k2=0;
  else counter_overlapk1k2=-1;

  if(counter_overlapk1j1 > 0) counter_overlapk1j1=0;
  else counter_overlapk1j1=1;

  if(counter_overlapk2j2 > 0) counter_overlapk2j2=0;
  else counter_overlapk2j2=1;

  if(overlap_j1j2 > 1) counter_overlapj1j2=0;
  else counter_overlapj1j2=-1;

  diff= counter_overlapk1k2+counter_overlapk1j1 + counter_overlapk2j2 + counter_overlapj1j2;

  return(diff);
}

void change_link(unsigned long long int i, unsigned long long int *J_slap, unsigned long long int *r_slap,
                 unsigned long long int j1, unsigned long long int j2){

  unsigned long long int k;

  //Andrea's comment: fprintf(stderr, "Rewiring %d to %d in %d\n", j1, j2, i);
  for(k = r_slap[i]; k< r_slap[i+1]; k++){
    if (J_slap[k] == j1){
      J_slap[k] = j2;
    }
  }
}





void write_edges_onfiles( unsigned long long  int  **J_slap, unsigned long long int **r_slap, unsigned int M, unsigned long long int N, double structural_overlap_value){

  unsigned long long int i, j,t=0;
  FILE *fp;
  char string[200];
  for (t=0;t<M;t++){
    sprintf(string,"./layer%llu.txt_%.2lf",t+1,structural_overlap_value);
    fp=fopen(string,"w+");
    for(i=0; i<N; i++){
      for (j=r_slap[t][i]; j<r_slap[t][i+1]; j++){
        // printf("%llu  %llu N: %llu max: %llu %llu\n", i,j,N, r_slap[t][i+1], J_slap[t][j]);
            if(i<J_slap[t][j])
            fprintf(fp, "%llu  %llu\n", i, J_slap[t][j]);
        }
      }
      fclose(fp);
  }
}

void write_edges_onfiles_top_to_bottom( unsigned long long int  **J_slap, unsigned long long int **r_slap, unsigned int M, unsigned long long int N, double structural_overlap_value){

  unsigned long long int i, j,t=0;
  FILE *fp;
  char string[200];
  for (t=0;t<M;t++){
    sprintf(string,"./layer%llu.txt_%.2lf_top_to_bottom",t+1,structural_overlap_value);
    fp=fopen(string,"w+");
    for(i=0; i<N; i++){
      for (j=r_slap[t][i]; j<r_slap[t][i+1]; j++){
            if(i<J_slap[t][j])
            fprintf(fp, "%llu  %llu\n", i, J_slap[t][j]);
        }
      }
      fclose(fp);
  }
}

void write_edges_onfiles_bottom_to_top( unsigned long long int  **J_slap, unsigned long long int **r_slap, unsigned int M, unsigned long long int N, double structural_overlap_value){

  unsigned long long int i, j,t=0;
  FILE *fp;
  char string[200];
  for (t=0;t<M;t++){
    sprintf(string,"./layer%llu.txt_%.2lf_bottom_to_top",t+1,structural_overlap_value);
    fp=fopen(string,"w+");
    for(i=0; i<N; i++){
      for (j=r_slap[t][i]; j<r_slap[t][i+1]; j++){
            if(i<J_slap[t][j])
            fprintf(fp, "%llu  %llu\n", i, J_slap[t][j]);
        }
      }
      fclose(fp);
  }
}




void decrease_overlap(unsigned long long int  **J_slap_layers, unsigned long long int **r_slap_layers, unsigned long long int *v_slap,
                unsigned long long int *r_slap_total, unsigned long long int N, unsigned int M, unsigned long long int *nodes_unique,
                struct layer *multiplex,unsigned long long int length_unique, unsigned int flag_decrease_increase,double TOLL)
{
  unsigned long long int j1, j2, j2_pos;  /* The first couple of nodes */
  unsigned long long int k1, k2, k2_pos;  /* The second couple of nodes */
  unsigned long long int i,s,t,z, alpha, temp_swap,overlap_j1j2,existence_flag;
  unsigned long long int numerator_structuraloverlap,denominator_structuraloverlap,flag,SIZE=100,dd=0,nedges_intersection=0,Ndistinct_edges=0;

  int t_step=0;
  double omega,KC,structural_overlap, structural_overlap_prev,delta_numerical=0.01;
  long long int MAX_LINKS=0;
  // If it is not possible to obtain the minimum value of overlap when rewiring due 
  // to degree sequences constraints, I stop after a certain number of iterations MAX_ITERATIONS
  long long int MAX_ITERATIONS=500000000;
  int delta_denom,count =1,flagga=0,finalflag=0;

   if (flag_decrease_increase == 0)
    delta_numerical=-delta_numerical;

  for (i = 0; i< M;i++)
  {
    if (MAX_LINKS < multiplex[i].nlinks)
      MAX_LINKS= multiplex[i].nlinks;
  }

  numerator_structuraloverlap=v_slap[M]/2;
  Ndistinct_edges=count_number_edge_aggregate(multiplex,M,N);
  denominator_structuraloverlap=Ndistinct_edges;

  structural_overlap=(1.0*(numerator_structuraloverlap))/(M*denominator_structuraloverlap);
  structural_overlap_prev = roundf(map_in_unitary(M,structural_overlap) * 100)/100 + delta_numerical;
  fprintf(stderr,"Number of links in the multiplex: %llu --- Number of links in the aggregate: %llu\n",numerator_structuraloverlap, denominator_structuraloverlap);
  fprintf(stderr,"structural overlap: %lf",map_in_unitary(M,structural_overlap));  
  fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);

 
  
  while (denominator_structuraloverlap!=numerator_structuraloverlap)
  {
    //This is the procedure to efficiently sample a link T from the multiplex
    s = rand() % v_slap[M]; //Select the first node of couple T// --> Random number between 0 and max number of links
    alpha = binary_search(v_slap, M, s); // Find the layer index such that v_slap[alpha] < s < v_slap[alpha+1]
    j2_pos= s-v_slap[alpha]; // This indicates the position of the second node of link T in the r-j-slap structure of the alpha layer
    // Find the first node of link T in the r_slap structure of the layer alpha, it will find the index such that it is verified
    // that r_slap_layers[j1] < j2_pos < r_slap_layers[j1+1]
    j1 = binary_search(r_slap_layers[alpha],r_slap_total[alpha],j2_pos);
    j2 = J_slap_layers[alpha][j2_pos];
    if(j2 < j1)
    {
      temp_swap=j1;
      j1=j2;
      j2=temp_swap;
    }
    //printf("The link: (%u,%u)\n",j1,j2);
    //This is for the selection of the first link -> Do a check for the overlap!//
    overlap_j1j2=check_overlap(j1,j2,alpha,r_slap_layers,J_slap_layers,r_slap_total,M);

    flag =0;
    if(overlap_j1j2 >= 1)
    {
      //printf("The link is: (%u, %u) --- overlap: %u\n",j1,j2,overlap_j1j2);
      //Selection for the second link

      for (t =0 ; t<100;t++)
      {
        z = rand() % r_slap_layers[alpha][r_slap_total[alpha]];
        k1 = binary_search(r_slap_layers[alpha],r_slap_total[alpha],z);
        k2 = J_slap_layers[alpha][z];
        //printf("The link is: (%u, %u)\n",k1,k2);

        if (is_neigh(J_slap_layers[alpha], r_slap_layers[alpha], r_slap_total[alpha], k2, j1) ||
            is_neigh(J_slap_layers[alpha], r_slap_layers[alpha], r_slap_total[alpha], k2, j2) ||
                k2 == j2 ||
                k2 == j1)
          continue;
        else
        {
          existence_flag=0;
          //check the existence of the link j1 -k1 and j2-k2 -> If the proposed links don't belong to the layer -> Make the swap
          existence_flag=check_existence(j1,k1,j2,k2,r_slap_layers,J_slap_layers, alpha);
          if(existence_flag==0)
          {
            flag =1;
            delta_denom=evaluate_Delta_denom(k1,k2,j1,j2,r_slap_layers,J_slap_layers,r_slap_total, M,overlap_j1j2);

            if(delta_denom >= 0)
            {

            //printf("First: make the swap at layer %u between link (%u,%u) and (%u,%u) -> Variation on the denominator: %d\n",alpha,j1,j2,k1,k2,delta_denom);
            change_link(j1, J_slap_layers[alpha], r_slap_layers[alpha], j2, k1); // The link j1-j2 becames j1-k1 replacing j2 -> k1
            change_link(k1, J_slap_layers[alpha], r_slap_layers[alpha], k2, j1); // The link k1-k2 becames k1-j1 replacing k2 -> j1
            change_link(j2, J_slap_layers[alpha], r_slap_layers[alpha], j1, k2); // The link j2-j1 becames j2-k2 replacing j1 -> k2
            change_link(k2, J_slap_layers[alpha], r_slap_layers[alpha], k1, j2); // The link k2-k1 becames k2-j2 replacing k1 -> k2
            dd =0;
            break;
            }
            else
              delta_denom=0;
          }
          else
          {
          //check the existence of the link j1 -k2 and j2-k1 -> If the proposed links don't belong to the layer -> Make the swap
          existence_flag=0;
          existence_flag=check_existence(j1,k2,j2,k1,r_slap_layers,J_slap_layers, alpha);
            if(existence_flag == 0)
            {
              flag =1;
              delta_denom=evaluate_Delta_denom(k1,k2,j2,j1,r_slap_layers,J_slap_layers,r_slap_total, M,overlap_j1j2);
              if(delta_denom >= 0)
              {
              //printf("Second: make the swap at layer %u between link (%u,%u) and (%u,%u) -> Variation on the denominator %d\n --Incrociato\n",alpha,j1,j2,k1,k2,delta_denom);
              change_link(j1, J_slap_layers[alpha], r_slap_layers[alpha], j2, k2); // The link j1-j2 becames j1-k2 replacing j2 -> k2
              change_link(k1, J_slap_layers[alpha], r_slap_layers[alpha], k2, j2); // The link k1-k2 becames k1-j2 replacing k2 -> j2
              change_link(j2, J_slap_layers[alpha], r_slap_layers[alpha], j1, k1); // The link j2-j1 becames j2-k1 replacing j1 -> k1
              change_link(k2, J_slap_layers[alpha], r_slap_layers[alpha], k1, j1); // The link k2-k1 becames k2-j1 replacing k1 -> j1
              dd =0;
              break;
              }
              else
                delta_denom=0;
            }
          }
        }
      }
    }
    dd++;

    if(flag==1)
    {

      denominator_structuraloverlap+=delta_denom;
      structural_overlap=(1.0*(numerator_structuraloverlap))/(M*denominator_structuraloverlap);
      //Count the number of steps trying to decrease the overlap
      t_step++;
      // If after a certain number of iterations, it cannot be possible to increase the overlap,
      // then stop
      if(t_step > MAX_ITERATIONS)
      {
        fprintf(stderr,"structural overlap: %lf ",map_in_unitary(M,structural_overlap));
        fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);
        fprintf(stderr,"Exit now after %llu iterations\n",MAX_ITERATIONS);
        exit(1);
      }
      if (fabs(structural_overlap_prev-map_in_unitary(M,structural_overlap)) < TOLL)
      {
        t_step=0;
        fprintf(stderr,"structural overlap: %lf ",map_in_unitary(M,structural_overlap));
        fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);
        if (fabs(structural_overlap_prev - 0.90) < TOLL || 
          fabs(structural_overlap_prev - 0.80) < TOLL || 
          fabs(structural_overlap_prev - 0.70) < TOLL || 
          fabs(structural_overlap_prev - 0.60) < TOLL || 
          fabs(structural_overlap_prev - 0.50) < TOLL || 
          fabs(structural_overlap_prev - 0.40) < TOLL || 
          fabs(structural_overlap_prev - 0.30) < TOLL || 
          fabs(structural_overlap_prev - 0.20) < TOLL || 
          fabs(structural_overlap_prev - 0.10) < TOLL )
        {
          //printf("structural_overlap_prev:%lf ---  current overlap: %lf \n",structural_overlap_prev,map_in_unitary(M,structural_overlap));
          write_edges_onfiles_top_to_bottom(J_slap_layers, r_slap_layers,M, N,structural_overlap_prev);
        }

        structural_overlap_prev += delta_numerical;
        count++;

      }
      if (count==SIZE)
        finalflag=1;
    }
  }
  printf("structural overlap: %lf ",map_in_unitary(M,structural_overlap));
  printf("-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);
  if (map_in_unitary(M,structural_overlap) == 0.00)
     write_edges_onfiles(J_slap_layers, r_slap_layers,M, N,0.0);
}


void increase_overlap(unsigned long long int  **J_slap_layers, unsigned long long int **r_slap_layers, unsigned long long int *v_slap,
                unsigned long long int *r_slap_total, unsigned long long int N, unsigned int M, unsigned long long int *nodes_unique,
                struct layer *multiplex,unsigned long long int length_unique, unsigned int flag_decrease_increase,double TOLL,double max_overlap)
{
   unsigned long long int j1, j2, j2_pos;  /* The first couple of nodes */
  unsigned long long int k1, k2, k2_pos;  /* The second couple of nodes */
  unsigned long long int i,s,t,z, alpha, temp_swap,overlap_j1j2,existence_flag;
  unsigned long long int numerator_structuraloverlap,denominator_structuraloverlap,flag,SIZE=100,dd=0,nedges_intersection=0,Ndistinct_edges=0;

  int t_step =0;
  double omega,KC,structural_overlap, structural_overlap_prev,delta_numerical=0.01;
  long long int MAX_LINKS=0;
  long long int MAX_ITERATIONS=500000000;
  int delta_denom,count =1,flagga=0,finalflag=0;

  if (flag_decrease_increase == 0)
    delta_numerical=-delta_numerical;

  for (i = 0; i< M;i++)
  {
    if (MAX_LINKS < multiplex[i].nlinks)
      MAX_LINKS= multiplex[i].nlinks;
  }

  numerator_structuraloverlap=v_slap[M]/2;
  Ndistinct_edges=count_number_edge_aggregate(multiplex,M,N);
  denominator_structuraloverlap=Ndistinct_edges;

  structural_overlap=(1.0*(numerator_structuraloverlap))/(M*denominator_structuraloverlap);
  structural_overlap_prev = roundf(map_in_unitary(M,structural_overlap) * 100)/100 + delta_numerical;
  fprintf(stderr,"Number of links in the multiplex: %llu --- Number of links in the aggregate: %llu\n",numerator_structuraloverlap, denominator_structuraloverlap);
  fprintf(stderr,"structural overlap: %lf",map_in_unitary(M,structural_overlap));  
  fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);

  while (denominator_structuraloverlap != MAX_LINKS )
  {
    //This is the procedure to efficiently sample a link T from the multiplex
    s = rand() % v_slap[M]; //Select the first node of couple T// --> Random number between 0 and max number of links
    alpha = binary_search(v_slap, M, s); // Find the layer index such that v_slap[alpha] < s < v_slap[alpha+1]
    j2_pos= s-v_slap[alpha]; // This indicates the position of the second node of link T in the r-j-slap structure of the alpha layer
    // Find the first node of link T in the r_slap structure of the layer alpha, it will find the index such that it is verified
    // that r_slap_layers[j1] < j2_pos < r_slap_layers[j1+1]
    j1 = binary_search(r_slap_layers[alpha],r_slap_total[alpha],j2_pos);
    j2 = J_slap_layers[alpha][j2_pos];
    if(j2 < j1)
    {
      temp_swap=j1;
      j1=j2;
      j2=temp_swap;
    }
    //printf("The link: (%u,%u)\n",j1,j2);
    // This is for the selection of the first link -> Do a check for the overlap!//
    overlap_j1j2=check_overlap(j1,j2,alpha,r_slap_layers,J_slap_layers,r_slap_total,M);
    //printf("HERE--- Layer:%u |||| %u -- %u --->%u ||| RANDOM :%u --- j2_pos: %u\n",alpha,j1,j2,overlap_j1j2,s,j2_pos);

    flag =0;
    if(overlap_j1j2 == 1)
    {
      //printf("The link is: (%u, %u) --- overlap: %u\n",j1,j2,overlap_j1j2);
      //Selection for the second link


      for (t =0 ; t<100;t++)
      {
        z = rand() % r_slap_layers[alpha][r_slap_total[alpha]];
        k1 = binary_search(r_slap_layers[alpha],r_slap_total[alpha],z);
        k2 = J_slap_layers[alpha][z];
        //printf("The link is: (%u, %u)\n",k1,k2);
        if (is_neigh(J_slap_layers[alpha], r_slap_layers[alpha], r_slap_total[alpha], k2, j1) ||
            is_neigh(J_slap_layers[alpha], r_slap_layers[alpha], r_slap_total[alpha], k2, j2) ||
                k2 == j2 ||
                k2 == j1)
          continue;
        else
        {
          existence_flag=0;
          //check the existence of the link j1 -k1 and j2-k2 -> If the proposed links don't belong to the layer -> Make the swap
          existence_flag=check_existence(j1,k1,j2,k2,r_slap_layers,J_slap_layers, alpha);
          if(existence_flag==0)
          {
            flag =1;
            delta_denom=evaluate_Delta_denom(k1,k2,j1,j2,r_slap_layers,J_slap_layers,r_slap_total, M,overlap_j1j2);

            if(delta_denom <= 0)
            {
             //printf("First: make the swap at layer %u between link (%u,%u) and (%u,%u) -> Variation on the denominator: %d\n",alpha,j1,j2,k1,k2,delta_denom);
             change_link(j1, J_slap_layers[alpha], r_slap_layers[alpha], j2, k1); // The link j1-j2 becames j1-k1 replacing j2 -> k1
             change_link(k1, J_slap_layers[alpha], r_slap_layers[alpha], k2, j1); // The link k1-k2 becames k1-j1 replacing k2 -> j1
             change_link(j2, J_slap_layers[alpha], r_slap_layers[alpha], j1, k2); // The link j2-j1 becames j2-k2 replacing j1 -> k2
             change_link(k2, J_slap_layers[alpha], r_slap_layers[alpha], k1, j2); // The link k2-k1 becames k2-j2 replacing k1 -> k2
             dd =0;
              break;
            }
            else
              delta_denom=0;
          }
          else
          {
          //check the existence of the link j1 -k2 and j2-k1 -> If the proposed links don't belong to the layer -> Make the swap
          existence_flag=0;
          existence_flag=check_existence(j1,k2,j2,k1,r_slap_layers,J_slap_layers, alpha);
            if(existence_flag == 0)
            {
              flag =1;
              delta_denom=evaluate_Delta_denom(k1,k2,j2,j1,r_slap_layers,J_slap_layers,r_slap_total, M,overlap_j1j2);
              if(delta_denom <=  0)
              {
              //printf("Second: make the swap at layer %u between link (%u,%u) and (%u,%u) -> Variation on the denominator %d\n --Incrociato\n",alpha,j1,j2,k1,k2,delta_denom);
              change_link(j1, J_slap_layers[alpha], r_slap_layers[alpha], j2, k2); // The link j1-j2 becames j1-k2 replacing j2 -> k2
              change_link(k1, J_slap_layers[alpha], r_slap_layers[alpha], k2, j2); // The link k1-k2 becames k1-j2 replacing k2 -> j2
              change_link(j2, J_slap_layers[alpha], r_slap_layers[alpha], j1, k1); // The link j2-j1 becames j2-k1 replacing j1 -> k1
              change_link(k2, J_slap_layers[alpha], r_slap_layers[alpha], k1, j1); // The link k2-k1 becames k2-j1 replacing k1 -> j1
              dd =0;
              break;
              }
              else
                delta_denom=0;
            }
          }
        }
      }
    }
    dd++;

    if(flag==1)
    {

      denominator_structuraloverlap+=delta_denom;
      structural_overlap=(1.0*(numerator_structuraloverlap))/(M*denominator_structuraloverlap);
      //Count the number of steps trying to increase the overlap
      t_step++;
      // If after a certain number of iterations, it cannot be possible to increase the overlap,
      // then stop
      if(t_step > MAX_ITERATIONS)
      {
        fprintf(stderr,"structural overlap: %lf ",map_in_unitary(M,structural_overlap));
        fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);
        fprintf(stderr,"Exit now after %llu iterations\n",MAX_ITERATIONS);
        exit(1);
      }

      if (fabs(structural_overlap_prev-map_in_unitary(M,structural_overlap)) < TOLL)
      {
        t_step=0;
      fprintf(stderr,"structural overlap: %lf ",map_in_unitary(M,structural_overlap));
      fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);
      if (fabs(structural_overlap_prev - 0.90) < TOLL || 
        fabs(structural_overlap_prev - 0.80) < TOLL || 
        fabs(structural_overlap_prev - 0.70) < TOLL || 
        fabs(structural_overlap_prev - 0.60) < TOLL || 
        fabs(structural_overlap_prev - 0.50) < TOLL || 
        fabs(structural_overlap_prev - 0.40) < TOLL || 
        fabs(structural_overlap_prev - 0.30) < TOLL || 
        fabs(structural_overlap_prev - 0.20) < TOLL || 
        fabs(structural_overlap_prev - 0.10) < TOLL)
      {
        write_edges_onfiles_bottom_to_top(J_slap_layers, r_slap_layers,M, N,structural_overlap_prev);
      
      }

      // Exit if I reach the max_overlap value inserted in input
      if (fabs(structural_overlap_prev - max_overlap) < TOLL)
        exit(1);

      structural_overlap_prev += delta_numerical;
        count++;

      }
      if (count==SIZE)
      {
        finalflag=1;
        //break;
      }

    }
  }
  fprintf(stderr,"structural overlap: %lf ",map_in_unitary(M,structural_overlap));
  fprintf(stderr,"-- Number of links in the aggregate: %llu\n",denominator_structuraloverlap);
  if (map_in_unitary(M,structural_overlap) == 1.00)
     write_edges_onfiles(J_slap_layers, r_slap_layers,M, N,1.0);
}
