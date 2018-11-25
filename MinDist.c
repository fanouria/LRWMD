/*
 ============================================================================
 Name        : MinDist.c
 Author      : fanath
 Version     : 1_3
 Copyright   :
 Description : Computes for every word of the vocabulary the euclidean
              distance to the closest word for a number of histograms
 ============================================================================

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#define R_W2V  				6423
#define C_W2V  				300
#define VocabularySize		6423
#define MAX_THREADS         300


unsigned int  *threads_num, *hists_num, *hist_len;

unsigned int *h_entries=NULL,*bows=NULL, *hists=NULL;
/*
 * h_entries    number of words that each document contains
 * bows         which words of our vocabulary exist in a certain document
 * hists		how many times each word appears in each document
 */
unsigned int * min_dpT;  //  minimum number of *hists_num per Thread
unsigned int * add_docT; //  how many threads are going to get an additional doc
float w2v[R_W2V][C_W2V], *mindist, *rwmd;

typedef struct
{
	unsigned int		hist;
	unsigned int		dpt; //*hists_num per thread
	float			   *mindist;
}               Mindiststr;

void * worker(void *arg);
void  MinDist(unsigned int hist, unsigned int wcount);



int main (int argc, char *argv[])
{

	/***************************************************************/
				/*parameters from terminal*/
	/***************************************************************/
	unsigned int init_arg[argc-1];
	unsigned int i = 0, hist1 = 0;
	char* s;
	char* value;

	if (argc!=5 )
	{
	   fprintf(stdout,"Wrong number of parameters.\n");
       fprintf(stdout,"parameters must be given as: 'threads=' value,'hists_num=' value,'hist_len=' value,'hist=' value");
       exit(1);
   	}

	printf("Number of init_arg: %d\n", argc);
	for (i=0; i < argc; i++)
	{
		printf("argument %d: (%s)\n", i, argv[i]);
	}

	i=1;

	if((strncmp(argv[i],"threads",strlen("threads"))))
	{
		fprintf(stdout,"TypeError in argument 'threads'.\n");
		fprintf(stdout,"parameters must be given as: 'threads=' value,'hists_num=' value,'hist_len=' value,'hist=' value");
		return -1;
	}
	else
	{
		s = strchr(argv[i], '=');
		value = s + 1;
		init_arg[i-1]=atoi(value);

		if((init_arg[i-1] > MAX_THREADS) || (init_arg[i-1] <= 0))
		{
			fprintf(stdout,"Error! Number of threads should be from 1 to %d:\t",MAX_THREADS);
			return -1;
		}

	}

	i++;
	if((strncmp(argv[i],"hists_num",strlen("hists_num"))))
	{
		fprintf(stdout,"TypeError in argument 'hists_num'\n");
		fprintf(stdout,"parameters must be given as: 'threads=' value,'hists_num=' value,'hist_len=' value,'hist=' value");
		return -1;
	}
	else
	{
		s = strchr(argv[i], '=');
		value = s+1;
		init_arg[i-1]=atoi(value);

		if(init_arg[i-1]<=0)
		{
			fprintf(stdout,"Number of histograms should be a positive number.\n");
			return -1;
		}
	}

	i++;
	if((strncmp(argv[i],"hist_len",strlen("hist_len"))))
	{
		fprintf(stdout,"TypeError in argument 'hist_len'\n");
		fprintf(stdout,"parameters must be given as: 'threads=' value,'hists_num=' value,'hist_len=' value,'hist=' value");
		return -1;
	}
	else
	{
		s = strchr(argv[i], '=');
		value = s+1;
		init_arg [i-1] = atoi(value);

		if(init_arg[i-1]<=0)
		{
			fprintf(stdout,"Length of histogram should be a positive number.\n");
			return -1;
		}
	}
	i++;
	if((strncmp(argv[i],"hist",strlen("hist"))))
	{
		fprintf(stdout,"TypeError in argument 'hist'\n");
		fprintf(stdout,"parameters must be given as: 'threads=' value,'hists_num=' value,'hist_len=' value,'hist=' value");
		return -1;
	}
	else
	{
		s = strchr(argv[i], '=');
		value = s+1;
		init_arg [i-1] = atoi(value);

		if(init_arg[i-1]<0)
		{
			fprintf(stdout,"Histogram should not be negative.\n");
			return -1;
		}
	}
	threads_num =malloc(sizeof(unsigned int));
	hists_num =malloc(sizeof(unsigned int));
	hist_len =malloc(sizeof(unsigned int));

	*threads_num = init_arg [0];
	*hists_num   = init_arg [1];
	*hist_len    = init_arg [2];
	 hist1		 = init_arg [3];

	fprintf(stdout,"threads_num= %u \t hists_num=%u \t hist_len=%u \t hist=%u \n", *threads_num, *hists_num, *hist_len,hist1);


	h_entries = malloc(*hists_num * 1 * sizeof(unsigned int));
	bows      = malloc(*hists_num * (*hist_len ) * sizeof(unsigned int));
	hists     = malloc(*hists_num * (*hist_len )  * sizeof(unsigned int));

	FILE *h_entries_file = NULL, *w2v_file = NULL, *hists_file = NULL, *bows_file  =NULL ;

	/***************************************************************/
			/*read data from the files*/
	/***************************************************************/

	/*open files*/
	h_entries_file = fopen("h_entries.bin","rb");
	if (h_entries_file == NULL) exit(1);

	w2v_file = fopen("w2v.bin","rb");
	if (w2v_file == NULL) exit(1);

	hists_file = fopen("hists.bin","rb");
	if (hists_file == NULL) exit(1);

	bows_file = fopen("bows.bin","rb");
	if (bows_file==NULL) exit(1);

	/*read data*/
	fread( h_entries, sizeof(int), *hists_num, h_entries_file );
	fread( w2v, sizeof(float), R_W2V  * C_W2V, w2v_file );
	fread( hists, sizeof(int), *hists_num * (*hist_len ), hists_file );
	fread( bows, sizeof(int), *hists_num * (*hist_len ), bows_file);


	/*close files*/
	fclose( h_entries_file );
	fclose( w2v_file );
	fclose( hists_file );
	fclose( bows_file );

	/***************************************************************/
				/*start threads and compute mindist*/
	/***************************************************************/
	min_dpT    =  malloc ( sizeof ( unsigned int ));
	add_docT   =  malloc ( sizeof ( unsigned int ));
	* min_dpT  =  *hists_num / *threads_num;
	* add_docT =  *hists_num % *threads_num;

/**************************** PHASE 1 ****************************/
	/*allocate memory for mindist,gdcounter,rwmd*/

	mindist    = malloc ( VocabularySize * (*hists_num) * sizeof( float ));
	rwmd       = malloc (( * hists_num )* sizeof( float ));


	pthread_t      *threads;
	pthread_attr_t  pthread_custom_attr;
	Mindiststr           *arg;

	threads = ( pthread_t *) malloc(( *threads_num ) *  sizeof( pthread_t ));
	pthread_attr_init( &pthread_custom_attr );
	arg = (Mindiststr *)malloc( sizeof( Mindiststr ) * ( *threads_num ));


/*  Time Computation in Linux 1/2
 *	struct timespec start, finish;
 *	double time_taken_RWMD;
 *	clock_gettime(CLOCK_REALTIME, &start);
 */
	clock_t el_time;
	el_time = clock();



	/* Spawn threads */
	 for (i = 0; i !=  (*threads_num); i++){

		if (i < * add_docT){
			arg[i].hist = i*(* min_dpT + 1);
			arg[i].dpt = * min_dpT + 1;
		}
		else{
			arg[i].hist = i * (* min_dpT) + * add_docT;
			arg[i].dpt = * min_dpT;
		}


		arg[i].mindist = mindist;
		pthread_create(&threads[i], &pthread_custom_attr, worker, (void *)(arg+i));

	 }


	/*wait for all threads to finish*/
	for (i = 0; i !=  (*threads_num); i++)
	{
		pthread_join(threads[i], NULL);

	}

/*  Time Computation in Linux 2/2
 * 	clock_gettime(CLOCK_REALTIME, &finish);
 * 	time_taken_RWMD = (finish.tv_sec - start.tv_sec);
 * 	time_taken_RWMD += (finish.tv_nsec - start.tv_nsec)/ 1000000000.0;
 */

	el_time = clock() - el_time;
	double time_taken_RWMD = ((double)el_time)/CLOCKS_PER_SEC; // in seconds
	printf("\nIt took %f seconds to compute RWMDs. \n", time_taken_RWMD);


	/***************************************************************/
			/*Print RWMDs into a txt file */
	/***************************************************************/

	FILE * fp;
	fp = fopen ("C:\\Users\\fanouriaath\\eclipse-workspace\\Rwmd\\mindist.txt","w");
	for (int hist = 0; hist != (*hists_num); hist++){
		for (int i = 0; i != VocabularySize ; i++)
		{

				  fprintf(fp, "%f ", *(mindist + VocabularySize*hist + i));

		}
		 fprintf(fp, "\n");

	}


	fclose(fp);

	free(mindist);
	return EXIT_SUCCESS;
}

void * worker(void *arg)
{
	Mindiststr   *p = (Mindiststr *) arg;
	int hist,i,imax;

	hist = p->hist;
	imax = p->dpt;

	for (i = 0; i != imax ; i++)
	{
		 MinDist( hist, *(h_entries + hist));
		 hist++;
	}

	return NULL;
}



void MinDist(unsigned int hist, unsigned int wcount){

	int w1j;
	float tmp, x;

	for (int i = 0; i != VocabularySize ; i++)
	{
		for (int j = 0; j != wcount ; j++)
		{
			w1j = *(bows +hist * (*hist_len ) + j);

			tmp = 0;

			for(int k = 0; k != C_W2V; k++)
			{
				x= w2v[i][k] - w2v[w1j][k];
				tmp = tmp + x * x;
			}

			if(j==0)
				*(mindist + VocabularySize*hist + i) = tmp;
			else
			{
				if (*(mindist + VocabularySize*hist + i) > tmp)
					*(mindist + VocabularySize*hist + i) = tmp;//min per row [x y z]
			}

		}
	}


}











