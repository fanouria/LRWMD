/*
 ============================================================================
 Name        : LRWMD.c
 Author      : fanath
 Version     : 1_3
 Copyright   :
 Description : Computation of RWMD for a given histogram in comparison with
  	  	  	  the rest histograms of the dataset using a precomputed mindist
  	  	  	  matrix computed with MinDist.c
 ============================================================================

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <pthread.h>
#include <math.h>


#define R_W2V  		    6423
#define C_W2V  		    300
#define VocabularySize	    6423
#define MAX_THREADS         300


unsigned int  *threads_num, *hists_num, *hist_len;

unsigned int *h_entries = NULL,*bows = NULL, *hists = NULL, *gdcounter;
/*
 * h_entries    number of words that each document contains
 * bows         which words of our vocabulary exist in a certain document
 * hists		how many times each word appears in each document
 */
unsigned int * min_dpT;  //  minimum number of *hists_num per Thread
unsigned int * add_docT; //  how many threads are going to get an additional doc
float w2v[R_W2V][C_W2V],*rwmd, * mindist = NULL;

typedef struct
{
	unsigned int		hist1;
	unsigned int		hist2;
	unsigned int		dpt; //*hists_num per thread
}               GDstr;

void * worker(void *arg);
void   RWMD(unsigned int hist1, unsigned int hist2, unsigned int wcount1, unsigned int wcount2);



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
	threads_num = malloc(sizeof(unsigned int));
	hists_num = malloc(sizeof(unsigned int));
	hist_len = malloc(sizeof(unsigned int));

	*threads_num = init_arg [0];
	*hists_num   = init_arg [1];
	*hist_len    = init_arg [2];
	 hist1       = init_arg [3];

	fprintf(stdout,"threads_num= %u \t hists_num=%u \t hist_len=%u \t hist=%u \n", *threads_num, *hists_num, *hist_len,hist1);


	h_entries = malloc(*hists_num * 1 * sizeof(unsigned int));
	bows      = malloc(*hists_num * (*hist_len )   * sizeof(unsigned int));
	hists     = malloc(*hists_num * (*hist_len )   * sizeof(unsigned int));
	mindist   = malloc(*hists_num * VocabularySize * sizeof(float));

	FILE *h_entries_file = NULL, *w2v_file = NULL, *hists_file = NULL, *bows_file = NULL, *mindist_file = NULL ;

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

	mindist_file = fopen("mindist.bin","rb");
	if (mindist_file==NULL) exit(1);

	/*read data*/
	fread( h_entries, sizeof(int), *hists_num, h_entries_file );
	fread( w2v, sizeof(float), R_W2V  * C_W2V, w2v_file );
	fread( hists, sizeof(int), *hists_num * (*hist_len ), hists_file );
	fread( bows, sizeof(int), *hists_num * (*hist_len ), bows_file);
	fread( mindist, sizeof(float), *hists_num * VocabularySize, mindist_file);


	/*close files*/
	fclose( h_entries_file );
	fclose( w2v_file );
	fclose( hists_file );
	fclose( bows_file );

	/***************************************************************/
	      /*start threads and compute rwmd*/
	/***************************************************************/
	min_dpT    =  malloc ( sizeof ( unsigned int ));
	add_docT   =  malloc ( sizeof ( unsigned int ));
	* min_dpT  =  *hists_num / *threads_num;
	* add_docT =  *hists_num % *threads_num;

	/* allocate memory for gd,gdcounter,rwmd */

	gdcounter  = malloc (( * hists_num + 1) * sizeof ( unsigned int ));
	rwmd       = malloc (( * hists_num )* sizeof( float ));

	unsigned int hist2, h_entries1, h_entries2;

	pthread_t      *threads;
	pthread_attr_t  pthread_custom_attr;
	GDstr           *arg;

	threads = ( pthread_t *) malloc(( *threads_num ) * sizeof( pthread_t ));
	pthread_attr_init( &pthread_custom_attr );
	arg = (GDstr *)malloc( sizeof( GDstr ) * ( *threads_num ));

	h_entries1 = *( h_entries + hist1 );

/*  Time Computation in Linux 1/2
 *	struct timespec start, finish;
 *	double time_taken_RWMD;
 *	clock_gettime(CLOCK_REALTIME, &start);
 */
	clock_t el_time;
	el_time = clock();

	for (hist2 = 0; hist2 != *hists_num; hist2++)
	{
		h_entries2 = *(h_entries+hist2);
		gdcounter[hist2+1] = h_entries1 * h_entries2 + *( gdcounter + hist2);
	}

	/* Spawn threads */
	 for (i = 0; i !=  (*threads_num); i++){

		arg[i].hist1 = hist1;

		if (i < * add_docT){
			arg[i].hist2 = i*(* min_dpT + 1);
			arg[i].dpt = * min_dpT + 1;
		}
		else{
			arg[i].hist2 = i * (* min_dpT) + * add_docT;
			arg[i].dpt = * min_dpT;
		}

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
	printf("\nIt took %f seconds to compute RWMDs for the given histogram. \n", time_taken_RWMD);


	/***************************************************************/
			/*Print RWMDs into a txt file */
	/***************************************************************/
	FILE * fp;
	fp = fopen ("C:\\Users\\fanouriaath\\eclipse-workspace\\Rwmd\\rwmd_threads.txt","w");
	for (hist2 = 0; hist2 != *hists_num; hist2++){
				fprintf(fp, "%f\n",*(rwmd+hist2));
	}
	fclose(fp);

	return EXIT_SUCCESS;
}

void * worker(void *arg)
{
	GDstr           *p = (GDstr *) arg;
	int hist1,hist2,i,imax;

	hist1 = p->hist1;
	hist2 = p->hist2;
	imax = p->dpt;

	for (i = 0; i !=  imax ; i++)
	{
		 RWMD(hist1, hist2, *(h_entries + hist1), *(h_entries + hist2));
		 hist2++;
	}

	return NULL;
}



void RWMD(unsigned int hist1, unsigned int hist2, unsigned int wcount1, unsigned int wcount2){

	int i, word;
	float d1 = 0, d2 = 0;

	//array multiplication

	for(i = 0; i != wcount2; i++)
	{
		word = *(bows +hist2 * (*hist_len ) + i);
		d1 += *(hists + ( (*hist_len ) * hist2) + i) * *(mindist + VocabularySize * hist1 + word);
	}

	for(i = 0; i != wcount1; i++)
	{
		word = *(bows +hist1 * (*hist_len ) + i);
		d2 += *(hists + ( (*hist_len ) * hist1) + i) * *(mindist + VocabularySize * hist2 + word);
	}

	if (d1 >= d2)
		*(rwmd + hist2) = d1;
	else{
		*(rwmd + hist2) = d2;
	}

}


