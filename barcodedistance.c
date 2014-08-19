#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int hamming_distance(char *s, char *t);
int levenshtein_distance(char *s, char *t);
int minimum(int a, int b, int c);

int main(int argc, char* argv[]) {
	FILE* bin_file; // file containing sequence and hits (from binning in perl script)
	FILE* lib_file; // file containing libsequence and sr_id
	FILE* out_file; // Output file
	
	int distance_threshold = 10000; // setting absurd high number by default, this can be finetuned through user input
	
	if (argc < 4 || argc > 5) {
		printf( "Error: Incorrect number of arguments.\nUsage: %s inputfile libraryfile outputfile [distance threshold]\n", argv[0] );
		return 1;
	}
	if (argc == 5)
		distance_threshold = atoi(argv[4]);
	
	/* Open files */
	bin_file = fopen(argv[1], "r");
	if (bin_file==NULL) {
		printf("Error: can't open binfile.\n");
		return 1;
	}
	lib_file = fopen(argv[2], "r");
	if (lib_file==NULL) {
		printf("Error: can't open libfile.\n");
		return 1;
	}
	out_file = fopen(argv[3], "w"); 
	if (out_file==NULL) {
		printf("Error: can't create outputfile.\n");
		return 1;
	}
	
	char inbuf[500]; // maybe a bit large
	char *item;
	
	/* Read bin file */
	printf("Reading binned file %s\n", argv[1]);
	int bin_file_count=0;
	int bin_lines = 0;
	while (fgets(inbuf, sizeof(inbuf), bin_file) != NULL)
		bin_lines++;
	rewind(bin_file);
	printf("Bin file length: %i\n", bin_lines);
	char **bin_file_seq_array = (char **)malloc(bin_lines * sizeof(char *));
	int *bin_file_hits_array = (int *)malloc(bin_lines * sizeof(int));
	char **barcode_array = (char **)malloc(bin_lines * sizeof(char *));
	while (fgets(inbuf, sizeof(inbuf), bin_file) != NULL) { // Assuming binned file has three columns (barcode, sequence, hits) delimited by \t and newlined by \n
		item = strtok(inbuf, "\t");
		barcode_array[bin_file_count] = strdup(item);
		item = strtok(NULL, "\t");
		bin_file_seq_array[bin_file_count] = strdup(item);
		item = strtok(NULL, "\n");
		bin_file_hits_array[bin_file_count] = atoi(item);
		bin_file_count++;
	}
	fclose(bin_file);
	
	/* Read lib file */
	printf("Reading library file %s\n", argv[2]);
	int lib_file_count=0;
	long lib_lines = 0;
	while (fgets(inbuf, sizeof(inbuf), lib_file) != NULL)
		lib_lines++;
	rewind(lib_file);
	printf("Lib file length: %li\n", lib_lines);
	char **lib_file_seq_array = (char **)malloc(lib_lines * sizeof(char *));
	char **lib_file_srid_array = (char **)malloc(lib_lines * sizeof(char *));
	while (fgets(inbuf, sizeof(inbuf), lib_file) != NULL) { // Assuming library file has two columns (sr_id, sequence) delimited by \t and newlined by \n
		item = strtok(inbuf, "\t");
		lib_file_srid_array[lib_file_count] = strdup(item);
//printf("Gelezen ID:%s\n",item);
		item = strtok(NULL, "\n");
		lib_file_seq_array[lib_file_count] = strdup(item);
//printf("Gelezen shRNA:%s\n",item);
		lib_file_count++;
	}
	fclose(lib_file);

	int i, j, k, distance;
	char tmpbuf[500];
	printf("Calculating distances ... ");
	fflush(stdout);
	for (i=0;i<lib_file_count;i++) { // Assume lib_file is smallest
		if((i % (lib_file_count / 100))==0) {
			printf("Analyzed %f percent of the library\n", 100*((double)i/(double)lib_file_count));
		}
		int out_file_count = 0;
		char **out_file_array = (char **)malloc(bin_file_count * sizeof(char *));

char *shrnasequence=lib_file_seq_array[i];			
//char shrnasequence[18];
//memcpy( shrnasequence, &pipo[0], 18 );
//shrnasequence[18] = '\0';

		for (j=0;j<bin_file_count;j++) {
			distance = hamming_distance(bin_file_seq_array[j], shrnasequence);
			//printf("Distance comparison of %s with %s resulted in a distance of %i\n", bin_file_seq_array[j],lib_file_seq_array[i],distance);					
			if (distance<=distance_threshold) {
				sprintf(tmpbuf, "%s\t%s\t%i\t%s\t%s\t%i", barcode_array[j], bin_file_seq_array[j], bin_file_hits_array[j], shrnasequence, lib_file_srid_array[i], distance);
				out_file_array[out_file_count] = strdup(tmpbuf);
				out_file_count++;
			}
		}
		for (k=0;k<out_file_count;k++) {
			fprintf(out_file, "%s\n", out_file_array[k]);
			free(out_file_array[k]);
		}
		free(out_file_array);
	}
	fclose(out_file);
	
	/* Free memory */
	free(bin_file_seq_array);
	free(bin_file_hits_array);
	free(lib_file_seq_array);
	free(lib_file_srid_array);
	
	return 0;
}

/* Computes Hamming */
int hamming_distance(char *s, char *t) {
	int i;
	int length = strlen(s);
	int distance = 0;
	for (i=0;i<length;i++)
		if (s[i] != t[i])
			distance++;
	return distance;
}

/* Computes Damerau-Levenshtein */
int levenshtein_distance(char *s,char *t) {
	//Step 1
	int k,i,j,n,m,cost,*d,distance;
	n=strlen(s); 
	m=strlen(t);
	if(n!=0&&m!=0) {
		d=(int*)malloc((sizeof(int))*(m+1)*(n+1));
		m++;
		n++;
		//Step 2	
		for(k=0;k<n;k++)
			d[k]=k;
		for(k=0;k<m;k++)
			d[k*n]=k;
		//Step 3 and 4	
		for(i=1;i<n;i++) {
			for(j=1;j<m;j++) {
				//Step 5
				if(s[i-1]==t[j-1])
					cost=0;
				else
					cost=1;
				//Step 6			 
				d[j*n+i]=minimum(d[(j-1)*n+i]+1,d[j*n+i-1]+1,d[(j-1)*n+i-1]+cost);
			}
		}
		distance=d[n*m-1];
		free(d);
		return distance;
	}
	else 
		return -1; // a negative return value means that one or both strings are empty.
}

/* Gets the minimum of three values */
int minimum(int a,int b,int c) {
	int min=a;
	if(b<min)
		min=b;
	if(c<min)
		min=c;
	return min;
}
