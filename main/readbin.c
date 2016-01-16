#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[]){
	int counter;
	FILE *ptr_myfile;
	int seed;
	float *point;
		point = malloc(sizeof(float));
	
	ptr_myfile=fopen(argv[1],"rb");

	fread(&seed, sizeof(unsigned int), 1, ptr_myfile);
	printf("%d\n", seed);
	for (counter=0; counter<atoi(argv[2]); counter++){
		fread(point, sizeof(point[0]), 3, ptr_myfile);
		printf("%lf\t%lf\t%lf\n", point[0], point[1], point[2]);
	}

	fclose(ptr_myfile);
	return 0;
	
}
