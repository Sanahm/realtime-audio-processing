#include "library.h"

int main(int argc, char * argv[]){
	//FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	FILE *fp = fopen("report.txt","r+");
	if(fp == NULL) {
		printf("can't open file \n");
		exit(0);
	}
	else {
		for(int i = 0; i <9; i++) {
			fseek(fp,0, SEEK_END);
			fprintf(fp,"%lf;",5.78);
		}
		printf("jhh");
		
	}
	fclose(fp);
	
	/*const char * filename = "impres"; int i = 0;
	double * data; unsigned  int size = 0,len; double * out,*yout;
	
	getImpulseResponse(filename, &data, &size );
	len = get_nextpow2(size);
	out = dgetmem(len);
	yout = dgetmem(len);
	memcpy(out, data, size);
	double * _data = fconv(data, size, data, size);
	//memcpy(_data, data, size);
	fftr(out, yout, len);
	printf("%d \n %f", len, _data[2*size-1]);
	
	/*fprintf(gnuplotPipe, "plot '-' title 'original' with lines \n");
	for(i = 0; i < len; i++) {
		fprintf(gnuplotPipe, "%f\n", out[i]);
	}
	*/
}
	