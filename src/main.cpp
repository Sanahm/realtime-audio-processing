#include "library.h"

int main(int argc, char * argv[]){
	//FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	
	
	
	const char * filename = "impres"; int i = 0;
	double * data; unsigned  int size = 0,len; double * out;
	
	getImpulseResponse(filename, &data, &size );
	len = get_nextpow2(2*size-1);
	out = dgetmem(len);
	double * _data = fconv(data, size, data, size);
	//memcpy(_data, data, size);
	//fft(_data, out, len);
	printf("%d \n %f", len, _data[2*size-1]);
	
	/*fprintf(gnuplotPipe, "plot '-' title 'original' with lines \n");
	for(i = 0; i < 2*size+5; i++) {
		fprintf(gnuplotPipe, "%f\n", _data[i]);
	}*/
	
}
	