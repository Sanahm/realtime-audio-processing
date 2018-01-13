/******************************************/
/*
  reverb.cpp
  by SANA M. and BERTRAND E. 2017-2018.

  This program opens a stream and add some reverberation to
  input before passing to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8


typedef signed short MY_TYPE;
#define FORMAT RTAUDIO_SINT16


typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32
*/
typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

#ifndef PI
#define PI 3.14159265359
#endif



#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


int getImpulseResponse(const char *filename, double ** data, unsigned int *size );//charge the impulse response in "filename"

double * fconv(double * x, unsigned int L, double * h, unsigned int M);//compute the convolution of x by h the return is of size L+M-1

int fft(double *x, double *y, const int m); //compute the fft of x and place the real part of the result in x and imaginary in y (initial x is overwrited)

int ifft(double *x, double *y, const int m);

int fftr(double *x, double *y, const int m);

int ifftr(double *x, double *y, const int l);

static int checkm(const int m);

int get_nextpow2(int n); //find the next power of 2 such that 2**m >= n and return 2**m

char *getmem(int leng, unsigned size);

double *dgetmem(int leng); //allocate memory

double get_process_time();


typedef struct {
	double * iR;
	unsigned int iRSize;
	double * tmp;
	unsigned int bufferBytes;
	FILE * report;
} MyData;

typedef struct {
	double threshold;
	double releaseFactor;
	double smoothingFactor;
	double * lookahead;
	double env = 0;
	double gain = 0;
	double sgain = 0;
	double env_tmp = 0;
	double sgain_tmp = 0;
} MyLimiterData;



double *_sintbl = 0;
int maxfftsize = 0;

int getImpulseResponse(const char * filename, double ** data, unsigned int *size ) {
	FILE * stream =  fopen(filename,"rb");
	if(stream == NULL) {
		printf("error while reading the file named %s \n", filename);
		return -1;
	}
	fseek(stream,0, SEEK_END);
	*size = ftell(stream)/sizeof(**data)/2;
	rewind(stream); // or you can use fseek(stream,0, SEEK_SET);to set back the curseur to the beginning of the file;

	*data = (double *) calloc(*size, sizeof(**data));
	return fread(*data, sizeof(**data), *size, stream);
}

double * fconv(double * x, unsigned int L, double * h, unsigned int M) {
	int i, N = get_nextpow2(L + M - 1);
	double * _rt = dgetmem(N), * _it = dgetmem(N), *y = dgetmem(L+M-1);
	double * _rx = dgetmem(N), * _ix = dgetmem(N);
	double * _rh = dgetmem(N), * _ih = dgetmem(N);
	
	memcpy(_rx, x, L*sizeof(double));
	memcpy(_rh, h, M*sizeof(double));
	fftr(_rx, _ix, N);
	fftr(_rh, _ih, N);
	for( i = 0; i < N; i++ ){
		_rt[i] = _rx[i]*_rh[i] - _ix[i]*_ih[i];
		_it[i] = _rx[i]*_ih[i] + _ix[i]*_rh[i];
	}
	ifftr(_rt, _it, N);
	memcpy(y, _rt, (L+M-1)*sizeof(double));
	free(_rt); free(_it); free(_rx); free(_ix); free(_rh); free(_ih);
	return y;	
}

///////////////////////////////
// FFT functions
int get_nextpow2(int n)
{
  int k = 1;
  while (k < n){
    k *= 2;
  }

  return k;
}

int fftr(double *x, double *y, const int m)
{
   int i, j;
   double *xp, *yp, *xq;
   double *yq;
   int mv2, n, tblsize;
   double xt, yt, *sinp, *cosp;
   double arg;

   mv2 = m / 2;

   /* separate even and odd  */
   xq = xp = x;
   yp = y;
   for (i = mv2; --i >= 0;) {
      *xp++ = *xq++;
      *yp++ = *xq++;
   }

   if (fft(x, y, mv2) == -1)    /* m / 2 point fft */
      return (-1);


   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }
   //printf("Debug: m=%i, maxfftsize=%i\n",m,maxfftsize);

   n = maxfftsize / m;
   sinp = _sintbl;
   cosp = _sintbl + maxfftsize / 4;

   xp = x;
   yp = y;
   xq = xp + m;
   yq = yp + m;
   *(xp + mv2) = *xp - *yp;
   *xp = *xp + *yp;
   *(yp + mv2) = *yp = 0;

   for (i = mv2, j = mv2 - 2; --i; j -= 2) {
      ++xp;
      ++yp;
      sinp += n;
      cosp += n;
      yt = *yp + *(yp + j);
      xt = *xp - *(xp + j);
      *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
      *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
   }

   xp = x + 1;
   yp = y + 1;
   xq = x + m;
   yq = y + m;

   for (i = mv2; --i;) {
      *xp++ = *(--xq);
      *yp++ = -(*(--yq));
   }

   return (0);
}

int ifftr(double *x, double *y, const int l)
{
   int i;
   double *xp, *yp;

   fftr(x, y, l);

   xp = x;
   yp = y;
   i = l;
   while (i--) {
      *xp++ /= l;
      *yp++ /= -l;
   }

   return (0);
}


static int checkm(const int m)
{
   int k;

   for (k = 4; k <= m; k <<= 1) {
      if (k == m)
         return (0);
   }
   fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n",m);

   return (-1);
}

int fft(double *x, double *y, const int m)
{
   int j, lmx, li;
   double *xp, *yp;
   double *sinp, *cosp;
   int lf, lix, tblsize;
   int mv2, mm1;
   double t1, t2;
   double arg;
   int checkm(const int);

   /**************
   * RADIX-2 FFT *
   **************/

   if (checkm(m))
      return (-1);

   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }

   lf = maxfftsize / m;
   lmx = m;

   for (;;) {
      lix = lmx;
      lmx /= 2;
      if (lmx <= 1)
         break;
      sinp = _sintbl;
      cosp = _sintbl + maxfftsize / 4;
      for (j = 0; j < lmx; j++) {
         xp = &x[j];
         yp = &y[j];
         for (li = lix; li <= m; li += lix) {
            t1 = *(xp) - *(xp + lmx);
            t2 = *(yp) - *(yp + lmx);
            *(xp) += *(xp + lmx);
            *(yp) += *(yp + lmx);
            *(xp + lmx) = *cosp * t1 + *sinp * t2;
            *(yp + lmx) = *cosp * t2 - *sinp * t1;
            xp += lix;
            yp += lix;
         }
         sinp += lf;
         cosp += lf;
      }
      lf += lf;
   }

   xp = x;
   yp = y;
   for (li = m / 2; li--; xp += 2, yp += 2) {
      t1 = *(xp) - *(xp + 1);
      t2 = *(yp) - *(yp + 1);
      *(xp) += *(xp + 1);
      *(yp) += *(yp + 1);
      *(xp + 1) = t1;
      *(yp + 1) = t2;
   }

   /***************
   * bit reversal *
   ***************/
   j = 0;
   xp = x;
   yp = y;
   mv2 = m / 2;
   mm1 = m - 1;
   for (lmx = 0; lmx < mm1; lmx++) {
      if ((li = lmx - j) < 0) {
         t1 = *(xp);
         t2 = *(yp);
         *(xp) = *(xp + li);
         *(yp) = *(yp + li);
         *(xp + li) = t1;
         *(yp + li) = t2;
      }
      li = mv2;
      while (li <= j) {
         j -= li;
         li /= 2;
      }
      j += li;
      xp = x + j;
      yp = y + j;
   }

   return (0);
}

int ifft(double *x, double *y, const int m)
{
   int i;

   if (fft(y, x, m) == -1)
      return (-1);

   for (i = m; --i >= 0; ++x, ++y) {
      *x /= m;
      *y /= m;
   }

   return (0);
}

double *dgetmem(int leng)
{
    return ( (double *)getmem(leng, sizeof(double)) );
}

char *getmem(int leng, unsigned size)
{
    char *p = NULL;

    if ((p = (char *)calloc(leng, size)) == NULL){
        fprintf(stderr, "Memory allocation error !\n");
        exit(3);
    }
    return (p);
}

double get_process_time() {
    struct rusage usage;
    if( 0 == getrusage(RUSAGE_SELF, &usage) ) {
        return (double)(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) +
               (double)(usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1.0e6;
    }
    return 0;
}

void usage( void ) {
	// Error function in case of incorrect command-line
	// argument specifications
	std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
	std::cout << "    where N = number of channels,\n";
	std::cout << "    fs = the sample rate,\n";
	std::cout << "    iDevice = optional input device to use (default = 0),\n";
	std::cout << "    oDevice = optional output device to use (default = 0),\n";
	std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
	std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
	exit( 0 );
}

/*Reverberation without fft method*/
int reverb1( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double /*streamTime*/, RtAudioStreamStatus status, void *myData )
{
	double dtime = get_process_time();
	
	MyData * pData = (MyData *) myData;
	double * h = (double *) pData->iR; //impulsionnel response
	double * x = (double *) inputBuffer; // signal x
	unsigned int M = pData->iRSize; //size of h
	unsigned int L = nBufferFrames; // size of x
	int i,k;


	double conv;
	int kmin,kmax;

	// Since the number of input and output channels is equal, we can do
	// a simple buffer copy operation here.
	
	if ( status ) std::cout << "Stream over/underflow detected." << std::endl;
	
	// Calcul of the convolution with fft
	for( i = 0; i < L+M-1; i++ ){
		conv=0;
		
		if(i>=M) kmin=i-M+1;
		else kmin=0;
		if(i<L) kmax=i;
		else kmax=L-1;
		
		for( k=kmin;k<=kmax;k++){
			conv += x[k]*h[i-k+1];
		}
		pData->tmp[i] += conv;
	}
	memcpy( outputBuffer, pData->tmp , (unsigned int) pData->bufferBytes);
	
	memmove(pData->tmp, pData->tmp + L, (M-1)*sizeof(double));
	
	bzero(pData->tmp + (M-1), L*sizeof(double));
	
	dtime = get_process_time() - dtime;
	
	// add to report
	if(pData->report != NULL){
		fseek(pData->report,0, SEEK_END);
		fprintf(pData->report,"%lf;",dtime);
	}
	printf("%lf\n", dtime);
	return 0;
}

/*Reverberation with fft method*/
int reverb2( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double /*streamTime*/, RtAudioStreamStatus status, void *myData )
{
	double dtime = get_process_time();
	MyData * pData = (MyData *) myData;
	double * h = (double *) pData->iR; //impulsionnel response
	double * x = (double *) inputBuffer; // signal x
	unsigned int M = pData->iRSize; //size of h
	unsigned int L = nBufferFrames; // size of x
	int i;
	
	double * conv = fconv(x, L, h, M);

	// Since the number of input and output channels is equal, we can do
	// a simple buffer copy operation here.
	
	if ( status ) std::cout << "Stream over/underflow detected." << std::endl;
	
	// Calcul of the convolution with fft
	for( i = 0; i < L+M-1; i++ ){
		pData->tmp[i] += conv[i];
	}
	
	memcpy( outputBuffer, pData->tmp , (unsigned int) pData->bufferBytes);
	
	memmove(pData->tmp, pData->tmp + L, (M-1)*sizeof(double));
	
	bzero(pData->tmp + (M-1), L*sizeof(double));
	
	dtime = get_process_time() - dtime;
	printf("%lf\n", dtime);	
	return 0;
}

/*Limiter implementation*/
int limiter( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double /*streamTime*/, RtAudioStreamStatus status, void *myLimiterData )
{
	MyLimiterData * pData = (MyLimiterData *) myLimiterData;
	double * x = (double *) inputBuffer; // signal x
	int i;
	

	for(i = 0; i < nBufferFrames; i++) {
		pData->env = MAX(fabs(x[i]), pData->env_tmp*pData->releaseFactor);
		pData->env_tmp = pData->env;
		
		if( pData->env < pData->threshold )
			pData->gain = 1;
		else
			pData->gain = 1 + pData->threshold - pData->env;

		pData->sgain = pData->sgain_tmp*pData->smoothingFactor + pData->gain*(1 - pData->smoothingFactor);
		pData->sgain_tmp = pData->sgain;
		
		//output
		x[i] = x[i]*pData->sgain;
	}

	memcpy( outputBuffer, x , nBufferFrames*sizeof(double));
	
	return 0;
}




int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  unsigned int *bytes = (unsigned int *) data;
  memcpy( outputBuffer, inputBuffer, *bytes );
  return 0;
}





int main( int argc, char *argv[] )
{
	MyData myData; MyLimiterData myLimiterData;
	unsigned int channels, fs, bufferBytes, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

	// Minimal command-line checking
	if (argc < 3 || argc > 7 ) usage();

	RtAudio adac;
	if ( adac.getDeviceCount() < 1 ) {
		std::cout << "\nNo audio devices found!\n";
		exit( 1 );
	}

	channels = (unsigned int) atoi(argv[1]);
	fs = (unsigned int) atoi(argv[2]);
	if ( argc > 3 )
		iDevice = (unsigned int) atoi(argv[3]);
	if ( argc > 4 )
		oDevice = (unsigned int) atoi(argv[4]);
	if ( argc > 5 )
		iOffset = (unsigned int) atoi(argv[5]);
	if ( argc > 6 )
		oOffset = (unsigned int) atoi(argv[6]);

	// Let RtAudio print messages to stderr.
	adac.showWarnings( true );

	// Set the same number of channels for both input and output.
	unsigned int bufferFrames = 512;
	RtAudio::StreamParameters iParams, oParams;
	iParams.deviceId = iDevice;
	iParams.nChannels = channels;
	iParams.firstChannel = iOffset;
	oParams.deviceId = oDevice;
	oParams.nChannels = channels;
	oParams.firstChannel = oOffset;

	if ( iDevice == 0 )
		iParams.deviceId = adac.getDefaultInputDevice();
	if ( oDevice == 0 )
		oParams.deviceId = adac.getDefaultOutputDevice();

	RtAudio::StreamOptions options;
	//options.flags |= RTAUDIO_NONINTERLEAVED;

	
	//Loading impulse response and myData
	getImpulseResponse("impres", &myData.iR, &myData.iRSize );
	myData.bufferBytes = bufferFrames * channels * sizeof( MY_TYPE );
	myData.tmp = (double *) calloc(bufferFrames + myData.iRSize - 1, sizeof( MY_TYPE ));
	myData.report = fopen("report.txt","r+");
	if(myData.report == NULL){
		printf("Unable to open the report file \n");
		exit(0);
	}
		
  
	// limiter data
  
	myLimiterData.threshold = 0.2;
	myLimiterData.releaseFactor = 0.35;
	myLimiterData.smoothingFactor = 0.7;
  
  
	//use reverb1 or reverb2
	try {
		adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &reverb2, /*(void *)&myData.bufferBytes */(void *)&myData, &options );
		
	//adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &limiter, /*(void *)&myData.bufferBytes */(void *)&myLimiterData, &options );		
	}
	catch ( RtAudioError& e ) {
		std::cout << '\n' << e.getMessage() << '\n' << std::endl;
		exit( 1 );
	}

	// Test RtAudio functionality for reporting latency.
	std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

	try {
		adac.startStream();

		char input;
		std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
		std::cin.get(input);

		// Stop the stream.
		adac.stopStream();
	}
	catch ( RtAudioError& e ) {
		std::cout << '\n' << e.getMessage() << '\n' << std::endl;
		goto cleanup;
	}

	cleanup:
	if ( adac.isStreamOpen() ) adac.closeStream();

	fclose(myData.report);
	return 0;
}
