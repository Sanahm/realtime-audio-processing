/**
* @file: include/library.h
* @author: M. SANA  <Mohamed.Sana@phelma.grenoble-inp.fr>
* @date: 24/12/2017 09:50:05
* @brief: 
*/

#ifndef _library_
#define _library_


#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

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
} MyData;

typedef struct {
	double threshold;
	double releaseFactor;
	double smoothingFactor;
	double * lookahead;
} MyLimiterData;

//g++ -Wall -D__LINUX_ALSA__ reverb.cpp library.cpp ../rtaudio-5.0.0/RtAudio.cpp -I "../include/" -I "../rtaudio-5.0.0/RtAudio.h" -o reverb -lasound -lpthread

#endif
