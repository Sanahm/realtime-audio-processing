/******************************************/
/*
  reverb.cpp
  by SANA M. and BERTRAND E. 2017-2018.

  This program opens a stream and add some reverberation to
  input before passing to the output.
*/
/******************************************/

#include "RtAudio.h"
#include "library.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

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
	//TODO
	//...
}

/*Reverberation with fft method*/
int reverb2( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double /*streamTime*/, RtAudioStreamStatus status, void *myData )
{
	MyData * pData = (MyData *) myData;
	double * h = (double *) pData->iR; //impulsionnel response
	double * x = (double *) inputBuffer; // signal x
	unsigned int M = pData->iRSize; //size of h
	unsigned int L = nBufferFrames; // size of x
	int i;
	
	double * conv = fconv(x, L, h, M);
	double * temp;
	memcpy(temp, pData->tmp, sizeof(double)*(L+M-1));
	// Since the number of input and output channels is equal, we can do
	// a simple buffer copy operation here.
	
	if ( status ) std::cout << "Stream over/underflow detected." << std::endl;
	
	// Calcul of the convolution with fft
	for( i = 0; i < L+M-1; i++ ){
		temp[i] += conv[i];
	}
	
	memcpy(pData->tmp, temp+L, (M-1)*sizeof(double));
	
	unsigned int bytes = (unsigned int ) pData->bufferBytes;
	memcpy( outputBuffer, temp, bytes );
	return 0;
}

/*Limiter implementation*/
int limiter( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double /*streamTime*/, RtAudioStreamStatus status, void *myLimiterData )
{
	MyLimiterData pData = (MyLimiterData *) myLimiterData;
	double * env = dgetmem(nBufferFrames);
	double * gain = dgetmem(nBufferFrames);
	double * sgain = dgetmem(nBufferFrames);
	int i, u = 0;
	
	env[0] = inputBuffer[0];
	for(i = 1; i < nBufferFrames; i++) {
		env[i] = MAX(inputBuffer[i], env[i-1]*pData->releaseFactor);
	}
	
	for(i = 0; i < nBufferFrames; i++) {
		if( env[i] < pData->threshold )
			gain[i] = 1;
		else
			gain[i] = 1 + pData->threshold - env[i];

		sgain[i] = u*pData->smoothingFactor + gain[i]*(1 - pData->smoothingFactor);
		u = sgain[i];
	}
	
	for(i = 0; i < nBufferFrames; i++) {
		outputBuffer[i] = inputBuffer[i]*sgain[i];
	}
	return 0;
}

int main( int argc, char *argv[] )
{
	MyData myData;
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
	myData.tmp = (double *) calloc( myData.iRSize - 1, sizeof( MY_TYPE ));
  
  
	//use reverb1 or reverb2
	try {
		adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &reverb2, (void *)&myData, &options );
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

	return 0;
}
