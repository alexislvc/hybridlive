#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include <portaudio.h>

#define _USE_MATH_DEFINES

#define SAMPLE_RATE     (48000)   /* Hertz */
#define MAXIMUM_FRAME_RATE  (60)   /* Maximum number of frames per second of the programme */
#define FRAMES_PER_BUFFER  (800)  /* Number of points acquired during the time window (FRAMES_PER_BUFFER = SAMPLE_RATE / MAXIMUM_FRAME_RATE) */
#define SL  (400)  /* Number of spectral lines in the fft --> number of usable frequencies (SL = FRAMES_PER_BUFFER/2) */ 
#define THRESHOLD_VALUE (450) /* Threshold value for the synchronisation scheme 
                               * The value has been chosen in such way that it's impossible to 
                               * obtain this value in a normal transmission */
#define LENGTH_OF_FRAME (28) /* Frame length added at the beginning of the payload for transmission */
#define LENGTH_OF_PAYLOAD (1536) /* Length in terms of bits for the payloads (32 bits for 16 points in 3d) */

#define PA_SAMPLE_TYPE  (paFloat32) /* In our application, we use the highest Portaudio resolution available */

double FREQUENCIES[SL - 1]; /* Table of usable frequencies depending on the parameters; f = 0 Hz is not used. */
double PHASES[16];   /* Array of phases that can be used according to the 16-PSK protocol */

char BITS[16][5] =  {   /* Array of bits to be transmitted */
    "0000", /* 0 */     /* The index of the String in this table can be translated into */
    "0001", /* 1 */     /* the index of the phase transmitted in the carrier signal in PHASES */
    "0010", /* 2 */
    "0011", /* 3 */     /* There may be a better option to implement this idea for transmission */
    "0100", /* 4 */     /* We can assign to each index its representation in binary */
    "0101", /* 5 */
    "0110", /* 6 */
    "0111", /* 7 */
    "1000", /* 8 */
    "1001", /* 9 */
    "1010", /* 10 */
    "1011", /* 11 */
    "1100", /* 12 */
    "1101", /* 13 */
    "1110", /* 14 */
    "1111"  /* 15 */
};

typedef union { /* Union useful for translating the floating point ieee754 to its float value */
	float f;
	struct
	{
		// Order is important.
		// Here the members of the union data structure
		// use the same memory (32 bits).
		// The ordering is taken
		// from the LSB to the MSB.
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} raw;
} myfloat;

int compute_array(double *tab, double start, double step, int len)
{   /* Fills the "tab" array with the parameters we will use */
    for (int i = 0; i < len ; i++)
        tab[i] = i * step + start;
    return 0;
}

void printarray(float *arr, int size)  
{   /* Prints the content of an array in parameters*/
    printf("\nElements of array are : \n");  
    for(int i = 0; i < size ;i++)  
        printf("%d : %f\n", i, arr[i]);
}  

void printarray_dec(int *arr, int size)  
{   /* Prints the content of an array in parameters*/
    printf("\nElements of array are : \n");  
    for(int i = 0; i < size ;i++)  
        printf("%d : %d\n", i, arr[i]);
}

float *generate_sine_wave(double frequencie, double phase)
{   /* Returns an array of size FRAMES_PER_BUFFER with sin values
     * at the desired frequency with the desired phase shift */
    static float res[FRAMES_PER_BUFFER];
    typedef unsigned long phase_t;
    double maxphase = (double)((phase_t)0-(phase_t)1)+1.0;

    phase_t delta = maxphase*frequencie/SAMPLE_RATE+0.5, iphase = - delta;

    for (long i = 0; i < FRAMES_PER_BUFFER; i++)
        res[i] = (float) cos((iphase += delta)/maxphase * (2 * M_PI) + phase);
    return res;
}

int *translate_message (int message[])
{   /* Takes a message and returns the index of the list of required PHASES */
    static int res[SL - 1];
    char temp[5];
    for (int iRes = 0; iRes < SL - 1; iRes++) 
    {
        int iMes = iRes * 4;
        sprintf(temp, "%d%d%d%d", message[iMes], message[iMes+1], message[iMes+2], message[iMes+3]);
        for (int ind = 0; ind < 16; ind++)
        {
            if (strcmp(temp, BITS[ind]) == 0)
                {
                    res[iRes] = ind;
                    break;
                }
        }
    }
    return res;
}

void encode_message(float *in, int list_index[])
{
    /* Generates the modulated signal to be transmitted */
    /* encode_message modifies the array given as parameters "in" without returning anything */
    for (int iFreq = 0; iFreq < SL - 1; iFreq++) /* f = 0 is not used */
    {
        float *sine = generate_sine_wave(FREQUENCIES[iFreq], PHASES[list_index[iFreq]]);
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {   /* It may be possible to optimise the process according to performance needs */
            if (iFreq == 0 )    
                in[i] = sine[i]; /* At the first iteration, we reset the signal sent */
            else    
                in[i] += sine[i];
        }
    }
}

int getBinary(int n, int i, int res[], int start)
{
	/* Get the binary representation of a number n up to i-bits */
    /* It adds it to the table given in the parameters. Starting at the index start.
     * This function has been realised in this way for the sake of 
     * the ieee754_converter function defined just below */
	int k;
	for (k = i - 1; k >= 0; k--) 
    {
		if ((n >> k) & 1)
			res[start + i - 1 - k] = 1;
		else
			res[start + i - 1 - k] = 0;
	}
    return 0;
}

int *ieee754_converter(myfloat var)
{   /* Returns an array with the ieee754 (32 bit) representation of the given myfloat /!\ */
    static int res[32];

    res[0] = var.raw.sign;
    getBinary(var.raw.exponent, 8, res, 1);
    getBinary(var.raw.mantissa, 23, res, 9);

    return res;
}

unsigned short get_checksum(int message[])
{   /* Take a binary message (without frame) as parameters and return the number
     * of "1's" it contains for the checksum at the receiver */
    unsigned short res = 0;
    for (int i = 0; i < LENGTH_OF_PAYLOAD; i++)
    {
        if (message[i] == 1)
            res += 1;
    }
    return res;
}

int *decimal_to_binary(unsigned short value)
{   /* Returns the 16-bit representations of the unsigned short given in the parameters */
    static int res[16];
    int temp[16];
    for (int i = 0; i < 16; i++)
    {
        temp[i] = value % 2;
        value = value / 2;
    }
    for (int j = 0; j < 16; j++) /* just flip temp to get the result in the right order */
        res[j] = temp[15 - j]; 
    return res;
}

int *generate_random_message()
{   /* Generates a message that could be created from mediapipe */
    /* 3 dimensions of 16 points coded in float (32 BITS) put one after the other */
    static int res[LENGTH_OF_PAYLOAD];
    myfloat var;
    int index_final; /* Index of the float of interest in the results array */
    int index_point; /* Index of the 3D point of interest in the results array */
    int index_coordinate; /* Index of the coordinate of interest in the 3D point of interest */

    for (int point = 0; point < 16; point++)
    {
        index_point = point * 32 * 3; /* At each new 3D point, we move the index by 32 bits 
                                       * for 3 floating point coordinates */
        for (int coord = 0; coord < 3; coord++)
        {
            index_coordinate = coord * 32; /* At each coordinate, the index is moved by 32 bits */
            /* Each coordinate being x, y, z */

            index_final = index_point + index_coordinate;
            var.f = (((float)rand()/(float)(RAND_MAX)) * 200) - (float)100; /* Generates random float number in [-100; 100]*/
            //printf("The chosen value is : %f \n", var.f); /* Debug purpose */
            int *listB = ieee754_converter(var);
            for (int i = 0; i < 32; i++)
                res[index_final + i] = listB[i];
        }
    }
    return res;
}

int *message_formatting(int message[])
{   /* Returns an array formatted according to the protocol and ready to be transmitted */
    /*      0 - 15 BITS : Frame number */
    /*     16 - 27 BITS : Checksum */
    /*   28 - 1563 BITS : message (Payload) */
    /* 1564 - 1596 BITS : Unused room (32 bits unused) */
    /* In this configuration, the length of the frame added to the message is 28.
     * If this length is to be changed in the future, the global parameter 
     * LENGTH_OF_FRAME must also be changed in the receiver program. */

    #if !defined(frameNum)
        static unsigned short frameNum = 0;
    #endif
    static int res[LENGTH_OF_FRAME + LENGTH_OF_PAYLOAD + 32]; /* The final binary message that will be sent */
    int index;
    
    frameNum = frameNum % USHRT_MAX;  /* Here we don't want to exceed the limit of unsigned short */
    int *binFrameNum = decimal_to_binary(frameNum);
    for (index = 0; index < 16; index++)    /* Adding the frame number */
        res[index] = binFrameNum[index];

    int *binChecksum = decimal_to_binary(get_checksum(message)); /* The 4 fist digit will always be 0 
                                                                  * because of the number of "ones" that can be sent */
    for (index = 16; index < 28; index++)   /* Adding the check sum */
        res[index] = binChecksum[index - 16 + 4]; /* The first 4 first digits of binCheckSum are always 0 */

    for (index = LENGTH_OF_FRAME; index < LENGTH_OF_FRAME + LENGTH_OF_PAYLOAD; index++) /* Adding the payload */
        res[index] = message[index - LENGTH_OF_FRAME];

    frameNum++;
    return res;
}

int init()
{   /* Initializes all the parameters of our program */
    compute_array(FREQUENCIES, MAXIMUM_FRAME_RATE, MAXIMUM_FRAME_RATE, SL - 1); /* Compute the array of usable frequencies */
    compute_array(PHASES, - M_PI, M_PI / 8, 16); /* Same but with the phases */
    return 0;
}

int main()
{
    init(); /* Initializes the FREQUENCIES and PHASES tables */
    srand(time(NULL)); /* Allows you to send a different random message each time */
    float in_message[FRAMES_PER_BUFFER]; /* Buffer where the message to sent will be written */

    /* Creates and processes the buffer_init which indicates the start of the transmission */
    float buffer_init[5]; for (int i = 0; i < 5; i++) buffer_init[i] = THRESHOLD_VALUE;

    /* --- PORTAUDIO PARAMETERS --- */
    PaStreamParameters outputParameters;
    PaStream *stream;
    PaError err;
    const PaDeviceInfo *outputDeviceInfo;

    /* Generates a random message */
    int *glo_message = message_formatting(generate_random_message());

    /* Processes "in_message" the signal to be transmitted */
    encode_message(in_message, translate_message(glo_message));

    /* Prepares for transmission via PortAudio */
    printf("PortAudio config : output sine wave. SR = %d, BufSize = %d\n", SAMPLE_RATE, FRAMES_PER_BUFFER);
    
    /* Initializes the portaudio library to be used 
     * This function MUST be called before using any other PortAudio API functions 
     */
    err = Pa_Initialize();
    if( err != paNoError ) goto error;

    /* -- setup output -- */
    /* Here we have defined the output device parameter */
    //outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    outputParameters.device = 0; /* Index of the output devices */
    printf( "Output device # %d.\n", outputParameters.device );
    outputDeviceInfo = Pa_GetDeviceInfo( outputParameters.device );
    printf( "    Name: %s\n", outputDeviceInfo->name );

    outputParameters.channelCount = 1;       /* mono output */
    outputParameters.sampleFormat = PA_SAMPLE_TYPE; /* 32 bit floating point output */
    outputParameters.suggestedLatency = outputDeviceInfo->defaultHighInputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    /* -- Open stream -- */
    err = Pa_OpenStream(
              &stream,
              NULL, /* no input */
              &outputParameters,
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              NULL, /* no callback, use blocking API */
              NULL ); /* no callback, so no callback userData */
    if( err != paNoError ) goto error;

    /* -- Start stream -- */
    err = Pa_StartStream( stream );
        if( err != paNoError ) goto error;
    
    /* These THRESHOLD_VALUEs mark the beginning of the transmission.
     * 5 THRESHOLD_VALUE are sent to initiate the transmission of a message */
    Pa_WriteStream( stream, buffer_init, 5 );

    /* The in_message is being transmitted */
    err = Pa_WriteStream( stream, in_message, FRAMES_PER_BUFFER );
            if( err != paNoError ) goto error;

    printf("Wire off.\n"); fflush(stdout);
    
    /* -- Now the stream stops -- */
    err = Pa_StopStream( stream );
    if( err != paNoError ) goto error;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto error;

    /* -- The Port audio library is disactivated -- */
    Pa_Terminate();

    /* Prints the transmitted message for debugging purposes */
    printarray_dec(translate_message(glo_message), SL - 1);

    return err;
    error:
        fprintf( stderr, "An error occured while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
	    // Print more information about the error.
	    if( err == paUnanticipatedHostError )
	    {
	    	const PaHostErrorInfo *hostErrorInfo = Pa_GetLastHostErrorInfo();
	    	fprintf( stderr, "Host API error = #%ld, hostApiType = %d\n", hostErrorInfo->errorCode, hostErrorInfo->hostApiType );
	    	fprintf( stderr, "Host API error = %s\n", hostErrorInfo->errorText );
	    }
        Pa_Terminate();
        return err;
}