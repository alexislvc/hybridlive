#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include <portaudio.h>

#define _USE_MATH_DEFINES

#define SAMPLE_RATE     (48000)   /* Hertz */
#define MAXIMUM_FRAME_RATE  (60)   /* Maximum frame per second of the program */
/*  The programme must transmit 60 frames per second at a constant rate.
 *  If at the time of tranmitting a data, mediapipe haven't produce it : 
 *  the programm will send a "0...0" message 
 */
#define FRAMES_PER_BUFFER  (800)  /* Number of acquired points during the time window (FRAMES_PER_BUFFER = SAMPLE_RATE / MAXIMUM_FRAME_RATE) */
#define SL  (400)  /* Number of spectral line in the fft --> number of usable frequencies (SL = FRAMES_PER_BUFFER/2) */

#define PA_SAMPLE_TYPE  (paFloat32) /* In our app we use the highest portaudio resolution available */

double FREQUENCIES[SL - 1]; /* Array of useable frequencies according to the parameters; f = 0 Hz is not used */
double PHASES[16];   /* Array of useable phases according to the protocole 16-PSK*/

char BITS[16][5] =  {   /* Array of bits to be transmitted */
    "0000", /* 0 */     /* The index of the String in this array can be translated to */
    "0001", /* 1 */     /* the index of the phase transmitted in the carrier signal in PHASES */
    "0010", /* 2 */
    "0011", /* 3 */     /* It may exist a better option to implement that idea for the transmission */
    "0100", /* 4 */     /* We can assign to each index his representation in binary */
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

typedef union { /* Union useful to the translation of float into ieee754 floating point */
	float f;    /* I don't get how the whole the translation work atm */
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
{   /* This function fill the array of the parameters we will use */
    for (int i = 0; i < len ; i++)
        tab[i] = i * step + start;
    return 0;
}

void printarray(float *arr, int size)  
{   /* Print the content of an array in parameters*/
    printf("\nElements of array are : \n");  
    for(int i = 0; i < size ;i++)  
        printf("%d : %f\n", i, arr[i]);
}  

void printarray_dec(int *arr, int size)  
{   /* Print the content of an array in parameters*/
    printf("\nElements of array are : \n");  
    for(int i = 0; i < size ;i++)  
        printf("%d : %d\n", i, arr[i]);
}

double round_to_decimal(float f) {
    /* Return a close double representation of the given float */
    return round(f * pow(10, 7)) / pow(10, 7);
}

double estimate_phase(double value)
{
    /* This function return the estimated phase given in parameters from the array PHASES*/
    /* It's currently using a decision tree algorithm but might be optimised */
    /* This algorithm still work with a little bit of noise (tested on O.1 rad noise magnitude) */
    if (value < - M_PI/16)
    {
        if (value < - (9 * M_PI)/16)
        {
            if (value < - (13 * M_PI)/16)
            {
                if (value < - (15 * M_PI)/16)   
                    return PHASES[0];
                else    
                    return PHASES[1]; 
            }
            else
            {
                if (value < - (11 * M_PI)/16)   
                    return PHASES[2];
                else    
                    return PHASES[3];
            }
        }
        else
        {
            if (value < - (5 * M_PI)/16)
            {
                if (value < - (7 * M_PI)/16)    
                    return PHASES[4];
                else    
                    return PHASES[5];
            }
            else 
            {
                if (value < - (3 * M_PI)/16)    
                    return PHASES[6];
                else    
                    return PHASES[7];
            }
        }
    }
    else
    {
        if (value < (7 * M_PI)/16)
        {
            if (value < (3 * M_PI)/16)
            {
                if (value < M_PI/16)    
                    return PHASES[8];
                else    
                    return PHASES[9];
            }
            else
            {
                if (value < (5 * M_PI)/16)  
                    return PHASES[10];
                else   
                    return PHASES[11];
            }
        }
        else
        {
            if (value < (11 * M_PI)/16)
            {
                if (value < (9 * M_PI)/16) 
                    return PHASES[12];
                else    
                    return PHASES[13];
            }
            else
            {
                if (value < (15 * M_PI)/16)
                {
                    if (value < (13 * M_PI)/16)  
                        return PHASES[14];
                    else    
                        return PHASES[15];
                }
                else
                    return PHASES[0];
            }
        }
    }
}

float *generate_sine_wave(double frequencie, double phase)
{   /* Return an array of size FRAMES_PER_BUFFER with sin values at desired frequencie */
    static float res[FRAMES_PER_BUFFER];
    typedef unsigned long phase_t;
    double maxphase = (double)((phase_t)0-(phase_t)1)+1.0;

    phase_t delta = maxphase*frequencie/SAMPLE_RATE+0.5, iphase = - delta;

    for (long i = 0; i < FRAMES_PER_BUFFER; i++)
        res[i] = (float) cos((iphase += delta)/maxphase * (2 * M_PI) + phase);
    return res;
}

int *translate_message (int message[])
{   /* Takes a message and return the index of the PHASES list needed */
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

int encode_message(float *in, int list_index[])
{
    /* Generate the modulate signal to be transmitted */
    /* encode_message modify the array given "in" parameters without returning anything */
    for (int iFreq = 0; iFreq < SL - 1; iFreq++) /* f = 0 is not used */
    {
        float *sine = generate_sine_wave(FREQUENCIES[iFreq], PHASES[list_index[iFreq]]);
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {   /* It's might be possible to optimise the process for performance need */
            if (iFreq == 0 )    
                in[i] = sine[i]; /* On the first iteration we reset the signal sent */
            else    
                in[i] += sine[i];
        }
    }
    return 0;
}

int getBinary(int n, int i, int res[], int start)
{
	/* Get the binary representation of a number n up to i-bits */
    /* it had it into the array given in parameters. Stating at the index start */
    /* It was made like this in the interest of the function ieee754_converter defined just bellow */
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
{   /* Return an array with the ieee754 (32-bits) representation of the given myfloat /!\ */
    static int res[32];

    res[0] = var.raw.sign;
    getBinary(var.raw.exponent, 8, res, 1);
    getBinary(var.raw.mantissa, 23, res, 9);

    return res;
}

unsigned short get_checksum(int message[])
{   /* Take a binary message in paramaters and return the number of "1" 
     * present in it for the checksum at the receiver */
    unsigned short res = 0;
    for (int i = 0; i < 1536; i++)
    {
        if (message[i] == 1)
            res += 1;
    }
    return res;
}

int *decimal_to_binary(unsigned short value)
{   /* Return the 16 Bits representations of the unsigned short given in parameters */
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

int *decimal_to_binary_4bits(int value)
{   /* Return the decimal value of the 4 Bits representations of the int given in parameters */
    static int res[4];
    int temp[4];
    for (int i = 0; i < 4; i++) /* Only used to decode 1 chunck into 4 bits */
    {
        temp[i] = value % 2;
        value = value / 2;
    }
    for (int j = 0; j < 4; j++) /* just flip temp to get the result in the right order */
        res[j] = temp[3 - j];
    return res;
}

int *generate_random_message()
{   /* Generate a message that could be transmitted from mediapipe */
    /* 3 dimensions of 16 points coded in float (32 BITS) put one after the other */
    static int res[1536];
    myfloat var;

    for (int point = 0; point < 16; point++)
    {
        for (int coord = 0; coord < 3; coord++)
        {
            int pos = point * 32 * 3 + coord * 32;
            var.f = (((float)rand()/(float)(RAND_MAX)) * 200) - (float)100; /* generate random number in [-100; 100]*/
            int *listB = ieee754_converter(var);
            for (int i = 0; i < 32; i++)
                res[pos + i] = listB[i];
        }
    }
    return res;
}

int *message_formatting(int message[])
{   /* The function return a array formatted according to the protocole and ready to be transmitted */
    /*    0 - 15 BITS : Frame number */
    /*   16 - 27 BITS : Checksum */
    /* 28 - 1563 BITS : message (Payload) */
    #if !defined(frameNum)
        static unsigned short frameNum = 0;
    #endif
    static int res[1596];
    int index;
    
    frameNum = frameNum % USHRT_MAX;  /* Here we don't want to exceed the limit of unsigned short */
    int *binFrameNum = decimal_to_binary(frameNum);
    for (index = 0; index < 16; index++)    /* Adding the frame number */
        res[index] = binFrameNum[index];

    int *binChecksum = decimal_to_binary(get_checksum(message)); /* The 4 fist digit will always be 0 */
    //printf ( "CheckSum = %hu", get_checksum(message));
    for (index = 16; index < 28; index++)   /* Adding the check sum */
        res[index] = binChecksum[index - 16 + 4]; /* The first 4 first digits of binCheckSum are always 0 */

    for (index = 28; index < 1564; index++) /* Adding the payload */
        res[index] = message[index - 28];

    frameNum++;
    return res;
}

int init()
{   /* Initialize all the parameters of our program */
    compute_array(FREQUENCIES, MAXIMUM_FRAME_RATE, MAXIMUM_FRAME_RATE, SL - 1); /* Compute the array of usable frequencies */
    compute_array(PHASES, - M_PI, M_PI / 8, 16); /* Same but with the phases */
    return 0;
}

int main()
{
    init(); /* Initiate the FREQUENCIES and PHASES arrays */
    srand(time(NULL)); /* Allows to have a different random message sent at each time */
    float in_message[FRAMES_PER_BUFFER]; /* Buffer where the message to sent will be written */

    /* Create and process the buffer_init that indicate the start of the transmission */
    float buffer_init[5]; for (int i = 0; i < 4; i++) buffer_init[i] = 1.0f; 
    buffer_init[4] = 2.0f;

    /* --- PORTAUDIO PARAMTERS --- */
    PaStreamParameters outputParameters;
    PaStream *stream;
    PaError err;
    const PaDeviceInfo *outputDeviceInfo;

    /* Generate a random message */
    int *glo_message = message_formatting(generate_random_message());

    /* Process "in" the signal to be transmitted */
    encode_message(in_message, translate_message(glo_message));

    /* Prepare the transmission through PortAudio */
    printf("PortAudio config : output sine wave. SR = %d, BufSize = %d\n", SAMPLE_RATE, FRAMES_PER_BUFFER);
    
    /* Initialize the portaudio library to be used 
     * This function MUST be called before using any other PortAudio API functions 
     */
    err = Pa_Initialize();
    if( err != paNoError ) goto error;

    /* -- setup output -- */
    /* Here we defined the parameter of the output device */
    //outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    outputParameters.device = 1;
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
    
    /* This one mark the begining of the transmission 
     * 4 ones (1.0f) are sent to initiate the transmission of a message
     * The transmission 
     */
    Pa_WriteStream( stream, buffer_init, 5 );

    /* The in_message is being transmitted */
    err = Pa_WriteStream( stream, in_message, FRAMES_PER_BUFFER );
            if( err != paNoError ) goto error;

    printf("Wire off.\n"); fflush(stdout);
    
    /* -- Now we stop the stream -- */
    err = Pa_StopStream( stream );
    if( err != paNoError ) goto error;

    /* -- We end the Port audio library -- */
    Pa_Terminate();

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