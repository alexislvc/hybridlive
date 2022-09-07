#include <math.h>

#include <fftw3.h>
#include <portaudio.h>

#define _USE_MATH_DEFINES

#define SAMPLE_RATE     (96000)    /* Hertz */
#define FRAMES_PER_BUFFER   (1600)  /* Number of points acquired during the time window */
#define SL  (800)  /* Number of spectral lines in the fft --> number of usable frequencies (SL = FRAMES_PER_BUFFER/2) */ 
#define THRESHOLD_VALUE (450) /* Threshold value for the synchronisation scheme 
                               * The value has been chosen in such way that it's impossible to 
                               * obtain this value in a normal transmission */
#define NUMBER_OF_TRANSMITTED_POINT (33) /* Number of 3D points transmited in the 96kHz configuration */
#define LENGTH_OF_FRAME (28) /* Frame length added at the beginning of the payload for transmission */
#define LENGTH_OF_PAYLOAD (3168) /* Length in terms of bits for the payloads (32 bits for 33 points in 3d) */

#define PA_SAMPLE_TYPE  (paFloat32) /* In our application, we use the highest Portaudio resolution available */

double PHASES[16];   /* Array of phases that can be used according to the 16-PSK protocol */

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

unsigned int convert_to_int(int* arr, int low, int high)
{   /* Converts a binary array to the corresponding integer.
     * This function has been realised in this way for the sake of 
     * the convert_to_float function defined just below */
	unsigned f = 0, i;
	for (i = high; i >= low; i--) {
		f = f + arr[i] * pow(2, high - i);
	}
	return f;
}

float convert_to_float(int *bin_float)
{   /* Converts the 32-bit ieee754 representation of a float to its 
     * float value */
    myfloat var; /* This (myFloat) var is useful to the conversion */
    
    /* Convert the least significant mantissa part (23 bits)
	 * to corresponding decimal integer */
	unsigned f = convert_to_int(bin_float, 9, 31);

	/* Assign integer representation of mantissa */
	var.raw.mantissa = f;

    /* Convert the exponent part (8 bits) to a corresponding decimal integer */
	f = convert_to_int(bin_float, 1, 8);

	/* Assign integer representation of the exponent */
	var.raw.exponent = f;

    /* Assign sign bit */
	var.raw.sign = bin_float[0];

    return var.f;
}

int compute_array(double *tab, double start, double step, int len)
{   /* Fills the "tab" array with the parameters we will use */
    for (int i = 0; i < len ; i++)
        tab[i] = i * step + start;
    return 0;
}

void printarray(float *arr, int size)  
{   /* Prints the content of an array in parameters */
    printf("\nElements of array are : \n");  
    for(int i = 0; i < size ;i++)  
        printf("%d : %f\n", i, arr[i]);
}  

void printarray_dec(int *arr, int size)  
{   /* Prints the content of an array in parameters */
    printf("\nElements of array are : \n");  
    for(int i = 0; i < size ;i++)  
        printf("%d : %d\n", i, arr[i]);
}

double round_to_double(float f) {
    /* Returns a close double representation of the given float */
    return round(f * pow(10, 7)) / pow(10, 7);
}

double estimate_phase(double value)
{
    /* Returns the estimated phase given in parameters, from the PHASES array. */
    /* It currently uses a decision tree algorithm, but it could be optimised. */
    /* This algorithm is designed to work with little noise. (tested on O.1 random noise magnitude) */
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

int *decode_message(fftw_complex *out)
{
    /* Takes the result of the fft as parameters and recreates the original index(phase) message */
    static int res[SL - 1];
    double temp_phase;

    for (int iFreq = 1; iFreq < SL; iFreq++) 
    {
        temp_phase = estimate_phase(atan2(out[iFreq][1], out[iFreq][0]));
        for (int ind = 0; ind < 16; ind++)
        {   /* It may be possible to optimise the process according to performance needs. */
            if (temp_phase == PHASES[ind]) /* Equality test on the double ... perhaps problematic in the future; */
            {
                res[iFreq - 1] = ind; /* -1 to get the right index in res */
                break;
            }
        }
    }
    return res;
}

int *decimal_to_binary_4bits(int value)
{   /* Returns the decimal value of the 4-bit representation of the value (int) given in the parameters */
    static int res[4];
    int temp[4];
    for (int i = 0; i < 4; i++) /* Used only to decode 1 chunk into 4 bits */
    {
        temp[i] = value % 2;
        value = value / 2;
    }
    for (int j = 0; j < 4; j++) /* Just flip temp to get the result in the right order */
        res[j] = temp[3 - j];
    return res;
}

int *back_to_binary(int encoded_message[])
{   /* Takes a phase message (from the fft phase analysis) and returns the original binary message */
    static int res[LENGTH_OF_FRAME + LENGTH_OF_PAYLOAD];
    for (int indMes = 0; indMes < SL - 1; indMes++)
    {
        int indRes = indMes * 4;
        int *temp = decimal_to_binary_4bits(encoded_message[indMes]);
        for (int i = 0; i < 4; i++)
            res[indRes + i] = temp[i]; 
    }
    return res;
}

int check_checksum(int *message)
{   /* Checks whether the checkSum calculated by the sender 
     * correspond to the payloads recieved. 
     * Return 1 if True, 0 if False */
    int checkSum = 0, calcSum = 0;
    for (int iCheck = 16; iCheck < LENGTH_OF_FRAME; iCheck++)
        checkSum += message[iCheck] * pow(2, (11 - iCheck) + 16);  /* Calculates the checkSum by 
                                                                    * the binary representation of int */
        /* The "+16" is there to compensate for the difference between iCheck and the desired exponent */
    printf("\n%d\n", checkSum); /* Here is a printout of the checksum sent with the message calculated by the sender */
    for (int iCalc = LENGTH_OF_FRAME; iCalc < LENGTH_OF_FRAME + LENGTH_OF_PAYLOAD; iCalc++)
        calcSum += message[iCalc];  /* re-count the number of ones */
    printf("\n%d\n", calcSum); /* Here is a printout of the checksum calculated by the receiver */
    if (checkSum == calcSum)
        return 1;
    else
        return 0;
}

int get_frame_number(int *message)
{   /* Returns the frame number of the given binary message (with the frame) */
    int frame_number_value = 0;

    for (int iFramNum = 0; iFramNum < 16; iFramNum++)
        frame_number_value += message[iFramNum] * pow(2, 15 - iFramNum);  /* calculate the frame number 
                                                                          * by the binary representation of int */
    return frame_number_value;
}

void back_to_float(float *res, int *message)
{   /* Takes the "float result array" and the result of the fft as parameters. 
     * Fills the float array "res" with the floats that have been transmitted
     * through the channel. */
    int temp[32]; /* Temporary buffer where the 32-bit representation of a float will be stored */
    int i;
    int index_in_message; /* Refers to the index of the payload in the "bin_message" array */

    for (i = 0; i < NUMBER_OF_TRANSMITTED_POINT * 3; i++)
    {
        index_in_message = i * 32 + LENGTH_OF_FRAME;
        for (int j = 0; j < 32; j++) 
        {   /* This loop is used to fill the temp buffer */
            /* This process might be optimized */
            temp[j] = message[index_in_message + j];
        }
        res[i] = convert_to_float(temp);
    }
}

int init()
{   /* Initializes all the parameters of our program */
    compute_array(PHASES, - M_PI, M_PI / 8, 16); /* Same but with the phases */
    return 0;
}

int main()
{
    init(); /* Initializes the PHASES table */
    float buffer[FRAMES_PER_BUFFER - 1]; /* Buffer for the actual payload message */
    /* The "-1" in the buffer's length is due to the synchronisation scheme */
    float buffer_init[1]; int recv_threshold = 0;
    float lowBound_threshold = THRESHOLD_VALUE - 0.05; /* The +/- 0.05 gives a margin */
    float upBound_threshold = THRESHOLD_VALUE + 0.05;  /* for possible noise on the channel */
    float lowBound_silent = - 0.05; /* This boundaries are here to defined the silent */
    float upBound_silent = 0.05;    /* phase of the transmission channel */
    int cond_wait = 1;

    /* --- PORTAUDIO PARAMTERS --- */
    PaStreamParameters inputParameters;
    PaStream *stream;
    PaError err;
    const PaDeviceInfo *inputDeviceInfo;

    /* --- FFTW PARAMETERS --- */
    fftw_complex *out;
    fftw_plan p;
    float res_fttw[NUMBER_OF_TRANSMITTED_POINT * 3]; /* Result from the ftt for the 16 points in 3d */

    /* Prepares for transmission via PortAudio */
    printf("PortAudio config : input sine wave. SR = %d, BufSize = %d\n", SAMPLE_RATE, FRAMES_PER_BUFFER);

    /* Initializes the portaudio library to be used 
     * This function MUST be called before using any other PortAudio API functions 
     */
    err = Pa_Initialize();
    if( err != paNoError ) goto error;

    /* -- setup input -- */
    /* Here we have defined the input device parameter */
    //inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    inputParameters.device = 0;
    printf( "Input device # %d.\n", inputParameters.device );
    inputDeviceInfo = Pa_GetDeviceInfo( inputParameters.device );
    printf( "    Name: %s\n", inputDeviceInfo->name );
    
    inputParameters.channelCount = 1; //mono input
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = inputDeviceInfo->defaultHighInputLatency ;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    /* -- Open stream -- */
    err = Pa_OpenStream(
              &stream,
              &inputParameters,
              NULL, /* no output */
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              NULL, /* no callback, use blocking API */
              NULL ); /* no callback, so no callback userData */
    if( err != paNoError ) goto error;

    /* -- Start stream -- */
    err = Pa_StartStream( stream );
    if( err != paNoError ) goto error;

    /* With this loop, the program looks for the start of the transmission :
     * Once at least 4 THRESHOLD_VALUE are detected, the program can proceed
     * to the decryption phase. */

    while(cond_wait)
    {
        Pa_ReadStream( stream, buffer_init, 1);

        if (lowBound_threshold < buffer_init[0] && buffer_init[0] < upBound_threshold )
        {   /* Nominal case : it waits for the THRESHOLD_VALUEs */
            recv_threshold++; /* Here it counts the number of THRESHOLD_VALUE received */
        } 
        else if (recv_threshold >= 4 && !(lowBound_silent < buffer_init[0] && buffer_init[0] < upBound_silent))
        {   /* Nominal case : It has received enough THRESHOLD_VALUE */
            cond_wait = 0; /* If it has received enough THRESHOLD_VALUE, the transmission can start */
        } 
        else if (recv_threshold < 4 && !(lowBound_silent < buffer_init[0] && buffer_init[0] < upBound_silent))
             /* In this case, we considered that it did not receive enough THRESHOLD_VALUE */
        {
            printf ("Unable to start the reception properly.\n");
            printf ("Not enough \"THRESHOLD_VALUE\" received.\n");
            printf ("A packet will be lost.\n\n");
            recv_threshold = 0;
            /* Here the program waits for a next eligible buffer */
            err = Pa_ReadStream( stream, buffer, FRAMES_PER_BUFFER - 1);
            if( err ) goto xrun;
        } 
    }

    /* -- Here is the loop where we get the data from the input -- */
    // You may get underruns or overruns if the output is not primed by PortAudio.
    /* Here we put the usable data in the "buffer" */
    err = Pa_ReadStream( stream, buffer, FRAMES_PER_BUFFER - 1);
    if( err ) goto xrun;
    /* The "-1" is there because we consider that we have already received the
     * beginning of the transmission in the buffer_init */

    printf("Wire off.\n"); fflush(stdout);
    
    /* -- Now the stream stops -- */
    err = Pa_StopStream( stream );
    if( err != paNoError ) goto error;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto error;

    /* -- The Port audio library is disactivated -- */
    Pa_Terminate();
    
    /* ------ FFTW PHASE -----*/
    /* Process out */
    int nc = ( FRAMES_PER_BUFFER / 2 ) + 1;
    out = fftw_malloc(sizeof(fftw_complex) * nc);

    /* The array of floats "buffer" is replaced by an array of doubles "in_d".
     * It is necessary for the argument of the fftw function during
     * the demodulation phase */
    double in_d[FRAMES_PER_BUFFER];
    in_d[0] = round_to_double(buffer_init[0]); /* The start of the transmission is in the buffer_init */
    for(int i = 1; i < FRAMES_PER_BUFFER; i++)
        in_d[i] = round_to_double(buffer[i - 1]);
    /* Defines the parameters of the plan p */
    p = fftw_plan_dft_r2c_1d (FRAMES_PER_BUFFER, in_d, out, FFTW_ESTIMATE);

    /* Excutes the fft on the plan p for the given parameters */
    fftw_execute(p);

    printarray_dec(decode_message(out), SL - 1);

    /* The received binary message is stored in the variable message */
    int *message = back_to_binary(decode_message(out));

    printf("\n CheckSUM = %d \n", check_checksum(message));
    printf("The frame number is : %d\n", get_frame_number(message));
    
    back_to_float(res_fttw, message);
    //printarray(res_fttw, NUMBER_OF_TRANSMITTED_POINT * 3); /* Debug purpose */

    fftw_destroy_plan(p);
    fftw_free(out);

    return err;
    xrun:
        printf("err = %d\n", err); fflush(stdout);
        if( stream ) {
           Pa_AbortStream( stream );
           Pa_CloseStream( stream );
        }
        Pa_Terminate();
        if( err & paInputOverflow )
           fprintf( stderr, "Input Overflow.\n" );
        if( err & paOutputUnderflow )
           fprintf( stderr, "Output Underflow.\n" );
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