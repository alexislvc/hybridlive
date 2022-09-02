#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <fftw3.h>
#include <portaudio.h>

#define _USE_MATH_DEFINES

#define SAMPLE_RATE     (48000)    /* Hertz */
#define MAXIMUM_FRAME_RATE  (60)   /* Maximum frame per second of the program */
#define FRAMES_PER_BUFFER   (800)  /* Number of acquired points during the time window (FRAMES_PER_BUFFER = SAMPLE_RATE / MAXIMUM_FRAME_RATE) */
#define SL  (400)  /* Number of spectral line in the fft --> number of usable frequencies (SL = FRAMES_PER_BUFFER/2) */ 

#define PA_SAMPLE_TYPE  (paFloat32) /* In our app we use the highest portaudio resolution available */

double PHASES[16];   /* Array of useable phases according to the protocole 16-PSK*/

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

int compute_array(double *tab, double start, double step, int len)
{   /* This function fill the array of the parameters we will use */
    for (int i = 0; i < len ; i++)
        tab[i] = i * step + start;
    return 0;
}

double round_to_double(float f) {
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

int *decode_message(fftw_complex *out)
{
    /* This function take the fft result in parameters and recreate the original index(phase) messsage */
    static int res[SL - 1];
    double temp_phase;

    for (int iFreq = 1; iFreq < SL; iFreq++) 
    {
        temp_phase = estimate_phase(atan2(out[iFreq][1], out[iFreq][0]));
        for (int ind = 0; ind < 16; ind++)
        {   /* It's might be possible to optimise the process for performance need */
            if (temp_phase == PHASES[ind]) /* Equality test on double ... maybe problematic in the future */
            {
                res[iFreq - 1] = ind; /* -1 to get the right postition in res */
                break;
            }
        }
    }
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

int *back_to_binary(int encoded_message[])
{   /* The function take an phase message (from the fft phase analyse) and return the original binary message */
    static int res[3168];
    for (int indMes = 0; indMes < SL - 1; indMes++)
    {
        int indRes = indMes * 4;
        int *temp = decimal_to_binary_4bits(encoded_message[indMes]);
        for (int i = 0; i < 4; i++)
            res[indRes + i] = temp[i]; 
    }
    return res;
}

int check_checksum(int message[])
{   /* The function verify if the check sum correspond to the payloads recieved. Return 1 if True, 0 if False */
    int checkSum = 0, calcSum = 0;
    for (int iCheck = 16; iCheck < 28; iCheck++)
        checkSum += message[iCheck] * pow(2, (11 - iCheck) + 16);  /* calculate the checkSum through the binary representation of int */
    printf("\n%d\n", checkSum);
    for (int iCalc = 28; iCalc < 1564; iCalc++)
        calcSum += message[iCalc];  /* re-count the number of ones */
    printf("\n%d\n", calcSum);
    if (checkSum == calcSum)
        return 1;
    else
        return 0;
}

int init()
{   /* Initialize all the parameters of our program */
    compute_array(PHASES, - M_PI, M_PI / 8, 16); /* Same but with the phases */
    return 0;
}

int main()
{
    init();
    float buffer[FRAMES_PER_BUFFER]; /* Buffer for the actual payload message */
    float buffer_init[1]; int recv_one = 0;
    int cond_wait = 1;

    /* --- PORTAUDIO PARAMTERS --- */
    PaStreamParameters inputParameters;
    PaStream *stream;
    PaError err;
    const PaDeviceInfo *inputDeviceInfo;

    /* --- FFTW PARAMETERS --- */
    fftw_complex *out;
    fftw_plan p;

    /* Prepare the transmission through PortAudio */
    printf("PortAudio config : input sine wave. SR = %d, BufSize = %d\n", SAMPLE_RATE, FRAMES_PER_BUFFER);

    /* Initialize the portaudio library to be used 
     * This function MUST be called before using any other PortAudio API functions 
     */
    err = Pa_Initialize();
    if( err != paNoError ) goto error;

    /* -- setup input -- */
    //inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    inputParameters.device = 1;
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

    /* With this loop the program looking for the start of the transmission
     * Once a "1.0" is detected, the program can move onto the decryption phase */
    while(cond_wait)
    {
        Pa_ReadStream( stream, buffer_init, 1);

        if (0.995 < buffer_init[0] && buffer_init[0] < 1.005 )
        {   /* Nominal case : we wait for ones (1.0f) */
            recv_one++; /* Here we count the number of one (1.0f received) */
        } 
        else if (recv_one >= 3 && (1.995 < buffer_init[0] && buffer_init[0] < 2.005 ))
        {   /* Nominal case : we've received enough ones and the two */
            cond_wait = 0; /* If we've received enough one and a two (2.0f) then the transmission can begin */
        } 
        else if ((!(recv_one >= 3) && (1.995 < buffer_init[0] && buffer_init[0] < 2.005 )))
             /* In this case we considered that we didn't received enough "one"
              * and received the required "two" to begin the transmission properly */
        {
            printf ("Unable to start the reception properly.\n");
            printf ("Not enough \"one\" received.\n");
            printf ("A packet will be lost.\n\n");
            recv_one = 0;
            printf(" the buffer is : %f\n", buffer_init[0]);
            /* Here we make the program wait until a next buffer will be admissible */
            err = Pa_ReadStream( stream, buffer, FRAMES_PER_BUFFER );
            if( err ) goto xrun;
        } 
        else if ((recv_one >= 3) && (!(1.995 < buffer_init[0] && buffer_init[0] < 2.005)))
        {
            /* In this case we considered that we have received enough "one"
              * but not received the required "two" to begin the transmission properly 
              * We can read the first point of the transmitted signal so the transmission
              * must be aborted.
              */
            printf ("Unable to start the reception properly.\n");
            printf ("The \"two\" have not been received.\n");
            printf ("A packet will be lost.\n\n");
            recv_one = 0;
            /* Here we make the program wait until a next buffer will be admissible */
            err = Pa_ReadStream( stream, buffer, FRAMES_PER_BUFFER );
            if( err ) goto xrun;
        }
    }

    /* -- Here's the loop where we get data from input -- */
    // You may get underruns or overruns if the output is not primed by PortAudio.
    err = Pa_ReadStream( stream, buffer, FRAMES_PER_BUFFER );
    if( err ) goto xrun;

    printf("Wire off.\n"); fflush(stdout);
    
    /* -- Now we stop the stream -- */
    err = Pa_StopStream( stream );
    if( err != paNoError ) goto error;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto error;

    Pa_Terminate();
    
    /* ------ FFTW PHASE -----*/
    /* Process out */
    int nc = ( FRAMES_PER_BUFFER / 2 ) + 1;
    out = fftw_malloc(sizeof(fftw_complex) * nc);

    /* We replace the float array "in" into a double array "in_d"
     * It is needed for the argument of the fftw function
     * during the demodulation phase
     */
    double in_d[FRAMES_PER_BUFFER];
    for(int i = 0; i < FRAMES_PER_BUFFER; i++)
        in_d[i] = round_to_double(buffer[i]);
    
    /* Defined the parameters of the plan p */
    p = fftw_plan_dft_r2c_1d (FRAMES_PER_BUFFER, in_d, out, FFTW_ESTIMATE);

    /* Excute the fft on the plan p for the given parameters */
    fftw_execute(p);

    printarray_dec(decode_message(out), SL - 1);

    //printf("\n CheckSUM = %d", check_checksum(back_to_binary(decode_message(out))));
    
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
        return -2;
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