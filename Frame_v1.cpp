#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <complex>
#include <iterator>
using namespace std;

//FFT implementation - works for DFT of lenghts of powers of 2
// bit inverse is set to True to get inverse FFT
unsigned int bitReverse(unsigned int x, int log2n){
	int n = 0;
	int mask = 0x1;
	for (int i=0; i < log2n; i++){
		n <<= 1;
		n |= (x & 1);
		x >>= 1;
		}
	return n;
}

const double PI = 3.1415926536;


template<class Iter_T>
void fft(Iter_T a, Iter_T b, int log2n,  bool inverse){
	typedef typename iterator_traits<Iter_T>::value_type complex;
	const complex J(0, 1);
	int n = 1 << log2n;

	for (unsigned int i=0; i < n; ++i) {
		b[bitReverse(i, log2n)] = a[i];
	}
 	complex sign;
 	if (inverse){sign = complex(1,0);}
	else {sign = complex(-1,0);}
	
	for (int s = 1; s <= log2n; ++s){
		int m = 1 << s;
		int m2 = m >> 1;
		complex w(1, 0);
		complex wm = exp(sign*J * (PI / m2));
		for (int j=0; j < m2; ++j){
			for (int k=j; k < n; k += m){
				complex t = w * b[k + m2];
				complex u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
				}
			w *= wm;
		}
	}
	
	if(inverse){
		complex divn(n,0);
		for (unsigned int i=0; i < n; ++i) {
			b[i]=b[i]/divn;
		}
	}	
}

//hilber transform - computes analytic signal - result found in c with signal kept in a
template<class Iter_T>
void hilbert (Iter_T a, Iter_T b, Iter_T c, int log2n){

	typedef typename iterator_traits<Iter_T>::value_type complex;

	int n = 1 << log2n;
	fft(a,b,log2n,false);
	
	complex mult2(2,0); 
	int nhalf = n>>1;
	
	//Taking only one sideband
	for(unsigned int i=1;i<nhalf;i++){
		b[i]=b[i]*mult2;
	}
	
	for(unsigned int i =nhalf+1; i<n; i++){
		b[i]=complex(0,0);
	}
	
	fft(b,c,log2n,true);
		
}


//Downsampler
// input - input array of samples, input_size - no. of input samples, decimate - integer factor
// by which you wish to downsample the input,  output - array of output samples
void downsample(const short int input[], int input_size, char decimate, double output[]){
	double ft= 2.0/decimate;	// cutoff frequency at twice the decimate sampling rate for anti-aliasing filter
	int taps = 31;					// no. of taps used in FIR - keep this odd
	int middle = (taps-1)/2;

	double lpf[taps];
	
	for(int i=0; i<taps; i++){
		if(i!=middle){
			lpf[i]=sin(2*PI*ft*(i-middle))/(PI*(i-middle))*(0.54-0.46*cos(PI*i/middle));
			//cout<<i<<" "<<lpf[i]<<endl;

		}
			else lpf[i] = 2*ft;
			//cout<<i<<" "<<lpf[i]<<endl;
		}

		int newlen = input_size/decimate;
		//cout<<"NEWLEN IS "<<newlen<<endl;
		//char ch; cin>>ch;
		
		//Applying anti-aliasing filter
		for(unsigned int i=0; i<newlen; i++){
			output[i] = 0;
			int k = i*decimate;	// we only need every decimate-th sample


			if(k<taps){
				for(unsigned int j=k; j>=0 && j<=k; j--){
					output[i] = output[i]+input[j]*lpf[k-j];
					}
				}
				else {
					for (unsigned int j=k; j>k-taps && j<=k; j--){
						output[i] = output[i]+input[j]*lpf[k-j];
						//cout<<"IF I: "<<i<<" K: "<<k<<" op "<<output[i]<<endl;

					}
				}				
			}		
}


// WAVE PCM soundfile format (you can find more in https://ccrma.stanford.edu/courses/422/projects/WaveFormat/ )
typedef struct header_file
{
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short int audio_format;
    short int num_channels;		// number of channels i.e 1-mono 2-stereo
    int sample_rate;			// sample_rate denotes the sampling rate.
    int byte_rate;
    short int block_align;
    short int bits_per_sample;		// number of bits per sample
    char subchunk2_id[4];
    int subchunk2_size;			// subchunk2_size denotes the number of samples.
} header;

typedef struct header_file* header_p;



int main()

{
	/* Input window length and hop length in msec */
	int window_time_ms = 1500;		// For window length of 1500ms, set window_time_ms to 30
	int hop_time_ms = 200;			// For hop length of 500ms, set hope_time_ms to 500

	/* Enter the name of wav file for processing and output file. File should be present in same folder as this code */
	FILE * infile = fopen("NST_Normal.wav","rb");		// Open wave file in read mode
	//FILE * infile = fopen("fhr-03_mono.wav","rb");		// Open wave file in read mode
	FILE * outfile = fopen("Rate.txt","wb");		// Create output file in write mode

	//FILE * outfile1 = fopen("Input.txt","wb");		// Create output file in write mode
	//FILE * outfile2 = fopen("Down.txt","wb");		// Create output file in write mode 
	//FILE * outfile3 = fopen("Env.txt", "wb");

	//cout << "Hop:" << infile << endl;	
	
	int count = 0;						// For counting number of frames in wave file.
	header_p meta = (header_p)malloc(sizeof(header));	// header_p points to a header struct that contains the wave file metadata fields
	int nb;							// variable storing number of byes returned
	if (infile)
	{
		/* Print header information i.e sampling rate, no. of samples, etc */

		fread(meta, 1, sizeof(header), infile);
		//fwrite(meta,1, sizeof(*meta), outfile);
		
		cout << " No. of bits per sample: "<<meta->bits_per_sample << endl;
		cout << " Num of channels: " << meta->num_channels << endl; 
		cout << " Size of Header file is "<<sizeof(*meta)<<" bytes" << endl;
		cout << " Sampling rate of the input wave file is "<< meta->sample_rate <<" Hz" << endl;
		cout << " Number of samples in wave file are " << meta->subchunk2_size << " samples" << endl;
		/* ----------------- Printing header info end -------------------------------------- */
	
		/* Convert time specifications to samples */
		int WINSIZE = ceil((meta->sample_rate)*(window_time_ms)*0.001);	// Round WINSIZE(samples) if WINSIZE is not an integer
		int HOPSIZE = ceil((meta->sample_rate)*(hop_time_ms)*0.001);	// Round HOPSIZE(samples) if HOPSIZE is not an integer
		short int buff16[WINSIZE];			// Defining buffer for holding samples that fall within the window.
		
		char decimate = 10;		//decimation by a factor of 10
		int WINLEN = WINSIZE/decimate;	//parameters in decimated samples
		int HOPLEN = HOPSIZE/decimate;
		int samp_rate = meta->sample_rate/decimate;	//sampling rate after downsampling
		double databuff[WINLEN];				// buffer to store downsampled signal		
		
		//short int sample_start,numb_samples;	// Variables for reading samples that fall within window for 1st frame		
		float hop_time = 0;			// Hop time for printing in output file
			
		for (int j = 0;j<WINSIZE ; j++)
			{
				buff16[j] = 0;
			}	
		cout << "Buffer size: "<< WINSIZE << endl;
		cout << "Hop size: " << HOPSIZE << endl;

		//parameters for ACF
		
		short int acf_jump = 10;		//difference between lags for which acf is calculated
		short int min_fhr = 60;	// Minimum and maximum FHR in beats per minute (bpm)
		short int max_fhr = 300;
		int lag_min = floor(samp_rate*60.0/max_fhr);
		int lag_max = floor(samp_rate*60.0/min_fhr);
 		int acf_len = floor((lag_max-lag_min)/acf_jump);	//length of ACF valarray
		vector<double> fhr;
		char d2;

		
		// Designing filter to remove noise with frequency cut-off = 50Hz
		int numtaps = 121;	//Number of FIR taps for Hamming window
		int middle = (numtaps-1)/2;
		double fir[numtaps];	//FIR LPF using Hamming window with cutoff freq = 50Hz
		double fc = 50.0/samp_rate;
		
		for(int i=0; i<numtaps; i++){
			if(i!=middle){
				fir[i]=sin(2*PI*fc*(i-middle))/(PI*(i-middle))*(0.54-0.46*cos(PI*i/middle));
			}
			else fir[i] = 2*fc;
			cout<<i<<" "<<fir[i]<<endl;
		} cin>>d2;
		
		
		char dummy;
		
		/* --------------- Converting time specifications to sample end-------------------------------------- */
		
		/* Start reading file frame by frame */
		while(hop_time<50)
		{
			if (count == 0){ 
				nb = fread(&buff16[0],2,WINSIZE,infile);// For 1st frame, read only samples that fall within window and reamining samples are 0.
				}
			else{ 
				for (int i=0;i<((WINLEN)-(HOPLEN));i++)                  // Shifting data in buffer so as to read new data equivalent only to Hop size 
					{
						buff16[i] = buff16[i+(HOPLEN)]; 
					}
				nb = fread(&buff16[(WINSIZE)-(HOPSIZE)],2,(HOPSIZE),infile); // Reading new data equivalent to HOPSIZE from 2nd frame onwards
				
			    }
			
			downsample(buff16, WINSIZE, decimate, databuff);
 			
			 //calculating envelope of signal
		 
		 	typedef complex<double> cx;
			int log2n = ceil(log2(WINLEN));
			int n = 1<<log2n;
			cx a[n];
			cx b[n];
			cx c[n];
			
			for (unsigned int i=0; i<n; i++){
				if(i<WINLEN){
					a[i] = cx(databuff[i],0);
				}
				else a[i] = cx(0,0);
			}

				
			hilbert (a,b,c,log2n);
			double env[WINLEN];
			for (unsigned int i=0; i<WINLEN; i++){
				env[i] = abs(c[i]);
				//cout<<"Env: "<<i<<" "<<env[i]<<endl;
			} 
			
			double data[WINLEN];
			
			
			//Filtering envelope using the previously defined FIR LPF
			for(unsigned int i=0; i<WINLEN; i++){
				data[i] = 0;

				if(i<numtaps){
					for(unsigned int j=i; j>=0 && j<=i; j--){
						data[i] = data[i]+env[j]*fir[i-j];
						//cout<<"Filt I: "<<i<<" "<<data[i]<<endl;
					}
				}
				else {
					for (unsigned int j=i; j>i-numtaps && j<=i; j--){
						data[i] = data[i]+env[j]*fir[i-j];
						//cout<<"Filt I: "<<i<<" "<<data[i]<<endl;
						//if(i==124) {cout<<"Data i: "<<data[i]<<" "<< "Env: "<<env[j]<<" at j = "<<j<<" "<<"FIR: "<<fir[i-j]<<endl;}


					}
				}
				
			} 
						 	
 			// Calculating norm
 			double norm = 0;		
			for(int i=0; i<WINLEN; i++){
				norm = norm + data[i]*data[i];
			}
				
			//Computing autocorrelation
			double acf[acf_len]; 
			double acfval;		
			
			for(int i=lag_min; i<lag_max; i=i+acf_jump){
					acfval = 0;
					for(int j = 0; j<WINLEN-i;j++) {
						acfval = acfval+data[j]*data[i+j];
					}
				acf[int((i-lag_min)/acf_jump)] = double(acfval/norm);		//assign this correctly
				}

			//finging position of max ACF value
			int maxpos = 0;
			for(int i=1; i<acf_len; i++){
				if(acf[i]>acf[maxpos]) {
					maxpos = i;
				}
			}

			FILE * acffile = fopen("ACF.txt","wb");		// Create output file in write mode
			
			for (int i=0; i<acf_len; i++){
				fprintf (acffile,"%d	%.05f \n", lag_min+i*acf_jump, acf[i]);
			}

			int lagmax = lag_min + maxpos*acf_jump;
			double curr_rate = samp_rate*60.0/lagmax;
			
			hop_time = hop_time + float(hop_time_ms*0.001);

			cout<<"Corr lag: "<<lagmax<<" at hop time: "<<hop_time<<endl;

			if(fhr.empty() || curr_rate-fhr.back()<=50){
				fhr.push_back(curr_rate);
				cout<<"Rate: "<<curr_rate<<endl;
					
				//Print time and FHR value to a file 	
				fprintf (outfile,"%.03f	%.02f \n", hop_time, curr_rate);	// Ouput file format: Time-Stamp(sec)   FHR Value
																			// Example          : 0.01		1234
				}
			 //------printing to file end --------------------------------- 
							
		
			count++;
		}
	
	cout << " Number of frames in the input wave file are " <<count << endl;  // Print total number of frames
	}



return 0;
}


