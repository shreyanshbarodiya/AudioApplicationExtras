import java.io.FileInputStream;
import java.util.*; 
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.Math.*;

public class FrameRead {
	static final double PI = 3.1415926536;
	
	public static int bitReverse(int x, int log2n){
		int n = 0;
		int mask = 0x1;
		for (int i=0; i < log2n; i++){
			n <<= 1;
			n |= (x & 1);
			x >>= 1;
			}
		return n;
	}
	
	public static void fft(ComplexNumber a[], ComplexNumber b[], int log2n,  Boolean inverse){
		//typedef typename iterator_traits<Iter_T>::value_type complex;
		ComplexNumber J = new ComplexNumber(0, 1);
		int n = 1 << log2n;

		for (int i=0; i < n; ++i) {
			b[bitReverse(i, log2n)] = a[i];
		}
		
	 	ComplexNumber sign;
	 	if (inverse){sign = new ComplexNumber(1,0);}
		else {sign = new ComplexNumber(-1,0);}
		
		for (int s = 1; s <= log2n; ++s){
			int m = 1 << s;
			int m2 = m >> 1;
			
			//verify
			ComplexNumber w = new ComplexNumber(1, 0);
			ComplexNumber copy = sign;
			copy.multiply(J);
			copy.multiply(new ComplexNumber(PI/m2, 0));
			ComplexNumber copy2 = new ComplexNumber(Math.exp(1.0), 0);
			copy2.exp(copy);
	
			ComplexNumber wm = copy2;

			//			ComplexNumber wm = exp(sign.multiply(J).multiply(PI / m2));
			
			for (int j=0; j < m2; ++j){
				for (int k=j; k < n; k += m){
					
					ComplexNumber t = w; 
					t.multiply(b[k+m2]);
					//complex t = w * b[k + m2];
					
					ComplexNumber u = b[k]; 
					//complex u = b[k];
					
					b[k].add(t);
					//b[k] = u + t;
					
					b[k+m2] = u;
					u.subtract(t);
					//b[k + m2] = u - t;
					
				}
				//w *= wm;
				w.multiply(wm);
			}
		}
		
		if(inverse){
			ComplexNumber divn = new ComplexNumber(n,0);
			for (int i=0; i < n; ++i) {
				b[i].divide(divn);
			}
		}
		
	}

	public static void hilbert(ComplexNumber a[], ComplexNumber b[], ComplexNumber c[], int log2n ){
		
		int n = 1 << log2n;
		fft(a,b,log2n,false);
		
		ComplexNumber mult2 = new ComplexNumber(2,0); 
		int nhalf = n>>1;
		
		//Taking only one sideband
		for(int i=1;i<nhalf;i++){
			b[i].multiply(mult2);
		}
		
		for(int i =nhalf+1; i<n; i++){
			b[i]=new ComplexNumber(0,0);
		}
		
		fft(b,c,log2n,true);
	}
	
	static void downsample(int input[], int input_size, char decimate, double output[]){
		double ft= 2.0/decimate;	// cutoff frequency at twice the decimate sampling rate for anti-aliasing filter
		int taps = 31;					// no. of taps used in FIR - keep this odd
		int middle = (taps-1)/2;

		double[] lpf = new double[taps];
		
		for(int i=0; i<taps; i++){
			if(i!=middle){
				lpf[i]= Math.sin(2*PI*ft*(i-middle))/(PI*(i-middle))*(0.54-0.46*Math.cos(PI*i/middle));
			}
				else lpf[i] = 2*ft;
			}

			int newlen = input_size/decimate;

			
			//Applying anti-aliasing filter
			for(int i=0; i<newlen; i++){
				output[i] = 0;
				int k = i*decimate;	// we only need every decimate-th sample


				if(k<taps){
					for(int j=k; j>=0 && j<=k; j--){
						output[i] = output[i]+input[j]*lpf[k-j];
						}
					}
					else {
						for (int j=k; j>k-taps && j<=k; j--){
							output[i] = output[i]+input[j]*lpf[k-j];
						}
					}				
				}		
	}
	
	public class header_file
	{
	    char chunk_id[] = new char[4];
	    int chunk_size;
	    char format[] = new char[4];
	    char subchunk1_id[] = new char[4];
	    int subchunk1_size;
	    int audio_format;
	    int num_channels;		// number of channels i.e 1-mono 2-stereo
	    int sample_rate;			// sample_rate denotes the sampling rate.
	    int byte_rate;
	    int block_align;
	    int bits_per_sample;		// number of bits per sample
	    char subchunk2_id[] = new char[4];
	    int subchunk2_size;			// subchunk2_size denotes the number of samples.
	};

	static FrameRead outer = new FrameRead();
	//TODO
	static header_file header = outer.new header_file();
	
	public static void main(String[] args) {
		/* Input window length and hop length in msec */
		int window_time_ms = 1500;		// For window length of 1500ms, set window_time_ms to 30
		int hop_time_ms = 200;			// For hop length of 500ms, set hope_time_ms to 500

		/* Enter the name of wav file for processing and output file. File should be present in same folder as this code */
		FileInputStream infile = null;
	    FileOutputStream outfile = null;

	      try {
	    	  infile = new FileInputStream("NST_Normal.wav");
	    	  outfile = new FileOutputStream("Rate.txt");
	         
/*	         int c;
	         while ((c = infile.read()) != -1) {
	        	 outfile.write(c);
	         }*/

	    	  int count = 0;						// For counting number of frames in wave file.
			
			
	    	  header_file meta = outer.new header_file();
			
			
			// header_p points to a header struct that contains the wave file metadata fields
	    	  int nb;							// variable storing number of byes returned
			
			
    		if (infile!=null)
    		{
    			/* Print header information i.e sampling rate, no. of samples, etc */

    			//fread(meta, 1, sizeof(header), infile);
    			
    			System.out.println(" No. of bits per sample: " + meta.bits_per_sample);
    			System.out.println(" Num of channels: " + meta.num_channels);
    			//System.out.println(" Size of Header file is " + (meta)<<" bytes");
    			System.out.println(" Sampling rate of the input wave file is " + meta.sample_rate + " Hz");
    			System.out.println(" Number of samples in wave file are " + meta.subchunk2_size + " samples");
    			
    			/* ----------------- Printing header info end -------------------------------------- */
    		
    			/* Convert time specifications to samples */
    			int WINSIZE = (int) Math.ceil((meta.sample_rate)*(window_time_ms)*0.001);	// Round WINSIZE(samples) if WINSIZE is not an integer
    			int HOPSIZE = (int) Math.ceil((meta.sample_rate)*(hop_time_ms)*0.001);	// Round HOPSIZE(samples) if HOPSIZE is not an integer
    			int buff16[] = new int[WINSIZE];			// Defining buffer for holding samples that fall within the window.
    			
    			char decimate = 10;		//decimation by a factor of 10
    			int WINLEN = WINSIZE/decimate;	//parameters in decimated samples
    			int HOPLEN = HOPSIZE/decimate;
    			int samp_rate = meta.sample_rate/decimate;	//sampling rate after downsampling
    			double databuff[] = new double[WINSIZE];				// buffer to store downsampled signal		
    			
    			//short int sample_start,numb_samples;	// Variables for reading samples that fall within window for 1st frame		
    			float hop_time = 0;			// Hop time for printing in output file
    				
    			for (int j = 0;j<WINSIZE ; j++)
    				{
    					buff16[j] = 0;
    				}	
    			
    			System.out.println(" Buffer size: " + WINSIZE);
    			System.out.println(" Hop size: " + HOPSIZE);

    			//parameters for ACF
    			
    			int acf_jump = 10;		//difference between lags for which acf is calculated
    			int min_fhr = 60;	// Minimum and maximum FHR in beats per minute (bpm)
    			int max_fhr = 300;
    			int lag_min = (int) Math.floor(samp_rate*60.0/max_fhr);
    			int lag_max = (int)Math.floor(samp_rate*60.0/min_fhr);
    	 		int acf_len = (int) Math.floor((lag_max-lag_min)/acf_jump);	//length of ACF valarray
    			//vector<double> fhr;
    	 		List<Double> fhr=new ArrayList<>(); 

    			
    			// Designing filter to remove noise with frequency cut-off = 50Hz
    			int numtaps = 121;	//Number of FIR taps for Hamming window
    			int middle = (numtaps-1)/2;
    			double fir[] = new double[numtaps];	//FIR LPF using Hamming window with cutoff freq = 50Hz
    			double fc = 50.0/samp_rate;
    			
    			for(int i=0; i<numtaps; i++){
    				if(i!=middle){
    					fir[i]=Math.sin(2*PI*fc*(i-middle))/(PI*(i-middle))*(0.54-0.46*Math.cos(PI*i/middle));
    				}
    				else fir[i] = 2*fc;
    				System.out.println(i + " " + fir[i]);
    			} 
    			
    			
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
    			 
    			 	//typedef complex<double> cx;
    				int log2n =(int) Math.ceil(Math.log(WINLEN)/Math.log(WINLEN));
    				int n = 1<<log2n;
    				ComplexNumber a[] = new ComplexNumber[n];
    				ComplexNumber b[] = new ComplexNumber[n];
    				ComplexNumber c[] = new ComplexNumber[n];
    				
    				for (int i=0; i<n; i++){
    					if(i<WINLEN){
    						a[i] = new ComplexNumber(databuff[i],0);
    					}
    					else a[i] = new ComplexNumber(0,0);
    				}

    					
    				hilbert (a,b,c,log2n);
    				double env[] = new double[WINLEN];
    				for (int i=0; i<WINLEN; i++){
    					env[i] = c[i].mod();
    				} 
    				
    				double data[] = new double[WINLEN];
    				
    				
    				//Filtering envelope using the previously defined FIR LPF
    				for(int i=0; i<WINLEN; i++){
    					data[i] = 0;

    					if(i<numtaps){
    						for(int j=i; j>=0 && j<=i; j--){
    							data[i] = data[i]+env[j]*fir[i-j];
    						}
    					}
    					else {
    						for (int j=i; j>i-numtaps && j<=i; j--){
    							data[i] = data[i]+env[j]*fir[i-j];
    						}
    					}
    					
    				} 
    							 	
    	 			// Calculating norm
    	 			double norm = 0;		
    				for(int i=0; i<WINLEN; i++){
    					norm = norm + data[i]*data[i];
    				}
    					
    				//Computing autocorrelation
    				double acf[] = new double[acf_len]; 
    				double acfval;		
    				
    				for(int i=lag_min; i<lag_max; i=i+acf_jump){
    						acfval = 0;
    						for(int j = 0; j<WINLEN-i;j++) {
    							acfval = acfval+data[j]*data[i+j];
    						}
    					acf[(int) ((i-lag_min)/acf_jump)] = (double) (acfval/norm);		//assign this correctly
    					}

    				//finging position of max ACF value
    				int maxpos = 0;
    				for(int i=1; i<acf_len; i++){
    					if(acf[i]>acf[maxpos]) {
    						maxpos = i;
    					}
    				}

//    				FILE * acffile = fopen("ACF.txt","wb");		// Create output file in write mode
    		    	FileOutputStream acffile = new FileOutputStream("ACF.txt");
    				
    		    	for (int i=0; i<acf_len; i++){
    					fprintf (acffile,"%d	%.05f \n", lag_min+i*acf_jump, acf[i]);
    				}

    				int lagmax = lag_min + maxpos*acf_jump;
    				double curr_rate = samp_rate*60.0/lagmax;
    				
    				hop_time = hop_time + (float) (hop_time_ms*0.001);

    				System.out.println("Corr lag: " + lagmax + " at hop time: " + hop_time);
    				
    				if(fhr.isEmpty() || curr_rate - fhr.get(fhr.size()-1) <=50){
    					fhr.add(curr_rate);
    					
        				System.out.println("Rate: " + curr_rate );    					
    						
    					fprintf (outfile,"%.03f	%.02f \n", hop_time, curr_rate);	// Ouput file format: Time-Stamp(sec)   FHR Value
    																				// Example          : 0.01		1234
    					}
    				 //------printing to file end --------------------------------- 
    								
    			
    				count++;
    			}
				System.out.println(" Number of frames in the input wave file are " + count );
    		
    		}
    		
    		infile.close();
    		outfile.close();
    	  
	      } catch (Exception e) {
			e.printStackTrace();
	      }
		
				

		
		
		
	}
}
