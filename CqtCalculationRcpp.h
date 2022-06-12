#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double pi = 3.141592653589793238463 ;

class CqtCalculationRcpp {
        public:
             const int SAMPLING_FREQUENCY = 2048;
             const int hop_length = 64;
             const bool center = true;
             const int bins_per_octave=12;
             const int fmin = 20;
             const int fmax = 1024;
             const int filter_scale=1;
             const int norm = 1;
             const String window = "hann";

            vec hanningWindow(int n);
            cx_vec applyWindowToSignal(cx_vec &signal);
            cx_vec signalProductWindowFft(cx_vec &signal);
            cx_vec calcProduct(cx_vec &v1,cx_vec &v2);
            cx_vec calcProduct2(cx_vec &v1,cx_vec &v2);
            cx_vec calcDiv(cx_vec &v1,vec &v2);
            vec energyOfSpectrum(cx_vec &signal);
            vec whiten(cx_vec &signal);
            double calcQ();
            double calcNBins();
            vec calcFreqs();
            double calcAlpha();
            vec calcLengths();
            int calcFFTLen();
            cx_mat calcKernels();
            cx_mat calcKernelWindows();
            cx_mat calcKernelSignals() ;
        private:
          double min, max;
};




//vec CqtCalculationRcpp::hanningWindow(int N)
//{
//  if (N == 1) {
//     return  {1};
//  } else {
//  	  vec w(N);
//  	if(N%2==0){
//  	  int half = N/2;
//  	  w.subvec(0,half - 1) = 0.54-0.46* cos( (regspace<vec>(0,1, half - 1 ) + 1 ) * 2 *  datum::pi / (N + 1) ) ;
//  	  w.subvec(half ,N - 1 ) = reverse(w.subvec(0,half - 1));
//  	} else {
//  	  int half = (N+1)/2;
//  	  w.subvec(0,half - 1) = 0.54-0.46* cos( (regspace<vec>(0,1, half - 1 ) + 1 ) * 2 *  datum::pi / (N + 1) ) ;
//  	  w.subvec(half ,N - 1 ) = reverse(w.subvec(0,half - 2));
//  	}
//
//    for(int i=N-1; i>=1; i--)
//        w[i] = w[i-1];
//    w[0] = 0.0;
//
//     return w;
//  }
//}


vec CqtCalculationRcpp::hanningWindow(int N)
{
  if (N == 1) {
     return  {1};
  } else {
     return 0.5-0.5*cos(2.0* datum::pi *  regspace<vec>(0,1, N -1 )  /N  );
  }
}




//def hamming(M, sym=True):
//    """The M-point Hamming window.
//    """
//    if M < 1:
//        return np.array([])
//    if M == 1:
//        return np.ones(1,'d')
//    odd = M % 2
//    if not sym and not odd:
//        M = M+1
//    n = np.arange(0,M)
//    w = 0.54-0.46*np.cos(2.0*np.pi*n/(M-1))
//    if not sym and not odd:
//        w = w[:-1]
//    return w


cx_vec CqtCalculationRcpp::applyWindowToSignal(cx_vec &signal) {
     cx_vec wSignal( signal.size() ) ;
     vec window = hanningWindow( signal.size());

    for (unsigned int i = 0 ; i < signal.size() ; i++) {
       wSignal[i] = window[i] * signal[i] ;
    }

    return wSignal;
}

//
cx_vec CqtCalculationRcpp::signalProductWindowFft(cx_vec &signal) {
     cx_vec windowSignal = applyWindowToSignal(signal);
    return fft(windowSignal);
}


cx_vec CqtCalculationRcpp::calcProduct(cx_vec &v1,cx_vec &v2) {
     cx_vec product(v1.size());
     for (unsigned int i = 0 ; i < v1.size() ; i++) {
        product[i] = v1[i] * v2[i];
     }
     return product;
}

cx_vec CqtCalculationRcpp::calcProduct2(cx_vec &v1,cx_vec &v2) {
  return v1 % v2;
}


cx_vec CqtCalculationRcpp::calcDiv(cx_vec &v1,vec &v2) {
     cx_vec div(v1.size());
     for (unsigned int i = 0 ; i < v1.size() ; i++) {
        div[i] = v1[i] / v2[i];
     }
     return div;
}

vec CqtCalculationRcpp::energyOfSpectrum(cx_vec &signal) {
     cx_vec windowSignalFft     = signalProductWindowFft(signal);
     cx_vec windowSignalFftConj = conj(windowSignalFft);
     cx_vec product(windowSignalFft.size());
     for (unsigned int i = 0 ; i < windowSignalFft.size() ; i++) {
        product[i] = windowSignalFft[i] * windowSignalFftConj[i] ;
     }

    return sqrt(real(product));
}


vec CqtCalculationRcpp::whiten(cx_vec &signal) {
     cx_vec spectrum     = signalProductWindowFft(signal);
     cx_vec spectrumConj = conj(spectrum);
     cx_vec product = calcProduct(spectrum,spectrumConj);
     vec mag = sqrt(real(product));
     return real(ifft(calcDiv(spectrum,mag))) * sqrt(signal.size()) * 0.5;
}


double CqtCalculationRcpp::calcQ() {
   return  ((double) filter_scale) /  ( pow(2.0, ( 1.0 / bins_per_octave ) ) - 1.0  );
}


double CqtCalculationRcpp::calcNBins() {
   return  ceil( bins_per_octave * log2(fmax / fmin) );
}


vec CqtCalculationRcpp::calcFreqs() {
//        return  self.fmin * 2.0 ** (np.r_[0:self.calcNBins()] / float(self.bins_per_octave))

     double nBins =  calcNBins();
     vec  range = regspace(0, nBins - 1);
     vec  freqs(range.size());
     for( unsigned int i = 0 ; i < range.size(); i ++ ) {
        freqs[i] = fmin * pow(2,range[i]/(double)bins_per_octave);
     }
     return  freqs;
}

double CqtCalculationRcpp::calcAlpha() {
     return pow(2, 1.0 / bins_per_octave) - 1.0;
}

vec CqtCalculationRcpp::calcLengths() {
     return ceil( calcQ() * SAMPLING_FREQUENCY  * calcAlpha() / calcFreqs() );
}


int CqtCalculationRcpp::calcFFTLen() {
     return (int)pow(2,ceil( log2((int)calcLengths().max())));
}

cx_mat CqtCalculationRcpp::calcKernels() {
   int fftLen = (int)calcFFTLen();
   int nBins = (int)calcNBins();
   vec lengths = calcLengths();
   vec freqs = calcFreqs();
   cx_mat kernals = zeros<cx_mat>(nBins,fftLen);

    for( int k = 0 ; k < nBins ; k++ ){
        int start = (int)ceil((fftLen-lengths[k])/2) - ( (int)lengths[k] % 2 == 1 ? 1 : 0);
        int high = (int)lengths[k]/2 - 1 ;
        int low = (int)lengths[k]/2 + ( (int)lengths[k] % 2 == 1 ? 1 : 0);

        cx_vec signalWindow =  hanningWindow(lengths[k]) % exp( regspace<cx_vec>(-low,1, high) * 2i * datum::pi * freqs[k] / SAMPLING_FREQUENCY ) / lengths[k];
        kernals.row(k).cols(start, start+lengths[k]-1) = (as<cx_rowvec>)(wrap(signalWindow / arma::norm(signalWindow,1)));
    }
    return kernals;
}



cx_mat CqtCalculationRcpp::calcKernelWindows() {
   int fftLen = (int)calcFFTLen();
   int nBins = (int)calcNBins();
   vec lengths = calcLengths();
   vec freqs = calcFreqs();
   cx_mat kernals = zeros<cx_mat>(nBins,fftLen);

    for( int k = 0 ; k < nBins ; k++ ){
        int start = (int)ceil((fftLen-lengths[k])/2) - ( (int)lengths[k] % 2 == 1 ? 1 : 0);

        kernals.row(k).cols(start, start+lengths[k]-1) = (as<cx_rowvec>)(wrap(hanningWindow(lengths[k]) ) );
    }
    return kernals;
}



cx_mat CqtCalculationRcpp::calcKernelSignals() {
   int fftLen = (int)calcFFTLen();
   int nBins = (int)calcNBins();
   vec lengths = calcLengths();
   vec freqs = calcFreqs();
   cx_mat kernals = zeros<cx_mat>(nBins,fftLen);

    for( int k = 0 ; k < nBins ; k++ ){
        int start = (int)ceil((fftLen-lengths[k])/2) - ( (int)lengths[k] % 2 == 1 ? 1 : 0);
        int high = (int)lengths[k]/2 - 1 ;
        int low = (int)lengths[k]/2 + ( (int)lengths[k] % 2 == 1 ? 1 : 0);

        cx_vec signalWindow =  exp( regspace<cx_vec>(-low,1, high) * 2i * datum::pi * freqs[k] / SAMPLING_FREQUENCY ) / lengths[k];
        kernals.row(k).cols(start, start+lengths[k]-1) = (as<cx_rowvec>)(wrap(signalWindow / arma::norm(signalWindow,1)));
    }
    return kernals;
}
