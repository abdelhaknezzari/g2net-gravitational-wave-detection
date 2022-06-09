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

            cx_vec hanningWindow(int n);
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

        private:
          double min, max;
};




cx_vec CqtCalculationRcpp::hanningWindow(int n)
{
  cx_vec w(n) ;
  if (n == 1) {
    w[0] = 1;
  } else {
    for (int i = 0 ; i < n ; i++) {
       w[i] = 0.5 - 0.5 * cos( (2 * pi * i) / n);
    }
  }
  return w;
}

cx_vec CqtCalculationRcpp::applyWindowToSignal(cx_vec &signal) {
     cx_vec wSignal( signal.size() ) ;
     cx_vec window = hanningWindow( signal.size());

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
     double nBins =  calcNBins();
     vec  range = regspace(0, nBins);
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
     vec freqs = calcFreqs() / calcAlpha();
     return ceil(SAMPLING_FREQUENCY * calcQ() / freqs);
}





//    def calcAlpha(self):
//            return  2.0 ** (1.0 / self.bins_per_octave) - 1.0
//
//    def calcLengths(self,freqs,alpha,Q):
//        return  np.ceil(Q * self.SAMPLING_FREQUENCY  / (freqs  / alpha))





//    def calcNBins(self):
//        return np.ceil( self.bins_per_octave * np.log2(self.fmax / self.fmin) )