#include "CqtCalculationRcpp.h"

CqtCalculationRcpp cqtCalculationRcpp;

// [[Rcpp::export]]
cx_vec hanningWindow(int n) {
   return cqtCalculationRcpp.hanningWindow(n);
}

// [[Rcpp::export]]
cx_vec applyWindowToSignal(cx_vec &signal) {
    return cqtCalculationRcpp.applyWindowToSignal(signal);
}

// [[Rcpp::export]]
cx_vec signalProductWindowFft(cx_vec &signal) {
    return cqtCalculationRcpp.signalProductWindowFft(signal);
}

// [[Rcpp::export]]
cx_vec calcProduct(cx_vec &v1,cx_vec &v2) {
    return cqtCalculationRcpp.calcProduct(v1,v2);
}


// [[Rcpp::export]]
cx_vec calcProduct2(cx_vec &v1,cx_vec &v2) {
    return cqtCalculationRcpp.calcProduct2(v1,v2);
}

// [[Rcpp::export]]
cx_vec calcDiv(cx_vec &v1,vec &v2) {
    return cqtCalculationRcpp.calcDiv(v1,v2);
}

// [[Rcpp::export]]
vec energyOfSpectrum(cx_vec &signal) {
    return cqtCalculationRcpp.energyOfSpectrum(signal);
}

// [[Rcpp::export]]
vec whiten(cx_vec &signal) {
     return cqtCalculationRcpp.whiten(signal);
}

// [[Rcpp::export]]
double calcQ() {
     return cqtCalculationRcpp.calcQ();
}

// [[Rcpp::export]]
double calcNBins() {
     return cqtCalculationRcpp.calcNBins();
}

// [[Rcpp::export]]
vec calcFreqs() {
     return cqtCalculationRcpp.calcFreqs();
}


// [[Rcpp::export]]
double calcAlpha() {
     return cqtCalculationRcpp.calcAlpha();
}


// [[Rcpp::export]]
vec calcLengths() {
     return cqtCalculationRcpp.calcLengths();
}





//    def calcAlpha(self):
//            return  2.0 ** (1.0 / self.bins_per_octave) - 1.0
//
//    def calcLengths(self,freqs,alpha,Q):
//        return  np.ceil(Q * self.SAMPLING_FREQUENCY  / (freqs  / alpha))


// [[Rcpp::export]]
cx_vec  waveSpectrum(cx_vec &waveform) {
  cx_vec spec = fft(conv(hanningWindow(waveform.size()),  waveform));
  return spec.t() * conj( spec) ;
//  double mag = sqrt( real(  spec.t() * conj( spec)  ) );
//
//  return fft(conv(hanningWindow(waveform.size()),  waveform));
}


//    spectrum = fft.fft(waveform * window)
//    mag = np.sqrt(np.real(spectrum*np.conj(spectrum)))
//    return np.real(fft.ifft(spectrum/mag)) * np.sqrt(len(waveform)/2)
