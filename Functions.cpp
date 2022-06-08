#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double pi = 3.141592653589793238463 ;

// [[Rcpp::export]]
cx_vec hanningWindow(int n) {
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



// [[Rcpp::export]]
cx_vec applyWindowToSignal(cx_vec &signal) {
     cx_vec wSignal( signal.size() ) ;
     cx_vec window = hanningWindow( signal.size());

    for (unsigned int i = 0 ; i < signal.size() ; i++) {
       wSignal[i] = window[i] * signal[i] ;
    }

    return wSignal;
}


// [[Rcpp::export]]
cx_vec signalProductWindowFft(cx_vec &signal) {
     cx_vec windowSignal = applyWindowToSignal(signal);
    return fft(windowSignal);
}

// [[Rcpp::export]]
vec energyOfSpectrum(cx_vec &signal) {
     cx_vec windowSignalFft     = signalProductWindowFft(signal);
     cx_vec windowSignalFftConj = conj(windowSignalFft);
     cx_vec product(windowSignalFft.size());
     for (unsigned int i = 0 ; i < windowSignalFft.size() ; i++) {
        product[i] = windowSignalFft[i] * windowSignalFftConj[i] ;
     }

    return sqrt(real(product));
}


// [[Rcpp::export]]
vec whiten(cx_vec &signal) {
     cx_vec spectrum = signalProductWindowFft(signal);
     vec mag = sqrt(real(dot(spectrum,conj(spectrum) ) ) );
     return mag;

//     cx_vec windowSignalFft     = signalProductWindowFft(signal);
//     cx_vec windowSignalFftConj = conj(windowSignalFft);
//     cx_vec product(windowSignalFft.size());
//     for (unsigned int i = 0 ; i < windowSignalFft.size() ; i++) {
//        product[i] = windowSignalFft[i] * windowSignalFftConj[i] ;
//     }
//
//    return sqrt(real(product));
}


//    def whiten(self,path,index):
//
//
//        spectrum = fft.fft(waveform * window)
//        mag = np.sqrt(np.real(spectrum*np.conj(spectrum)))
//        return np.real(fft.ifft(spectrum/mag)) * np.sqrt(len(waveform)/2)



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
