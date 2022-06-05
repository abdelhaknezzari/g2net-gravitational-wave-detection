#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double pi = 3.141592653589793238463 ;

// [[Rcpp::export]]
cx_vec hanningWindow(int n) {
  cx_vec w( n == 1 ? 1 : n  ) ;
  int N = n ;

  if (n == 1) {
    w[0] = 1;
  } else {
    n = n ;
    for (int i = 0 ; i < n ; i++) {
       w[i] = 0.5 - 0.5 * cos( (2 * pi * i) / N);
    }
  }
  return w;
}


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
