#include "CqtCalculationRcpp.h"

CqtCalculationRcpp cqtCalculationRcpp;

// [[Rcpp::export]]
vec hanningWindow(int n) {
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


// [[Rcpp::export]]
int calcFFTLen() {
     return cqtCalculationRcpp.calcFFTLen();
}

// [[Rcpp::export]]
cx_mat calcKernels() {
    return cqtCalculationRcpp.calcKernels();
}



// [[Rcpp::export]]
cx_mat calcKernelWindows() {
   return cqtCalculationRcpp.calcKernelWindows();
}


// [[Rcpp::export]]
cx_mat calcKernelSignals() {
   return cqtCalculationRcpp.calcKernelSignals();
}


// [[Rcpp::export]]
cx_mat calcKernels3( cx_mat signals, cx_mat windows) {
   return windows % signals;
}




