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
cx_vec whiten(cx_vec &signal) {
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
cx_vec padding( cx_vec &signal ) {
   return cqtCalculationRcpp.padding(signal);
}

// [[Rcpp::export]]
vec paddingParam( vec &signal, int pad ) {
   return cqtCalculationRcpp.padding(signal,pad);
}

// [[Rcpp::export]]
cx_mat slideSignal( cx_vec &signal , int stride, int numberOfRows ) {
   return cqtCalculationRcpp.slideSignal(signal,stride,numberOfRows);
}

// [[Rcpp::export]]
cx_mat calcConv( cx_mat &matr, cx_vec &signal , int stride) {
   return cqtCalculationRcpp.calcConv(matr,signal,stride);
}


//////////////////////////////////

// [[Rcpp::export]]
cx_vec generateWaveWhitenNormalise( cx_vec &signal ) {
   return cqtCalculationRcpp.generateWaveWhitenNormalise(signal);
}






// [[Rcpp::export]]
vec generateWaveWhiten( vec &signal ) {
   return cqtCalculationRcpp.whiten(signal);
}


// [[Rcpp::export]]
double getMaxWhitenSignal( vec &signal ) {
   vec signalWhiten = cqtCalculationRcpp.whiten(signal);
   return signalWhiten.max();
}

// [[Rcpp::export]]
vec generateWhitenPad( vec &signal, int pad) {
  vec signalWhiten =  cqtCalculationRcpp.whiten(signal);
  return cqtCalculationRcpp.padding(signalWhiten,pad);
}


// [[Rcpp::export]]
mat calcCqt( cx_vec &signal) {
   cx_mat cqt = cqtCalculationRcpp.calcKernels();
   cx_vec signalPadded = cqtCalculationRcpp.padding(signal);
   cx_mat convMat = cqtCalculationRcpp.calcConv(cqt,signalPadded,cqtCalculationRcpp.hop_length);
   return sqrt(pow(real(convMat),2) + pow(imag(convMat),2));
}



// [[Rcpp::export]]
cx_vec signalPaddingWhitenNormalise( cx_vec &signal, int pad ) {
   cx_vec witenSignal = cqtCalculationRcpp.generateWaveWhitenNormalise(signal);
   return cqtCalculationRcpp.padding(witenSignal,pad);
}


// [[Rcpp::export]]
vec signalPadding( vec &signal, int pad ) {
   return cqtCalculationRcpp.padding(signal,pad);
}



// [[Rcpp::export]]
mat generateCQTConvWaveReal( vec &signal, int stride, int pad) {
  cx_mat cqt = cqtCalculationRcpp.calcKernels();
  mat cqtReal = real(cqt);

  vec signalWhiten =  cqtCalculationRcpp.whiten(signal);
  signalWhiten = signalWhiten / signalWhiten.max();
  vec signalPadded = cqtCalculationRcpp.padding(signalWhiten,pad);

//  cx_vec signalPadded = signalPaddingWhitenNormalise(signal,pad);


//  cx_vec witenSignal = cqtCalculationRcpp.generateWaveWhitenNormalise(signal);
//  cx_vec signalPadded = cqtCalculationRcpp.padding(witenSignal,pad);

  mat signalMat = cqtCalculationRcpp.slideSignal(signalPadded,stride,cqtReal.n_cols);
  return real(cqtReal*signalMat);
}


// [[Rcpp::export]]
mat generateCQTConvWaveImag( vec &signal, int stride, int pad) {
  cx_mat cqt = cqtCalculationRcpp.calcKernels();
  mat cqtImag = imag(cqt);

  vec signalWhiten =  cqtCalculationRcpp.whiten(signal);
  signalWhiten = signalWhiten / signalWhiten.max();

  vec signalPadded = cqtCalculationRcpp.padding(signalWhiten,pad);
  mat signalMat = cqtCalculationRcpp.slideSignal(signalPadded,stride,cqtImag.n_cols);
  return cqtImag*signalMat;
}


// [[Rcpp::export]]
mat generateCQTConvWave( vec &signal) {
  cx_mat cqt = cqtCalculationRcpp.calcKernels();
  mat cqtImag = imag(cqt);
  mat cqtReal = real(cqt);

  vec signalWhiten =  cqtCalculationRcpp.whiten(signal);
  signalWhiten = signalWhiten / signalWhiten.max();

  vec signalPadded = cqtCalculationRcpp.padding(signalWhiten);
  mat signalMat = cqtCalculationRcpp.slideSignal(signalPadded,cqtCalculationRcpp.hop_length,cqtImag.n_cols);

  mat cqtCalc = sqrt( pow(cqtReal*signalMat,2) + pow(cqtImag*signalMat,2) );
  return cqtCalc.each_col() % sqrt(calcLengths());
// return cqtCalc;
}



// [[Rcpp::export]]
mat calcMatDistance( mat &a,mat &b) {
  return sqrt(pow(a-b,2));
}


// [[Rcpp::export]]
mat matMultVec(mat &matA, vec &vecB) {
    mat result = matA.each_row() % vecB.t();
    return  result;
}




// [[Rcpp::export]]
vec calcMatDistanceVec( vec &a,vec &b) {
  return sqrt(pow(a-b,2));
}

