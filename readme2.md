## Deep Learning:
Detect signals corresponding to gravitational waves

### Transform waves:
The signals can not feed the neural network directly, they need first to be filtered to remove unnecessary information, like noise
### Whitening signal
minimize the effects of leakage
     When doing operational noise and vibration measurements, the Hanning window is commonly used.
     Random data has spectral leakage due to the abrupt cutoff at the beginning and end of the time block
     random signal is composed of many different frequencies
     Hanning windows are often used with random data because they have moderate impact on the frequency resolution and amplitude accuracy of the resulting frequency spectrum
     The Hanning window starts at a value of zero and ends at a value of zero
     This gradual transition between 0 and 1 ensures a smooth change in amplitudes when multiplying the measured signal by the window, which helps reduce the spectral leakage.
     The benefit is not that the captured signal is perfectly replicated. The main benefit is that the leakage is now confined over a smaller frequency range, instead of affecting the entire frequency bandwidth of the measurement.

[Window Types: Hanning, Flattop, Uniform, Tukey, and Exponential (siemens.com)](https://community.sw.siemens.com/s/article/window-types-hanning-flattop-uniform-tukey-and-exponential)
[Windows and Spectral Leakage (siemens.com)](https://community.sw.siemens.com/s/article/windows-and-spectral-leakage)

Most data recorded from a gravitational-wave interferometer carry information across a wide band of frequencies, typically up to a few kiloHertz, but often it is the case that the low-frequency amplitude dwarfs that of the high-frequency content, making discerning high-frequency features difficult.

We employ a technique called ‘whitening’ to normalize the power at all frequencies so that excess power at any frequency is more obvious.
We demonstrate below with an auxiliary signal recording transmitted power in one of the interferometer arms, which recorded two large glitches with a frequency of around 5-50Hz.
[3. Whitening a TimeSeries — GWpy 2.0.1.dev185+gb5f79110 documentation](https://gwpy.github.io/docs/latest/examples/timeseries/whiten.html)


Spectral Whitening is the process of making the Magnitude spectrum Uniform.
For an image it makes the Magnitude Spectrum more continuous rather than having few frequencies jumping around here and there. Basically the word "Whitening" comes from White Process whose spectrum is just a constant at all frequencies. But if you do that to an image it'll make no sense. So in effect you'd want a rather jumpy and jittery Spectrum to look more smooth without overly inducing noise.



### CQT: Constant Q transform
In many cases, such as that of musical signals, a constant Q transform gives a better representation of spectral
data than the computationally efficient fast Fourier transform. Various solutions to this problem using constant Q
filterbanks or a "bounded Q" Fourier transform have been
proposed (Harris, 1976; Schwede, 1983; Mont-Reynaud,
1985). The music group at Marseilles has proposed a"wavelet transform" (Kronland-Martinet, 1988). Brown (1991)
describes results applied to musical signals based on a direct
evaluation of the DFT for the desired components.
We have calculated a constant Q transform based on
transforming a fast Fourier transform into the log frequency
domain. The FFT is calculated using a standard FFT program, and the entire calculation takes only slightly longer to
run than the FFT since there are few operations involved in
the computation of the transformation. The transformation
is based upon the following. A constant Q transform can be
calculated directly (Brown, 1991 ) by evaluating

[http://academics.wellesley.edu/Physics/brown/pubs/cq1stPaper.pdf](http://academics.wellesley.edu/Physics/brown/pubs/cq1stPaper.pdf)
transform signal to frequency domain.
it is a type of wavelet transform.



Constant-Q transform (CQT) here refers to a technique
that transforms a time-domain signal x(n) into the timefrequency domain so that the center frequencies of the frequency bins are geometrically spaced and their Q-factors
are all equal. In effect, this means that the frequency resolution is better for low frequencies and the time resolution is better for high frequencies. The CQT is essentially a wavelet transform, but here the term CQT is preferred since it underlines the fact that we are considering transforms with relatively high Q-factors, equivalent
to 12–96 bins per octave. This renders many of the conventional wavelet transform techniques inadequate; for example methods based on iterated filterbanks would require
filtering the input signal hundreds of times.
The CQT is well-motivated from both musical and perceptual viewpoints. The fundamental frequencies (F0s) of
the tones in Western music are geometrically spaced: in
the standard 12-tone equal temperament, for example, the
F0s obey Fk = 440Hz × 2
k/12, where k ∈ [−50, 40] is an
[https://iem.kug.ac.at/fileadmin/media/iem/projects/2010/smc10_schoerkhuber.pdf](https://iem.kug.ac.at/fileadmin/media/iem/projects/2010/smc10_schoerkhuber.pdf)



The CQT-based time-frequency analysis provides variable spectro-temporal resolution with higher frequency resolution at lower frequencies. 
Since lower-frequency regions of speech signal contain more emotion-related information than higher-frequency regions, 
the increased low-frequency resolution of CQT makes it more promising for SER than standard short-time Fourier transform (STFT). 
We present a comparative analysis of short-term acoustic features based on STFT and CQT for SER with deep neural network (DNN) as a back-end classifier. We optimize different parameters for both features. The CQT-based features outperform the STFT-based spectral features for SER experiments. Further experiments with cross-corpora evaluation demonstrate that the CQT-based systems provide better generalization with out-of-domain training data.


CQT refers to a time-frequency representation where the frequency bins are geometrically spaced and the 
Q-factors (ratios of the center frequencies to bandwidths) of all bins are equal. 
An inverse transform is proposed which enables a reasonable-quality (around 55dB signal-to-noise ratio) reconstruction of the original signal from its CQT coefficients. 
Here CQTs with high Q-factors, equivalent to 12–96 bins per octave, are of particular interest. 
The proposed method is flexible with regard to the number of bins per octave, the applied window function, 
and the Q-factor, and is particularly suitable for the analysis of music signals. 
A reference implementation of the proposed methods is published as a Matlab toolbox. 
The toolbox includes user-interface tools that facilitate spectral data visualization 
and the indexing and working with the data structure produced by the CQT
[https://www.semanticscholar.org/paper/CONSTANT-Q-TRANSFORM-TOOLBOX-FOR-MUSIC-PROCESSING-Sch%C3%B6rkhuber/4cef10ea66e40ad03f434c70d745a4959cea96dd](https://www.semanticscholar.org/paper/CONSTANT-Q-TRANSFORM-TOOLBOX-FOR-MUSIC-PROCESSING-Sch%C3%B6rkhuber/4cef10ea66e40ad03f434c70d745a4959cea96dd)


### Why CQT:





### How to calculate
### Simplify python:

### Accelerate CQT calculation Using python:
### Use RCPP
Translate python to CPP (Armadillo)



