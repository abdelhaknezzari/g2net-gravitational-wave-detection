# from nnAudio.Spectrogram import CQT1992v2
#
# import torch
# from scipy.signal import butter, lfilter
#
#
# import torch.nn as nn
# import numpy as np
#
# from scipy.signal import get_window
#
# from torch.nn.functional import conv1d, fold
#
#
#
# from scipy import signal
# from scipy.fftpack import fft
#
#
# from nnAudio.librosa_functions import *

from scipy.signal import get_window
import sys
from scipy.signal import butter, lfilter
from scipy import signal, fft

import torch.nn as nn
import torch
from torch.nn.functional import conv1d, fold

import numpy as np

from scipy.signal import get_window
from scipy import signal



from nnAudio.librosa_functions import *



lowcut=50
highcut=500
id= "fac5791f7b"
BATCH_SIZE = 256
EPOCHS = 1
EXAMPLE_IDENTIFIER_1 = "00000e74ad"
EXAMPLE_IDENTIFIER_0 = "00001f4945"
RANDOM_SAMPLE_SIZE = 1
PERFORM_FITTING = True
SAMPLING_FREQUENCY = 2048
SAMPLES_1 = 4096
SAMPLES_3 = 3 * SAMPLES_1
USE_TRAIN_SUBSET = False
USE_TEST_SUBSET = False
SUBSET_SIZE = 1024
LEARNING_RATE = 0.001
TRAIN_TEST_SPLIT = 0.95



def get_array(identifier):
    path = f"{identifier}.npy"
    return np.load(path)

def getWaves(path):
    return np.load(path)

def test1(index):
    return int(index)


def whiten0(path,index):
    waves = np.load(path)
    return waves[int(index)]

def whiten1(path):
    waveform = whiten0(path,0)
    return signal.hann(waveform.size)

def whiten2(path):
    waveform = whiten0(path,0)
    window = whiten1(path)
    return waveform * window

def whiten3(path):
    windowWave = whiten2(path)
    return fft.fft(windowWave)

def whiten4(path):
    windowWaveSpec = whiten3(path)
    return windowWaveSpec*np.conj(windowWaveSpec)

def whiten5(path):
    spectrum = whiten4(path)
    return np.sqrt(np.real(spectrum))

def whiten6(path):
    waveform = whiten0(path,0)
    window = signal.hann(waveform.size)
    spectrum = fft.fft(waveform * window)
    mag = np.sqrt(np.real(spectrum*np.conj(spectrum)))
    return np.real(fft.ifft(spectrum/mag)) * np.sqrt(len(waveform)/2)

def getWave(path):
    waveform = whiten6(path)
    return waveform / np.max(waveform)

def get_cqt_spectrogram_of_data2(path):
    cqt = CQT1992v2(sr=SAMPLING_FREQUENCY, hop_length=64, fmin=20, fmax=1024, bins_per_octave=12, norm=1, window='hann', center=True, pad_mode='reflect', trainable=False, output_format='Magnitude', verbose=False)
    waveform = torch.from_numpy(getWave(path)).float()
    cqt_image = cqt(waveform)
    cqt_image = np.array(cqt_image)
    cqt_image = np.transpose(cqt_image, (1,2,0))
    return cqt_image

def whiten(waveform):
    window = signal.hann(waveform.size)
    spectrum = fft.fft(waveform * window)
    mag = np.sqrt(np.real(spectrum*np.conj(spectrum)))
    return np.real(fft.ifft(spectrum/mag)) * np.sqrt(len(waveform)/2)

def get_whitened_data(id, detector):
    return whiten(get_array(id)[detector])


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def generateCQTNpy(id,lowcut, highcut):
    whitened_data0 = get_whitened_data(id, 0)
    whitened_data1 = get_whitened_data(id, 1)
    whitened_data2 = get_whitened_data(id, 2)

    bandpassed_data0 = butter_bandpass_filter(whitened_data0, lowcut, highcut, whitened_data0.size)
    bandpassed_data1 = butter_bandpass_filter(whitened_data1, lowcut, highcut, whitened_data1.size)
    bandpassed_data2 = butter_bandpass_filter(whitened_data2, lowcut, highcut, whitened_data2.size)

    image0 = get_cqt_spectrogram_of_data(whitened_data0)
    image1 = get_cqt_spectrogram_of_data(whitened_data1)
    image2 = get_cqt_spectrogram_of_data(whitened_data2)

    np.save("image.npy", np.concatenate( (image0.reshape(69*65 )  , image1.reshape(69*65 )  , image2.reshape(69*65 ) ) , axis = 0)  )


def generateCQTNpy2(id,lowcut, highcut):
    whitened_data0 = get_whitened_data(id, 0)
    whitened_data1 = get_whitened_data(id, 1)
    whitened_data2 = get_whitened_data(id, 2)
    image0 = get_cqt_spectrogram_of_data(whitened_data0)
    image1 = get_cqt_spectrogram_of_data(whitened_data1)
    image2 = get_cqt_spectrogram_of_data(whitened_data2)

    return np.concatenate( (image0.reshape(69*65 )  , image1.reshape(69*65 )  , image2.reshape(69*65 ) ))


def broadcast_dim(x):
    if x.dim() == 2:
        x = x[:, None, :]
    elif x.dim() == 1:
        x = x[None, None, :]
    elif x.dim() == 3:
        pass
    return x





class CQT2022v2(nn.Module):
    def __init__(
            self,
            sr=22050,
            hop_length=512,
            fmin=32.70,
            fmax=None,
            bins_per_octave=12,
            filter_scale=1):

        super().__init__()

        self.hop_length = hop_length
        self.center = True
        self.pad_mode = "reflect"
        self.output_format = 'Magnitude'

        Q = float(filter_scale) / (2 ** (1 / bins_per_octave) - 1)

        cqt_kernels, self.kernel_width, lenghts, freqs = create_cqt_kernels( Q, sr, fmin, bins_per_octave, 1, "hann", fmax )

        self.register_buffer("lenghts", lenghts)
        self.frequencies = freqs

        cqt_kernels_real = torch.tensor(cqt_kernels.real).unsqueeze(1)
        cqt_kernels_imag = torch.tensor(cqt_kernels.imag).unsqueeze(1)

        self.register_buffer("cqt_kernels_real", cqt_kernels_real)
        self.register_buffer("cqt_kernels_imag", cqt_kernels_imag)


    def forward(self, x):
        padding = nn.ReflectionPad1d(self.kernel_width // 2)

        x = padding(broadcast_dim(x))
        CQT_real = conv1d(x, self.cqt_kernels_real, stride=self.hop_length)
        CQT_imag = -conv1d(x, self.cqt_kernels_imag, stride=self.hop_length)

        CQT_real *= torch.sqrt(self.lenghts.view(-1, 1))
        CQT_imag *= torch.sqrt(self.lenghts.view(-1, 1))

        return torch.sqrt(CQT_real.pow(2) + CQT_imag.pow(2) + 1e-8)



    def forward_manual(self, x):
        padding = nn.ReflectionPad1d(self.kernel_width // 2)
        x = padding( broadcast_dim(x))

        CQT_real = conv1d(x, self.cqt_kernels_real, stride=self.hop_length)
        CQT_imag = conv1d(x, self.cqt_kernels_imag, stride=self.hop_length)

        CQT = torch.sqrt(CQT_real.pow(2) + CQT_imag.pow(2))
        return CQT * torch.sqrt(self.lenghts.view(-1, 1))

class simpleCQT():

    def __init__(self):
        self.hannWindows = []
        self.SAMPLING_FREQUENCY = 2048
        self.hop_length = 64
        self.center = True
        self.bins_per_octave=12
        self.fmin = 20
        self.fmax = 1024
        self.filter_scale=1
        self.norm = 1
        self.window = "hann"

    def broadcast_dim(self,x):
        if x.dim() == 2:
            x = x[:, None, :]
        elif x.dim() == 1:
            x = x[None, None, :]
        elif x.dim() == 3:
            pass
        return x

    def calcQ(self):
        return float(self.filter_scale) / (2 ** (1 / self.bins_per_octave) - 1)

    def calcNBins(self):
        return np.ceil( self.bins_per_octave * np.log2(self.fmax / self.fmin) )

    def calcFreqs(self,n_bins):
        return  self.fmin * 2.0 ** (np.r_[0:n_bins] / float(self.bins_per_octave))

    def calcAlpha(self):
            return  2.0 ** (1.0 / self.bins_per_octave) - 1.0

    def calcLengths(self,freqs,alpha,Q):
        return  np.ceil(Q * self.SAMPLING_FREQUENCY  / (freqs  / alpha))

    def calcFFTLen(self,lengths):
        max_len = int(max(lengths))
        return int(2 ** (np.ceil(np.log2(max_len))))

    def calcKernels(self,n_bins,freqs,lengths,fftLen):
        tempKernel = np.zeros((int(n_bins), int(fftLen)), dtype=np.complex64)
        for k in range(0, int(n_bins)):
            freq = freqs[k]
            l = lengths[k]

            if l % 2 == 1:
                start = int(np.ceil(fftLen / 2.0 - l / 2.0)) - 1
            else:
                start = int(np.ceil(fftLen / 2.0 - l / 2.0))

            window_dispatch = get_window(self.window, int(l), fftbins=True)
            sig = window_dispatch * np.exp(np.r_[-l // 2 : l // 2] * 1j * 2 * np.pi * freq / self.SAMPLING_FREQUENCY ) / l
            tempKernel[k, start : start + int(l)] = sig / np.linalg.norm(sig, self.norm)
        return tempKernel


    def create_cqt_kernels(self):
        Q = self.calcQ()
        n_bins = self.calcNBins()
        freqs = self.calcFreqs(n_bins)
        alpha = self.calcAlpha()
        lengths = self.calcLengths(freqs  , alpha,Q)
        fftLen = self.calcFFTLen(lengths)
        tempKernel = self.calcKernels(n_bins,freqs,lengths,fftLen)
        return tempKernel, fftLen, torch.tensor(lengths).float(), freqs


    def getHannWindow(self, size):
        if len(self.hannWindows) == 0:
            self.hannWindows = signal.hann(int(size))
        return self.hannWindows

    def whiten(self,waveform):
        window = self.getHannWindow(waveform.size)
        spectrum = fft.fft(waveform * window)
        mag = np.sqrt(np.real(spectrum*np.conj(spectrum)))
        return np.real(fft.ifft(spectrum/mag)) * np.sqrt(len(waveform)/2)

    def get_cqt_spectrogram_of_data(self,data):
        waveform = torch.from_numpy(data / np.max(data)).float()
        cqt_image = self.calcCQT(waveform)
        cqt_image = np.array(cqt_image)
        cqt_image = np.transpose(cqt_image, (1,2,0))
        return cqt_image


    def generateCQTNpyFromFile(self,path):
        waves = np.load(path)
        image0 = self.get_cqt_spectrogram_of_data(self.whiten(waves[0]))
        image1 = self.get_cqt_spectrogram_of_data(self.whiten(waves[1]))
        image2 = self.get_cqt_spectrogram_of_data(self.whiten(waves[2]))
        return np.concatenate( (image0.reshape(69*65 )  , image1.reshape(69*65 )  , image2.reshape(69*65 ) ))


    def calcCQT(self,x):


        self.cqt_kernels, self.kernel_width, self.lenghts, self.freqs = self.create_cqt_kernels( )

        self.cqt_kernels_real = torch.tensor(self.cqt_kernels.real).unsqueeze(1)
        self.cqt_kernels_imag = torch.tensor(self.cqt_kernels.imag).unsqueeze(1)

        padding = nn.ReflectionPad1d(self.kernel_width // 2)
        x = padding( self.broadcast_dim(x))
        CQT_real = conv1d(x, self.cqt_kernels_real, stride=self.hop_length)
        CQT_imag = conv1d(x, self.cqt_kernels_imag, stride=self.hop_length)
        CQT = torch.sqrt(CQT_real.pow(2) + CQT_imag.pow(2))

        return CQT * torch.sqrt(self.lenghts.view(-1, 1))

def create_cqt_kernels(
        Q,
        fs,
        fmin,
        bins_per_octave=12,
        norm=1,
        window="hann",
        fmax=None
):
    n_bins = np.ceil( bins_per_octave * np.log2(fmax / fmin))
    freqs = fmin * 2.0 ** (np.r_[0:n_bins] / float(bins_per_octave))

    alpha = 2.0 ** (1.0 / bins_per_octave) - 1.0
    lengths = np.ceil(Q * fs / (freqs  / alpha))
    max_len = int(max(lengths))
    fftLen = int(2 ** (np.ceil(np.log2(max_len))))

    tempKernel = np.zeros((int(n_bins), int(fftLen)), dtype=np.complex64)

    for k in range(0, int(n_bins)):
        freq = freqs[k]
        l = lengths[k]
        if l % 2 == 1:
            start = int(np.ceil(fftLen / 2.0 - l / 2.0)) - 1
        else:
            start = int(np.ceil(fftLen / 2.0 - l / 2.0))

        window_dispatch = get_window(window, int(l), fftbins=True)
        sig = window_dispatch * np.exp(np.r_[-l // 2 : l // 2] * 1j * 2 * np.pi * freq / fs) / l


        tempKernel[k, start : start + int(l)] = sig / np.linalg.norm(sig, norm)

    return tempKernel, fftLen, torch.tensor(lengths).float(), freqs

def get_cqt_spectrogram_of_data(data):
    cqt = CQT2022v2(sr=SAMPLING_FREQUENCY, hop_length=64, fmin=20, fmax=1024, bins_per_octave=12, filter_scale = 1)
    waveform = torch.from_numpy(data / np.max(data)).float()
    cqt_image = cqt(waveform)
    cqt_image = np.array(cqt_image)
    cqt_image = np.transpose(cqt_image, (1,2,0))
    return cqt_image


def generateCQTNpyFromFile(path):
    waves = np.load(path)
    image0 = get_cqt_spectrogram_of_data(whiten(waves[0]))
    image1 = get_cqt_spectrogram_of_data(whiten(waves[1]))
    image2 = get_cqt_spectrogram_of_data(whiten(waves[2]))
    return np.concatenate( (image0.reshape(69*65 )  , image1.reshape(69*65 )  , image2.reshape(69*65 ) ))


#generateCQTNpy(id,lowcut, highcut)





