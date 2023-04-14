# Python SPRiNT
# Author: Luc Wilson (2023)

from math import floor
from statistics import median
from itertools import compress
import numpy as np
import os
import csv
import matplotlib
import fooof
from fooof.data import FOOOFResults
from fooof.sim.gen import gen_periodic
from fooof.sim.gen import gen_aperiodic
from fooof.objs.utils import combine_fooofs
from copy import deepcopy

# F should be a nxm numpy array
# where n is the number of channels, m is the number of samples
# opt should be a dictionary

def SPRiNT_stft_py(F, opt):
    ''' SPRiNT_stft_py: Compute a locally averaged short-time Fourier transform
    (for use in SPRiNT)

    Inputs
    F - Time series (nxm numpy array),
    where n is number of channels, m is number of samples
    opt - Model settings/hyperparameters (dict)

    Segments of this function were adapted from the Brainstorm software package:
    https://neuroimage.usc.edu/brainstorm
    Tadel et al. (2011)

    Author: Luc Wilson (2023)
    '''
    # Get sampling frequency
    n_chan = F.shape
    # print(n_chan)
    if not len(n_chan) > 1:
        n_sample = n_chan[0]
    else:
        n_sample = n_chan[1]

    sfreq = opt['sfreq']
    avgWin = opt['WinAverage']

    # Initialize returned values
    ind_good = 1  # index for kept data

    # WINDOWING
    Lwin = round(opt['WinLength'] * opt['sfreq'])  # n data points per window
    Loverlap = round(Lwin * opt["WinOverlap"] / 100)  # n data points in overlap

    # If window is too small
    Messages = []
    if Lwin < 50:
        print("Time window too small, please increase and run process again.")
        return
    # If window is bigger than the data
    elif Lwin > n_sample:
        Lwin = len(F[0])
        Lwin = Lwin - (Lwin % 2)  # Make sure the number of samples is even
        Loverlap = 0
        Nwin = 1
        print("Time window too large, using entire recording for spectrum.")
    # Else: there is at least one full time window
    else:
        Lwin = Lwin - (Lwin % 2)  # Make sure the number of samples is even
        Nwin = (n_sample - Loverlap) // (Lwin - Loverlap)

    # Positive frequency bins spanned by FFT
    FreqVector = sfreq / 2 * np.linspace(0, 1, round(Lwin / 2 + 1))
    # Determine hamming window shape/power
    Win = np.hanning(Lwin)
    WinNoisePowerGain = sum(Win**2)
    # Initialize STFT, time matrices
    ts = np.full((Nwin-(avgWin-1)),np.nan)
    if len(n_chan) > 1:
        TF = np.full((len(F), Nwin - (avgWin - 1), len(FreqVector)), np.nan)
        TFtmp = np.full((len(F), avgWin, len(FreqVector)), np.nan)
    else:
        TF = np.full((1, Nwin - (avgWin - 1), len(FreqVector)), np.nan)
        TFtmp = np.full((1, avgWin, len(FreqVector)), np.nan)
    # Calculate FFT for each window
    TFfull = np.zeros((len(F), Nwin, len(FreqVector)))
    for iWin in range(Nwin):
        # print(iWin)
        # Build indices
        iTimes = list(range((iWin)*(Lwin-Loverlap),Lwin+(iWin)*(Lwin-Loverlap)))
        center_time = floor(median(np.add(np.array(iTimes),1))-\
            (avgWin-1)/2*(Lwin-Loverlap))/200
        if len(n_chan) > 1:
            Fwin = F[:, iTimes]
            Fwin = Fwin- Fwin.mean(axis=1, keepdims=True)
        else:
            # print(iTimes)
            Fwin = F[iTimes]
            # print(Fwin)
            Fwin = Fwin - Fwin.mean()
        # Apply a Hann window to signal
        Fw = np.multiply(Fwin,Win)
        # Compute FFT
        Ffft = np.fft.fft(Fw, Lwin)
        if len(n_chan) > 1:
            TFwin = Ffft[:, :Lwin//2+1] * \
                np.sqrt(2 / (sfreq * WinNoisePowerGain))
            TFwin[:, [0, -1]] = TFwin[:, [0, -1]] / np.sqrt(2)
            # print(TFwin)
            TFwin = np.abs(TFwin)**2
            TFfull[:,iWin,:] = TFwin
            TFtmp[:, iWin % avgWin, :] = TFwin
            # print(iWin%avgWin)
        else:
            TFwin = Ffft[:Lwin//2+1] * np.sqrt(2 / (sfreq * WinNoisePowerGain))
            TFwin[[0, -1]] = TFwin[[0, -1]] / np.sqrt(2)
            # print(TFwin)
            TFwin = np.abs(TFwin)**2
            TFfull[:,iWin,:] = TFwin
            TFtmp[:, iWin % avgWin, :] = TFwin
        if np.isnan(TFtmp[-1, -1, -1]):
            pass
            # continue  # Do not record anything until transient is gone
        else:
            # Save STFTs for window
            TF[:, iWin - (avgWin - 1), :] = np.mean(TFtmp, axis=1)
            ts[iWin - (avgWin - 1)] = center_time
            # print(center_time)

    output = {
        "TF": TF,
        "freqs": FreqVector,
        "ts": ts,
        "options": opt
    }

    return output


def SPRiNT_remove_outliers(fooof_chan, ts, opt):
    ''' SPRiNT_remove_outliers: helper function to remove outlier peaks
    according to user-defined specifications

    Input
    fooof_chan - fooof model for a given channel (FOOOFObject)
    ts - sampled times (numpy array)
    opt - Model settings/hyperparameters (dict)

    Author: Luc Wilson (2023)
    '''
    n_peaks_changed = True
    peaks = fooof_chan.get_params('gaussian_params')
    n_peaks = len(peaks)
    npeak_bytime = fooof_chan.n_peaks_

    peaks_tmp = peaks
    n_peaks_tmp = n_peaks
    npeak_bytime_tmp = deepcopy(npeak_bytime)

    while n_peaks_changed:
        n_peaks_changed = False
        remove = [False for _ in range(n_peaks_tmp)]
        n_peaks_tmp2 = n_peaks_tmp

        for p in range(n_peaks_tmp):
            n_close = 0
            close_t = [(np.abs(peaks_tmp[p,3] - peaks_tmp[r,3])\
                <= opt['maxTime']) for r in range(n_peaks_tmp2)]
            close_f = [(np.abs(peaks_tmp[p,0] - peaks_tmp[r,0])\
                <= opt['maxFreq']) for r in range(n_peaks_tmp2)]

            for r in range(n_peaks_tmp2):
                if close_t[r] and close_f[r]:
                    n_close +=1 # had to flesh out to avoid bugs

            if n_close < (opt['minNear']+1): # fewer than min neighbors
                remove[p] = True # remove this peak
                npeak_bytime_tmp[int(peaks_tmp[p,3])] -= 1 # one less peak here
                n_peaks_tmp -= 1 # one less peak overall
                n_peaks_changed = True

        peaks_tmp = peaks_tmp[[not bln for bln in remove]]

    fg = list([])
    for t in range(len(ts)):
        tmp = fooof_chan.get_fooof(t)
        if npeak_bytime_tmp[t] != npeak_bytime[t]:
            if npeak_bytime_tmp[t] == 0:
                # all peaks at this time were removed
                gaus_pars = []
                pk_fit = gen_periodic(tmp.freqs, [])
                ap_pars = tmp._simple_ap_fit(tmp.freqs,tmp.power_spectrum)
                ap_fit = gen_aperiodic(tmp.freqs, ap_pars)
            else:
                # some peaks at this time were removed
                # finds the indices where time = t
                gaus_pars = peaks_tmp[np.where(peaks_tmp[:,3] == t)[0],:3]
                pk_fit = gen_periodic(tmp.freqs, np.ndarray.flatten(gaus_pars))
                ap_pars = tmp._simple_ap_fit(tmp.freqs,tmp.power_spectrum-pk_fit)
                ap_fit = gen_aperiodic(tmp.freqs, ap_pars)
            error = np.abs(tmp.power_spectrum - ap_fit - pk_fit).mean()
            r_squared = np.corrcoef(tmp.power_spectrum, ap_fit + pk_fit)[0][1]**2
            tmp.add_results(FOOOFResults(ap_pars,\
                tmp._create_peak_params(gaus_pars), r_squared, error, gaus_pars))
        fg.append(tmp)
    fg = combine_fooofs(fg)
    return fg



# Begin user-specific code
os.chdir('/Users/lucwilson/Desktop')
# Get data (others may not need, depending on data format)
with open('sample_ts.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    n = []
    for row in csv_reader:
        n.append(row)

# Set hyperparameters
opt = {
    "sfreq": 200,  # Input sampling rate
    "WinLength": 1,  # STFT window length
    "WinOverlap": 50,  # Overlap between sliding windows (in %)
    "WinAverage": 5, # Number of overlapping windows being averaged
    "rmoutliers": 1, # Apply peak post-processing
    "maxTime": 6, # Maximum distance of nearby peaks in time (in n windows)
    "maxFreq": 2.5, # Maximum distance of nearby peaks in frequency (in Hz)
    "minNear": 3, # Minimum number of similar peaks nearby (using above bounds)
    }

# generate numpy array for data
if len(n) == 1:
   # convert F into a numpy array
    n = n[0]
    n = [float(x) for x in n]
    F = np.array(n)
else:
    n = [[float(x) for x in n[i]] for i in range(len(n))]
    F = np.array(n)

# print(F)
# F = F[0,0:2400]

# expects a numpy array
output = SPRiNT_stft_py(F,opt)

# For architecture
# print(output['ts'])
# print(output['TF'].shape)
# print(output['TF'][0,0,:])

# run fooof across channels and time
# Only issue: Does not use previous window's exponent estimate, haven't seen
# discrepancies yet (...)
fg = fooof.FOOOFGroup(peak_width_limits=[2, 6],\
    min_peak_height=0.5, max_n_peaks = 3)
fgs = fooof.fit_fooof_3d(fg, output['freqs'], output['TF'], freq_range=[1, 40])

# before removing outliers
print(fgs[0].get_params('peak_params'))

fgs = [SPRiNT_remove_outliers(fgs[i], output['ts'], opt) for i in range(2)]

# after removing outliers
print(fgs[0].get_params('peak_params'))
