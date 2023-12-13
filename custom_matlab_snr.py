"""
matlab_snr.py
20230831

python equivalent to matlab snr(x, fs, n) function with custom
fundamental frequency: snr(x, fs, n, tf)

fork from https://github.com/hrtlacek/SNR
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.stats import norm
import scipy


def custom_mean(targetFrequency, snr):
    """
    count the mean snr among the snr results where target frequency is not 0.
    Args:
        targetFrequency (np.ndarray):  each row contains a set of target
            frequency(s). 0 in it means it is just a placeholder.
        snr (np.ndarray): snr results whose shape is same as targetFrequency.

    Returns:
        mean_snr (np.ndarray): mean snr results, column vector, whose length
            is equal to the row number of snr.
    """
    valid_mask = targetFrequency > 1e-8
    valid_count = np.sum(valid_mask, axis=-1)
    mean_snr = np.sum(snr * valid_mask, axis=-1) / valid_count
    return mean_snr


def snr_helper(waveform, sampleRate, targetFrequency, numberHarmonics=1,
               isPlot=False, isLogScale=True):
    """
    count snr for each frequency in targetFrequency
    Args:
        waveform (numpy.ndarray): waveforms, each row contains one piece of
            signal in time domain.
        sampleRate (int|float): sample rate in Hz.
        targetFrequency (numpy.ndarray|None): each row contains a set of target
            frequency(s). If it is None, target frequency will be found
            at the maximum of PSD.
        numberHarmonics (int): number of harmonic to exclude from the SNR
            computation, specified as a positive integer scalar. Default: 1.
        isPlot (bool): whether to plot the SNR result. Default: False.
        isLogScale (bool): whether to show frequency in log scale while
            plotting snr results. Default: True.

    Returns:
        (numpy.ndarray): b by n matrix contains the snr results, where b
            is the number of waveforms, and n is the number of target
            frequencies
    """
    if targetFrequency is None:
        snr_results, _ = snr(
            waveform=waveform,
            sampleRate=sampleRate,
            numberHarmonics=numberHarmonics,
            targetFrequency=targetFrequency,
            isPlot=isPlot,
            isLogScale=isLogScale
        )
    else:
        b, n = targetFrequency.shape
        snr_results = np.zeros([b, n], dtype=float)
        if n == 1:
            snr_results[:, 0], _ = snr(
                waveform=waveform,
                sampleRate=sampleRate,
                numberHarmonics=numberHarmonics,
                targetFrequency=targetFrequency,
                isPlot=isPlot,
                isLogScale=isLogScale
            )
        else:
            # count the snr for each target frequency
            for i in range(n):
                targetFrequencyRolled = np.roll(
                    targetFrequency, -1*i, axis=-1
                )
                snr_results[:, i], _ = snr(
                    waveform=waveform,
                    sampleRate=sampleRate,
                    numberHarmonics=numberHarmonics,
                    targetFrequency=targetFrequencyRolled,
                    isPlot=isPlot,
                    isLogScale=isLogScale
                )
    return snr_results


def snr(waveform, sampleRate, numberHarmonics, targetFrequency=None,
        alias=False, isPlot=False, isLogScale=True):
    """
    Count SNR with the same procedure as matlab's snr(x, fs, n).
    Different from matlab's snr function, whose fundamental frequency
    is found by find the frequency whose corresponding PSD is the
    maximum among all the frequencies, the fundamental frequency could
    be set in here. If two or more target frequencies are set, this
    function will treat the first one as the "fundamental" frequency,
    and treat others as "harmonics" frequency.

    A designed usage is:
        snr(x, fs, 1, target_frequency)

    this usage will count the signal to noise and distortion (SINAD) indeed,
    and its result will be closed to the snr result counting from time domain.

    Args:
        waveform (numpy.ndarray): waveforms, each row contains one piece of
            signal in time domain.
        sampleRate (int|float): sample rate in Hz
        numberHarmonics (int): number of harmonic to exclude from the SNR
            computation, specified as a positive integer scalar.
        targetFrequency (numpy.ndarray|None): each row contains a set of target
            frequency(s). In each row, the first frequency will be seen as
            the fundamental frequency whose power will be treated as signal
            power, and the others will be treated as harmonics whose power
            will be excluded from counting the noise power. If it is None,
            fundamental frequency will be found at the maximum of PSD.
            Default: None
        isPlot (bool): whether to plot the SNR result. Default: False.
        isLogScale (bool): whether to show frequency in log scale while
            plotting snr results. Default: True.

    Returns:
        (numpy.ndarray): a vector contains the snr results

    """
    # input check
    if targetFrequency is not None:
        assert waveform.shape[0] == targetFrequency.shape[0], \
            f"{waveform.shape[0]} waveform(s) are input, " \
            f"but {targetFrequency.shape[0]} set(s) of target frequency(s) " \
            f"are" \
            f"provided, NOT EQUAL!!!"
    assert numberHarmonics > 0, f'number of harmonics({numberHarmonics}) ' \
                                 f'should > 0'

    b, n = waveform.shape

    # get periodogram, parametrized like in matlab
    window = ('kaiser', 38)
    freq, psd = sig.periodogram(
        waveform, fs=sampleRate, window=window)
    # freq [nfft, ]
    # psd [b, nfft]
    # freq_expand = np.broadcast_to(freq, [b, freq.size])  # [b, n]

    # save a copy of the original PSD estimates
    origPsd = psd.copy()

    if targetFrequency is None:
        nGivenHarm = 0
    else:
        nGivenHarm = targetFrequency.shape[-1] - 1
    # pre-allocate harmonic table
    psdHarmPow = np.zeros([b, numberHarmonics + nGivenHarm])
    psdHarmFreq = np.zeros([b, numberHarmonics + nGivenHarm])
    harmIdx = np.zeros([b, numberHarmonics + nGivenHarm, 2], dtype=int)
    dcIdx = np.zeros([b, 2], dtype=int)

    rbw = equivalentNoiseBandwidth(window, n, sampleRate)

    # bump DC component by 3dB and remove it.
    psd[:, 0] *= 2.
    psd_bump_dc = psd.copy()
    _, _, _, iLeft, iRight = getToneFromPSD(psd_bump_dc, freq, rbw, 0)
    for i in range(b):
        if iLeft[i] >= 0 and iRight[i] >= 0:
            if targetFrequency is not None and \
                    freq[iRight[i]] >= targetFrequency[i, 0]:
                # target frequency has been recognized as DC
                # only remove 0 Hz psd
                psd[i, 0] = 0
                dcIdx[i] = 0, 0
            else:
                # remove whole dc component
                psd[i, iLeft[i]:iRight[i]+1] = 0
                dcIdx[i] = iLeft[i], iRight[i]

    # get an estimate of the actual frequency / amplitude
    if targetFrequency is None:
        Pfund, Ffund, iFund, iLeft, iRight = getToneFromPSD(
            psd_bump_dc, freq, rbw)
    else:
        Pfund, Ffund, iFund, iLeft, iRight = getToneFromPSD(
            psd_bump_dc, freq, rbw, toneFreq=targetFrequency[:, 0])

    # fill into 1st harmonic
    psdHarmPow[:, 0] = 10 * np.log10(psd[range(b), iFund]+np.finfo(float).eps)  # advanced indexing
    psdHarmFreq[:, 0] = freq[iFund]

    # remove fundamental psd
    harmIdx[:, 0, 0] = iLeft[:]
    harmIdx[:, 0, 1] = iRight[:]
    for i in range(b):
        psd[i, iLeft[i]:iRight[i]+1] = 0.

    # remove harmonic content
    for h in range(2, numberHarmonics+1):
        toneFreq = Ffund * h
        if alias:
            toneFreq = aliasToNyquist(toneFreq, sampleRate)
        harmPow, _, iHarm, iLeft, iRight = getToneFromPSD(
            psd_bump_dc, freq, rbw, toneFreq)
        psdHarmPow[:, h-1] = 10 * np.log10(
            psd[range(b), iHarm] + np.finfo(float).eps)  # advanced indexing
        psdHarmFreq[:, h-1] = freq[iHarm]
        # remove the power of this tone
        harmIdx[:, h-1, 0] = iLeft[:]
        harmIdx[:, h-1, 1] = iRight[:]
        for i in range(b):
            if harmPow[i] > 0:
                psd[i, iLeft[i]:iRight[i]+1] = 0.

    # remove given harmonic content
    if nGivenHarm:
        for idxGivenHarm in range(nGivenHarm):
            toneFreq = targetFrequency[:, idxGivenHarm+1]
            if alias:
                toneFreq = aliasToNyquist(toneFreq, sampleRate)
            harmPow, _, iHarm, iLeft, iRight = getToneFromPSD(
                psd_bump_dc, freq, rbw, toneFreq)
            psdHarmPow[:, numberHarmonics+idxGivenHarm] = 10 * np.log10(
                psd[range(b), iHarm] + np.finfo(float).eps)  # advanced indexing
            psdHarmFreq[:, numberHarmonics+idxGivenHarm] = freq[iHarm]
            # remove the power of this tone
            harmIdx[:, numberHarmonics+idxGivenHarm, 0] = iLeft[:]
            harmIdx[:, numberHarmonics+idxGivenHarm, 1] = iRight[:]
            for i in range(b):
                if harmPow[i] > 0:
                    psd[i, iLeft[i]:iRight[i] + 1] = 0.

    # get an estimate of the noise floor by computing the median
    # noise power of the non-harmonic region
    estimatedNoiseDensity = np.median(psd[psd > 0], axis=-1)

    # extrapolate estimated noise density into dc/signal/harmonic regions
    psd[psd < np.finfo(float).eps] = estimatedNoiseDensity

    # prevent estimate from obscuring low peaks
    psd = np.minimum(psd, origPsd)

    # compute the noise distortion.
    totalNoise = bandpower(psd, freq)  # [b,]

    snr_result = 10 * np.log10(Pfund / totalNoise + np.finfo(float).eps)  # [b,]
    noisePower = 10 * np.log10(totalNoise + np.finfo(float).eps)  # [b,]

    if isPlot:
        index = 0
        plotSNR(
            psd=origPsd[index],
            freq=freq,
            dcIndizes=dcIdx[index],
            fullHarmonicBin=harmIdx[index],
            harmPeakFreqs=psdHarmFreq[index],
            harmPeakPow=psdHarmPow[index],
            snr=snr_result[index],
            rbw=rbw,
            isLogScale=isLogScale
        )

    return snr_result, noisePower


def aliasToNyquist(f, fs):
    """

    Args:
        f (int|float|numpy.ndarray): frequency with shape [b,]
        fs (int|float|numpy.ndarray): sample rate with shape [b,]

    Returns:
        (int|float|numpy.ndarray): aliased frequency.
    """
    f = np.remainder(f, fs)
    if f > fs/2:
        f = fs - f
    return f


def equivalentNoiseBandwidth(window, Nx, sample_rate):
    """
    returns the two-sided equivalent noise bandwidth
    (in Hz) for a uniformly sampled window whose coefficients are specified
    in the vector WINDOW, where Fs is the sampling rate of the window.
    MATLAB enbw(window, fs) equivalent.
    Args:
        window (str | tuple | array):
            Desired window to use. If `window` is a string or tuple, it is
            passed to `scipy.get_window` to generate the window values,
            which are
            DFT-even by default. See `get_window` for a list of windows and
            required parameters. If `window` is array_like it will be used
            directly as the window and its length must be nperseg. Defaults
            to 'boxcar'.\
        Nx (int): The number of samples in the window.
        sample_rate (int): sample rate in Hz.

    Returns:

    """
    window = sig.get_window(window, Nx)
    # compute normalized ENBW
    bw_t = np.mean(window**2) / np.square(np.mean(window))
    bw = bw_t * sample_rate / window.size
    return bw


def getToneFromPSD(psd, freq, rbw=None, toneFreq=None):
    """
    Equivalent to MATLAB signal.internal.getToneFromPSD(Pxx, F, rbw)
    Retrieve the power and frequency of a windowed sinusoid
    Args:
        psd (numpy.ndarray): power spectrogram density with shape [b, nfft+1]
        freq (numpy.ndarray): frequency with shape [nfft+1,]
        rbw (float|None):the resolution bandwidth.
        toneFreq (float|ndarray|None): tone frequency to be located with shape
            [b,]. A float value means all the samples share this tone
            frequency. None means finding the tone with maximum psd.

    Returns:
        tuple (power, toneFreqTrue, idxTone, idxLeft, idxRight), where

        - power(numpy.ndarray): tone powers with shape [b,];

        - toneFreqTrue(numpy.ndarray): tone central frequencies with shape [b,];

        - idxTone(numpy.ndarray): indices of tones with shape [b,].

        - idxLeft(numpy.ndarray): freq[idxLeft[i]] is the left border
          frequency of tone in i-th psd.

        - idxRight(numpy.ndarray): freq[idxRigt[i]] is the right border
          boarder frequency of tone in i-th psd.

    """
    b, n = psd.shape
    if toneFreq is None:
        idxTone = np.argmax(psd, axis=-1)  # [b,]
    else:
        if isinstance(toneFreq, float) or isinstance(toneFreq, int):
            idxTone = np.argmin(np.abs(freq - toneFreq), axis=-1)  # scalar
            idxTone *= np.ones(psd.shape[0], dtype=int)  # [b,]
        else:  # ndarray
            idxTone = np.argmin(
                np.abs(np.expand_dims(freq, axis=0) -
                       np.reshape(toneFreq, [-1, 1])), axis=-1)  # [b,]
        # look for local peak in vicinity of tone
        for i in range(b):
            iLeftBin = max(0, idxTone[i] - 1)
            iRightBin = min(idxTone[i] + 1, n-1)
            idxMax = np.argmax(psd[i, iLeftBin:iRightBin+1])
            idxTone[i] = iLeftBin + idxMax

    power = np.ones(b) * -1.
    toneFreqTrue = np.ones(b) * -1.
    allIdxLeft = np.ones(b, dtype=int) * -1
    allIdxRight = np.ones(b, dtype=int) * -1

    for i in range(b):
        # sidelobes treated as noise
        idxLeft = idxTone[i] - 1
        idxRight = idxTone[i] + 1
        # roll down slope to left
        while idxLeft >= 0 and psd[i, idxLeft] <= psd[i, idxLeft+1]:
            idxLeft -= 1
        # roll down slope to right
        while idxRight < n and psd[i, idxRight-1] >= psd[i, idxRight]:
            idxRight += 1
        # provide indices to the tone border (inclusive)
        idxLeft += 1
        idxRight -= 1
        # compute the central moment in the neighborhood of the peak
        freqFund = freq[idxLeft:idxRight+1]
        psdFund = psd[i, idxLeft:idxRight+1]
        toneFreqTrue[i] = np.sum(freqFund * psdFund / np.sum(psdFund))

        # report back the integrated power in this band
        if idxLeft < idxRight:
            # more than one bin
            power[i] = bandpower(psd[i, idxLeft:idxRight+1],
                              freq[idxLeft:idxRight+1])
        elif 0 < idxRight < n-1:
            # otherwise just use the current bin
            power[i] = psd[i, idxRight] * \
                    (freq[idxRight + 1] - freq[idxRight - 1]) / 2
        else:
            # otherwise just use the average bin width
            power[i] = psd[i, idxRight] * np.mean(np.diff(freq))

        # protect against nearby tone invading the window kernel
        if rbw is not None and power[i] < rbw * psd[i, idxTone[i]]:
            power[i] = rbw * psd[i, idxTone[i]]
            toneFreqTrue[i] = freq[idxTone[i]]
        # update
        allIdxLeft[i] = idxLeft
        allIdxRight[i] = idxRight
    idxLeft = allIdxLeft
    idxRight = allIdxRight
    return power, toneFreqTrue, idxTone, idxLeft, idxRight


def bandpower(psd, freq):
    """
    estimate bandpower,
    see https://de.mathworks.com/help/signal/ref/bandpower.html
    """
    return scipy.integrate.trapz(psd, freq)


def plotSNR(psd, freq, dcIndizes, fullHarmonicBin, harmPeakFreqs,
            harmPeakPow, snr, rbw=None, isLogScale=True):
    """
    plot SNR, equivalent to MATLAB's plotSNR
    Args:
        psd (numpy.ndarray): power spectrogram density vector.
        freq (numpy.ndarray): frequency vector
        dcIndizes (numpy.ndarray|None): start and end indices of DC component.
        fullHarmonicBin (numpy.ndarray|None): n harmonic (include
            fundamental) components start and end indices of frequency
            with shape [n, 2].
        harmPeakFreqs (numpy.ndarray|None): n harmonic (include
            fundamental) components peak frequency with shape [n,].
        harmPeakPow (numpy.ndarray|None): n harmonic (include
            fundamental) components peak power with shape [n,].
        rbw (float|None):the resolution bandwidth.
        snr (float): snr result.

    Returns:

    """
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 12}

    matplotlib.rc('font', **font)
    # tnrfont = {'fontname': 'Times New Roman'}
    fig, ax = plt.subplots()
    arrowprops = dict(
        arrowstyle="->",
        connectionstyle="angle,angleA=0,angleB=90,rad=10")
    arrowprops1 = dict(
        arrowstyle="simple",
        connectionstyle="arc3")
    bbox = dict(boxstyle="round", fc="0.8", alpha=0.4)
    offset = 10
    if rbw is not None:
        psd *= rbw
        harmPeakPow += 10 * np.log10(rbw + np.finfo(float).eps)
    psd = 10 * np.log10(psd + np.finfo(float).eps)
    start_f = fullHarmonicBin[0, 0]
    end_f = fullHarmonicBin[0, 1]+1

    if isLogScale:
        # plot noise
        plt.semilogx(freq, psd, c='k', label='noise')

        # plot fundamental
        plt.semilogx(
            freq[start_f:end_f], psd[start_f:end_f], c='r',
            linewidth=3, label='fundamental')
    else:
        # plot noise
        plt.plot(freq, psd, c='k', label='noise')

        # plot fundamental
        plt.plot(
            freq[start_f:end_f], psd[start_f:end_f], c='r',
            linewidth=3, label='fundamental')

    ax.annotate(f"F",
                (harmPeakFreqs[0], harmPeakPow[0]),
                xytext=(0, 4), textcoords='offset points',
                bbox=bbox, arrowprops=None)

    # plot DC component
    if dcIndizes is not None:
        plt.plot(
            freq[dcIndizes[0]:dcIndizes[1] + 1],
            psd[dcIndizes[0]:dcIndizes[1] + 1], c='g', linewidth=2,
            label='$DC$' if dcIndizes[1]>dcIndizes[0] else None)

    # plot harmonics
    if fullHarmonicBin is not None:
        for i in range(1, fullHarmonicBin.shape[0]):
            start = fullHarmonicBin[i, 0]
            end = fullHarmonicBin[i, 1]+1
            # whether overlay with fundamental
            if start_f <= start <= end_f:
                start = end_f

            if start_f <= end <= end_f:
                end = start_f

            if start_f > start and end_f < end:
                pass  # just plot all the interval
            if isLogScale:
                plt.semilogx(freq[start:end], psd[start:end], 'b',
                         label='harmonics / other signals' if i == 1 else
                         None)
            else:
                plt.plot(freq[start:end], psd[start:end], 'b',
                         label='harmonics / other signals' if i == 1 else
                         None)

            ax.annotate(f"f{i+1}",
                        (harmPeakFreqs[i], harmPeakPow[i]),
                        xytext=(0, 30), textcoords='offset points',
                        bbox=bbox, arrowprops=arrowprops)

    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'PSD (dB re $\rm{rad}^2 /\rm{Hz}$)')
    plt.title(f'SNR={snr :.4f} dB')
    plt.legend()
    plt.grid()
    plt.show()
    return


if __name__ == '__main__':
    print('generating mono frequency signal...')
    fs = 4e4
    t = np.arange(1/fs, 0.05, 1/fs)  # seconds
    f = 200

    clean = np.sin(2 * np.pi * f * t)
    np.random.seed(111)  # reproducability

    def awgn(s, snr):
        # https://stackoverflow.com/a/53688043
        # https://www.cnblogs.com/skykill/p/7474136.html
        # Convert to linear Watt units
        ps_vs_pn = 10 ** (snr / 10)
        signal_power = np.sum(s ** 2) / len(s)
        nosie_power = signal_power / ps_vs_pn
        # Generate noise samples
        mean_noise = 0
        noise = np.random.normal(
            mean_noise, np.sqrt(nosie_power), len(s))
        return s + noise

    target_snr = 10  # dB
    signal = awgn(clean, target_snr)

    plt.plot(t, signal, c='k', label='signal')
    plt.plot(t, clean, c='r', label='clean')
    plt.xlabel('Time (s)')
    plt.ylabel(r'Value')
    plt.title(f'signal SNR={target_snr :.4f} dB')
    plt.legend()
    plt.grid()
    plt.show()

    all_snr = snr_helper(
        waveform=np.reshape(signal, [1, -1]),
        sampleRate=fs,
        numberHarmonics=1,
        targetFrequency=np.ones([1, 1]) * f,
        isPlot=True,
        isLogScale=True
    )
    print(f'snr of 1 harmonic: {all_snr[0, 0]} dB')

    # snr from time domain
    snr_t = 10 * np.log10(
        (np.sum(signal ** 2) / len(signal)) /
        (np.sum((signal-clean) ** 2) / len(signal))
    )
    print(f'SNR from time domain: {snr_t:.4f} dB')

    # dual frequency
    print('generating dual frequency signal...')
    f2 = 500
    clean2 = np.sin(2 * np.pi * f2 * t) + clean
    signal2 = awgn(clean2, target_snr)

    plt.plot(t, signal2, c='k', label='signal')
    plt.plot(t, clean2, c='r', label='clean')
    plt.xlabel('Time (s)')
    plt.ylabel(r'Value')
    plt.title(f'signal SNR={target_snr :.4f} dB')
    plt.legend()
    plt.grid()
    plt.show()

    all_snr = snr_helper(
        waveform=np.reshape(signal2, [1, -1]),
        sampleRate=fs,
        numberHarmonics=1,
        targetFrequency=np.array([[f, f2]]),
        isPlot=True,
        isLogScale=True
    )
    mean_snr = np.mean(all_snr)
    print(f'SNR:\n\toverall (mean): {mean_snr} dB'
          f'\n\t{f} Hz: {all_snr[0, 0]} dB'
          f'\n\t{f2} Hz: {all_snr[0, 1]} dB')

    sqrt_sum_square_snr = np.sqrt(np.sum(np.square(all_snr)))
    print(f'\toverall (sqrt sum square): {sqrt_sum_square_snr} dB')

    # snr from time domain
    snr_t = 10 * np.log10(
        (np.sum(signal2 ** 2) / len(signal2)) /
        (np.sum((signal2 - clean2) ** 2) / len(signal2))
    )
    print(f'SNR from time domain: {snr_t:.4f} dB')

    pass
