"""Compute time-frequency decomposition base on Morlet or Multitaper methods.

This script contains the function:

1. _tf_decomp used to decompose the sinal in tf domains using Morlet or
    Multitaper
1. _create_kernel: Create a kernel to smooth the spectra (either boxcar or
    hanning)
2. _smooth_kernel: Perform the smoothing operation on the spectra based on the
   convolution theorem
"""
# Authors : Vinicius Lima <vinicius.lima.cordeiro@gmail.com >
#           Etienne Combrisson <e.combrisson@gmail.com>
#
# License : BSD (3-clause)

import xarray as xr
import numpy as np

from mne.time_frequency import tfr_array_morlet, tfr_array_multitaper
from scipy.signal import fftconvolve

from frites.io import set_log_level, logger
from frites.conn import conn_io


def _tf_decomp(data, sf, freqs, mode='morlet', n_cycles=7.0, mt_bandwidth=None,
               decim=1, kw_cwt={}, kw_mt={}, n_jobs=1):
    """Time-frequency decomposition using MNE-Python.
    Parameters
    ----------
    data : array_like
        Electrophysiological data of shape (n_trials, n_chans, n_times)
    sf : float
        Sampling frequency
    freqs : array_like
        Central frequency vector.
    mode : {'morlet', 'multitaper'}
        Spectrum estimation mode can be either: 'multitaper' or 'morlet'.
    n_cycles : array_like | 7.
        Number of cycles to use for each frequency. If a float or an integer is
        used, the same number of cycles is going to be used for all frequencies
    mt_bandwidth : int | float | array_like | None
        The bandwidth of the multitaper windowing function in Hz. Only used in
        'multitaper' mode.
    decim : int | 1
        To reduce memory usage, decimation factor after time-frequency
        decomposition. default 1 If int, returns tfr[…, ::decim]. If slice,
        returns tfr[…, decim].
    kw_cwt : dict | {}
        Additional arguments sent to the mne-function
        :py:`mne.time_frequency.tfr_array_morlet`
    kw_mt : dict | {}
        Additional arguments sent to the mne-function
        :py:`mne.time_frequency.tfr_array_multitaper`
    Returns
    -------
    out : array_like
        Time-frequency transform of shape (n_epochs, n_chans, n_freqs, n_times)
    """
    if mode == 'morlet':
        out = tfr_array_morlet(
            data, sf, freqs, n_cycles=n_cycles, output='complex', decim=decim,
            n_jobs=n_jobs, **kw_cwt)
    elif mode == 'multitaper':
        # In case multiple values are provided for mt_bandwidth
        # the MT decomposition is done separatedly for each
        # Frequency center
        if isinstance(mt_bandwidth, (list, tuple, np.ndarray)):
            # Arrays freqs, n_cycles, mt_bandwidth should have the same size
            assert len(freqs) == len(n_cycles) == len(mt_bandwidth)
            out = []
            for f_c, n_c, mt in zip(freqs, n_cycles, mt_bandwidth):
                out += [tfr_array_multitaper(
                    data, sf, [f_c], n_cycles=float(n_c), time_bandwidth=mt,
                    output='complex', decim=decim, n_jobs=n_jobs, **kw_mt)]
            out = np.stack(out, axis=2).squeeze()
        elif isinstance(mt_bandwidth, (type(None), int, float)):
            out = tfr_array_multitaper(
                data, sf, freqs, n_cycles=n_cycles,
                time_bandwidth=mt_bandwidth, output='complex', decim=decim,
                n_jobs=n_jobs, **kw_mt)
    else:
        raise ValueError('Method should be either "morlet" or "multitaper"')

    return out

###############################################################################
###############################################################################
#                             SPECTRA ESTIMATION
###############################################################################
###############################################################################


def wavelet_spec(
        data, freqs=None, roi=None, times=None, sfreq=None,
        foi=None, sm_times=.5, sm_freqs=1, sm_kernel='hanning', mode='morlet',
        n_cycles=7., mt_bandwidth=None, decim=1, kw_cwt={},
        kw_mt={}, block_size=None, n_jobs=-1, verbose=None):
    """Wavelet-based single-trial auto-spectra.

    Parameters
    ----------
    data : array_like
        Electrophysiological data. Several input types are supported :

            * Standard NumPy arrays of shape (n_epochs, n_roi, n_times)
            * mne.Epochs
            * xarray.DataArray of shape (n_epochs, n_roi, n_times)

    freqs : array_like
        Array of frequencies.
    roi : array_like | None
        ROI names of a single subject. If the input is an xarray, the
        name of the ROI dimension can be provided
    times : array_like | None
        Time vector array of shape (n_times,). If the input is an xarray, the
        name of the time dimension can be provided
    sfreq : float | None
        Sampling frequency
    foi : array_like | None
        Extract frequencies of interest. This parameters should be an array of
        shapes (n_freqs, 2) defining where each band of interest start and
        finish.
    sm_times : float | .5
        Number of points to consider for the temporal smoothing in seconds. By
        default, a 500ms smoothing is used
    sm_freqs : int | 1
        Number of points for frequency smoothing. By default, 1 is used which
        is equivalent to no smoothing
    kernel : {'square', 'hanning'}
        Kernel type to use. Choose either 'square' or 'hanning'
    mode : {'morlet', 'multitaper'}
        Spectrum estimation mode can be either: 'multitaper' or 'morlet'.
    n_cycles : array_like | 7.
        Number of cycles to use for each frequency. If a float or an integer is
        used, the same number of cycles is going to be used for all frequencies
    mt_bandwidth : array_like | None
        The bandwidth of the multitaper windowing function in Hz. Only used in
        'multitaper' mode.
    decim : int | 1
        To reduce memory usage, decimation factor after time-frequency
        decomposition. default 1 If int, returns tfr[…, ::decim]. If slice,
        returns tfr[…, decim].
    kw_cwt : dict | {}
        Additional arguments sent to the mne-function
        :py:`mne.time_frequency.tfr_array_morlet`
    kw_mt : dict | {}
        Additional arguments sent to the mne-function
        :py:`mne.time_frequency.tfr_array_multitaper`
    block_size : int | None
        Number of blocks of trials to process at once. This parameter can be
        use in order to decrease memory load. If None, all trials are used. If
        for example block_size=2, the number of trials are subdivided into two
        groups and each group is process one after the other.
    n_jobs : int | 1
        Number of jobs to use for parallel computing (use -1 to use all
        jobs). The parallel loop is set at the pair level.

    Returns
    -------
    sxx : xarray.DataArray
        Auto-spectra, DataArray of shape (n_trials, n_roi, n_freqs, n_times)
    """
    set_log_level(verbose)

    # _________________________________ INPUTS ________________________________
    # inputs conversion
    data, cfg = conn_io(
        data, times=times, roi=roi, agg_ch=False, win_sample=None,
        pairs=None, sort=True, block_size=block_size, sfreq=sfreq,
        freqs=freqs, foi=foi, sm_times=sm_times, sm_freqs=sm_freqs,
        name=f'Spectra (mode={mode})', verbose=verbose,
    )

    # extract variables
    x, trials, attrs = data.data, data['y'].data, cfg['attrs']
    times, n_trials = data['times'].data, len(trials)
    indices, sfreq = cfg['blocks'], cfg['sfreq']
    freqs, need_foi, foi_idx = cfg['freqs'], cfg['need_foi'], cfg['foi_idx']
    f_vec, sm_times, sm_freqs = cfg['f_vec'], cfg['sm_times'], cfg['sm_freqs']
    roi, n_rois, n_freqs = data["roi"].data, data.sizes["roi"], len(freqs)

    # temporal decimation
    if isinstance(decim, int):
        times = times[::decim]
        sm_times = int(np.round(sm_times / decim))
        sm_times = max(sm_times, 1)

    # kernel smoothing definition
    kernel = _create_kernel(sm_times, sm_freqs, kernel=sm_kernel)

    # show info
    logger.info(f"Compute auto- and cross-spectra ("
                f"n_freqs={n_freqs}, decim={decim}, sm_times={sm_times}, "
                f"sm_freqs={sm_freqs})")

    # _______________________________ SPECTRA _______________________________
    # compute auto-spectra on blocks of trials
    sxx = np.zeros((n_trials, n_rois, len(f_vec), len(times)))

    for tr in indices:
        # --------------------------- TIME-FREQUENCY --------------------------
        # time-frequency decomposition
        w = _tf_decomp(
            x[tr, ...], sfreq, freqs, n_cycles=n_cycles, decim=decim,
            mode=mode, mt_bandwidth=mt_bandwidth, kw_cwt=kw_cwt, kw_mt=kw_mt,
            n_jobs=n_jobs)

        # ----------------------------- COHERENCE -----------------------------
        # auto spectra (faster that w * w.conj())
        s_auto = w.real ** 2 + w.imag ** 2
        # smooth the auto spectra
        s_auto = _smooth_spectra(s_auto, kernel)

        # mean inside frequency sliding window (if needed)
        if need_foi:
            sxx[tr, ...] = _foi_average(s_auto, foi_idx)
        else:
            sxx[tr, ...] = s_auto

        # Free some memory
        del s_auto

    # _________________________________ OUTPUTS _______________________________
    # configuration
    cfg = dict(
        sfreq=sfreq, sm_times=sm_times, sm_freqs=sm_freqs, sm_kernel=sm_kernel,
        mode=mode, n_cycles=n_cycles, mt_bandwidth=mt_bandwidth, decim=decim
    )

    # conversion
    sxx = xr.DataArray(sxx, dims=('trials', 'roi', 'freqs', 'times'),
                       name='sxx', coords=(trials, roi, f_vec, times),
                       attrs={**attrs, **cfg})

    return sxx

###############################################################################
###############################################################################
#                         SPECTRA SMOOTHING METHODS
###############################################################################
###############################################################################


def _create_kernel(sm_times, sm_freqs, kernel='hanning'):
    """2D (freqs, time) smoothing kernel.
    Parameters
    ----------
    sm_times : int, array_like
        Number of points to consider for the temporal smoothing,
        if it is an array it will be considered that the kernel
        if frequence dependent.
    sm_freqs : int
        Number of points to consider for the frequency smoothing
    kernel : {'square', 'hanning'}
        Kernel type to use. Choose either 'square' or 'hanning'
    Returns
    -------
    kernel : array_like
        Smoothing kernel of shape (sm_freqs, sm_times)
    """
    scale = isinstance(sm_times, np.ndarray)

    if scale:
        # I know this piece of code is terrible ='D
        logger.info("For frequency dependent kernel sm_freqs is not used"
                    "")
        # Number of kernels
        n_kernel = len(sm_times)
        # Get the size of the biggest kernel
        max_size = sm_times.max()
        # Container for the padded kernel
        s_pad = np.zeros((n_kernel, max_size), dtype=np.float32)
        # Store kernel for each frequency
        s = []

        def __pad_kernel(s):
            for i in range(n_kernel):
                #  print(f"{s[i]}")
                pad_size = int(max_size-len(s[i]))
                # The len(s[i])%2 corrects in case the len is odd
                s_pad[i, :] = np.pad(
                    s[i], (pad_size//2, pad_size//2+pad_size % 2))
            return s_pad

    if kernel == 'square':
        if not scale:
            return np.full((sm_freqs, sm_times), 1. / (sm_times * sm_freqs))
        else:
            for i in range(n_kernel):
                s += [np.ones(sm_times[i])/sm_times[i]]
            # Pad with zeros
            return __pad_kernel(s)
    elif kernel == 'hanning':
        if not scale:
            hann_t, hann_f = np.hanning(sm_times), np.hanning(sm_freqs)
            hann = hann_f.reshape(-1, 1) * hann_t.reshape(1, -1)
            return hann / np.sum(hann)
        else:
            for i in range(n_kernel):
                hann = np.hanning(sm_times[i])
                s += [hann/np.sum(hann)]
            return __pad_kernel(s)
    else:
        raise ValueError(f"No kernel {kernel}")


def _smooth_spectra(spectra, kernel, scale=False, decim=1):
    """Smoothing spectra.
    This function assumes that the frequency and time axis are respectively
    located at positions (..., freqs, times).
    Parameters
    ----------
    spectra : array_like
        Spectra of shape (..., n_freqs, n_times)
    kernel : array_like
        Smoothing kernel of shape (sm_freqs, sm_times)
    decim : int | 1
        Decimation factor to apply after the kernel smoothing
    Returns
    -------
    sm_spectra : array_like
        Smoothed spectra of shape (..., n_freqs, n_times)
    """
    # fill potentially missing dimensions
    while kernel.ndim != spectra.ndim:
        kernel = kernel[np.newaxis, ...]
    # smooth the spectra
    if not scale:
        axes = (-2, -1)
    else:
        axes = -1

    spectra = fftconvolve(spectra, kernel, mode='same', axes=axes)
    # return decimated spectra
    return spectra[..., ::decim]


def _foi_average(conn, foi_idx):
    """Average inside frequency bands.

    The frequency dimension should be located at -2.

    Parameters
    ----------
    conn : np.ndarray
        Array of shape (..., n_freqs, n_times)
    foi_idx : array_like
        Array of indices describing frequency bounds of shape (n_foi, 2)

    Returns
    -------
    conn_f : np.ndarray
        Array of shape (..., n_foi, n_times)
    """
    # get the number of foi
    n_foi = foi_idx.shape[0]

    # get input shape and replace n_freqs with the number of foi
    sh = list(conn.shape)
    sh[-2] = n_foi

    # compute average
    conn_f = np.zeros(sh, dtype=conn.dtype)
    for n_f, (f_s, f_e) in enumerate(foi_idx):
        conn_f[..., n_f, :] = conn[..., f_s:f_e, :].mean(-2)
    return conn_f
