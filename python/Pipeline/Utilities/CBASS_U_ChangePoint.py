#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 09:36:28 2024

Change point algorithm to detect locomotion onset and offset from traces of running speed

@author: Quentin Perrenoud
"""

import numpy as np
import math
from scipy.signal import convolve
import matplotlib.pyplot as plt


def CBASS_U_ChangePoint(db1_trace, in_win_len, bl_rem_low=False, bl_plot=False):
    """
    Finds the indices of the onset and offset of periods of high variance in a 1D time series.

    Parameters:
        db1_trace (array): A 1D time series.
        in_win_len (int): The length of the window for the analysis in number of samples.
        bl_rem_low (bool, optional): If True, removes epochs where variance is high but average is low. Default is False.
        bl_plot (bool, optional): If True, plots the result for visual inspection. Default is False.

    Returns:
        tuple: Indices of the onset (in1_on_idx) and offset (in1_off_idx) of high variance epochs.
    """
    # Ensure window length is even
    if in_win_len % 2 != 0:
        in_win_len += 1

    # Calculate some useful variables
    bl1_numeric = np.isfinite(db1_trace).astype(float)
    db_trace_SD = np.std(db1_trace)

    # Moving window calculations
    in1_win = np.ones(in_win_len)
    db1_mov_av = convolve(db1_trace, in1_win, mode='same')
    in1_dof = convolve(bl1_numeric, in1_win, mode='same')
    db1_mov_av_sq = convolve(db1_trace**2, in1_win, mode='same')
    db1_mov_SD = np.sqrt((db1_mov_av_sq - (db1_mov_av**2 / in1_dof)) / (in1_dof - 1))

    # Threshold the moving SD
    bl1_sig = db1_mov_SD > db_trace_SD
    bl1_sig[0] = bl1_sig[-1] = 0
    in1_on_idx = np.where(np.diff(bl1_sig.astype(int)) == 1)[0] + 1
    in1_off_idx = np.where(np.diff(bl1_sig.astype(int)) == -1)[0] + 1

    # Remove close triggers
    in1_rem_idx = np.where((in1_on_idx[1:] - in1_off_idx[:-1]) <= in_win_len)[0]
    in1_on_idx = np.delete(in1_on_idx, in1_rem_idx + 1)
    in1_off_idx = np.delete(in1_off_idx, in1_rem_idx)

    # Renormalize moving average
    db1_mov_av /= in1_dof

    # Forward and backward z-scores
    db1_z_forward = np.concatenate([
        [np.nan] * in_win_len,
        db1_mov_av[in_win_len + in_win_len//2:-in_win_len//2] / db1_mov_SD[in_win_len//2:-in_win_len//2 - in_win_len],
        [np.nan] * in_win_len
    ])
    db1_z_backward = np.concatenate([
        [np.nan] * in_win_len,
        db1_mov_av[in_win_len//2:-in_win_len - in_win_len//2] / db1_mov_SD[in_win_len + in_win_len//2:-in_win_len//2],
        [np.nan] * in_win_len
    ])

    # Refine indices
    for ii in range(len(in1_on_idx)):
        if in1_on_idx[ii] == 0:
            continue
        in_beg_win = max(in1_on_idx[ii] - in_win_len//2, 0)
        in_end_win = min(in1_on_idx[ii] + in_win_len//2, len(db1_trace))
        in1_on_idx[ii] = in_beg_win + np.argmax(db1_z_forward[in_beg_win:in_end_win])

    for ii in range(len(in1_off_idx)):
        if in1_off_idx[ii] == len(db1_trace) - 1:
            continue
        in_beg_win = max(in1_off_idx[ii] - in_win_len//2, 0)
        in_end_win = min(in1_off_idx[ii] + in_win_len//2, len(db1_trace))
        in1_off_idx[ii] = in_beg_win + np.argmax(db1_z_backward[in_beg_win:in_end_win])

    # Remove short periods
    in1_rem_idx = np.where((in1_off_idx - in1_on_idx) <= in_win_len//2)[0]
    in1_on_idx = np.delete(in1_on_idx, in1_rem_idx)
    in1_off_idx = np.delete(in1_off_idx, in1_rem_idx)

    # Remove low average epochs
    if bl_rem_low:
        db1_av_epoch = np.array([np.mean(np.abs(db1_trace[start:end]))
                                 for start, end in zip(in1_on_idx, in1_off_idx)])
        in1_rem_idx = np.where(db1_av_epoch <= db_trace_SD)[0]
        in1_on_idx = np.delete(in1_on_idx, in1_rem_idx)
        in1_off_idx = np.delete(in1_off_idx, in1_rem_idx)

    # Optional plotting
    if bl_plot:
        plt.figure()
        plt.plot(db1_trace, 'k')
        for idx in in1_on_idx:
            plt.axvline(x=idx, color='g', linestyle='--')
        for idx in in1_off_idx:
            plt.axvline(x=idx, color='r', linestyle='--')
        plt.show()

    return in1_on_idx, in1_off_idx