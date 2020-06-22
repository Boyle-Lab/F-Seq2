#!urs/bin/env python3

"""F-Seq Version 2.0

"""

__author__ = 'Nanxiang(Samuel) Zhao'
__email__ = 'samzhao@umich.edu'
__version__ = 2.0

import argparse
import sys
import functools
import pandas as pd
import pybedtools
import numpy as np
import KDEpy
import multiprocessing as mp
import time
import math
from scipy.signal import find_peaks
from scipy import stats
from numba import jit
from statsmodels.stats.multitest import multipletests
import os


def read_in_file(file_name):
    """ Read in file(s).
    """
    if file_name[0][-3:] == 'bam':
        bed = pd.concat([pybedtools.BedTool(file).bam_to_bed().to_dataframe() for file in file_name],
                              axis=0)
    elif file_name[0][-3:] == 'bed':
        bed = pd.concat(map(functools.partial(pd.read_csv, sep='\t', header=None), file_name),
                              axis=0).rename(columns={0:'chrom', 1:'start', 2:'end', 3: 'name', 4: 'score', 5:'strand'})
    else:
        sys.exit('Error: Input file format should be either ended with .bam or .bed, please convert files and try again.')
    bed = bed.drop(['name', 'score'], axis=1)
    bed['cuts'] = bed['start'] * (bed['strand'] == '+') + \
                        bed['end'] * (bed['strand'] == '-')

    return bed


def calculate_fragment_size(cuts_df, verbose=False):
    """Estimate fragment size using the chr with most cuts.
    This estimation is done by conducting Wilcoxon-Mann-Whitney rank sum test. These fragments are usually < 200bp
    so we can assume that the average is < 500.

    Input:
        cuts_df: (pd.Dataframe) the cuts dataframe with most cuts, must have cols =
        ['chrom', 'start', 'end', 'strand', 'cuts']

    Return:
        fragment_size: (int) estimated fragment size
    """
    max_range = 500 # max search window for each + cut
    prior = 150
    max_cut_to_use = 1000000 #1 million is a good estimate
    if max_cut_to_use > cuts_df.shape[0] - max_range:
        max_cut_to_use = cuts_df.shape[0] - max_range

    cuts_df_mod = cuts_df[:max_cut_to_use].copy()
    pos_cuts_ls = cuts_df_mod[cuts_df_mod['strand'] == '+']['cuts'].values
    cuts_ls = cuts_df_mod['cuts'].values
    neg_strand_ls = (cuts_df_mod['strand'] == '-').values

    if verbose:
        print('Calculating fragment size (specify -f 0 if using DNase HS data) ...',
              flush=True)

    fragment_size = wmw_rank_searchsorted(cuts_ls=cuts_ls, pos_cuts_ls=pos_cuts_ls, max_range=max_range, prior=prior,
                                          neg_strand_ls=neg_strand_ls)
    #fragment_size = int(np.median(np.concatenate([cuts_ls[(cuts_ls > (current_cuts - max_range + prior)) &
                                                          #(cuts_ls < (current_cuts + max_range + prior)) &
                                                          #neg_strand_ls] - current_cuts
                                                  #for current_cuts in pos_cuts_ls], axis=0)))

    if verbose:
        print('Done!', flush=True)

    return fragment_size


def wmw_rank_searchsorted(cuts_ls, pos_cuts_ls, max_range, prior, neg_strand_ls):
    """ Wilcon-Mann-Whitney rank sum test using sorted array.
    """
    cuts_ls_index = np.argsort(cuts_ls)
    cuts_ls = cuts_ls[cuts_ls_index]
    neg_strand_ls = neg_strand_ls[cuts_ls_index]

    min_values = np.array([current_cuts - max_range + prior for current_cuts in pos_cuts_ls])
    max_values = min_values + 2 * max_range

    i0_ls = np.searchsorted(cuts_ls, min_values, 'left')
    i1_ls = np.searchsorted(cuts_ls, max_values, 'right')

    fragment_size = int(np.median(np.concatenate(
        [cuts_ls[i0:i1][neg_strand_ls[i0:i1]] - current_cuts for current_cuts, i0, i1 in
         zip(pos_cuts_ls, i0_ls, i1_ls)])))

    return fragment_size


def calculate_bandwidth(feature_length):
    """Calculate bandwidth.
    """
    bandwidth = float(feature_length / 2 / 3) # 3 standard deviations

    return bandwidth


def calculate_threshold(size, ncuts, std, bandwidth):
    """Calculate threshold for peak calling.

    Input:
        size: sum range of cuts across all chroms in bp
        ncuts: total number of cuts
        std: standard deviation for calculating threshold
        window_size: window size in bp, window_size = int(2860/2) when bandwidth = 100

    Return:
        threshold: threshold for peak calling
    """
    # calculate window_size
    window_size = 1400.0; java_min_value = 2 ** (-149)
    while True:
        window_size += 1
        x = window_size / bandwidth
        v = math.exp(-(x * x) / 2) / math.sqrt(math.pi * 2)
        if v < java_min_value:
            break
    window_size = int(window_size - 1)

    total_window = 1 + window_size * 2
    cut_density = int((ncuts / size) * total_window)
    threshold_iterations = 10000

    np.random.seed(934) # may not need it
    simulated_results = [KDEpy.FFTKDE(bw=bandwidth).fit(np.random.randint(total_window, size=cut_density)).evaluate(
        np.arange(-1, total_window + 1))[window_size + 1] for _ in range(threshold_iterations)]
    # the following has more readability
    #simulated_results = []
    #for _ in range(threshold_iterations):
        #simulated_cuts = np.random.randint(total_window, size=cut_density)
        #kdepy_kde = KDEpy.FFTKDE(bw=bandwidth)
        #kdepy_kde.fit(simulated_cuts)
        #simulated_results.append(kdepy_kde.evaluate(np.arange(-1, total_window + 1))[window_size + 1])
    kdepy_result = np.array(simulated_results) * cut_density

    threshold = np.mean(kdepy_result) + std * np.std(kdepy_result)

    return threshold


def init_child(lock_):
    global lock
    lock = lock_


def run_kde(cuts_df, cuts_df_control, chrom,
            fragment_offset, bandwidth, ncuts,
            fragment_offset_control, bandwidth_control, ncuts_control,
            verbose=True):
    """Run kde with control.
    """
    start_time = time.time()
    weights_control = np.ones(cuts_df_control.shape[0])
    weights = np.ones(cuts_df.shape[0])


    ### 1. run kde for input ###
    kdepy_result, first_cut, num_cuts = calculate_kde(cuts_df=cuts_df, weights=weights, fragment_offset=fragment_offset,
                                                      bandwidth=bandwidth, num_cuts=False)
    last_cut = cuts_df['cuts'].max()
    del cuts_df, weights


    ### 2. run kde for control ###
    kdepy_result_control, first_cut_control, _ = calculate_kde(cuts_df=cuts_df_control, weights=weights_control,
                                                               fragment_offset=fragment_offset_control,
                                                               bandwidth=bandwidth_control, num_cuts=num_cuts)
    del cuts_df_control, weights_control

    # align kde_result_control with kde_result
    kdepy_result_control = align_array(kde_result_control=kdepy_result_control, kde_result_shape=kdepy_result.shape[0],
                                       first_cut_control=first_cut_control, first_cut=first_cut)
    #maybe we can run kde for control first => max one kdepy_result at run time
    #first_cut = cuts_df['cuts'].min()
    #kde_result_shape needs some effort
    #so maybe not?

    # unload memory
    with lock:
        control_np_tmp_name = f'{args.o}/{args.name}_{chrom}_control.dat'
        control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='w+', shape=kdepy_result_control.shape)
        control_np_tmp[:] = kdepy_result_control[:]
        del control_np_tmp, kdepy_result_control, first_cut_control
    # if without control: control_np_tmp_name = False


    ### 3. call peaks using input ###
    result_df = call_peaks(chrom=chrom, first_cut=first_cut, kdepy_result=kdepy_result,
                           min_height=threshold, min_distance=min_distance, min_prominence=min_prominence)


    ### 4. calculate q value ###
    summit_abs_pos_array = result_df['summit'].values - first_cut

    # calculate query value array and lambda background
    lambda_bg, query_value = prepare_stats_test(result_df=result_df, kdepy_result=kdepy_result,
                                                control_np_tmp_name=control_np_tmp_name,
                                                summit_abs_pos_array=summit_abs_pos_array,
                                                window_size=window_size, use_max=False)

    # calculate lambda local
    if control_np_tmp_name is False:
        control_np_tmp = kdepy_result
    else:
        control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='r', shape=kdepy_result.shape)
    #with lock: #intermediate output
        #np.savez_compressed(f'{tmp_dir}/{chrom}.npz', arr_0=kdepy_result)
    del kdepy_result
    with lock:
        lambda_local = find_local_lambda(control_np_tmp=control_np_tmp, control_np_tmp_name=control_np_tmp_name,
                                         summit_abs_pos_array=summit_abs_pos_array, lambda_bg=lambda_bg,
                                         window_size=window_size, use_max=False)
        # memory intensive, only allow one process
    del control_np_tmp, summit_abs_pos_array

    # calculate q value
    result_df = calculate_q_value_local(query_value=query_value, lambda_local=lambda_local, result_df=result_df,
                                        distribution_name=distribution_name, threshold_q=0.05)


    end_time = time.time()
    if verbose:
        print(f'{chrom}: first={first_cut}, last={last_cut} \n'
              f'{chrom}: Completed in {end_time - start_time:.3f} seconds. \n'
              f'------------------------------------- \n', flush=True)


    with lock: #move things to storage
        os.remove(control_np_tmp_name)
        #shutil.move(control_np_tmp_name, f'{storage_dir}/{chrom}_control.dat')
        #shutil.move(f'{tmp_dir}/{chrom}.npz', f'{storage_dir}/{chrom}.npz')


    return result_df


def run_kde_wo_control(cuts_df, chrom,
                       fragment_offset, bandwidth, ncuts,
                       verbose=True):
    """Run kde without control.
    In this case, control is the same as input(treatment).
    """
    start_time = time.time()
    #weights_control = np.ones(cuts_df_control.shape[0])
    weights = np.ones(cuts_df.shape[0])


    ### 1. run kde for input ###
    kdepy_result, first_cut, num_cuts = calculate_kde(cuts_df=cuts_df, weights=weights, fragment_offset=fragment_offset,
                                                      bandwidth=bandwidth, num_cuts=False)
    last_cut = cuts_df['cuts'].max()
    del cuts_df, weights


    ### 2. run kde for control ###
    #kdepy_result_control, first_cut_control, _ = calculate_kde(cuts_df=cuts_df_control, weights=weights_control,
                                                               #fragment_offset=fragment_offset_control,
                                                               #bandwidth=bandwidth_control, num_cuts=num_cuts)
    #del cuts_df_control, weights_control

    # align kde_result_control with kde_result
    #kdepy_result_control = align_array(kde_result_control=kdepy_result_control, kde_result_shape=kdepy_result.shape[0],
                                       #first_cut_control=first_cut_control, first_cut=first_cut)
    #maybe we can run kde for control first => max one kdepy_result at run time
    #first_cut = cuts_df['cuts'].min()
    #kde_result_shape needs some effort
    #so maybe not?

    # unload memory
    #with lock:
        #control_np_tmp_name = f'{args.o}/{args.name}_{chrom}_control.dat'
        #control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='w+', shape=kdepy_result_control.shape)
        #control_np_tmp[:] = kdepy_result_control[:]
        #del control_np_tmp, kdepy_result_control, first_cut_control
    # if without control: control_np_tmp_name = False


    ### 3. call peaks using input ###
    result_df = call_peaks(chrom=chrom, first_cut=first_cut, kdepy_result=kdepy_result,
                           min_height=threshold, min_distance=min_distance, min_prominence=min_prominence)


    ### 4. calculate q value ###
    summit_abs_pos_array = result_df['summit'].values - first_cut

    # calculate query value array and lambda background
    lambda_bg, query_value = prepare_stats_test(result_df=result_df, kdepy_result=kdepy_result,
                                                control_np_tmp_name=False,
                                                summit_abs_pos_array=summit_abs_pos_array,
                                                window_size=window_size, use_max=False)

    # calculate lambda local
    #if control_np_tmp_name is False:
        #control_np_tmp = kdepy_result
    #else:
        #control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='r', shape=kdepy_result.shape)
    #with lock: #intermediate output
        #np.savez_compressed(f'{tmp_dir}/{chrom}.npz', arr_0=kdepy_result)
    #del kdepy_result
    with lock:
        lambda_local = find_local_lambda(control_np_tmp=kdepy_result, control_np_tmp_name=False,
                                         summit_abs_pos_array=summit_abs_pos_array, lambda_bg=lambda_bg,
                                         window_size=window_size, use_max=False)
    #del control_np_tmp, summit_abs_pos_array

    # calculate q value
    result_df = calculate_q_value_local(query_value=query_value, lambda_local=lambda_local, result_df=result_df,
                                        distribution_name=distribution_name, threshold_q=0.05)


    end_time = time.time()
    if verbose:
        print(f'{chrom}: first={first_cut}, last={last_cut} \n'
              f'{chrom}: Completed in {end_time - start_time:.3f} seconds. \n'
              f'------------------------------------- \n', flush=True)


    #with lock: #move things to storage
        #os.remove(control_np_tmp_name)
        #shutil.move(control_np_tmp_name, f'{storage_dir}/{chrom}_control.dat')
        #shutil.move(f'{tmp_dir}/{chrom}.npz', f'{storage_dir}/{chrom}.npz')


    return result_df


def calculate_kde(cuts_df, weights, fragment_offset, bandwidth, num_cuts=False):
    """Calculate kde array without all chrs scaling.
    Individual chrom scaling and num cuts scaling.

    """
    cuts = cuts_df['cuts'].values
    first_cut = cuts.min()
    last_cut = cuts.max()
    if fragment_offset == 0:
        kdepy_kde = KDEpy.FFTKDE(bw=bandwidth).fit(cuts, weights=weights)
        kdepy_result = kdepy_kde.evaluate(np.arange(first_cut - 1, last_cut + 2))[1:-2]
        # kdepy_result = (kdepy_result * 20000000) / ncuts * (weights.sum())
    else:
        cuts_w_offset = ((cuts_df['start'] + fragment_offset) * (cuts_df['strand'] == '+') +
                         (cuts_df['end'] - fragment_offset) * (cuts_df['strand'] == '-')).astype('int').values
        kde_evaluation_start = min(cuts_w_offset.min(), first_cut);
        kde_evaluation_end = max(cuts_w_offset.max(), last_cut)
        # assert (first_cut - kde_evaluation_start) >= 0
        kdepy_kde = KDEpy.FFTKDE(bw=bandwidth).fit(cuts_w_offset, weights=weights)
        try:
            kdepy_result = kdepy_kde.evaluate(np.arange(kde_evaluation_start - 1, kde_evaluation_end + 2))[1:-2]
        except:
            np.fft.restore_all()  # if encounter the bug for MLK https://github.com/IntelPython/mkl_fft/issues/24
            kdepy_result = kdepy_kde.evaluate(np.arange(kde_evaluation_start - 1, kde_evaluation_end + 2))[1:-2]
        kdepy_result = kdepy_result[(first_cut - kde_evaluation_start):(last_cut - kde_evaluation_start)]
        # kdepy_result = (kdepy_result * 20000000) / ncuts * (weights.sum()) / 2
        # this scaling is simulating results with offset (left window positive cuts and right window negative cuts) by
        # multiplying all signal with 0.5

    if not num_cuts:
        num_cuts = cuts.shape[0]
    kdepy_result = (kdepy_result * num_cuts * scaling_factor).astype(np.float32) # lets use all reads, more related to threshold

    return kdepy_result, first_cut, num_cuts


def align_array(kde_result_control, kde_result_shape, first_cut_control, first_cut):
    """ Extract corresponding kde_result_control given kde_result's length and first cut.
    Fill nan if kde_result_control is shorter, then interpolate based on neighboring value.
    Directly overwrite kde_result_control to avoid extra copy.

    Return:
        kde_result_control: (np.array) same length and same first_cut as kde_result
    """
    if first_cut > first_cut_control:
        # align
        abs_start = first_cut - first_cut_control
        # extract or concatenate depends on aligned length
        sig_shape = kde_result_shape
        control_shape = kde_result_control.shape[0] - abs_start
        if sig_shape > control_shape:
            # concatenate
            kde_result_control = kde_result_control[abs_start:]
            kde_result_control = np.concatenate(
                [kde_result_control, np.full(sig_shape - control_shape, fill_value=np.nan)])
            #kde_result_control = np.concatenate(
                #[kde_result_control, np.zeros(sig_shape - control_shape, dtype=np.float32)])
        else:
            # extract
            kde_result_control = kde_result_control[abs_start:abs_start + sig_shape]
    else:
        # align
        abs_start = first_cut_control - first_cut
        # extract or concatenate depends on aligned length
        sig_shape = kde_result_shape - abs_start
        control_shape = kde_result_control.shape[0]
        if sig_shape > control_shape:
            # concatenate
            kde_result_control = np.concatenate([np.full(abs_start, fill_value=np.nan), kde_result_control,
                                                 np.full(sig_shape - control_shape, fill_value=np.nan)])
            #kde_result_control = np.concatenate([np.zeros(abs_start, dtype=np.float32), kde_result_control,
                                                 #np.zeros(sig_shape - control_shape, dtype=np.float32)])
        else:
            # extract
            kde_result_control = kde_result_control[:sig_shape]
            kde_result_control = np.concatenate([np.full(abs_start, fill_value=np.nan), kde_result_control])
            #kde_result_control = np.concatenate([np.zeros(abs_start, dtype=np.float32), kde_result_control])

    mask = np.isnan(kde_result_control)
    kde_result_control[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), kde_result_control[~mask])

    return kde_result_control


def call_peaks(chrom, first_cut, kdepy_result, min_height, min_prominence=None, min_distance=None, min_peak_len=0,
               merge_peak_distance=1, continous_regions_stats=False):
    """Call peaks satisfying certeria.

    Input:
        chrom: (str) chrom name
        first_cut: (int) first cut on chrom
        kde_result: (np.array)
        min_height: (float) peak calling threshold
        min_prominence: (float)
        min_distance: (float) min distance between peak calls
        min_peak_len: (int) min length of a region considered to be a peak
        merge_peak_distance: (float) distance to allow merge between peaks

    Return:
        peak_regions: (pd.DataFrame) [chr, start, end, summit, signal]
    """
    peak_regions = find_contiguous_regions(kdepy_result >= min_height)
    peak_indexes, properties = find_peaks(kdepy_result, height=min_height, distance=min_distance,
                                          prominence=min_prominence)
    peak_regions = pd.concat([pd.Series([chrom] * peak_regions.shape[0], name='chrom'),
                              pd.DataFrame(peak_regions, columns=['start', 'end'])], axis=1)
    peak_indexes = pd.concat([pd.Series([chrom] * peak_indexes.shape[0], name='chrom'),
                              pd.Series(peak_indexes, name='start'),
                              pd.Series(peak_indexes + 1, name='end'),
                              pd.Series(properties['peak_heights'], name='score'),
                              pd.Series(properties['prominences'], name='strand')], axis=1)  # change
    peak_regions['start'] += first_cut
    peak_regions['end'] += first_cut
    peak_indexes['start'] += first_cut
    peak_indexes['end'] += first_cut
    peak_regions = pybedtools.BedTool.from_dataframe(peak_regions)
    peak_indexes = pybedtools.BedTool.from_dataframe(peak_indexes)
    peak_regions = peak_regions.intersect(peak_indexes, wa=True, wb=True).merge(c=[5, 7, 8], o='collapse',
                                                                                d=merge_peak_distance).to_dataframe()  # change

    peak_regions = peak_regions[peak_regions['end'] - peak_regions['start'] >= min_peak_len].copy()
    # print(peak_regions)

    peak_regions.rename(columns={'name': 'summit'}, inplace=True)

    if continous_regions_stats:
        print((peak_regions.iloc[:, 2] - peak_regions.iloc[:, 1]).describe())

    # focus on summits
    peak_regions = peak_regions.astype({'summit': 'str', 'score': 'str', 'strand': 'str'})
    peak_regions['summit'] = peak_regions['summit'].str.split(',')
    peak_regions['score'] = peak_regions['score'].str.split(',')
    peak_regions['strand'] = peak_regions['strand'].str.split(',')  # change
    peak_regions = explode_multiindex(df=peak_regions, lst_cols=['summit', 'score', 'strand'])  # change
    peak_regions = peak_regions.astype(
        {'summit': 'int32', 'score': 'float32', 'strand': 'float32'})  # memory: may delete strand finally?
    peak_regions['summit_end'] = peak_regions['summit'] + 1

    return peak_regions


def find_contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index. Notice this idx is half closed [start, end).
    From https://stackoverflow.com/a/4495197/7674821
    """

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]  # Edit

    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx


def explode_multiindex(df, lst_cols, fill_value='', preserve_index=False):
    """Adapt from https://stackoverflow.com/a/40449726/7674821
    Example:
    >>> df
           aaa  myid        num          text
        0   10     1  [1, 2, 3]  [aa, bb, cc]
        1   11     2         []            []
        2   12     3     [1, 2]      [cc, dd]
        3   13     4         []            []

    >>> explode_multiindex(df, ['num','text'], fill_value='')
           aaa  myid num text
        0   10     1   1   aa
        1   10     1   2   bb
        2   10     1   3   cc
        3   11     2
        4   12     3   1   cc
        5   12     3   2   dd
        6   13     4

    """
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
            and len(lst_cols) > 0
            and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
        col: np.repeat(df[col].values, lens)
        for col in idx_cols},
        index=idx).assign(**{col: np.concatenate(df.loc[lens > 0, col].values)
                             for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens == 0, idx_cols], sort=False)
               .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res


def prepare_stats_test(result_df, kdepy_result, control_np_tmp_name, summit_abs_pos_array, window_size, use_max=False):
    """Prepare for stats tests. Focus on memory.

    Notice big change here: no lower limit for average!!!


    Return:
        lambda_bg: (np.float32) background lambda
        query_value: (np.array) query values on summits with kdepy_result

    """
    minimum_value = 0.0001  # max or average value should not smaller than this

    if control_np_tmp_name is False:
        control_np_tmp = kdepy_result
    else:
        control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='r', shape=kdepy_result.shape)


    ### 1. find max (or average) signal for each non-overlapping window, then estimate lambda_bg ###
    num_windows = int(control_np_tmp.shape[0] / window_size)
    if use_max:
        value_windows = np.max(control_np_tmp[:num_windows * window_size].reshape(num_windows, window_size), axis=1)

    else:
        value_windows = np.average(control_np_tmp[:num_windows * window_size].reshape(num_windows, window_size), axis=1)
        #value_windows = control_np_tmp[:]
    value_windows[value_windows < minimum_value] = minimum_value

    lambda_bg = max(lambda_bg_lower_bound, np.average(value_windows)) # lower bound for lambda_bg #change
    if control_np_tmp_name:
        del control_np_tmp # remove reference or flush out to disk


    ### 2. pre-process ###
    if use_max:
        query_value = result_df['strand'].values
    else:
        query_value = extract_value_around_window(kdepy_result, summit_abs_pos_array, window_length=window_size,
                                                  operation='average')

    return lambda_bg, query_value


def extract_value_around_window(result_sig, summit_abs_pos_array, window_length, operation='average'):
    """ Extract average (or max or var) value around center of summit in the region of window_length.
    Input:
        results_sig: (np.array) 1d
        summit_abs_pos_array: (np.array) 1d and is absolute position in result_sig
        window_length: (int) window size to calculate average or find max
    Return:
        window_length_array: (np.array) 1d contains the average (or max or var) value centered at each summit with window_length
    """
    if operation == 'average':
        window_length_array = np.average(
            result_sig[np.clip(np.arange(-window_length / 2, window_length / 2, dtype='int32')
                               + summit_abs_pos_array[:, np.newaxis], a_min=0, a_max=result_sig.shape[0] - 1)], axis=1)
    elif operation == 'max':
        window_length_array = np.max(result_sig[np.clip(np.arange(-window_length / 2, window_length / 2, dtype='int32')
                                                        + summit_abs_pos_array[:, np.newaxis], a_min=0,
                                                        a_max=result_sig.shape[0] - 1)], axis=1)
    elif operation == 'var':
        window_length_array = np.var(result_sig[np.clip(np.arange(-window_length / 2, window_length / 2, dtype='int32')
                                                        + summit_abs_pos_array[:, np.newaxis], a_min=0,
                                                        a_max=result_sig.shape[0] - 1)], axis=1)

    return window_length_array.astype(np.float32)  # memory


def find_local_lambda(control_np_tmp, control_np_tmp_name, summit_abs_pos_array, lambda_bg, window_size, use_max=False):
    """Find local lambda.
    Memory intensive.

    Return:
        lambda_local: (np.array) max lambda for each summit
    """
    # 3. find all local lambda
    oneK_length = 1000; fiveK_length = 5000; tenK_length = 10000

    # extract average (or max) value around each summit window
    if use_max:
        tenK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=tenK_length, operation='max')
        fiveK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=fiveK_length, operation='max')
        oneK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=oneK_length, operation='max') #change
        if control_np_tmp_name:
            #region_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=window_size, operation='max') #change
            oneK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=oneK_length, operation='max')
    else:
        tenK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=tenK_length, operation='average')
        fiveK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=fiveK_length, operation='average')
        oneK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=oneK_length, operation='average') #change
        if control_np_tmp_name:
            #region_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=window_size, operation='average') #change
            oneK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=oneK_length, operation='average')

    # lambda_local = max(lambda_bg, lambda_region, lambda_1k, lambda_5k, lambda_10k)
    if control_np_tmp_name:
        lambda_local = np.max(np.vstack([oneK_array, fiveK_array, tenK_array]), axis=0)  # change
        del oneK_array, fiveK_array, tenK_array, summit_abs_pos_array #change
    else:
        lambda_local = np.max(np.vstack([oneK_array, fiveK_array, tenK_array]), axis=0) #change
        del fiveK_array, tenK_array, summit_abs_pos_array, oneK_array
    np.clip(lambda_local, a_min=lambda_bg, a_max=None, out=lambda_local)  # memory

    return lambda_local

    # TODO:
    #   1. For memory efficiency, we could only create one array which is tenK_array and slice it along the way


def calculate_q_value_local(query_value, lambda_local, result_df,
                            distribution_name='poisson', threshold_p=False, negative_log=True, threshold_q=0.05):
    """Calculate q-value of summits with BH FDR method.
    Using prominence if max or using signal if average. Default is using average.
    Using local poisson or local nbinom.

    Input:
        result_sig: (np.array) kde signal in bp
        result_df: (pd.DataFrame) required columns ['start', 'end', 'summit', 'summit_end', 'strand'], 'strand' contains actually prominence value
        first_cut: (int) first cut of result_sig and control_sig if have control
        control_sig: (np.array) control kde signal in bp. Same length with result_sig
        window_size: (int) window size to extract maximum signal across all
        distribution: (str) ['poisson', 'nbinom']
        threshold_q: (float) filter out summit whose q-value is larger than the threshold, not in negative log scale. Default is 0.05.
        negative_log: (boolean) whether -log10(q_value), save as separate column
        use_max: (boolean) use maximum for each window; otherwise, use average.

    Return:
        result_df_significance: (pd.DataFrame) result_df and corresponding p-value. Notice 'start' and 'end' could be replicated, but not 'summit' and 'summit_end'


    """
    # 4. fit a distribution and calculate p-value for each summit using prominence
    if distribution_name == 'poisson':
        pvalue_array = [stats.poisson.sf(value, mu=lambda_value) for value, lambda_value in
                        zip(query_value, lambda_local)]
    elif distribution_name == 'nbinom':
        # assume p=0.5
        var_local = lambda_local * 2

        param_ls = [calculate_nbinom_param(mu=mu_value, var=var_value) for mu_value, var_value in
                    zip(lambda_local, var_local)]
        pvalue_array = [stats.nbinom.sf(value, p=param[0], n=param[1]) for value, param in
                        zip(query_value, param_ls)]
    elif distribution_name == 'binom':
        var_local = lambda_local / 2
        param_ls = [calculate_binom_param(mu=mu_value, var=var_value) for mu_value, var_value in
                    zip(lambda_local, var_local)]
        pvalue_array = [stats.binom.sf(value, p=param[0], n=param[1]) for value, param in
                        zip(query_value, param_ls)]

    result_df['p_value'] = pvalue_array
    result_df['query_value'] = query_value # change
    result_df['lambda_local'] = lambda_local # change
    result_df.sort_values(['p_value'], inplace=True)
    if threshold_p:
        result_df = result_df[result_df['p_value'] < threshold_p].copy()


    # 5. convert p value to q value
    #(reject, adjusted_pvalue_array, _, _) = multipletests(pvals=result_df['p_value'], alpha=threshold_q,
                                                          #method='fdr_bh', is_sorted=True)

    # 6. post process
    #result_df['q_value'] = adjusted_pvalue_array
    #result_df = result_df[reject].copy()
    #if negative_log:
        #result_df['-log(q_value)'] = -np.log10(result_df['q_value'])


    return result_df

def interpolate_poisson_p_value(result_df):
    """
    Input:
        result_df: (pd.DataFrame) contain columns ['query_value', 'lambda_local']

    Returns:
        result_df: (pd.DataFrame) with new column ['log10_p_value_interpolated']

    """
    result_df['query_int'] = result_df['query_value'].astype(int)
    result_df = result_df.round({'lambda_local': 2}) #change
    for query_int, lambda_local in result_df.loc[:, ['query_int', 'lambda_local']].drop_duplicates().values:
        query_int = int(query_int)
        rows_to_select = (result_df['query_int'] == query_int) & (result_df['lambda_local'] == lambda_local)
        poisson_func = stats.poisson(mu=lambda_local)
        result_df.loc[rows_to_select, '-log10_p_value_interpolated'] = interpolate_log_p_value(
            result_df.loc[rows_to_select, 'query_value'].values,
            query_int, poisson_func.sf(query_int), query_int + 1, poisson_func.sf(query_int + 1))

    return result_df


def interpolate_log_p_value(query_value_array, xa, ya, xb, yb):
    """ Interpolate p value based on query value, and two points a = (xa, ya), b = (xb, yb)
    Interpolation is in log space.
    Returns:
        interpolated_p_value: (np.array) same size as query_value_array, in -log10 scale.
    """
    interpolated_p_value = (query_value_array - xa) * (np.log10(yb) - np.log10(ya)) / (xb - xa) + np.log10(ya)
    interpolated_p_value[np.isnan(interpolated_p_value)] = 0.0
    return -interpolated_p_value


@jit
def calculate_nbinom_param(mu, var):
    """Input (mu, var), output (p, r)
    """
    p = mu / var
    return p, p * mu / (1 - p)


@jit
def calculate_binom_param(mu, var):
    """Input (mu, var), output (p, n)
    """
    p = var / mu
    return p, mu / p


def calculate_lambda_bg_lower_bound(bandwidth_control, ncuts_treatment, ncuts_control, pseudo_count=1):
    """ Calculate lambda background lower bound.
    Deciding feature is bandwidth_control, and same scaling as control.
    """
    return stats.norm(0, bandwidth_control).pdf(0) * pseudo_count / ncuts_control * ncuts_treatment


if __name__ == '__main__':
    ### 1. user setting ###
    parser = argparse.ArgumentParser(
        description='This command executes fseq which is a feature density estimator for high-throughput sequence tags')

    parser.add_argument('input_file', help='treatment file(s) in bam or bed format. If multiple files, separate '
                                           'by space; if one file, it can contain all chromosomes', nargs='+')
    parser.add_argument('-f', help="fragment size of treatment (default=estimated from data). For DNase HS data (5' ends)"
                                   " set -f 0", default=False, type=int)
    parser.add_argument('-l', help='feature length of treatment (default=600 for DNase, '
                                   'recommend 50 for narrow ChIP-seq, 1000 for broad ChIP-seq)', default=600, type=int)

    parser.add_argument('-control_file', help='control file(s) in bam or bed format. If multiple files, separate '
                                           'by space; if one file, it can contain all chromosomes', nargs='+')
    parser.add_argument('-fc', help="fragment size of control (default=estimated from data). For DNase HS data (5' ends)"
                                   " set -f 0", default=False, type=int)
    #parser.add_argument('-lc', help='feature length of control (default=12000)', default=12000, type=int)

    parser.add_argument('-t', help='threshold (standard deviations) (default=4.0 for broad peaks,'
                                   'recommend 8.0 for narrow peaks)', default=4.0, type=float)
    parser.add_argument('-v', help='verbose output', default=False, action='store_true')
    parser.add_argument('-cpus', help='number of cpus to use (default=min(unique_chrs, all-2))', default=False,
                        type=int)
    parser.add_argument('-o', help='output directory (default=current directory)')
    parser.add_argument('-name', help='prefix for output files (default=fseq_run)', default='fseq_run')
    args = parser.parse_args()


    ### 2. deal with input (treatment) ###
    input_bed = read_in_file(args.input_file)

    # parameters for input
    feature_length = args.l
    if args.f is False:
        fragment_size = calculate_fragment_size(
            input_bed[input_bed['chrom'] == input_bed.groupby('chrom')['chrom'].count().idxmax()], verbose=args.v)
    else:
        fragment_size = args.f
    fragment_offset = int(fragment_size / 2)
    bandwidth = calculate_bandwidth(feature_length)
    ncuts = input_bed.shape[0]

    threshold = calculate_threshold(size=input_bed.groupby('chrom')['cuts'].agg(np.ptp).sum(),
                                    ncuts=ncuts, std=args.t, bandwidth=bandwidth)
    sequence_length = input_bed.iloc[0, 2] - input_bed.iloc[0, 1]


    ### 3. deal with control ###
    if args.control_file:
        control_bed = read_in_file(args.control_file)

        # parameters for control
        feature_length_control = feature_length * 20
        if args.fc is False:
            fragment_size_control = calculate_fragment_size(
                control_bed[control_bed['chrom'] == control_bed.groupby('chrom')['chrom'].count().idxmax()], verbose=args.v)
        else:
            fragment_size_control = args.fc
        fragment_offset_control = int(fragment_size_control / 2)
        bandwidth_control = calculate_bandwidth(feature_length_control)
        ncuts_control = control_bed.shape[0]

        #threshold_control = calculate_threshold(size=control_bed.groupby('chrom')['cuts'].agg(np.ptp).sum(), ncuts=ncuts_control,
                                        #std=args.t, bandwidth=bandwidth_control)
        sequence_length_control = control_bed.iloc[0, 2] - control_bed.iloc[0, 1]

        if (args.control_file[0][-3:] == 'bam') and (args.input_file[0][-3:] == 'bam'):
            assert control_bed.iloc[-1, 0] == input_bed.iloc[-1, 0], "pybedtools incomplete write-in error, " \
                                                                     "please free memory and try again."

    # general parameters
    scaling_factor = feature_length * 20 / 50
    threshold *= scaling_factor
    #threshold_control *= scaling_factor
    #threshold = 3.116 # 4std = 1.657, 8std = 3.116  # speed up for simulated datasets
    #threshold_control = 0.805 # 4std = 0.500, 8std = 0.805  # speed up for simulated dataset
    if args.control_file:
        lambda_bg_lower_bound = calculate_lambda_bg_lower_bound(bandwidth_control=bandwidth_control,
                                                                ncuts_control=ncuts_control, ncuts_treatment=ncuts,
                                                                pseudo_count=1)
    else:
        lambda_bg_lower_bound = 0
        #lambda_bg_lower_bound = calculate_lambda_bg_lower_bound(bandwidth_control=bandwidth,
                                                                #ncuts_control=ncuts, ncuts_treatment=ncuts,
                                                                #pseudo_count=1)
    lambda_bg_lower_bound *= scaling_factor
    min_prominence = threshold # 0.001
    min_distance = fragment_size
    window_size = fragment_size
    distribution_name = 'poisson'


    if args.v:
        print('=====================================', flush=True)
        print(f'F-Seq Version {__version__}', flush=True)
        print('=====================================', flush=True)
        print('Settings: ', flush=True)
        print(f'\tbandwidth = {bandwidth}', flush=True)
        print(f'\tthreshold = {threshold:.3f}', flush=True)
        print(f'\tlambda bg lower bound = {lambda_bg_lower_bound:.3f}', flush=True)
        print(f'\test. fragment size = {fragment_size}', flush=True)
        print(f'\tsequence length = {sequence_length}', flush=True)
        print(f'\ttotal cuts = {ncuts}', flush=True)
        print('-------------------------------------', flush=True)
        if args.control_file:
            print(f'\tcontrol bandwidth = {bandwidth_control}', flush=True)
            print(f'\tcontrol est. fragment size = {fragment_size_control}', flush=True)
            print(f'\tcontrol sequence length = {sequence_length_control}', flush=True)
            print(f'\tcontrol total cuts = {ncuts_control}', flush=True)
            print('-------------------------------------', flush=True)


    ### 4. run kde ###
    chrom_ls = input_bed['chrom'].unique()
    if args.control_file:
        input_param_ls = [[input_bed[input_bed['chrom'] == chrom]] + [control_bed[control_bed['chrom'] == chrom]] +
                          [chrom, fragment_offset, bandwidth, ncuts,
                           fragment_offset_control, bandwidth_control, ncuts_control, args.v] for chrom in chrom_ls]
    else:
        input_param_ls = [[input_bed[input_bed['chrom'] == chrom]] +
                          [chrom, fragment_offset, bandwidth, ncuts,
                           args.v] for chrom in chrom_ls]

    # cpus = min(chrom_ls.shape[0], fseq.mp.cpu_count() - 2)
    cpus = args.cpus
    lock = mp.Lock()

    with mp.Pool(processes=cpus, initializer=init_child, initargs=(lock,)) as pool:
        if args.control_file:
            results = pool.starmap(run_kde, input_param_ls)
        else:
            results = pool.starmap(run_kde_wo_control, input_param_ls)
    result_df = pd.concat(results)

    ### 5. prepare and output ###
    result_df = interpolate_poisson_p_value(result_df)
    result_df.to_hdf(path_or_buf=f'{args.o}/{args.name}.h5', key='result_df', mode='a')

    print('Thanks for using F-seq.')


    # TODO:
    #   1. interpolate_poisson_p_value() needs optimization
    #   2. consider no calculation of p-value inside run_kde(), directly interpolate all together
