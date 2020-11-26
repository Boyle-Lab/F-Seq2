#!urs/bin/env python3

"""F-Seq Version 2.0 module.

"""

__author__ = 'Nanxiang(Samuel) Zhao'
__email__ = 'samzhao@umich.edu'
__version__ = '2.0.1'

import functools
import math
import os
import sys
import time

import KDEpy
import h5py
import numpy as np
import pandas as pd
import pyBigWig
import pybedtools
from scipy import stats
from scipy.ndimage import filters
from scipy.signal import find_peaks, fftconvolve
from statsmodels.stats.multitest import multipletests


def read_in_file(file_name, pe, pe_fragment_size_range, nucleosome_size, chrom_ls_to_process):
    """ Read in file(s).
    """
    if not pe:
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
        fragment_size_before_filter = np.NaN
        fragment_size_after_filter = np.NaN
        filter_out_rate = np.NaN

    else:
        if file_name[0][-3:] == 'bam':
            bed = pd.concat([pybedtools.BedTool(file).bam_to_bed(bedpe=True, mate1=True).to_dataframe() for file in file_name], axis=0)
        elif file_name[0][-3:] == 'bed':
            bed = pd.concat(map(functools.partial(pd.read_csv, sep='\t', header=None), file_name), axis=0).rename(
                columns={0:'chrom', 1:'start', 2:'end', 4: 'score', 5:'strand'})
        else:
            sys.exit('Error: Input file format should be either ended with .bam or .bed, please convert files and try again.')
        bed['start_new'] = (bed[['start', 'score']]).min(axis=1)
        bed['end_new'] = (bed[['end', 'strand']]).max(axis=1)
        try:
            bed['cuts'] = ((bed['start_new'] + bed['end_new']) / 2).astype(np.int32)
        except TypeError:
            bed = bed.loc[(bed['strand'] != -1) & (bed['start'] != -1), :].astype(
                {'start_new': np.int32, 'end_new': np.int32})
            bed['cuts'] = ((bed['start_new'] + bed['end_new']) / 2).astype(np.int32)
        bed['strand'] = '.'
        bed = bed.loc[:, ['chrom', 'start_new', 'end_new', 'strand', 'cuts']]
        bed.rename(columns={'start_new': 'start', 'end_new': 'end'}, inplace=True)
        fragment_size_before_filter = (bed['end'] - bed['start']).mean()
        fragment_count_before_filter = bed.shape[0]
        if pe_fragment_size_range is not False:
            if pe_fragment_size_range[0] != 'auto':
                pe_fragment_size_range[0] = int(pe_fragment_size_range[0])
                if pe_fragment_size_range[1] != 'inf':
                    pe_fragment_size_range[1] = int(pe_fragment_size_range[1])
                    bed = bed.loc[((bed['end'] - bed['start']) >= pe_fragment_size_range[0]) &
                                  ((bed['end'] - bed['start']) <= pe_fragment_size_range[1]), :]
                else:
                    bed = bed.loc[((bed['end'] - bed['start']) >= pe_fragment_size_range[0]), :]
            else:
                mononucleosome_range = [nucleosome_size[0], nucleosome_size[1]] # [180, 247]
                dinucleosome_range = [nucleosome_size[2], nucleosome_size[3]] # [315, 473]
                trinucleosome_range = [nucleosome_size[4], nucleosome_size[5]] # [558, 615]
                fragment_length_series = (bed['end'] - bed['start'])
                mask = ((fragment_length_series >= mononucleosome_range[0]) & (fragment_length_series <= mononucleosome_range[1])) | (
                    (fragment_length_series >= dinucleosome_range[0]) & (fragment_length_series <= dinucleosome_range[1])) | (
                    (fragment_length_series >= trinucleosome_range[0]))# & (fragment_length_series <= trinucleosome_range[1]))
                bed = bed.loc[~mask, :]
        fragment_size_after_filter = (bed['end'] - bed['start']).mean()
        fragment_count_after_filter = bed.shape[0]
        filter_out_rate = float(fragment_count_before_filter - fragment_count_after_filter) / fragment_count_before_filter

    if chrom_ls_to_process:
         bed = bed.loc[bed['chrom'].isin(chrom_ls_to_process)]

    bed.sort_values(['chrom', 'start'], inplace=True)

    return bed, [fragment_size_before_filter, fragment_size_after_filter, filter_out_rate]



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
    if cuts_df.shape[0] < max_range:
        max_cut_to_use = cuts_df.shape[0]
    elif max_cut_to_use > cuts_df.shape[0] - max_range:
        max_cut_to_use = cuts_df.shape[0] - max_range

    cuts_df_mod = cuts_df[:max_cut_to_use].copy()
    pos_cuts_ls = cuts_df_mod[cuts_df_mod['strand'] == '+']['cuts'].values
    cuts_ls = cuts_df_mod['cuts'].values
    neg_strand_ls = (cuts_df_mod['strand'] == '-').values

    fragment_size = wmw_rank_searchsorted(cuts_ls=cuts_ls, pos_cuts_ls=pos_cuts_ls, max_range=max_range, prior=prior,
                                          neg_strand_ls=neg_strand_ls)

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


def calculate_threshold(size, ncuts, std, std_p, bandwidth):
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
    kdepy_result = np.array(simulated_results) * cut_density

    threshold = np.mean(kdepy_result) + std * np.std(kdepy_result)

    peak_region_threshold = np.mean(kdepy_result) + std_p * np.std(kdepy_result)

    return threshold, peak_region_threshold


def init_child(lock_):
    global lock
    lock = lock_


def run_kde(cuts_array, start_array, end_array, strand_array,
            cuts_array_control, start_array_control, end_array_control, strand_array_control, chrom, params):
    """Run kde with control.
    """
    start_time = time.time()


    ### 1. run kde for input ###
    kdepy_result, first_cut, num_cuts = calculate_kde(cuts_array=cuts_array, start_array=start_array,
                                                      end_array=end_array,
                                                      strand_array=strand_array,
                                                      fragment_offset=params.fragment_offset,
                                                      bandwidth=params.bandwidth, scaling_factor=params.scaling_factor,
                                                      num_cuts=False)
    last_cut = np.max(cuts_array)
    del cuts_array, start_array, end_array, strand_array


    ### 2. run kde for control ###
    kdepy_result_control, first_cut_control, _ = calculate_kde(cuts_array=cuts_array_control,
                                                               start_array=start_array_control,
                                                               end_array=end_array_control,
                                                               strand_array=strand_array_control,
                                                               fragment_offset=params.fragment_offset_control,
                                                               bandwidth=params.bandwidth_control,
                                                               scaling_factor=params.scaling_factor,
                                                               num_cuts=num_cuts)
    lambda_bg_lower_bound_ind = calculate_lambda_bg_lower_bound(bandwidth_control=params.bandwidth_control,
                                                                ncuts_control=cuts_array_control.shape[0],
                                                                ncuts_treatment=num_cuts,
                                                                pseudo_count=1) #6.30
    lambda_bg_lower_bound_ind *= params.scaling_factor
    del cuts_array_control, start_array_control, end_array_control, strand_array_control

    # align kde_result_control with kde_result
    kdepy_result_control = align_array(kde_result_control=kdepy_result_control, kde_result_shape=kdepy_result.shape[0],
                                       first_cut_control=first_cut_control, first_cut=first_cut)


    ### 2.5 estimate total number of tests ###
    # unload memory
    with lock: # 7.23
        # 1. if output sig, weight and output
        if params.sig_format:
            with h5py.File(params.treatment_np_tmp_name, mode='a', libver='latest') as sig_file:
                try:
                    sig_file.create_dataset(chrom,
                                            data=np.round(np.divide(kdepy_result,
                                                           kdepy_result_control, out=np.zeros_like(kdepy_result),
                                                           where=kdepy_result_control>lambda_bg_lower_bound_ind),
                                                          params.sig_float_precision).astype(np.float16),
                                            dtype='f2')#, compression='gzip')
                    sig_file.attrs[chrom] = first_cut
                except RuntimeError:
                    del sig_file[chrom]
                    sig_file.create_dataset(chrom,
                                            data=np.round(np.divide(kdepy_result,
                                                           kdepy_result_control, out=np.zeros_like(kdepy_result),
                                                           where=kdepy_result_control>lambda_bg_lower_bound_ind),
                                                          params.sig_float_precision).astype(np.float16),
                                            dtype='f2')#, compression='gzip')
                    sig_file.attrs[chrom] = first_cut

        # 2. calculate lambda_bg which needs all control signal in memory
        lambda_bg = calculate_lambda_bg(kdepy_result_control, window_size=params.window_size,
                                        lambda_bg_lower_bound=params.lambda_bg_lower_bound)

        # 3. unload memory of control sig
        control_np_tmp_name = f'{params.temp_dir_name}/{params.name}_{chrom}_control.dat'
        control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='w+', shape=kdepy_result_control.shape)
        control_np_tmp[:] = kdepy_result_control[:]
        del control_np_tmp, kdepy_result_control, first_cut_control
        # if without control: control_np_tmp_name = False


    ### 3. call peaks using input ###
    result_df = call_peaks(chrom=chrom, first_cut=first_cut, kdepy_result=kdepy_result,
                           min_height=params.threshold, min_distance=params.min_distance, min_prominence=params.min_prominence,
                           peak_region_threshold=params.peak_region_threshold) #7.27


    ### 4. prepare to calculate q value ###
    try:
        summit_abs_pos_array = result_df['summit'].values - first_cut
    except KeyError:
        return pd.DataFrame() # no peaks called

    # calculate query value array
    # change prepare_stats_test to calculate_query_value, if need control_value, return it in find_local_lambda
    query_value = calculate_query_value(result_df=result_df, kdepy_result=kdepy_result,
                                        summit_abs_pos_array=summit_abs_pos_array, window_size=params.window_size)

    # calculate lambda local
    control_np_tmp = np.memmap(control_np_tmp_name, dtype=np.float32, mode='r', shape=kdepy_result.shape)
    del kdepy_result
    with lock:
        lambda_local = find_local_lambda(control_np_tmp=control_np_tmp, control_np_tmp_name=control_np_tmp_name,
                                         summit_abs_pos_array=summit_abs_pos_array, lambda_bg=lambda_bg,
                                         window_size=params.window_size, sparse_data=params.sparse_data, use_max=False)
        # memory intensive, only allow one process
    del control_np_tmp, summit_abs_pos_array

    # prepare for interpolation of p value and q value calculation
    result_df['query_value'] = query_value
    result_df['lambda_local'] = lambda_local


    end_time = time.time()
    if params.v:
        print(f'\t{chrom}: first={first_cut}, last={last_cut}, completed in {end_time - start_time:.3f} seconds.',
              flush=True)


    # cleanup
    with lock: #move things to storage
        os.remove(control_np_tmp_name)

    return result_df


def run_kde_wo_control(cuts_array, start_array, end_array, strand_array, chrom, params):
    """Run kde without control.
    In this case, control is the same as input(treatment).
    """
    start_time = time.time()


    ### 1. run kde for input ###
    kdepy_result, first_cut, num_cuts = calculate_kde(cuts_array=cuts_array, start_array=start_array,
                                                      end_array=end_array, strand_array=strand_array,
                                                      fragment_offset=params.fragment_offset,
                                                      bandwidth=params.bandwidth, scaling_factor=params.scaling_factor,
                                                      num_cuts=False)
    last_cut = np.max(cuts_array)
    del cuts_array, start_array, end_array, strand_array


    ### 3. call peaks using input ###
    result_df = call_peaks(chrom=chrom, first_cut=first_cut, kdepy_result=kdepy_result,
                           min_height=params.threshold, min_distance=params.min_distance, min_prominence=params.min_prominence,
                           peak_region_threshold=params.peak_region_threshold) #7.27


    ### 4. prepare to calculate q value ###
    try:
        summit_abs_pos_array = result_df['summit'].values - first_cut
    except KeyError:
        return pd.DataFrame() # no peaks called

    # calculate query value array and lambda background
    lambda_bg = calculate_lambda_bg(kdepy_result_control=kdepy_result, window_size=params.window_size,
                                    lambda_bg_lower_bound=params.lambda_bg_lower_bound)
    query_value = calculate_query_value(result_df=result_df, kdepy_result=kdepy_result,
                                        summit_abs_pos_array=summit_abs_pos_array, window_size=params.window_size)

    # calculate lambda local
    with lock:
        lambda_local = find_local_lambda(control_np_tmp=kdepy_result, control_np_tmp_name=False,
                                         summit_abs_pos_array=summit_abs_pos_array, lambda_bg=lambda_bg,
                                         window_size=params.window_size, sparse_data=params.sparse_data, use_max=False)

    # prepare for interpolation of p value and q value calculation
    result_df['query_value'] = query_value
    result_df['lambda_local'] = lambda_local

    with lock: # 7.23
        # 1. if output sig, weight and output
        if params.sig_format:
            with h5py.File(params.treatment_np_tmp_name, mode='a', libver='latest') as sig_file:
                try:
                    sig_file.create_dataset(chrom,
                                            data=np.round(kdepy_result, params.sig_float_precision).astype(np.float16),
                                            dtype='f2')#, compression='gzip')
                    sig_file.attrs[chrom] = first_cut
                except RuntimeError:
                    del sig_file[chrom]
                    sig_file.create_dataset(chrom,
                                            data=np.round(kdepy_result, params.sig_float_precision).astype(np.float16),
                                            dtype='f2')#, compression='gzip')
                    sig_file.attrs[chrom] = first_cut


    end_time = time.time()
    if params.v:
        print(f'\t{chrom}: first={first_cut}, last={last_cut}, completed in {end_time - start_time:.3f} seconds.',
              flush=True)

    return result_df


def calculate_kde(cuts_array, start_array, end_array, strand_array, fragment_offset, bandwidth, scaling_factor, num_cuts=False):
    """Calculate kde array without all chrs scaling.
    Individual chrom scaling and num cuts scaling.

    """
    #cuts = cuts_df['cuts'].values
    first_cut = np.min(cuts_array)
    last_cut = np.max(cuts_array)
    if fragment_offset == 0:
        kdepy_kde = KDEpy.FFTKDE(bw=bandwidth).fit(cuts_array)
        try:
            kdepy_result = kdepy_kde.evaluate(np.arange(first_cut - 1, last_cut + 2))[1:-2]
        except ValueError:
            np.fft.restore_all()  # if encounter the bug for MLK https://github.com/IntelPython/mkl_fft/issues/24
            kdepy_result = kdepy_kde.evaluate(np.arange(first_cut - 1, last_cut + 2))[1:-2]
        # kdepy_result = (kdepy_result * 20000000) / ncuts * (weights.sum())
    else:
        cuts_w_offset = ((start_array + fragment_offset) * (strand_array) +
                         (end_array - fragment_offset) * (~strand_array)).astype(np.int32)
        # cuts_w_offset = ((start_array + fragment_offset) * (strand_array == '+') +
        #                  (end_array - fragment_offset) * (strand_array == '-')).astype('int').values
        kde_evaluation_start = min(cuts_w_offset.min(), first_cut);
        kde_evaluation_end = max(cuts_w_offset.max(), last_cut)
        # assert (first_cut - kde_evaluation_start) >= 0
        kdepy_kde = KDEpy.FFTKDE(bw=bandwidth).fit(cuts_w_offset)
        try:
            kdepy_result = kdepy_kde.evaluate(np.arange(kde_evaluation_start - 1, kde_evaluation_end + 2))[1:-2]
        except ValueError:
            np.fft.restore_all()  # if encounter the bug for MLK https://github.com/IntelPython/mkl_fft/issues/24
            kdepy_result = kdepy_kde.evaluate(np.arange(kde_evaluation_start - 1, kde_evaluation_end + 2))[1:-2]
        kdepy_result = kdepy_result[(first_cut - kde_evaluation_start):(last_cut - kde_evaluation_start)]
        # kdepy_result = (kdepy_result * 20000000) / ncuts * (weights.sum()) / 2
        # this scaling is simulating results with offset (left window positive cuts and right window negative cuts) by
        # multiplying all signal with 0.5

    if not num_cuts:
        num_cuts = cuts_array.shape[0]
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


def est_num_total_tests(chrom, kdepy_result, kdepy_result_control, lambda_bg_lower_bound_ind, min_distance, window_size):
    """
    output control_value for inspection?
    """
    # 1. call all summits on kdepy_result
    peak_indexes, properties = find_peaks(kdepy_result, height=0, distance=min_distance,
                                          prominence=0)

    # 2. calculate lambda_region on kdepy_result_control
    lambda_region = extract_value_around_window(kdepy_result_control, peak_indexes, window_length=window_size,
                                                operation='average')

    # 3. count num of summits above lambda_bg_lower_bound_ind
    num_tests_after_filter = np.sum(lambda_region > lambda_bg_lower_bound_ind)
    num_tests_before_filter = lambda_region.shape[0]

    np.savetxt(f'./Fseq2_output/est_num_total_{chrom}.txt', lambda_region)

    return num_tests_before_filter, num_tests_after_filter


def call_peaks(chrom, first_cut, kdepy_result, min_height, peak_region_threshold, min_prominence=None, min_distance=None,
               min_peak_len=0, merge_peak_distance=1, continous_regions_stats=False):
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
    try:
        peak_regions = find_contiguous_regions(kdepy_result >= peak_region_threshold)
    except IndexError:
        print(f'no peaks find on {chrom}')
        return pd.DataFrame()
    peak_indexes, properties = find_peaks(kdepy_result, height=min_height, distance=min_distance,
                                          prominence=min_prominence)
    peak_regions = pd.concat([pd.Series([chrom] * peak_regions.shape[0], name='chrom'),
                              pd.DataFrame(peak_regions, columns=['start', 'end'])], axis=1)
    peak_indexes = pd.concat([pd.Series([chrom] * peak_indexes.shape[0], name='chrom'),
                              pd.Series(peak_indexes, name='start'),
                              pd.Series(peak_indexes + 1, name='end'),
                              pd.Series(properties['peak_heights'], name='score'),
                              pd.Series(properties['prominences'], name='strand')], axis=1)
    peak_regions['start'] += first_cut
    peak_regions['end'] += first_cut
    peak_indexes['start'] += first_cut
    peak_indexes['end'] += first_cut
    peak_regions = pybedtools.BedTool.from_dataframe(peak_regions)
    peak_indexes = pybedtools.BedTool.from_dataframe(peak_indexes)
    try:
        peak_regions = peak_regions.intersect(peak_indexes, wa=True, wb=True).merge(c=[5, 7, 8], o='collapse',
                                                                                    d=merge_peak_distance).to_dataframe()
    except pybedtools.helpers.BEDToolsError:
        print(f'no peaks find on {chrom}')
        return pd.DataFrame()

    peak_regions = peak_regions[peak_regions['end'] - peak_regions['start'] >= min_peak_len].copy()
    # print(peak_regions)

    peak_regions.rename(columns={'name': 'summit'}, inplace=True)

    if continous_regions_stats:
        print((peak_regions.iloc[:, 2] - peak_regions.iloc[:, 1]).describe())

    # focus on summits
    peak_regions = peak_regions.astype({'summit': 'str', 'score': 'str', 'strand': 'str'})
    peak_regions['summit'] = peak_regions['summit'].str.split(',')
    peak_regions['score'] = peak_regions['score'].str.split(',')
    peak_regions['strand'] = peak_regions['strand'].str.split(',')
    peak_regions = explode_multiindex(df=peak_regions, lst_cols=['summit', 'score', 'strand'])
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


def prepare_stats_test(result_df, kdepy_result, control_np_tmp_name, summit_abs_pos_array, window_size,
                       lambda_bg_lower_bound, use_max=False, return_control_value=False):
    """Prepare for stats tests. Focus on memory.

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

    lambda_bg = max(lambda_bg_lower_bound, np.average(value_windows)) # lower bound for lambda_bg

    if return_control_value:
        control_value = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=window_size,
                                                    operation='average') #6.30

    if control_np_tmp_name:
        del control_np_tmp # remove reference or flush out to disk


    ### 2. pre-process ###
    if use_max:
        query_value = result_df['strand'].values
    else:
        query_value = extract_value_around_window(kdepy_result, summit_abs_pos_array, window_length=window_size,
                                                  operation='average')
    if return_control_value: #6.30
        return lambda_bg, query_value, control_value
    else:
        return lambda_bg, query_value


def calculate_query_value(result_df, kdepy_result, summit_abs_pos_array, window_size, use_max=False):
    ### 2. pre-process ###
    if use_max:
        query_value = result_df['strand'].values
    else:
        query_value = extract_value_around_window(kdepy_result, summit_abs_pos_array, window_length=window_size,
                                                  operation='average')
    return query_value


def calculate_lambda_bg(kdepy_result_control, window_size, lambda_bg_lower_bound, use_max=False):

    minimum_value = 0.0001  # max or average value should not smaller than this

    ### 1. find max (or average) signal for each non-overlapping window, then estimate lambda_bg ###
    num_windows = int(kdepy_result_control.shape[0] / window_size)
    if use_max:
        value_windows = np.max(kdepy_result_control[:num_windows * window_size].reshape(num_windows, window_size), axis=1)

    else:
        value_windows = np.average(kdepy_result_control[:num_windows * window_size].reshape(num_windows, window_size), axis=1)
    value_windows[value_windows < minimum_value] = minimum_value

    lambda_bg = max(lambda_bg_lower_bound, np.average(value_windows)) # lower bound for lambda_bg

    return lambda_bg


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


def find_local_lambda(control_np_tmp, control_np_tmp_name, summit_abs_pos_array, lambda_bg, window_size,
                      sparse_data=False, use_max=False):
    """Find local lambda.
    Memory heavy.

    Return:
        lambda_local: (np.array) max or average lambda for each summit
    """
    # 3. find all local lambda
    oneK_length = 1000; fiveK_length = 5000; tenK_length = 10000; tenfiveK_length = 15000

    # extract average (or max) value around each summit window
    if use_max:
        if control_np_tmp_name or sparse_data:
            # region_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=window_size, operation='max')
            oneK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=oneK_length, operation='max')
            # tenfiveK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=tenfiveK_length, operation='max')
        fiveK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=fiveK_length, operation='max')
        tenK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=tenK_length, operation='max')
    else:
        if control_np_tmp_name or sparse_data:
            # region_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=window_size, operation='average')
            oneK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=oneK_length, operation='average')
            # tenfiveK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=tenfiveK_length, operation='average')
        fiveK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=fiveK_length, operation='average')
        tenK_array = extract_value_around_window(control_np_tmp, summit_abs_pos_array, window_length=tenK_length, operation='average')


    if control_np_tmp_name or sparse_data:
        lambda_local = np.max(np.vstack([oneK_array, fiveK_array, tenK_array]), axis=0)
        del oneK_array, fiveK_array, tenK_array, summit_abs_pos_array
    else:
        lambda_local = np.max(np.vstack([fiveK_array, tenK_array]), axis=0)
        del fiveK_array, tenK_array, summit_abs_pos_array
    np.clip(lambda_local, a_min=lambda_bg, a_max=None, out=lambda_local)  # memory

    return lambda_local


def calculate_q_value_local(query_value, lambda_local, result_df,
                            lambda_bg_lower_bound_ind, control_value,
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
    result_df['query_value'] = query_value
    result_df['lambda_local'] = lambda_local
    if control_value:
        result_df['control_value'] = control_value
    result_df['lambda_bg_lower_bound_ind'] = lambda_bg_lower_bound_ind
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
    result_df = result_df.round({'lambda_local': 2})
    with np.errstate(divide='ignore', invalid='ignore'):
        for query_int, lambda_local in result_df.loc[:, ['query_int', 'lambda_local']].drop_duplicates().values:
            query_int = int(query_int)
            rows_to_select = (result_df['query_int'] == query_int) & (result_df['lambda_local'] == lambda_local)
            poisson_func = stats.poisson(mu=lambda_local)
            result_df.loc[rows_to_select, '-log10_p_value_interpolated'] = interpolate_log_p_value(
                result_df.loc[rows_to_select, 'query_value'].values,
                query_int, poisson_func.sf(query_int), query_int + 1, poisson_func.sf(query_int + 1))

    # change inf to max value besides inf
    mask_not_inf = result_df['-log10_p_value_interpolated'] != np.inf
    result_df.loc[~mask_not_inf, '-log10_p_value_interpolated'] = result_df.loc[mask_not_inf, '-log10_p_value_interpolated'].max()

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


def calculate_q_value(result_df, p_thr, q_thr, num_peaks):
    """

    :param result_df: (pd.DataFrame) requires ['-log10_p_value_interpolated']
    :param p_value_cutoff:
    :param q_value_cutoff:
    :param num_peak:
    :return:
    """
    result_df['p_value'] = 10**(-result_df['-log10_p_value_interpolated'])
    result_df.sort_values(['p_value'], ascending=True, inplace=True)
    result_df['q_value'] = -np.log10(multipletests(pvals=result_df['p_value'], method='fdr_bh', is_sorted=True)[1])

    if num_peaks is not False:
        result_df = result_df.head(num_peaks)
    elif q_thr is not False:
        result_df = result_df.loc[(result_df['q_value'] > -np.log10(q_thr)), :]
    elif (p_thr is not False) and (p_thr != 'False'):
        result_df = result_df.loc[(result_df['-log10_p_value_interpolated'] > -np.log10(float(p_thr))), :]

    return result_df


# @jit
def calculate_nbinom_param(mu, var):
    """Input (mu, var), output (p, r)
    """
    p = mu / var
    return p, p * mu / (1 - p)


# @jit
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


def read_chrom_file(file_name):
    with open(file_name, 'r') as file:
        chrom_ls_to_process = [line.strip() for line in file]

    return chrom_ls_to_process


def read_chrom_size_file(file_name):
    with open(file_name, 'r') as file:
        chr_size_dic = {}
        for line in file:
            line = line.strip().split()
            chr_size_dic[line[0]] = int(line[1])

    return chr_size_dic


def gaussian_kernel_fft(sig, sigma):
    lw = int(4 * sigma + 0.5)

    return fftconvolve(sig, filters._gaussian_kernel1d(sigma, 0, lw), mode='same')


def output_sig(sig_format, treatment_np_tmp_name, out_dir, out_name, chr_size_dic, sig_float_precision, gaussian_smooth_sigma=False):

    #start_time = time.time()

    with h5py.File(treatment_np_tmp_name, mode='r', libver='latest') as sig_file:
        chrom_line_ls = list(sig_file.keys())
        if sig_format == 'wig':
            with open(f'{out_dir}/{out_name}.wig', 'w') as output_wig:
                output_wig.write(f'track type=wiggle_0 name={out_name}\n')
                if gaussian_smooth_sigma:
                    for chrom_line in chrom_line_ls:
                        kdepy_result = np.round(gaussian_kernel_fft(sig_file[chrom_line][:], sigma=gaussian_smooth_sigma),
                                                 sig_float_precision).astype(np.float16)
                        output_wig.write(f'fixedStep chrom={chrom_line} start={sig_file.attrs[chrom_line]} step=1 span=1\n')
                        output_wig.write('\n'.join(kdepy_result.astype('str')) + '\n')
                else:
                    for chrom_line in chrom_line_ls:
                        kdepy_result = np.round(sig_file[chrom_line][:], sig_float_precision).astype(np.float16)
                        output_wig.write(f'fixedStep chrom={chrom_line} start={sig_file.attrs[chrom_line]} step=1 span=1\n')
                        output_wig.write('\n'.join(kdepy_result.astype('str')) + '\n')
        elif sig_format == 'bigwig':
            with pyBigWig.open(f'{out_dir}/{out_name}.bw', 'w') as output_bw:
                output_bw.addHeader([(chrom_line, chr_size_dic[chrom_line]) for chrom_line in chrom_line_ls])
                if gaussian_smooth_sigma:
                    for chrom_line in chrom_line_ls:
                        # kdepy_result = gaussian_kernel_fft(sig_file[chrom_line][:],
                        #                                    sigma=gaussian_smooth_sigma).astype(np.float64)
                        kdepy_result = np.round(gaussian_kernel_fft(sig_file[chrom_line][:], sigma=gaussian_smooth_sigma),
                                                 sig_float_precision).astype(np.float16)
                        #kdepy_result[kdepy_result == 0] = np.NaN
                        output_bw.addEntries(chrom_line, int(sig_file.attrs[chrom_line]), values=kdepy_result, span=1,
                                             step=1)
                else:
                    for chrom_line in chrom_line_ls:
                        #kdepy_result = sig_file[chrom_line][:].astype(np.float64)
                        kdepy_result = np.round(sig_file[chrom_line][:], sig_float_precision).astype(np.float16)
                        #kdepy_result[kdepy_result == 0] = np.NaN
                        output_bw.addEntries(chrom_line, int(sig_file.attrs[chrom_line]), values=kdepy_result, span=1,
                                             step=1)

    #end_time = time.time()
    #print(f'time = {end_time-start_time:.3f} seconds.')

    return


def narrowPeak_writer(result_df, peak_type, name, out_dir, prior_pad_summit=0, score_method_for_peaks='max',
                      sort_by='pValue'):
    """Write summits or peaks in narrowPeak Format.
    Format Info: https://genome.ucsc.edu/FAQ/FAQformat.html

    Note: signalValue is not the value accounted for local bias.
        Instead use score which is -log10(p_value_interpolated) * 10.

    Input:
        result_df: (pd.Dataframe) required columns ['chrom', 'summit', 'summit_end', 'query_value', '-log10_p_value_interpolated', 'q_value']

    """
    result_df = result_df.round({'query_value': 3, '-log10_p_value_interpolated': 3, 'q_value': 3})

    if peak_type == 'summit':

        result_df = result_df.loc[:,
                    ['chrom', 'summit', 'summit_end', 'query_value', '-log10_p_value_interpolated', 'q_value']]

        result_df['summit'] = result_df['summit'] - int(prior_pad_summit / 2)
        result_df['summit_end'] = result_df['summit_end'] + int(prior_pad_summit / 2)
        result_df['peak'] = int(prior_pad_summit / 2)

        result_df.rename(columns={'summit': 'chromStart', 'summit_end': 'chromEnd',
                                  'query_value': 'signalValue', '-log10_p_value_interpolated': 'pValue',
                                  'q_value': 'qValue'}, inplace=True)

        if sort_by == 'pValue':
            result_df.sort_values(['pValue'], ascending=False, inplace=True)
        elif sort_by == 'chromAndStart':
            result_df.sort_values(['chrom', 'chromStart'], inplace=True)

        result_df['name'] = [name + '_summit_' + str(i) for i in range(result_df.shape[0])]
        result_df['score'] = (result_df['pValue'] * 10).astype(int)
        result_df['strand'] = '.'

        result_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                   'signalValue', 'pValue', 'qValue', 'peak']].to_csv(f'{out_dir}/{name}_summits.narrowPeak', sep='\t',
                                                                      header=None, index=None)

    elif peak_type == 'peak':

        result_df = result_df.loc[:,
                    ['chrom', 'start', 'end', 'query_value', '-log10_p_value_interpolated', 'q_value', 'summit']]

        result_df.sort_values(['chrom', 'start'], inplace=True)
        result_df['summit'] = result_df['summit'] - result_df['start']
        try:
            result_df = pybedtools.BedTool.from_dataframe(result_df).merge(c=[5, 4, 5, 6, 7],
                                                                           o=[score_method_for_peaks, 'collapse',
                                                                              'collapse', 'collapse',
                                                                              'collapse']).to_dataframe()
        except:
            if result_df.empty:
                result_df.to_csv(f'{out_dir}/{name}_peaks.narrowPeak', sep='\t', header=None, index=None)
                return

        if sort_by == 'pValue':
            result_df.sort_values(['strand'], ascending=False, inplace=True)
        elif sort_by == 'chromAndStart':
            result_df.sort_values(['chrom', 'start'], inplace=True)

        result_df['assigned_name'] = [name + '_peak_' + str(i) for i in range(result_df.shape[0])]
        result_df['assigned_strand'] = '.'
        result_df['name'] = (result_df['name'] * 10).astype(int)

        result_df[['chrom', 'start', 'end', 'assigned_name', 'name', 'assigned_strand',
                   'score', 'strand', 'thickStart', 'thickEnd']].to_csv(f'{out_dir}/{name}_peaks.narrowPeak', sep='\t',
                                                                        header=None, index=None)

    return


def bed_chunking(df, df_control, chrom_ls, params):
    """Chunk bed into a list
    [[cuts_array, start_array, end_array, strand_array,
    cuts_array_control, start_array_control, end_array_control, strand_array_control, chrom], ...]
    """
    chunked_ls = []
    if df_control is not False:
        for chrom in chrom_ls:
            mask = (df['chrom'] == chrom)
            mask_control = (df_control['chrom'] == chrom)
            chunked_ls.append([df.loc[mask, 'cuts'].to_numpy(dtype=np.int32, copy=True),
                               df.loc[mask, 'start'].to_numpy(dtype=np.int32, copy=True),
                               df.loc[mask, 'end'].to_numpy(dtype=np.int32, copy=True),
                               df.loc[mask, 'strand'].to_numpy(copy=True) == '+',
                               df_control.loc[mask_control, 'cuts'].to_numpy(dtype=np.int32, copy=True),
                               df_control.loc[mask_control, 'start'].to_numpy(dtype=np.int32, copy=True),
                               df_control.loc[mask_control, 'end'].to_numpy(dtype=np.int32, copy=True),
                               df_control.loc[mask_control, 'strand'].to_numpy(copy=True) == '+',
                               chrom, params])
            df.drop(df.index[mask], axis=0, inplace=True)
            df_control.drop(df_control.index[mask_control], axis=0, inplace=True)
    else:
        for chrom in chrom_ls:
            mask = (df['chrom'] == chrom)
            chunked_ls.append([df.loc[mask, 'cuts'].to_numpy(dtype=np.int32, copy=True),
                               df.loc[mask, 'start'].to_numpy(dtype=np.int32, copy=True),
                               df.loc[mask, 'end'].to_numpy(dtype=np.int32, copy=True),
                               df.loc[mask, 'strand'].to_numpy(copy=True) == '+',
                               chrom, params])
            df.drop(df.index[mask], axis=0, inplace=True)

    return chunked_ls

