#!urs/bin/env python3

"""F-Seq Version 2.0 Prototype

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
import pyBigWig
from itertools import chain
import multiprocessing as mp
import time


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
        print('Calculating fragment size (this may take a long time, specify -f 0 if using DNase HS data) ...',
              flush=True)

    fragment_size = int(np.median(np.concatenate([cuts_ls[(cuts_ls > (current_cuts - max_range + prior)) &
                                                          (cuts_ls < (current_cuts + max_range + prior)) &
                                                          neg_strand_ls] - current_cuts
                                                  for current_cuts in pos_cuts_ls], axis=0)))

    if verbose:
        print('Done!', flush=True)

    return fragment_size


def calculate_bandwidth(feature_length):
    """Calculate bandwidth.
    """
    bandwidth = float(feature_length / 2 / 3) # 3 standard deviations

    return bandwidth


def calculate_threshold(size, ncuts, std, window_size=int(2860/2)):
    """Calculate threshold for peak calling.

    Input:
        size: sum range of cuts across all chroms in bp
        ncuts: total number of cuts
        std: standard deviation for calculating threshold
        window_size: window size in bp

    Return:
        threshold: threshold for peak calling
    """
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
    kdepy_result = (np.array(simulated_results) * 20000000) / ncuts * cut_density

    threshold = np.mean(kdepy_result) + std * np.std(kdepy_result)

    return threshold


def run_kde(cuts_df, chrom, bf, pf, fragment_offset, bandwidth, verbose=False):
    """Estimate KDE.

    Return:
        [chrom, first_cut, kde_results]: a list contains chrom name, first cut position, and kde results from first cut
            to last cut in np.array
    """
    start_time = time.time()
    if bf:
        with pyBigWig.open(bf) as background_bw:
            background_bw_weights = np.array(list(chain.from_iterable(
                [background_bw.values(chrom, start, start + 1) * duplicate
                 for start, duplicate in zip(np.unique(cuts_df['start'].values),
                                             cuts_df.groupby('start')['start'].count().values)])))
            # bigwig is 0-index, wig is 1-index
        background_bw_weights[background_bw_weights == 0] = 1 / 5
        background_bw_weights = 1 / background_bw_weights
    if pf:
        with pyBigWig.open(pf) as ploidy_bw:
            ploidy_bw_weights = np.array(list(chain.from_iterable(
                [ploidy_bw.values(chrom, start, start + 1) * duplicate
                 for start, duplicate in zip(np.unique(cuts_df['start'].values),
                                             cuts_df.groupby('start')['start'].count().values)])))
        ploidy_bw_weights = ploidy_bw_weights / 1000

    if bf and pf:
        weights = background_bw_weights / ploidy_bw_weights
    elif bf:
        weights = background_bw_weights
    elif pf:
        weights = ploidy_bw_weights
    else:
        weights = np.ones(cuts_df.shape[0])

    cuts = cuts_df['cuts'].values
    first_cut = cuts[0]
    last_cut = cuts[-1]
    if fragment_offset == 0:
        kdepy_kde = KDEpy.FFTKDE(bw=bandwidth).fit(cuts, weights=weights)
        kdepy_result = kdepy_kde.evaluate(np.arange(first_cut - 1, last_cut + 2))[1:-2]
        kdepy_result = (kdepy_result * 20000000) / ncuts * (weights.sum())
    else:
        cuts_w_offset = ((cuts_df['start'] + fragment_offset) * (cuts_df['strand'] == '+') +
                         (cuts_df['end'] - fragment_offset) * (cuts_df['strand'] == '-')).astype('int').values
        kde_evaluation_start = min(cuts_w_offset[0], first_cut); kde_evaluation_end = max(cuts_w_offset[-1], last_cut)
        assert (first_cut - kde_evaluation_start) >= 0
        kdepy_kde = KDEpy.FFTKDE(bw=bandwidth).fit(cuts_w_offset, weights=weights)
        kdepy_result = kdepy_kde.evaluate(np.arange(kde_evaluation_start - 1, kde_evaluation_end + 2))[1:-2]
        kdepy_result = kdepy_result[(first_cut - kde_evaluation_start):(last_cut - kde_evaluation_start)]
        kdepy_result = (kdepy_result * 20000000) / ncuts * (weights.sum()) / 2
        # this scaling is simulating results with offset (left window positive cuts and right window negative cuts) by
        # multiplying all signal with 0.5
    end_time = time.time()

    if verbose:
        print(f'{chrom}: first={first_cut}, last={last_cut}', flush=True)
        print(f'{chrom}: Completed in {end_time - start_time:.3f} seconds.', flush=True)
        print('-------------------------------------')

    # todo:
    #   1. F-seq 1.85 timing includes output wig. Output bw takes longer time and we do not measure it here.
    #   We also need to measure the time for calling peaks.

    return [chrom, first_cut, kdepy_result]


if __name__ == '__main__':
    # command line input settings
    parser = argparse.ArgumentParser(
        description='This command executes fseq which is a feature density estimator for high-throughput sequence tags')
    parser.add_argument('input_file', help='input file(s) in bam or bed format. If multiple input files, separate '
                                           'by space; if one input file, it can contain all chromosomes', nargs='+')
    parser.add_argument('-b', help='background directory (default=none)')
    parser.add_argument('-bf', help='background bigwig file for all chromosomes (default=none)')
    parser.add_argument('-f', help="fragment size (default=estimated from data). For DNase HS data (5' ends)"
                                   " set -f 0", default=False, type=int)
    parser.add_argument('-l', help='feature length (default=600)', default=600, type=int)
    parser.add_argument('-o', help='output directory (default=current directory)')
    parser.add_argument('-of', help='output format bigwig, wig, bed, or npf (default bigwig)')
    parser.add_argument('-p', help='ploidy/input directory (default=none)')
    parser.add_argument('-pf', help='ploidy/input bigwig file for all chromosome (default=none)')
    parser.add_argument('-s', help='wiggle track step (default=1)', default=1, type=int)
    parser.add_argument('-t', help='threshold (standard deviations) (default=4.0)', default=4.0, type=float)
    parser.add_argument('-v', help='verbose output', default=False, action='store_true')

    args = parser.parse_args()

    CHR_SIZE_DIC = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
                    'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
                    'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
                    'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
                    'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}  # hg19

    ### 1. Read in input files ###
    if args.input_file[0][-3:] == 'bam':
        input_bed = pd.concat([pybedtools.BedTool(file).bam_to_bed().to_dataframe() for file in args.input_file],
                              axis=0)
    elif args.input_file[0][-3:] == 'bed':
        input_bed = pd.concat(map(functools.partial(pd.read_csv, sep='\t', header=None), args.input_file),
                              axis=0).rename(columns={0:'chrom', 1:'start', 2:'end', 3: 'name', 4: 'score', 5:'strand'})
    else:
        sys.exit('Error: Input file format should be either .bam or .bed, please convert files and try again.')
    input_bed = input_bed.drop(['name', 'score'], axis=1)
    input_bed['cuts'] = input_bed['start'] * (input_bed['strand'] == '+') + \
                        input_bed['end'] * (input_bed['strand'] == '-')

    # TODO:
    #  1. check file type?
    #  2. check strand implicitly later
    #  3. maybe use dask to speed up calculations on dataframe, and solve if dataframe not fit in memory
    #  4. chrom size input for other genome build; this is only needed for bigwig output


    ### 2. Set up all parameters ###
    if args.f is False:
        fragment_size = calculate_fragment_size(
            input_bed[input_bed['chrom'] == input_bed.groupby('chrom')['chrom'].count().idxmax()], verbose=args.v)
    else:
        fragment_size = args.f
    fragment_offset = int(fragment_size / 2)
    feature_length = args.l
    bandwidth = calculate_bandwidth(feature_length)
    sequence_length = input_bed.iloc[0, 2] - input_bed.iloc[0, 1]
    ncuts = input_bed.shape[0]
    threshold = calculate_threshold(size=input_bed.groupby('chrom')['cuts'].agg(np.ptp).sum(), ncuts=ncuts, std=args.t)

    if args.v:
        print('=====================================', flush=True)
        print(f'F-Seq Version {__version__} Prototype', flush=True)
        print('=====================================', flush=True)
        print('Settings: ', flush=True)
        print(f'\tbandwidth = {bandwidth}', flush=True)
        print(f'\tthreshold = {threshold:.3f}', flush=True)
        print(f'\test. fragment size = {fragment_size}', flush=True)
        print(f'\tsequence length = {sequence_length}', flush=True)
        print('-------------------------------------', flush=True)

    # TODO:
    #   1. check if sequence(read)_length consistent


    ### 3. Parallel kde
    chrom_ls = input_bed['chrom'].unique()
    input_param_ls = [[input_bed[input_bed['chrom'] == chrom]] + [chrom, args.bf, args.pf, fragment_offset, bandwidth,
                                                                  args.v] for chrom in chrom_ls]
    cpus = min(chrom_ls.shape[0], mp.cpu_count() - 2)

    with mp.Pool(processes=cpus) as pool:
        results = pool.starmap(run_kde, input_param_ls)


    ### 4. Output results
    for result in results:
        with pyBigWig.open(args.o+'/'+result[0]+'.bw', 'w') as output_bw:
            output_bw.addHeader([(result[0], CHR_SIZE_DIC[result[0]])])
            output_bw.addEntries(result[0], result[1], values=result[2], span=1, step=1)
    if args.v:
        print('Thank you for using F-seq.')

    # TODO:
    #   1. specify step
    #   2. avoid overwrite files
