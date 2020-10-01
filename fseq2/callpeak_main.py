#!/usr/bin/env python

"""F-Seq Version2 call peak main script.

"""
import multiprocessing as mp
import sys

from numpy import ptp
from pandas import concat, DataFrame

from fseq2 import fseq2


def main(args):
    ### 1.0 chrom to process ###
    if args.v:
        print('=====================================', flush=True)
        print(f'F-Seq Version {fseq2.__version__}', flush=True)
        print('=====================================', flush=True)
        print('#1: Read in files and calculate parameters', flush=True)

    if args.chrom_size_file:
        chr_size_dic = fseq2.read_chrom_size_file(args.chrom_size_file)
        chrom_ls_to_process = list(chr_size_dic)
    else:
        chrom_ls_to_process = False
        if (args.sig_format == 'bigwig'):
            sys.exit('Error: please specify chrom_size_file to output bigwig signal.')


    ### 1.1 deal with input (treatment) ###
    input_bed, fragment_size_ls = fseq2.read_in_file(args.treatment_file, args.pe, args.pe_fragment_size_range,
                                                     args.nucleosome_size, chrom_ls_to_process)

    # parameters for input
    feature_length = args.l
    if args.pe:
        fragment_size = 0
    elif args.f is False:
        fragment_size = fseq2.calculate_fragment_size(
            input_bed[input_bed['chrom'] == input_bed.groupby('chrom')['chrom'].count().idxmax()], verbose=False)
    else:
        fragment_size = args.f
    fragment_offset = int(fragment_size / 2)
    bandwidth = fseq2.calculate_bandwidth(feature_length)
    ncuts = input_bed.shape[0]

    threshold, peak_region_threshold = fseq2.calculate_threshold(size=input_bed.groupby('chrom')['cuts'].agg(ptp).sum(),
                                                           ncuts=ncuts, std=args.t, std_p=args.tp, bandwidth=bandwidth)  # 7.27
    sequence_length = input_bed.iloc[0, 2] - input_bed.iloc[0, 1]


    ### 1.2 deal with control ###
    if args.control_file:
        control_bed, fragment_size_ls_control = fseq2.read_in_file(args.control_file, args.pe, args.pe_fragment_size_range,
                                                                   args.nucleosome_size, chrom_ls_to_process)

        # parameters for control
        feature_length_control = feature_length * 20  # 7.21
        if args.pe:
            fragment_size_control = 0
        elif args.fc is False:
            fragment_size_control = fseq2.calculate_fragment_size(
                control_bed[control_bed['chrom'] == control_bed.groupby('chrom')['chrom'].count().idxmax()],
                verbose=args.v)
        else:
            fragment_size_control = args.fc
        fragment_offset_control = int(fragment_size_control / 2)
        bandwidth_control = fseq2.calculate_bandwidth(feature_length_control)
        ncuts_control = control_bed.shape[0]

        # threshold_control = calculate_threshold(size=control_bed.groupby('chrom')['cuts'].agg(np.ptp).sum(), ncuts=ncuts_control,
        # std=args.t, bandwidth=bandwidth_control)
        sequence_length_control = control_bed.iloc[0, 2] - control_bed.iloc[0, 1]

        if (args.control_file[0][-3:] == 'bam') and (args.treatment_file[0][-3:] == 'bam'):
            assert input_bed.iloc[-1, 0] in control_bed.iloc[:, 0].values, "pybedtools incomplete write-in error, " \
                                                                           "please free memory and try again."
            # assert control_bed.iloc[-1, 0] == input_bed.iloc[-1, 0], "pybedtools incomplete write-in error, " \
            # "please free memory and try again."


    ### 1.3 general params ###
    scaling_factor = feature_length * 20 / 50
    threshold *= scaling_factor
    peak_region_threshold *= scaling_factor  # 7.27
    # threshold_control *= scaling_factor
    if args.control_file:
        lambda_bg_lower_bound = fseq2.calculate_lambda_bg_lower_bound(bandwidth_control=bandwidth_control,
                                                                ncuts_control=ncuts_control, ncuts_treatment=ncuts,
                                                                pseudo_count=1)
    else:
        lambda_bg_lower_bound = 0
        # lambda_bg_lower_bound = calculate_lambda_bg_lower_bound(bandwidth_control=bandwidth,
        # ncuts_control=ncuts, ncuts_treatment=ncuts,
        # pseudo_count=1)
    lambda_bg_lower_bound *= scaling_factor
    min_prominence = threshold
    if fragment_size == 0:
        min_distance = args.nfr_upper_limit
        window_size = args.nfr_upper_limit
    else:
        min_distance = fragment_size
        window_size = fragment_size
    distribution_name = 'poisson'
    if args.sig_format:
        sig_float_precision = 2
        gaussian_smooth_sigma = int(1000 / 6)
        if args.sig_format == 'np_array':
            treatment_np_tmp_name = f'{args.o}/{args.name}_sig.h5'
        elif (args.sig_format == 'wig') or (args.sig_format == 'bigwig'):
            treatment_np_tmp_name = f'{args.temp_dir_name}/{args.name}_sig.h5'


    if args.v:
        print('#1: Done', flush=True)
        print('#1: Settings:\n', flush=True)
        print(f'\tbandwidth = {bandwidth:.3f}', flush=True)
        print(f'\tthreshold = {threshold:.3f}', flush=True)
        print(f'\tlambda bg lower bound = {lambda_bg_lower_bound:.3f}', flush=True)
        print(f'\test. fragment size = {fragment_size}', flush=True)
        if args.pe:
            print(f'\tavg. pe fragment size before filter = {fragment_size_ls[0]:.3f}', flush=True)
            print(f'\tavg. pe fragment size after filter = {fragment_size_ls[1]:.3f}', flush=True)
            print(f'\tfilter out rate  = {fragment_size_ls[2]:.3f}', flush=True)
        print(f'\tsequence length = {sequence_length}', flush=True)
        print(f'\ttotal cuts = {ncuts}', flush=True)
        if args.control_file:
            print(f'\n\tcontrol bandwidth = {bandwidth_control:.3f}', flush=True)
            print(f'\tcontrol est. fragment size = {fragment_size_control}', flush=True)
            if args.pe:
                print(f'\tcontrol avg. pe fragment size before filter = {fragment_size_ls_control[0]:.3f}', flush=True)
                print(f'\tcontrol avg. pe fragment size after filter = {fragment_size_ls_control[1]:.3f}', flush=True)
                print(f'\tcontrol filter out rate  = {fragment_size_ls_control[2]:.3f}', flush=True)
            print(f'\tcontrol sequence length = {sequence_length_control}', flush=True)
            print(f'\tcontrol total cuts = {ncuts_control}', flush=True)
        print('\n-------------------------------------', flush=True)

    # params into run_kde() or run_kde_wo_control()
    args.fragment_offset = fragment_offset
    args.bandwidth = bandwidth

    if args.control_file:
        args.fragment_offset_control = fragment_offset_control
        args.bandwidth_control = bandwidth_control

    args.scaling_factor = scaling_factor
    args.threshold = threshold
    args.min_distance = min_distance
    args.min_prominence = min_prominence
    args.peak_region_threshold = peak_region_threshold
    args.window_size = window_size
    args.distribution_name = distribution_name
    args.lambda_bg_lower_bound = lambda_bg_lower_bound

    if args.sig_format:
        args.treatment_np_tmp_name = treatment_np_tmp_name
        args.sig_float_precision = sig_float_precision


    ### 2.1 run kde ###
    if args.v:
        print('#2: Reconstruct signal and call peaks\n', flush=True)

    chrom_ls = input_bed['chrom'].unique()
    if args.control_file:
        input_param_ls = fseq2.bed_chunking(df=input_bed, df_control=control_bed, chrom_ls=chrom_ls, params=args)
        input_bed = DataFrame()
        control_bed = DataFrame()
    else:
        input_param_ls = fseq2.bed_chunking(df=input_bed, df_control=False, chrom_ls=chrom_ls, params=args)
        input_bed = DataFrame()

    cpus = args.cpus  # cpus = min(chrom_ls.shape[0], fseq.mp.cpu_count() - 2)
    ctx = mp.get_context('spawn')
    lock = ctx.Lock()

    with ctx.Pool(processes=cpus, initializer=fseq2.init_child, initargs=(lock,)) as pool:
        if args.control_file:
            results = pool.starmap(fseq2.run_kde, input_param_ls)
        else:
            results = pool.starmap(fseq2.run_kde_wo_control, input_param_ls)
    result_df = concat(results)

    ### 2.2 calculate p-value and q-value ###
    if args.v:
        print('\n#2: Compute p-value and q-value', flush=True)

    result_df = fseq2.interpolate_poisson_p_value(result_df)
    result_df = fseq2.calculate_q_value(result_df=result_df, p_thr=args.p_thr, q_thr=args.q_thr,
                                        num_peaks=args.num_peaks)
    if args.v:
        print('#2: Done', flush=True)
        print('-------------------------------------', flush=True)


    ### 3. prepare and output ###
    if args.v:
        print('#3: Write output', flush=True)

    #result_df.to_hdf(path_or_buf=f'{args.o}/{args.name}.h5', key='result_df', mode='w')  # for inspection
    if args.prior_pad_summit == 'fragment_size':
        args.prior_pad_summit = fragment_size
    fseq2.narrowPeak_writer(result_df=result_df, peak_type='summit', name=args.name, out_dir=args.o,
                      prior_pad_summit=args.prior_pad_summit, sort_by=args.sort_by)
    fseq2.narrowPeak_writer(result_df=result_df, peak_type='peak', name=args.name, out_dir=args.o, sort_by=args.sort_by)

    if args.sig_format:
        fseq2.output_sig(sig_format=args.sig_format, treatment_np_tmp_name=treatment_np_tmp_name,
                   out_dir=args.o, out_name=args.name, chr_size_dic=chr_size_dic,
                   sig_float_precision=sig_float_precision, gaussian_smooth_sigma=gaussian_smooth_sigma) #note if output np_array, then no gaussian smooth

    if args.v:
        print('#3: Done', flush=True)
        print('-------------------------------------', flush=True)
        print(f'Thanks for using F-seq{fseq2.__version__}!\n', flush=True)


    return
