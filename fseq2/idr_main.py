#!/usr/bin/env python

"""F-Seq Version2 idr main script.

"""
from fseq2.idr_2_0_3.idr import *


def main(args):

    args = parse_args(args)

    # load and merge peaks
    merged_peaks, signal_type = load_samples(args)
    s1 = numpy.array([pk.signals[0] for pk in merged_peaks])
    s2 = numpy.array([pk.signals[1] for pk in merged_peaks])

    # build the ranks vector
    idr.log("Ranking peaks", 'VERBOSE')
    r1, r2 = build_rank_vectors(merged_peaks)

    if args.only_merge_peaks:
        localIDRs, IDRs = None, None
    else:
        if len(merged_peaks) < 20:
            error_msg = "Peak files must contain at least 20 peaks post-merge"
            error_msg += "\nHint: Merged peaks were written to the output file"
            write_results_to_file(
                merged_peaks, args.output_file,
                args.output_file_type, signal_type)
            raise ValueError(error_msg)

        localIDRs = fit_model_and_calc_local_idr(
            r1, r2,
            starting_point=(
                args.initial_mu, args.initial_sigma,
                args.initial_rho, args.initial_mix_param),
            max_iter=args.max_iter,
            convergence_eps=args.convergence_eps,
            fix_mu=args.fix_mu, fix_sigma=args.fix_sigma)

        if args.use_best_multisummit_IDR:
            localIDRs = correct_multi_summit_peak_IDR_values(
                localIDRs, merged_peaks)
        IDRs = calc_global_IDR(localIDRs)

        if args.plot:
            assert len(args.samples) == 2
            plot(args, [s1, s2], [r1, r2], IDRs)

    num_peaks_passing_thresh = write_results_to_file(
        merged_peaks,
        args.output_file,
        args.output_file_type,
        signal_type,
        localIDRs=localIDRs,
        IDRs=IDRs,
        max_allowed_idr=args.idr_threshold,
        soft_max_allowed_idr=args.soft_idr_threshold,
        useBackwardsCompatibleOutput=args.use_old_output_format)

    args.output_file.close()


def parse_args(args):

    args.output_file = open(f"{args.o}/{args.output_file}", "w")
    idr.log_ofp = args.log_output_file

    if args.output_file_type is None:
        if args.input_file_type in ('narrowPeak', 'broadPeak', 'bed'):
            args.output_file_type = args.input_file_type
        else:
            args.output_file_type = 'bed'

    if args.verbose:
        idr.VERBOSE = True

    global QUIET
    if args.quiet:
        idr.QUIET = True
        idr.VERBOSE = False

    if args.dont_filter_peaks_below_noise_mean is True:
        idr.FILTER_PEAKS_BELOW_NOISE_MEAN = False

    if args.allow_negative_scores is True:
        idr.ONLY_ALLOW_NON_NEGATIVE_VALUES = False

    assert idr.DEFAULT_IDR_THRESH == 1.0
    if args.idr_threshold == None and args.soft_idr_threshold == None:
        args.idr_threshold = idr.DEFAULT_IDR_THRESH
        args.soft_idr_threshold = idr.DEFAULT_SOFT_IDR_THRESH
    elif args.soft_idr_threshold == None:
        assert args.idr_threshold != None
        args.soft_idr_threshold = args.idr_threshold
    elif args.idr_threshold == None:
        assert args.soft_idr_threshold != None
        args.idr_threshold = idr.DEFAULT_IDR_THRESH

    numpy.random.seed(args.random_seed)

    if args.plot:
        try:
            import matplotlib
        except ImportError:
            idr.log("WARNING: matplotlib does not appear to be installed and " \
                    + "is required for plotting - turning plotting off.",
                    level="WARNING")
            args.plot = False

    return args
