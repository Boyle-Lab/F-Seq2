#!/usr/bin/env python

"""F-Seq Version2 call peak and idr main script.

"""
import sys

from fseq2 import idr_2_0_3 as idr
from fseq2.callpeak_main import main as callpeak_main
from fseq2.fseq2 import __version__ as fseq_version
from fseq2.idr_2_0_3 import __version__ as idr_version
from fseq2.idr_2_0_3.idr import PossiblyGzippedFile
from fseq2.idr_main import main as idr_main


def main(args):

    ### 1. sample 1 setup and run fseq2 ###
    if args.v:
        print(f'Workflow#1: fseq{fseq_version} callpeak with treatment_file_1, control_file_0', flush=True)
    args.treatment_file = [args.treatment_file_1]
    args.control_file = [args.control_file_1, args.control_file_2]
    args.name = args.name_1
    args.sort_by = 'pValue'
    args.num_peaks = 300000
    args.p_thr = False
    args.q_thr = False

    callpeak_main(args)


    ### 2. sample 2 setup and run fseq2 ###
    if args.v:
        print(f'Workflow#2: fseq{fseq_version} callpeak with treatment_file_2, control_file_0', flush=True)
    args.treatment_file = [args.treatment_file_2]
    args.name = args.name_2

    callpeak_main(args)


    ### 3. sample 0 setup and run fseq2 ###
    if args.v:
        print(f'Workflow#3: fseq{fseq_version} callpeak with treatment_file_0, control_file_0', flush=True)
    args.treatment_file = [args.treatment_file_1, args.treatment_file_2]
    args.name = args.name_1 + '_' + args.name_2

    callpeak_main(args)


    ### 4. idr setup and run idr ###
    if args.v:
        print(f'Workflow#4: idr{idr_version} with samples(1,2), oracle peak list(3)', flush=True)
    args = prepare_idr_arg(args)

    idr_main(args)


    return


def prepare_idr_arg(args):
    """Recommended settings for IDR with F-Seq2 results.
    """

    args.samples = [f'{args.o}/{args.name_1}_summits.narrowPeak',
                    f'{args.o}/{args.name_2}_summits.narrowPeak']
    args.samples = [PossiblyGzippedFile(filename) for filename in args.samples]
    args.peak_list = PossiblyGzippedFile(f'{args.o}/{args.name_1}_{args.name_2}_summits.narrowPeak')
    args.input_file_type = 'narrowPeak'
    args.rank = 'p.value'
    args.output_file = f'{args.name_1}_{args.name_2}_conservative_IDR_thresholded_peaks.narrowPeak'
    args.output_file_type = 'narrowPeak'
    if args.plot:
        args.plot = f'{args.name_1}_{args.name_2}.png'
    args.verbose = False
    args.quiet = not args.v

    # default
    args.log_output_file = sys.stderr
    args.use_old_output_format = False
    args.use_nonoverlapping_peaks = False
    args.peak_merge_method = None
    args.initial_mu = float(idr.DEFAULT_MU)
    args.initial_sigma = float(idr.DEFAULT_SIGMA)
    args.initial_rho = float(idr.DEFAULT_RHO)
    args.initial_mix_param = float(idr.DEFAULT_MIX_PARAM)
    args.fix_mu = False
    args.fix_sigma = False
    args.dont_filter_peaks_below_noise_mean = False
    args.use_best_multisummit_IDR = False
    args.allow_negative_scores = False
    args.random_seed = 0
    args.max_iter = int(idr.MAX_ITER_DEFAULT)
    args.convergence_eps = float(idr.CONVERGENCE_EPS_DEFAULT)
    args.only_merge_peaks = False
    args.version = f"IDR {idr.__version__}"

    return args
