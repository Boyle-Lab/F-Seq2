#!/usr/bin/env python

"""Tests for `fseq2` package mains."""

import filecmp
import shutil
import subprocess


def test_callpeak_idr():
    treatment_file_1 = './tests/integration/fixtures/treatment_1.bam'
    treatment_file_2 = './tests/integration/fixtures/treatment_2.bam'
    control_file_1 = './tests/integration/fixtures/control_1.bam'
    control_file_2 = './tests/integration/fixtures/control_2.bam'
    chrom_size_file = './tests/integration/fixtures/hg19.chrom.sizes'
    output_dir = './tests/integration/output'

    subprocess.run(
        f"fseq2 callpeak_idr {treatment_file_1} {treatment_file_2} -control_file_1 {control_file_1} -control_file_2 {control_file_2} "
        f"-sig_format wig -chrom_size_file {chrom_size_file} -l 50 -t 4.0 -o {output_dir} -v -cpus 2 "
        f"-name_1 result_1 -name_2 result_2", shell=True)

    # fseq2 callpeak_idr ./tests/integration/fixtures/treatment_1.bam ./tests/integration/fixtures/treatment_2.bam
    # -control_file_1 ./tests/integration/fixtures/control_1.bam -control_file_2 ./tests/integration/fixtures/control_2.bam
    # -sig_format wig -chrom_size_file ./tests/integration/fixtures/hg19.chrom.sizes -l 50 -t 4.0 -o ./tests/integration/output
    # -v -cpus 2 -name_1 result_1 -name_2 result_2

    assert 1


def test_compare_results():
    sig_file = './tests/integration/output/result_1_result_2.wig'
    bed_file = './tests/integration/output/result_1_result_2_conservative_IDR_thresholded_peaks.narrowPeak'


    # compare files
    assert filecmp.cmp(f1=sig_file, f2='./tests/integration/fixtures/result_1_result_2.wig', shallow=False) == True
    assert filecmp.cmp(f1=bed_file,
                       f2='./tests/integration/fixtures/result_1_result_2_conservative_IDR_thresholded_peaks.narrowPeak',
                       shallow=False) == True


def test_rmtree():
    output_dir = './tests/integration/output'
    shutil.rmtree(output_dir)

    assert 1

# @pytest.fixture
# def example_args(tmp_path):
#     parser = argparse.ArgumentParser()
#     args = parser.parse_args()
#
#     args.treatment_file_1 = './tests/integration/fixtures/treatment_1.bed'
#     args.treatment_file_2 = './tests/integration/fixtures/treatment_2.bed'
#     args.control_file_1 = './tests/integration/fixtures/control_1.bed'
#     args.control_file_2 = './tests/integration/fixtures/control_2.bed'
#
#     args.pe = False
#     args.chrom_size_file = './tests/integration/fixtures/hg19.chrom.sizes'
#     args.o = './tests/unit' # str(tmp_path.absolute())
#     args.name_1 = 'result_1'
#     args.name_2 = 'result_2'
#     args.sig_format = 'wig'
#     args.v = False
#
#     args.f = False
#     args.l = 50
#     args.fc = False
#     args.t = 4.0
#     args.cpus = 1
#
#     args.prior_pad_summit = 'fragment_size'
#     args.idr_threshold = None
#     args.soft_idr_threshold = None
#     args.plot = True
#
#     temp_dir = tempfile.TemporaryDirectory()
#     temp_dir_name = temp_dir.name
#     args.temp_dir_name = temp_dir_name
#
#     #fseq2 callpeak_idr ./BamFiles/treatment_1.bed ./BamFiles/treatment_2.bed -control_file_1
#     # ./BamFiles/control_1.bed -control_file_2 ./BamFiles/control_2.bed
#     # -l 50 -o ./Fseq2_output/ -name_1 result_1 -name_2 result_2 -sig_format wig
#     # -v -cpus 2 -chrom_size_file ./BamFiles/hg19.chrom.sizes
#
#     # args = argparse()
#     # args.o = str(tmp_path)
#     # args.sig_format = 'wig'
#     #args = '/Users/samuel/PycharmProjects/fseq2/fseq2/tests/integration/fixtures/my_file.txt'
#     #print(type(tmp_path))
#     # args = tmp_path.joinpath("my_file_1.txt")
#     # args.write_text('Hi!This is Sam!')
#     # args_str = args.as_posix()
#
#     return args

# def test_callpeak_idr_main():
#
#    parser = argparse.ArgumentParser()
#    args = parser.parse_args()
#
#    args.treatment_file_1 = './tests/integration/fixtures/treatment_1.bed'
#    args.treatment_file_2 = './tests/integration/fixtures/treatment_2.bed'
#    args.control_file_1 = './tests/integration/fixtures/control_1.bed'
#    args.control_file_2 = './tests/integration/fixtures/control_2.bed'
#
#    args.pe = False
#    args.chrom_size_file = './tests/integration/fixtures/hg19.chrom.sizes'
#    args.o = './tests/unit' # str(tmp_path.absolute())
#    args.name_1 = 'result_1'
#    args.name_2 = 'result_2'
#    args.sig_format = 'wig'
#    args.v = False
#
#    args.f = False
#    args.l = 50
#    args.fc = False
#    args.t = 4.0
#    args.cpus = 1
#
#    args.prior_pad_summit = 'fragment_size'
#    args.idr_threshold = None
#    args.soft_idr_threshold = None
#    args.plot = True
#
#    temp_dir = tempfile.TemporaryDirectory()
#    temp_dir_name = temp_dir.name
#    args.temp_dir_name = args.o
#
#    # run main
#    submain(args)
#    #with pytest.raises(ValueError, match="Peak files must contain at least 20 peaks post-merge"):
#
#
#    assert 1
