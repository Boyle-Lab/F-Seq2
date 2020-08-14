# F-Seq2: a feature density estimator for high-throughput sequence tags

Tag sequencing using high-throughput sequencing technologies are employed to identify specific sequence features such as 
DNase-seq, ATAC-seq, ChIP-seq, and FAIRE-seq. To intuitively summarize and display individual sequence data as an 
accurate and interpretable signal, we have developed the original [F-Seq](http://fureylab.web.unc.edu/software/fseq/) 
[GitHub](https://github.com/aboyle/F-seq), a software package that generates a continuous tag sequence density 
estimation allowing identification of biologically meaningful sites whose output can be displayed directly in the UCSC 
Genome Browser. 

F-Seq2 is a complete rewrite and improvement of the original version in Python. F-Seq2 implements a dynamic parameter to 
conduct statistical analysis in “continuous” Poisson distribution. Combining the dynamic test with the 
statistically rigorous read distribution modeled by KDE allows us to more robustly account for local bias and solve ties 
when ranking candidate summits. F-Seq2 reconstructs a more biologically accurate signal by weighting the reads with 
user-input control data. A higher quality signal leads to better sanity-check data, compare and combine signals 
from different datasets.



## Table of contents

1. [Installation](./INSTALL.md)
2. [Usage](#usage)
    - [`callpeak`](#callpeak)
    - [`callpeak_idr`](#callpeak_idr)
    - [`idr`](#idr)
3. [Output files and formats](#output-files-and-formats)
4. [Examples](#examples)
5. [Reference](#reference)
6. [Troubleshooting](#troubleshooting)



## Installation
See [here](./INSTALL.md).



## Usage

```
fseq2 [-h] [--version]
    {callpeak, callpeak_idr, idr}
```

Available subcommands

Subcommand | Description
-----------|----------
`callpeak` | F-Seq2 main function to call peaks from alignment results.
`callpeak_idr` | Call peaks and follow by IDR framework with recommended settings.
`idr` | A wrapper for [IDR package](https://github.com/nboley/idr) for customized IDR analysis.



## `callpeak`
#### Command line input:
##### `-treatment_file`
REQUIRED argument for fseq2. Treatment file(s) in bam or bed format. If specifiy multiple files (separated by space), 
they are considered as one treatment experiment. See [here](./INPUT_FORMAT.md) for more details about input format.

##### `-control_file`
Control file(s) corresponding to treatment file(s).

##### `-pe`
Paired-end mode. If this flag on, treatment (and control) file(s) are paired-end data, either in format of BAMPE or BEDPE. 
Default is False to treat all data as single-end. See [here](./INPUT_FORMAT.md) for more details about paired-end mode.

##### `-chrom_size_file` 
A file specify chrom sizes, where each line has one chrom and its size. This is required if output signal format is bigwig. 
Note if this file is specified, fseq2 only process the chroms in this file. Default is False to process all and cannot output bigwig.

##### `-o`
Output directory. Default is current directory.

##### `-name`
Prefix for all output files. This overrides exisiting files. Default is `fseq2_result`.

##### `-sig_format`
Signal format for reconstructed signal. Available format `wig`, `bigwig`, `np_array`. Note if choose `np_array`, arrays 
for each chrom are stored in [`NAME_sig.h5`](#name_sigh5) with `chrom` as key, and no gaussian smooth applied. Default is False, without output signal.

##### `-sort_by`
Sort peaks and summits by `pValue` or `chromAndStart`. Default is `chromAndStart`.

##### `-v`
Verbose output. Default is False.  

##### `-f`
Fragment size of treatment data. Default is to estimate from data. This determines shift size where `offset = fragment_size/2`. 
For DNase-seq and ATAC-seq data, set `-f 0`. 

##### `-l`
Feature length for treatment data. Default is 600. Recommend 50 for TF ChIP-seq, 600 for DNase-seq and ATAC-seq, 
1000 for histone ChIP-seq.

##### `-fc`
Fragment size of control data.

##### `-t`
Threshold (standard deviations) to call candidate summits. Default is 4.0. Recommend 4.0 for broad peaks, 
8.0 for sharp peaks.

##### `-p_thr`
P value threshold. Default is 0.01.

##### `-q_thr`
Q value (FDR) threshold. Default is not set and use `p_thr`. If set, only use `q_thr`.

##### `-cpus`
Number of cpus to use. Default is 1.

##### `-nfr_upper_limit`
Nucleosome free region upper limit. Default is 150. Used as window_size and min_distance when `-f 0`.

##### `-pe_fragment_size_range`
Effective only if `-pe` on. Keep PE fragments whose size within the range. Default is False, without any selection. Useful for ATAC-seq data:  
(1) to call peaks on nucleosome free regions, specify: 0 150  
(2) to call peaks on nucleosome centers, specify: 150 inf

##### `-prior_pad_summit`
Prior knowledge about peak length which only padded into `NAME_summits.narrowPeak`. Default is 0. 
Useful for IDR analysis: in `callpeak_idr`, we set it to est. fragment size.

##### `-num_peaks`
Maximum number of peaks called. Default is not set. If set, overrides `p_thr` and `q_thr`.



## `callpeak_idr`
#### Command line input:
Most arguments are shared between `callpeak` and `callpeak_idr`. Here are the unique ones.  
> Notice if it is `-` or `--` ahead of arguments. `--` arguments are from IDR package. `-` are from fseq2.
##### `-treatment_file_1`
Treatment file in bam or bed format as replicate 1.

##### `-treatment_file_2`
Treatment file in bam or bed format as replicate 2.

##### `-control_file_1`
Control file in bam or bed format, paired with replicate 1 treament file.

##### `-control_file_2`
Control file in bam or bed format, paired with replicate 2 treament file.

##### `-name_1`
Prefix for output files for replicate 1 (default=`fseq2_result_1`).

##### `-name_2`
Prefix for output files for replicate 2 (default=`fseq2_result_2`).

##### `-prior_pad_summit`
Prior knowledge about peak length which only padded into `NAME_summits.narrowPeak`. Default is est. fragment size.

##### `--idr_threshold`
Only return peaks with a global idr threshold below this value. Default: report all peaks.

##### `--soft_idr_threshold`
Report statistics for peaks with a global idr below this value but return all peaks with an idr below --idr Default: 0.05.

##### `--plot`
Plot IDR results. Specify False if no plot. Default is to plot to `NAME_1_NAME_2.png`. Can specify other name here. 
Notice this is different from original IDR package which is only a flag.



## `idr`
#### Command line input and output:
See original [IDR documentation](https://github.com/nboley/idr#usage).  
> Notice all single letter arguments are removed to avoid conflict with fseq2, e.g. no `-s`, use `--samples`



## Output files and formats
#### `NAME_summits.narrowPeak` 
BED6+4 format
1. chrom
2. chromStart 
3. chromEnd 
4. name - `NAME_summit_num`, num is sorted by either `Pvalue` or `chromAndStart`.
5. score - `int(10*-log10(pValue))`.
6. strand - `.`
7. signalValue - Average treatment signal value given window size.
8. pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
9. qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
10. peak - 0 if no specification of `-prior_pad_summit`.
  
#### `NAME_peaks.narrowPeak` 
Similar to summit file except that it can contain multiple summits information. 
For 7-10 columns, if multiple summits in a peak, output a comma separated list for each column.
1. chrom
2. chromStart 
3. chromEnd 
4. name - `NAME_peak_num`, num is sorted by either `Pvalue` or `chromAndStart`.
5. score - Max `int(10*-log10(pValue))` of all summits.
6. strand - `.`
7. signalValue
8. pValue
9. qValue
10. peak - Relative summit position(s) to peak start.

#### `NAME.bw` and `NAME.wig`
Reconstructed signal files which can be displayed directly in the UCSC Genome Browser. 
Recommend `bw` for efficient indexing in the browser.  

#### `NAME_sig.h5`
Reconstructed signal file without any smoothing. Signal is stored for each chrom in `np.array` and accessed by key `chrom`.   
For example:
```
>>> with h5py.File(NAME_sig.h5, mode='r') as sig_file:
...     signal = sig_file['chr1'][:] # read in all signal on chr1
```

#### `NAME_1_NAME_2_conservative_IDR_thresholded_peaks.narrowPeak` and `NAME_1_NAME_2.png`
Generated by `fseq2 callpeak_idr`.    Detailed format information is [here](https://github.com/nboley/idr#output).



## Examples

#### DNase-seq data
```
$ fseq2 callpeak treatment_file.bam -f 0 -l 600 -t 4.0 -v -cpus 5
```

#### ATAC-seq data
Paired-end ATAC-seq data, and call peaks on open chromatin regions, without calling on nucleosomes
```
$ fseq2 callpeak treatment_file.bam -f 0 -l 600 -t 4.0 -pe -nfr_upper_limit 150 -pe_fragment_size_range 0 150
```

#### ChIP-seq data
TF ChIP-seq data
```
$ fseq2 callpeak treatment_file.bed -control_file control_file.bed -l 50 -t 8.0 -sig_format bigwig -chrom_size_file ./tests/integration/fixtures/hg19.chrom.sizes  
```



## Troubleshooting

##### 1. Install error on mac Mojave: 
```
fatal error: 'ios' file not found 
#include "ios"
```
Solution:  
add `CFLAGS='-stdlib=libc++'` in front of `pip install`
```
$ CFLAGS='-stdlib=libc++' pip install fseq2
```

##### 2. Memory error

Solution:  
try with less CPUs


##### 3. `NotImplementedError: "bamToBed" does not appear to be installed`

Solution:  
update bedtools >= 2.29.0


##### 4. Warnings when `-pe`
Mostly likely bam file is not sorted by name.  
Solution:  
see [here](./INPUT_FORMAT.md)

##### 5. Too few peaks after multi-test correction
This may indicate poor data quality.  
Solution:  
use `-p_thr` instead of `-q_thr`