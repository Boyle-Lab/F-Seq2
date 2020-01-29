# Fseq2
F-seq version 2

## Usage
```
fseq.py [-h] [-b B] [-bf BF] [-f F] [-l L] [-o O] [-of OF] [-p P]
               [-pf PF] [-s S] [-t T] [-v]
               input_file [input_file ...]

This command executes fseq which is a feature density estimator for high-
throughput sequence tags

positional arguments:
  input_file  input file(s) in bam or bed format. If multiple input files,
              separate by space; if one input file, it can contain all
              chromosomes

optional arguments:
  -h, --help  show this help message and exit
  -b B        background directory (default=none)
  -bf BF      background bigwig file for all chromosomes (default=none)
  -f F        fragment size (default=estimated from data). For DNase HS data
              (5' ends) set -f 0
  -l L        feature length (default=600)
  -o O        output directory (default=current directory)
  -of OF      output format bigwig, wig, bed, or npf (default bigwig)
  -p P        ploidy/input directory (default=none)
  -pf PF      ploidy/input bigwig file for all chromosome (default=none)
  -s S        wiggle track step (default=1)
  -t T        threshold (standard deviations) (default=4.0)
  -v          verbose output
  ```
Example  
```
python3 fseq.py Duke_GM12878_Chr21.bed Duke_GM12878_Chr22.bed -bf wgEncodeDukeMapabilityUniqueness20bp.bigWig -f 0 -v
```
example input data is available [here](https://www.dropbox.com/sh/1g12t8qkxgs0psx/AABKzAY-EcSXFwfo6Q0cSa1la?dl=0)
