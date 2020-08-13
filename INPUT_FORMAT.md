## Input formats

### `BED`

The BED format doc can be found at [UCSC genome browser
website](http://genome.ucsc.edu/FAQ/FAQformat#format1).

The required columns as input for fseq2:
- 1.chrom
- 2.start
- 3.end 
- 6.strand


### `BAMPE` and `BEDPE`
The BAMPE format is a BAM format containing paired-end information.  
The BEDPE format doc can be found at [BEDTools general usage](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).

By default, fseq2 cannot detect if BAM or BED is paired-end data, and treat them as single-end data. Only with flag (`-pe`) on, 
BAM or BED is recognized as paired-end data.  
Paired-end mode overrides `-f` and set it to 0.
 
> Note on SE vs. PE  
SE: fseq2 use 5' end to model distribution  
PE: fseq2 use middle point to model distribution


If format is BAMPE, fseq2 requires BAM sorted by name
```
$ samtools sort -n in.bam -o out.bam
```

PE data has advantage of knowing fragment length. This can be useful for filtering for desired fragments with certain length, 
e.g. in ATAC-seq, filter out fragments with length >= 150 to avoid calling peaks on nucleosomes.



