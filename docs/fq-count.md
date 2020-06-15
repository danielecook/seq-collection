# fq-count

Count the number of reads in a FASTQ and calculate additional metrics.

```
Counts lines in a FASTQ

Usage:
  fq-count [options] [fastq ...]

Arguments:
  [fastq ...]      Input FASTQ

Options:
  -t, --header               Output the header
  -b, --basename             Add basename column
  -a, --absolute             Add column for absolute path
  -h, --help                 Show this help
```

#### Example

```bash
sc fq-count --header -b *.fq
```

|   reads |   gc_content |   gc_bases |   n_bases |   bases | basename              |
|--------:|-------------:|-----------:|----------:|--------:|:----------------------|
|       8 |     0.53125  |         17 |         0 |      32 | dup.fq                |
|       8 |     0.53125  |         17 |         0 |      32 | dup.fq.gz             |
|       1 |     0.35     |         21 |         0 |      60 | illumina_1.fq         |
|       1 |     0.35     |         21 |         0 |      60 | illumina_2.fq         |
|       1 |     1        |        101 |         0 |     101 | illumina_2000_2500.fq |
|       6 |     0.35     |        126 |         0 |     360 | illumina_3.fq         |
|       1 |     1        |        101 |         0 |     101 | illumina_3000_4000.fq |
|       1 |     0.35     |         21 |         0 |      60 | illumina_4.fq         |
|       1 |     0.35     |         21 |         0 |      60 | illumina_6.fq         |
|       1 |     0.35     |         21 |         0 |      60 | illumina_7.fq         |
|       2 |     0.333333 |         14 |         0 |      42 | illumina_8.fq         |
|       1 |     0.35     |         21 |         0 |      60 | illumina_hiseq_x.fq   |
|       4 |     0.5      |          8 |         0 |      16 | nodup.fq              |
|       9 |     0        |          0 |         0 |       9 | novaseq.fq            |
|       2 |     0.430556 |         62 |         0 |     144 | sra.fq                |
