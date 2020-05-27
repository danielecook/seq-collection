# rand

```
rand

Generate random genomic positions and ranges

Usage:
  rand [options] input

Arguments:
  input            Input FASTA, BAM, or VCF or BAM

Options:
  -n, --sites=SITES          Number of sites (default: 10)
  -b, --bed=BED              BED (0-based) of regions to restrict to
  -s, --seq=SEQ              Output additional information
  -1, --one                  Output 1-based coordinates
  -h, --help                 Show this help
```

### 0 vs. 1 based coordinates

`rand` automatically outputs coordinates that correspond to the input format. So:

* fasta - 0-based
* bam - 0-based
* vcf - 1-based

You can force fasta / bam output to be 1-based using the `-1` flag.

### Options

* `n`, `--sites`
* `b`, `--bed`
* `-1` - Output 1-based coordinates instead of 0-based coordinates