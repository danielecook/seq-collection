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
  -d, --dist=DIST            Output regions following a distribution ex: N(1,5) [see docs] (default: 0)
  -p, --pattern=PATTERN      A regular expression to use for chromosomes to keep
  -1, --one                  Output 1-based coordinates
  -h, --help                 Show this help

```

The `rand` command is useful for generating random sites and ranges. You can apply the `rand` command to `fasta`, `bam`, or `vcf` files.

If applied to a `fasta` file, the `rand` command will output the sequences.

### Options

#### `--pattern`

Use `--pattern` to restrict chromosomes in output. For example, the human genome contains additional contigs for alternative haplotypes, or unplaced contigs (e.g. `chrUn_gl000216`). See this [UCSC FAQ](http://genome.ucsc.edu/FAQ/FAQdownloads#download10) for more information.

```bash
sc rand --pattern "chr[0-9MXY]+$" hg19.fasta # only outputs chr1-22,X,Y,M
```

### Examples

Sample regions with a normal distribution of \\(N(\mu=5,\sigma=3)\\).

```bash
sc rand --dist=5,3 hg19.fasta
```

Sample regions with a uniform distribution of \\(U(min=2,max=10))

```bash
sc rand --dist=2-10 hg19.fasta
```

```
chr4	43446646	43446653	ATCTTTTA
chr2	146664347	146664354	aagattga
chr14	29047935	29047937	tca
...
```

Sample

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