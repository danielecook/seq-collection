![Build](https://github.com/danielecook/seq-collection/workflows/Build/badge.svg)

# seq-collection

A collection of useful sequencing utilities.

```
Sequence data utilities (Version 0.0.2)

Usage:
  sc [options] COMMAND

Commands:

  VCF:

    json             Convert a VCF to JSON
    tsv              Convert a VCF to TSV

  FASTQ:

    fq-meta          Output metadata for FASTQ
    fq-count         Counts lines in a FASTQ
    fq-dedup         Removes exact duplicates from FASTQ Files

  MULTI:

    iter             Iterate genomic ranges from a BAM or VCF for parallel execution

  BAM:

    insert-size      Calculate insert-size metrics

Options:
  --debug                    Debug
  -h, --help                 Show this help
  ```

# Command Documentation


## Usage

Install using `nim c sc.nim`

## Development

As I have done in the past, I intend to develop `seq-kit` over a period of months or years. In general, this means I'll add a new tool or utility when it becomes apparent that I need it.

However, I am open to PRs, requests, and feedback. Please let me know what you think.

I intend to port some commands over from [VCF-kit](https://github.com/AndersenLab/VCF-kit), but will broaden the applicability to FASTQs, BAMs, and perhaps even other items.
 
## Tools

### FASTA

#### fa-gc

Use `fa-gc` to calculate the GC content surrounding a given genomic position.

The window sizes correspond to the distance upstream and downstream used to calculate GC content for. For example:

* 50 --> "100bp" window, but a 101bp width is examined.


Example output:

```
> sc fa-gc --pos pos.tsv hg19.fa.gz 50 3200 500000
```

| chrom   |   pos |   gc_100 |   gc_6400 |   gc_1000000 |
|:--------|------:|---------:|----------:|-------------:|
| chr1    | 50000 |   0.3663 |  0.373535 |     0.430292 |
| chr1    | 50020 |   0.3663 |  0.373535 |     0.430288 |
| chr1    | 50050 |   0.3564 |  0.373067 |     0.430283 |

__Notes__

* GC content calculations omit `N` bases.
* Output is sorted by chromosome, position.

### FASTQ

#### fq-dedup

The `fq-dedup` command de-duplicates a FASTQ by read ID (e.g. `@@D00446:1:140101_HWI-D00446_0001_C8HN4ANXX:8:2210:1238:2018`). Ideally, this should never happen!

The command uses a [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filter) to identify duplicates, and has to read through the file twice, and output the original FASTQ.

```bash
sc fq-dedup myfastq.fq.gz 2> dup.err | gzip > dedupped.fq.gz
```

fq-dedup can read both `.fq.gz` and `.fq` files. It sends the deduplicated FASTQ to stdout.

__Output__

Once complete, the following is sent to stderr:

```
total_reads: 2500000
duplicates 1086043
false-positive: 0
false-positive-rate: 0.0
```

The false-positive values are for diagnostics only based on reads initially labeled as duplicates by the bloom filter that were later found not to be.

__Benchmark__

```
2.5M Reads; 1M+ duplicates; 2015 MacBook Pro
0m58.738s
```

#### fq-count

Count the number of reads in a FASTQ and other metrics.

#### fq-meta

__Scenario:__ You are given an old dusty hard drive packed with sequence data. Your collaborator says "We have some great sequencing data here, if only someone could analyze it." You peek at the filesystem and discover that FASTQs have been renamed, removing crucial information about how they were generated. Your collaborator, however, recalls certain details about which data was sequenced on which sequencer and he has a list of sequencing barcodes and associated samples that you can match on." If only there was a way to determine the barcodes, sequencer, or other essential metadata for each FASTQ...

If this scenerio sounds familiar, `fq-meta`™ is for you. It samples the first few sequences of a FASTQ and outputs a summary for every FASTQ including:

* A best guess as to the type of sequencer based on instrument and flowcell information
* A best guess for the quality score format
* The barcode/index used to multiplex the sample which can be useful for identifying samples
* The machine, run, and other metadata about the FASTQ
* Information derived from the filename
* The file location

Data is output in a tab-delimited format so that this tool can be used to assemble a database of FASTQs and their associated information. You can do so like this:

```bash
sc fq-meta --header > my_fq_database.txt # Use this to output just the variable names
sc fq-meta *.fq.gz >> my_fq_database.txt
```

The resulting dataset can be combined with other metadata and filtered to select samples for processing in a pipeline.

You can also parallelize the operation with [GNU-parallel](https://www.gnu.org/software/parallel/).

```bash
sc fq-meta --header > my_fq_database.txt # Use this to output just the variable names
parallel -j 8 sc fq-meta ::: find . -name '*.fq.gz' >> my_fq_database.txt
```

__Example output__

| machine | sequencer      | prob_sequencer        | flowcell  | flowcell_description              | run | lane | sequence_id | index1   | index2 | qual_format          | qual_phred | qual_multiple | min_qual | max_qual | n_lines | basename              | absolute_path |
|---------|----------------|-----------------------|-----------|-----------------------------------|-----|------|-------------|----------|--------|----------------------|------------|---------------|----------|----------|---------|-----------------------|---------------|
| D00446  | HiSeq2000/2500 | high:machine+flowcell | C8HN4ANXX | High Output (8-lane) v4 flow cell | 1   | 8    |             | GCTCGGTA |        | Sanger;Illumina 1.8+ | Phred+33   | TRUE          | 14       | 14       | 1       | illumina_2000_2500.fq | …             |
| K00100  | HiSeq3000/4000 | high:machine+flowcell | H300JBBXX | (8-lane) v1 flow cell             | 33  | 6    |             | GCCAAT   |        | Sanger;Illumina 1.8+ | Phred+33   | TRUE          | 14       | 14       | 1       | illumina_3000_4000.fq | …             |
| D00209  | HiSeq2000/2500 | high:machine+flowcell | CACDKANXX | High Output (8-lane) v4 flow cell | 258 | 6    |             | CGCAGTT  |        | Sanger;Illumina 1.8+ | Phred+33   | TRUE          | 0        | 37       | 1       | illumina_6.fq         | …             |
| D00209  | HiSeq2000/2500 | high:machine+flowcell | CACDKANXX | High Output (8-lane) v4 flow cell | 258 | 6    |             | GAGCAAG  |        | Sanger;Illumina 1.8+ | Phred+33   | TRUE          | 0        | 37       | 1       | illumina_7.fq         | …             |

__Input____

`fq-meta` accepts both gzipped FASTQs (`.fq.gz`, `.fastq.gz` ~ inferred from `.gz` extension) and raw text FASTQs.

##### Assembling an FQ-Database

```
sc fq-meta --header > fastq_db.tsv # Prints just the header
sc fq-meta sample_1_R1.fq.gz sample_1_R2.fq.gz sample_2_R1.fq.gz sample_2_R2.fq.gz >> fastq_db.tsv
```

### BAM

#### insert-size

Calculate the insert-size of a bam or a set of bams. Bams are estimated by evaluating up to the 99.5th percentile of read insert-sizes. This gives numbers that are very close to Picard but a lot faster. 

```bash
sc insert-size --header input.bam # One bam
sc insert-size --header *.bam # Multiple bams
```

__Options__

* `--verbose` Output information about progress.
* `--dist` Output the frequency distribution of insert sizes.

__Output__

|   median |   mean |   std_dev |   min |   percentile_99.5 |   max_all |   n_reads |   n_accept |   n_use | sample   |
|---------:|-------:|----------:|------:|------------------:|----------:|----------:|-----------:|--------:|:---------|
|      179 |  176.5 |    63.954 |    38 |               358 |       359 |       237 |        101 |     100 | AB1      |


Calculate insert-size metrics on a set of bams.

### VCF

#### json (VCF to JSON conversion)

Convert a VCF to JSON. This is most useful in the context of a web-service, where you can serve variant data for the purposes of visualization, browsing, etc. For an example of what this might look like, see the [elegans variation variant browser](https://elegansvariation.org/data/browser/).

The `json` command can be used to parse `ANN` columns (effect annotations) by specifying `--annotation`. Additionally, FORMAT fields (GT, DP, etc) can be zipped with sample names.

```bash
Usage:
  sc json [options] vcf region

Arguments:
  vcf              VCF to convert to JSON
  region           Region

Options:
  -i, --info=INFO            comma-delimited INFO fields; Use 'ALL' for everything
  -f, --format=FORMAT        comma-delimited FORMAT fields; Use 'ALL' for everything
  -s, --samples=SAMPLES      Set Samples (default: ALL)
  -p, --pretty               Prettify result
  -a, --array                Output as a JSON array instead of ind. JSON lines
  -z, --zip                  Zip sample names with FORMAT fields (e.g. {'sample1': 25, 'sample2': 34})
  -n, --annotation           Parse ANN Fields
  --pass                     Only output variants where FILTER=PASS
  --debug                    Debug
  -h, --help                 Show this help
```

__Example__

```
sc json --format=GT --zip --pretty tests/data/test.vcf.gz
```

Outputs:

```json
{
  "CHROM": "I",
  "POS": 41947,
  "ID": ".",
  "REF": "A",
  "ALT": [
    "T"
  ],
  "QUAL": 999.0,
  "FILTER": [
    "PASS"
  ],
  "FORMAT": {
    "GT": {
      "AB1": [
        2,
        2
      ],
      "...": [
        0,
        0
      ]
    }
  }
}
```

You can also specify custom `SGT` and `TGT` output formats which transform `GT` fields from their
typical output.

* `GT` - Outputs genotypes as [[0, 0], [0, 1], [1, 1], ...
* `SGT` - Outputs genotypes as `0/0`, `0/1`, `1/1`, ...
* `TGT` - Outputs genotypes as `T/T`, `A/T`, `A/A`, ...

### Multi

#### iter

The `iter` command operates on BAM/CRAM and VCF/BCF files, and is used to generate genomic ranges that can be used to process genomic data in chunks. It works well with tools such as `xargs` or [gnu-parallel](https://www.gnu.org/software/parallel/).

__Example__

```bash
sc iter test.bam 100,000 # Iterate on bins of 100k base pairs
sc iter test.bam 100000 # Also valid
sc iter test.bam 1e6 # Also valid

# Outputs
> I:0-999999
> I:1000000-1999999
> I:2000000-2999999
> I:3000000-3999999
> I:4000000-4999999
```

This list of genomic ranges can be used to process a BAM or VCF in parallel:

```bash

function process_chunk {
  # Code to process chunk
  vcf=$1
  region=$2
  # e.g. bcftools call -m --region 
  echo bcftools call --region $region $vcf # ...
}

# Export the function to make it available to GNU parallel
export -f process_chunk

parallel --verbose process_chunk ::: test.bam ::: $(sc iter test.bam)

```

You can also set the `[width]` option to 0 to generate a list of chromosomes.

