[![docs](https://img.shields.io/badge/seq--collection-Documentation-blue)](https://www.danielecook.com/seq-collection)
![Build](https://github.com/danielecook/seq-collection/workflows/Build/badge.svg)

# seq-collection

A collection of useful sequencing utilities.

__[Documentation](https://www.danielecook.com/seq-collection/)__

## Development

As I have done in the past, I intend to develop `seq-kit` over a period of months or years. In general, this means I'll add a new tool or utility when it becomes apparent that I need it.

However, I am open to PRs, requests, and feedback. Please let me know what you think.

I intend to port some commands over from [VCF-kit](https://github.com/AndersenLab/VCF-kit), but will broaden the applicability to FASTQs, BAMs, and perhaps even other items.

## Commands

### FASTA

* __fa-gc__ - Calculates GC content surrounding a given genomic position at specified window sizes.

### FASTQ

* __fq-dedup__ - The `fq-dedup` command de-duplicates a FASTQ by read ID. 
* __fq-count__ - Count the number of reads in a FASTQ and other metrics.
* __fq-meta__ - Provides basic information and a best guess as to sequencer from FASTQs.

### BAM

* __insert-size__ - Calculates insert-size quickly.
* __json__ - Convert VCF output to JSON.
* __sample__ - Samples variants from a VCF.

### Multi

* __iter__ - generates genomic ranges for parallelizing commands.
* __rand__ - Generates random sites and regions from FASTA, BAM, and VCF
