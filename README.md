# seq-collection

[![Build Status](https://travis-ci.org/danielecook/seq-collection.svg?branch=development)](https://travis-ci.org/danielecook/seq-collection)

A collection of useful sequencing utilities.

## Usage

Install using `nim c sc.nim`

You will need to install a few dependencies. I will eventually create a nimble installer.

## Development

As I have done in the past, I intend to develop `seq-kit` over a period of months or years. In general, this means I'll add a new tool or utility when it becomes apparent that I need it.

However, I am open to PRs, requests, and feedback. Please let me know what you think.

I intend to port some commands over from [VCF-kit](https://github.com/AndersenLab/VCF-kit), but will broaden the applicability to FASTQs, BAMs, and perhaps even other items.
 
## Tools

### fq-meta

__Scenario:__ - You are given an old dusty hard drive packed with sequence data. Your collaborator says "We have some great sequencing data here, if only someone could analyze it." You peek at the filesystem and discover that FASTQs have been renamed, removing crucial information about how they were generated. Your collaborator, however, recalls certain details about which data was sequenced on which sequencer and he has a list of sequencing barcodes and associated samples that you can match on." If only there was a way to determine the barcodes, sequencer, or other essential metadata for each FASTQ...

If this scenerio sounds familiar, `fq-meta`™ is for you. It samples the first few sequences of a FASTQ and outputs a summary including:

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

You can also parallelize the operation with (GNU-parallel)[https://www.gnu.org/software/parallel/].

```bash
sc fq-meta --header > my_fq_database.txt # Use this to output just the variable names
parallel -j 8 sc fq-meta ::: find . -name '*.fq.gz' >> my_fq_database.txt
```

#### Example Output



#### Assembling an FQ-Database

```
sc fq-meta --header > fastq_db.tsv # Prints just the header
sc fq-meta sample_1_R1.fq.gz sample_1_R2.fq.gz sample_2_R1.fq.gz sample_2_R2.fq.gz >> fastq_db.tsv
```

### json (VCF to JSON conversion)

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
  --debug                    Debug
  -h, --help                 Show this help
```

