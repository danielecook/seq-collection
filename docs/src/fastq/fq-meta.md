# fq-meta


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