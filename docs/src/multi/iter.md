
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

# iter
