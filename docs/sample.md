# Sample

Randomly sample sites from a VCF

```
sample

Randomly sample a VCF

Usage:
  sample [options] vcf

Arguments:
  vcf              Variant file

Options:
  --bed=BED                  A set of bed regions to restrict sampling to
  -t, --types=TYPES          Variant types to sample (all,snps,mnps,indels (default: all)
  -n, --sites=SITES          Number of sites to sample (default: 10)
  -h, --help                 Show this help
```

The `sample` command selects random variants. The algorithm works as follows.

1. Chromosomes are weighted by length and randomly selected by their weights.
2. A position is randomly selected.
3. 1kb upstream is searched for variants and returned. If not, the process is repeated until the desired number of sites is returned.
4. A bloom filter is used to reduce the likelihood of duplicates. The key used is the chromosome + position of each variant which prevents sampling of additional alleles located at the same position in subsequent draws.

## Options

### `--types`

Select only sites that are `snps`, `mnps` (Multi-nucleotide polymorphism), or `indels`.