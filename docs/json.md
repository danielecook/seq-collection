# JSON

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
