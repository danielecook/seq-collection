# seq-collection

A collection of useful sequencing utilities.

## Usage

Install using `nim c sk.nim`

You will need to install a few dependencies. I will eventually create a nimble installer.

## Development

As I have done in the past, I intend to develop `seq-kit` over a period of months or years. In general, this means I'll add a new tool or utility when it becomes apparent that I need it.

However, I am open to PRs, requests, and feedback. Please let me know what you think.

I intend to port some commands over from [VCF-kit](https://github.com/AndersenLab/VCF-kit), but will broaden the applicability to FASTQs, BAMs, and perhaps even other items.

## Tools

### json (VCF to JSON conversion)

Convert a VCF to JSON. This is most useful in the context of a web-service, where you can serve variant data for the purposes of visualization, browsing, etc. For an example of what this might look like, see the [elegans variation variant browser](https://elegansvariation.org/data/browser/).

The `json` command can be used to parse `ANN` columns (effect annotations) by specifying `--annotation`. Additionally, FORMAT fields (GT, DP, etc) can be zipped with sample names.

```
Usage:
  sk json [options] vcf region

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
