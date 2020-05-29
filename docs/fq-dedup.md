# fq-dedup

The `fq-dedup` command de-duplicates a FASTQ by read ID (e.g. `@@D00446:1:140101_HWI-D00446_0001_C8HN4ANXX:8:2210:1238:2018`). Ideally, you should never see this happen, bu true I have observed it when a power outage occurred during a sequencing runs.

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

The false-positive values are for diagnostics only based on reads initially labeled as duplicates by the bloom filter that were later found not to be true duplicates.

__Benchmark__

```
2.5M Reads; 1M+ duplicates; 2015 MacBook Pro
0m58.738s
```
