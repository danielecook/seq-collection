# Insert-Size

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

__Benchmark__

Below is a comparison of seq-collection with picard:

![](https://www.danielecook.com/insert-size-benchmark.png)

The results are also very very close:

![](https://www.danielecook.com/insert_size_compare.png)