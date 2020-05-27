# fa-gc

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