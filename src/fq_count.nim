import streams
import zip/gzipfiles
import utils/helpers
import strutils
import os

const fq_count_header* = ["reads",
                          "gc_content",
                          "gc_bases",
                          "n_bases",
                          "bases"].join("\t")


proc fq_count*(fastq: string, basename: bool, absolute: bool) =
    #[
        Deduplicates reads by ID in FASTQs

        Based on Brent Pedersens implementation:
        https://gist.github.com/brentp/640806
    ]#

    var 
        i = 0
        gc_cnt: int64
        n_cnt: int64
        total_len: int64
        n_reads: int
        #n_dups: int

    let stream: Stream =
        if fastq[^3 .. ^1] == ".gz":
            newGZFileStream(fastq)
        else:
            newFileStream(fastq, fmRead)
    if stream == nil:
        quit_error("Unable to open file: " & fastq, 2)

    for line in lines(stream):
        i.inc()
        if (i mod 4) == 1:
            n_reads.inc()
        if (i mod 4) == 2:
            gc_cnt.inc(line.count("G") + line.count("C"))
            n_cnt.inc(line.count("N"))
            total_len.inc(line.len)
    
    var output = [$n_reads,
                $(gc_cnt.float / (total_len - n_cnt).float),
                $gc_cnt,
                $n_cnt,
                $total_len]

    echo output_w_fnames(output.join("\t"), fastq, basename, absolute)