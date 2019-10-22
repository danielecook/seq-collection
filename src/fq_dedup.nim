import tables
import sequtils
import strutils
import utils/gz
import os
import re
import sets
import streams
import zip/gzipfiles
import utils/helpers
import bitvector
import bloom


type
    H = uint64

proc fq_dedup*(fastq: string) =
    #[
        Deduplicates reads by ID in FASTQs

        Based on Brent Pedersens implementation:
        https://gist.github.com/brentp/640806
    ]#

    var 
        i = 0
        n_reads: int
        n_dups: int
        check = initCountTable[string]()
        putative_false_positives = initCountTable[string]()

    var bloom = newBloomFilter[uint32](1e8.int, 0.0001, 0, true)

    let stream: Stream =
        if fastq[^3 .. ^1] == ".gz":
            newGZFileStream(fastq)
        else:
            newFileStream(fastq, fmRead)
    if stream == nil:
        quit_error("Unable to open file: " & fastq, 2)
    
    # Iterate once to place dups in.
    for record in stream.lines:
        if i mod 4 == 0:
            if bloom.lookup($record):
                check.inc(record, 1)
            bloom.insert($record)
        i.inc()

    n_reads = i div 4.int

    if check.len == 0:
        stderr.writeLine("No Duplicates Found")
        stderr.writeLine("Copying fq to stdout")

    stream.setPosition(0)
    var write_ln = true
    i = 0
    for record in stream.lines:
        i.inc()
        if (i-1) mod 4 == 0:
            if record in check == false:
                echo record
                write_ln = true
                continue
            else:
                putative_false_positives.inc(record, 1)
                if putative_false_positives[record] > 1:
                    write_ln = false
                    n_dups.inc(1)
                    continue
                echo record
                write_ln = true
        elif write_ln:
            echo record
    
    stderr.writeLine("total_reads: ", n_reads)
    stderr.writeLine("duplicates ", n_dups)
    var fp = 0
    for v in putative_false_positives.values():
        if v == 1:
            fp.inc(1)
    stderr.writeLine("false-positive: ", fp)
    stderr.writeLine("false-positive-rate: ", fp.float / n_dups.float)


    stream.close()