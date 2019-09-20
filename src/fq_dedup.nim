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
import utils/bloom


proc fq_string(s: Stream): string =
    var record = s.readLine()
    return record & "\n" & s.readLine() & "\n" & s.readLine() & "\n" & s.readLine() 

proc fq_dedup*(fastq: string) =

    var 
        i = 0
        record: string
        bloom: BloomFilter
        false_positives: int
        check = initCountTable[string]()
        putative_false_positives = initCountTable[string]()

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
            if bloom.contains(record):
                check.inc(record, 1)
            bloom.incl(record)
        i.inc()
        if i mod 10000 == 0:
            stderr.writeLine $i

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
                    continue
                echo record
                write_ln = true
        elif write_ln:
            echo record


    stream.close()