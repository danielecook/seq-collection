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
    while not stream.atEnd():
        record = stream.readLine()
        #record = stream.fq_string()
        
        if bloom.contains(record):
            check.inc(record, 1)
        bloom.incl(record)
        discard stream.readLine()
        discard stream.readLine()
        discard stream.readLine()
        i.inc()
        if i mod 1000 == 0:
            stderr.writeLine($i)

    stream.setPosition(0)
    while not stream.atEnd():
        record = stream.readLine()
        if record in check == false:
            echo record
            echo stream.readLine()
            echo stream.readLine()
            echo stream.readLine()
            continue
        
        putative_false_positives.inc(record, 1)
        if putative_false_positives[record] > 1:
            continue
        echo record
        echo stream.readLine()
        echo stream.readLine()
        echo stream.readLine()

    stream.close()