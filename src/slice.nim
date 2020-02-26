import hts
import streams
import zip/gzipfiles
import strutils
import deques
import strformat
import tables
import sequtils
import logging
import math
import utils/helpers

proc cmd_slice*(bamfile: string, pos_file: string) =
    #[
        Slice a bam cleanly, returning only reads found within a set of bounds
    ]#
    var 
        b: Bam
        bed_ln: region_t
       
    # Parse genomic positions file
    let stream: Stream =
        if pos_file[^3 .. ^1] == ".gz":
            newGZFileStream(pos_file)
        else:
            newFileStream(pos_file, fmRead)

    open(b, bamfile, index=true)
    for line in lines(stream):
        bed_ln = bed_line_to_region(line)
        for record in b.query(bed_ln.chrom, bed_ln.start, bed_ln.stop):
            echo record.tostring()
